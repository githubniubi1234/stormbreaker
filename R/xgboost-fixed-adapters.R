sb_xgboost_binary_learner_ids <- function() {
  c(
    "xgboost_binary",
    "xgboost_binary_dart",
    "xgboost_binary_linear"
  )
}

sb_xgboost_cox_learner_ids <- function() {
  c(
    "xgboost_cox",
    "xgboost_cox_dart",
    "xgboost_cox_linear"
  )
}

sb_xgboost_aft_learner_ids <- function() {
  c(
    "xgboost_aft",
    "xgboost_aft_dart",
    "xgboost_aft_linear"
  )
}

sb_xgboost_survival_learner_ids <- function() {
  c(sb_xgboost_cox_learner_ids(), sb_xgboost_aft_learner_ids())
}

sb_xgboost_learner_ids <- function() {
  c(sb_xgboost_binary_learner_ids(), sb_xgboost_survival_learner_ids())
}

sb_xgboost_booster <- function(learner) {
  if (learner %in% c(
    "xgboost_binary_dart",
    "xgboost_cox_dart",
    "xgboost_aft_dart"
  )) {
    return("dart")
  }
  if (learner %in% c(
    "xgboost_binary_linear",
    "xgboost_cox_linear",
    "xgboost_aft_linear"
  )) {
    return("gblinear")
  }
  "gbtree"
}

sb_xgboost_available <- function() {
  requireNamespace("xgboost", quietly = TRUE)
}

sb_xgboost_selector_gain <- function(model) {
  importance <- xgboost::xgb.importance(feature_names = model$features, model = model$booster)
  score <- numeric(length(model$features))
  if (nrow(importance)) score[match(importance$Feature, model$features)] <- importance$Gain
  score
}

sb_xgboost_params <- function(learner) {
  if (!learner %in% sb_xgboost_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed XGBoost learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  objective <- if (learner %in% sb_xgboost_binary_learner_ids()) {
    "binary:logistic"
  } else if (learner %in% sb_xgboost_cox_learner_ids()) {
    "survival:cox"
  } else {
    "survival:aft"
  }
  booster <- sb_xgboost_booster(learner)
  if (identical(booster, "gblinear")) {
    params <- list(
      nrounds = 100L,
      booster = booster,
      eta = 0.1,
      lambda = 1,
      alpha = 0,
      updater = "coord_descent",
      feature_selector = "cyclic",
      nthread = 1L,
      seed = 1L,
      objective = objective
    )
    if (learner %in% sb_xgboost_aft_learner_ids()) {
      params <- c(params, list(
        aft_loss_distribution = "normal",
        aft_loss_distribution_scale = 1
      ))
    }
    return(params)
  }
  tree_params <- list(
    nrounds = 100L,
    eta = 0.1,
    max_depth = 3L,
    min_child_weight = 1,
    subsample = 1,
    colsample_bytree = 1
  )
  if (identical(booster, "dart")) {
    params <- c(
      tree_params["nrounds"],
      list(booster = booster),
      tree_params[-1L],
      list(
        sample_type = "uniform",
        normalize_type = "tree",
        rate_drop = 0.1,
        skip_drop = 0,
        one_drop = FALSE,
        nthread = 1L,
        seed = 1L,
        objective = objective
      )
    )
    if (learner %in% sb_xgboost_aft_learner_ids()) {
      params <- c(params, list(
        aft_loss_distribution = "normal",
        aft_loss_distribution_scale = 1
      ))
    }
    return(params)
  }
  params <- c(
    tree_params,
    list(
      nthread = 1L,
      seed = 1L,
      objective = objective
    )
  )
  if (learner %in% sb_xgboost_aft_learner_ids()) {
    params <- c(params, list(
      aft_loss_distribution = "normal",
      aft_loss_distribution_scale = 1
    ))
  }
  params
}

sb_validate_xgboost_params <- function(learner, params) {
  expected <- names(sb_xgboost_params(learner))
  valid <- is.list(params) && identical(names(params), expected) &&
    is.numeric(params$nrounds) && length(params$nrounds) == 1L &&
    !is.na(params$nrounds) && params$nrounds >= 1 &&
    identical(params$nrounds, as.integer(params$nrounds)) &&
    is.numeric(params$nthread) && length(params$nthread) == 1L &&
    !is.na(params$nthread) && params$nthread >= 1 &&
    identical(params$nthread, as.integer(params$nthread)) &&
    is.numeric(params$seed) && length(params$seed) == 1L &&
    !is.na(params$seed) && params$seed >= 1 &&
    identical(params$seed, as.integer(params$seed)) &&
    is.character(params$objective) && length(params$objective) == 1L &&
    identical(
      params$objective,
      sb_xgboost_params(learner)$objective
    )
  booster <- sb_xgboost_booster(learner)
  if (booster %in% c("gbtree", "dart")) {
    numeric_controls <- unlist(params[c(
      "eta", "max_depth", "min_child_weight", "subsample",
      "colsample_bytree"
    )], use.names = FALSE)
    valid <- valid && !anyNA(numeric_controls) &&
      all(is.finite(numeric_controls)) && params$eta > 0 &&
      params$max_depth >= 1 && params$min_child_weight >= 0 &&
      params$subsample > 0 && params$subsample <= 1 &&
      params$colsample_bytree > 0 && params$colsample_bytree <= 1
  }
  if (identical(booster, "dart")) {
    dropout_controls <- unlist(
      params[c("rate_drop", "skip_drop")],
      use.names = FALSE
    )
    valid <- valid && identical(params$booster, booster) &&
      !anyNA(dropout_controls) && all(is.finite(dropout_controls)) &&
      all(dropout_controls >= 0 & dropout_controls <= 1) &&
      identical(params$sample_type, "uniform") &&
      identical(params$normalize_type, "tree") &&
      is.logical(params$one_drop) && length(params$one_drop) == 1L &&
      !is.na(params$one_drop)
  }
  if (identical(booster, "gblinear")) {
    linear_controls <- unlist(
      params[c("eta", "lambda", "alpha")],
      use.names = FALSE
    )
    valid <- valid && identical(params$booster, booster) &&
      !anyNA(linear_controls) && all(is.finite(linear_controls)) &&
      params$eta > 0 && params$lambda >= 0 && params$alpha >= 0 &&
      identical(params$updater, "coord_descent") &&
      identical(params$feature_selector, "cyclic")
  }
  if (learner %in% sb_xgboost_aft_learner_ids()) {
    valid <- valid &&
      identical(params$aft_loss_distribution, "normal") &&
      is.numeric(params$aft_loss_distribution_scale) &&
      length(params$aft_loss_distribution_scale) == 1L &&
      !is.na(params$aft_loss_distribution_scale) &&
      is.finite(params$aft_loss_distribution_scale) &&
      params$aft_loss_distribution_scale > 0
  }
  if (!valid) {
    sb_abort(
      "The fixed XGBoost parameter specification is invalid.",
      "stormbreaker_bad_xgboost_parameters"
    )
  }
  params$nrounds <- as.integer(params$nrounds)
  if (booster %in% c("gbtree", "dart")) {
    params$max_depth <- as.integer(params$max_depth)
  }
  params$nthread <- as.integer(params$nthread)
  params$seed <- as.integer(params$seed)
  params
}

sb_require_xgboost <- function() {
  if (!sb_xgboost_available()) {
    sb_abort(
      "Package `xgboost` is required for a fixed XGBoost learner.",
      "stormbreaker_missing_xgboost"
    )
  }
  invisible(TRUE)
}

sb_validate_xgboost_survival_inputs <- function(
  x, time, event, learner = "xgboost_cox"
) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  is_aft <- learner %in% sb_xgboost_aft_learner_ids()
  if (any(time <= 0)) {
    sb_abort(
      paste0(
        "A fixed XGBoost-", if (is_aft) "AFT" else "Cox",
        " analysis fold requires strictly positive times."
      ),
      if (is_aft) {
        "stormbreaker_xgboost_aft_nonpositive_time"
      } else {
        "stormbreaker_xgboost_cox_nonpositive_time"
      }
    )
  }
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }
  if (!any(event == 1L)) {
    sb_abort(
      paste0(
        "A fixed XGBoost-", if (is_aft) "AFT" else "Cox",
        " analysis fold must contain at least one event."
      ),
      if (is_aft) {
        "stormbreaker_xgboost_aft_no_event"
      } else {
        "stormbreaker_xgboost_cox_no_event"
      }
    )
  }
  list(x = x, time = time, event = as.integer(event))
}

sb_fit_xgboost <- function(x, endpoint, learner, params) {
  params <- sb_validate_xgboost_params(learner, params)
  if (learner %in% sb_xgboost_binary_learner_ids()) {
    x <- sb_validate_predictor_frame(x, "x")
    endpoint <- sb_validate_binary_outcome_vector(endpoint, nrow(x))
    if (!all(c(0L, 1L) %in% endpoint)) {
      sb_abort(
        "A fixed XGBoost binary analysis fold must contain both classes.",
        "stormbreaker_xgboost_binary_single_class"
      )
    }
    label <- as.numeric(endpoint)
  } else {
    inputs <- sb_validate_xgboost_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      learner
    )
    x <- inputs$x
    label <- if (learner %in% sb_xgboost_cox_learner_ids()) {
      ifelse(inputs$event == 1L, inputs$time, -inputs$time)
    } else {
      NULL
    }
  }
  sb_require_xgboost()

  dtrain <- if (learner %in% sb_xgboost_aft_learner_ids()) {
    matrix <- xgboost::xgb.DMatrix(
      as.matrix(x),
      nthread = params$nthread
    )
    xgboost::setinfo(matrix, "label_lower_bound", inputs$time)
    xgboost::setinfo(
      matrix,
      "label_upper_bound",
      ifelse(inputs$event == 1L, inputs$time, Inf)
    )
    matrix
  } else {
    xgboost::xgb.DMatrix(
      as.matrix(x),
      label = label,
      nthread = params$nthread
    )
  }
  booster_type <- sb_xgboost_booster(learner)
  engine_params <- if (identical(booster_type, "gblinear")) {
    list(
      objective = params$objective,
      booster = booster_type,
      eta = params$eta,
      lambda = params$lambda,
      alpha = params$alpha,
      updater = params$updater,
      feature_selector = params$feature_selector,
      nthread = params$nthread
    )
  } else {
    tree_params <- list(
      objective = params$objective,
      booster = booster_type,
      eta = params$eta,
      max_depth = params$max_depth,
      min_child_weight = params$min_child_weight,
      subsample = params$subsample,
      colsample_bytree = params$colsample_bytree,
      nthread = params$nthread
    )
    if (identical(booster_type, "dart")) {
      tree_params$sample_type <- params$sample_type
      tree_params$normalize_type <- params$normalize_type
      tree_params$rate_drop <- params$rate_drop
      tree_params$skip_drop <- params$skip_drop
      tree_params$one_drop <- params$one_drop
    }
    tree_params
  }
  if (learner %in% sb_xgboost_aft_learner_ids()) {
    engine_params$aft_loss_distribution <- params$aft_loss_distribution
    engine_params$aft_loss_distribution_scale <-
      params$aft_loss_distribution_scale
  }
  booster <- withr::with_seed(
    params$seed,
    xgboost::xgb.train(
      params = engine_params,
      data = dtrain,
      nrounds = params$nrounds,
      verbose = 0
    )
  )
  if (!inherits(booster, "xgb.Booster")) {
    sb_abort(
      "Fixed XGBoost did not return a usable booster.",
      "stormbreaker_incomplete_xgboost_fit"
    )
  }
  structure(
    list(
      learner = learner,
      features = names(x),
      params = params,
      booster = booster
    ),
    class = "stormbreaker_xgboost_fixed"
  )
}

sb_xgboost_prediction <- function(model, x) {
  if (!inherits(model, "stormbreaker_xgboost_fixed")) {
    sb_abort(
      "The supplied model is not a fitted fixed XGBoost learner.",
      "stormbreaker_bad_xgboost_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fixed XGBoost schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  sb_require_xgboost()
  dnew <- xgboost::xgb.DMatrix(
    as.matrix(x),
    nthread = model$params$nthread
  )
  value <- as.numeric(stats::predict(
    model$booster,
    dnew,
    outputmargin = model$learner %in% sb_xgboost_survival_learner_ids()
  ))
  if (model$learner %in% sb_xgboost_aft_learner_ids()) {
    value <- -value
  }
  valid <- length(value) == nrow(x) && !anyNA(value) && all(is.finite(value))
  if (model$learner %in% sb_xgboost_binary_learner_ids()) {
    valid <- valid && all(value >= 0 & value <= 1)
  }
  if (!valid) {
    sb_abort(
      "Fixed XGBoost produced invalid predictions.",
      "stormbreaker_invalid_xgboost_prediction"
    )
  }
  unname(value)
}

sb_predict_xgboost_cox_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  if (!model$learner %in% sb_xgboost_cox_learner_ids()) {
    sb_abort(
      "Survival probability requires the fixed XGBoost-Cox learner.",
      "stormbreaker_bad_xgboost_fit"
    )
  }
  inputs <- sb_validate_xgboost_survival_inputs(
    x_analysis, time, event, model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fixed XGBoost-Cox schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  analysis_risk <- sb_xgboost_prediction(model, inputs$x)
  assessment_risk <- sb_xgboost_prediction(model, x_assessment)
  baseline <- sb_breslow_baseline_hazard(
    inputs$time,
    inputs$event,
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    baseline,
    source = "xgboost_cox_breslow_calibration",
    ties = "breslow",
    risk = assessment_risk
  )
}

sb_xgboost_aft_survival_probability <- function(
  risk, prediction_horizon, scale
) {
  valid <- is.numeric(risk) && length(risk) > 0L && !anyNA(risk) &&
    all(is.finite(risk)) && is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L && !is.na(prediction_horizon) &&
    is.finite(prediction_horizon) && prediction_horizon > 0 &&
    is.numeric(scale) && length(scale) == 1L && !is.na(scale) &&
    is.finite(scale) && scale > 0
  if (!valid) {
    sb_abort(
      "Fixed XGBoost-AFT received invalid prediction inputs.",
      "stormbreaker_invalid_xgboost_aft_prediction"
    )
  }
  margin <- -risk
  probability <- stats::pnorm(
    (margin - log(prediction_horizon)) / scale
  )
  if (length(probability) != length(risk) || anyNA(probability) ||
    any(!is.finite(probability)) ||
    any(probability < 0 | probability > 1)) {
    sb_abort(
      "Fixed XGBoost-AFT produced invalid survival probabilities.",
      "stormbreaker_invalid_xgboost_aft_prediction"
    )
  }
  unname(probability)
}

sb_predict_xgboost_aft_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  if (!model$learner %in% sb_xgboost_aft_learner_ids()) {
    sb_abort(
      "Survival probability requires a fixed XGBoost-AFT learner.",
      "stormbreaker_bad_xgboost_fit"
    )
  }
  inputs <- sb_validate_xgboost_survival_inputs(
    x_analysis, time, event, model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fixed XGBoost-AFT schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  risk <- sb_xgboost_prediction(model, x_assessment)
  list(
    prediction_horizon = prediction_horizon,
    source = "xgboost_aft_normal_native",
    ties = "not_applicable",
    probability = sb_xgboost_aft_survival_probability(
      risk,
      prediction_horizon,
      model$params$aft_loss_distribution_scale
    )
  )
}

sb_xgboost_model_raw <- function(model) {
  sb_require_xgboost()
  raw <- xgboost::xgb.save.raw(model$booster, raw_format = "json")
  if (!is.raw(raw) || length(raw) == 0L) {
    sb_bad_learner_state("The fitted XGBoost booster could not be saved.")
  }
  raw
}

sb_extract_xgboost_learner_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_xgboost_fixed") ||
    !model$learner %in% sb_xgboost_learner_ids()) {
    sb_bad_learner_state("The fitted XGBoost learner state is invalid.")
  }
  baseline_hazard <- NULL
  if (model$learner %in% sb_xgboost_cox_learner_ids() &&
    "prediction_horizon" %in% names(context)) {
    probability <- sb_predict_xgboost_cox_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    baseline_hazard <- list(
      prediction_horizon = probability$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      cumulative_baseline_hazard = probability$cumulative_baseline_hazard,
      baseline_survival = probability$baseline_survival,
      source = probability$source,
      ties = probability$ties
    )
  }
  survival_probability <- NULL
  if (model$learner %in% sb_xgboost_aft_learner_ids() &&
    "prediction_horizon" %in% names(context)) {
    probability <- sb_predict_xgboost_aft_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x,
      context$prediction_horizon
    )
    survival_probability <- list(
      prediction_horizon = probability$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      source = probability$source,
      ties = probability$ties
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      model_raw = sb_xgboost_model_raw(model),
      risk_reference = if (model$learner %in% sb_xgboost_cox_learner_ids()) {
        "xgboost_native"
      } else if (model$learner %in% sb_xgboost_aft_learner_ids()) {
        "xgboost_aft_normal_native"
      } else {
        NULL
      },
      risk_formula = if (model$learner %in% sb_xgboost_cox_learner_ids()) {
        "xgboost_raw_margin"
      } else if (model$learner %in% sb_xgboost_aft_learner_ids()) {
        "negative_xgboost_aft_log_time_margin"
      } else {
        NULL
      },
      baseline_hazard = baseline_hazard,
      survival_probability = survival_probability
    ),
    class = c(
      "stormbreaker_xgboost_fixed_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_load_xgboost_state_model <- function(state, learner) {
  x <- state[["model_raw"]]
  if (!is.raw(x) || length(x) == 0L || !is.list(state[["params"]])) {
    sb_bad_learner_state("The saved XGBoost learner state is invalid.")
  }
  sb_require_xgboost()
  booster <- tryCatch(
    if ("as_booster" %in% names(formals(xgboost::xgb.load.raw))) {
      xgboost::xgb.load.raw(x, as_booster = TRUE)
    } else {
      xgboost::xgb.load.raw(x)
    },
    error = function(error) {
      sb_abort(
        paste0(
          "The saved XGBoost model could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_xgboost_state",
        parent = error
      )
    }
  )
  if (!inherits(booster, "xgb.Booster")) {
    sb_bad_learner_state("The saved XGBoost booster is invalid.")
  }
  structure(
    list(
      learner = learner,
      features = state[["features"]],
      params = state[["params"]],
      booster = booster
    ),
    class = "stormbreaker_xgboost_fixed"
  )
}

sb_predict_xgboost_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_xgboost_learner_ids()) {
    sb_bad_learner_state("The saved XGBoost learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_load_xgboost_state_model(state, learner)
  prediction <- sb_xgboost_prediction(model, x)
  if (learner %in% sb_xgboost_binary_learner_ids()) {
    return(list(probability = prediction))
  }

  result <- list(risk = prediction)
  baseline <- state[["baseline_hazard"]]
  if (!is.null(baseline)) {
    valid <- is.list(baseline) &&
      is.numeric(baseline$prediction_horizon) &&
      length(baseline$prediction_horizon) == 1L &&
      !is.na(baseline$prediction_horizon) &&
      is.finite(baseline$prediction_horizon) &&
      baseline$prediction_horizon > 0 &&
      is.character(baseline$time_unit) && length(baseline$time_unit) == 1L &&
      !is.na(baseline$time_unit) && nzchar(baseline$time_unit)
    if (!valid) {
      sb_bad_learner_state(
        "The saved XGBoost-Cox baseline-hazard state is invalid."
      )
    }
    result$prediction_horizon <- baseline$prediction_horizon
    result$time_unit <- baseline$time_unit
    result$survival_probability <- sb_survival_probability_from_hazard(
      prediction,
      baseline$cumulative_baseline_hazard
    )
  } else if (!is.null(state[["survival_probability"]])) {
    probability <- state[["survival_probability"]]
    valid <- learner %in% sb_xgboost_aft_learner_ids() &&
      is.list(probability) &&
      is.numeric(probability$prediction_horizon) &&
      length(probability$prediction_horizon) == 1L &&
      !is.na(probability$prediction_horizon) &&
      is.finite(probability$prediction_horizon) &&
      probability$prediction_horizon > 0 &&
      is.character(probability$time_unit) &&
      length(probability$time_unit) == 1L &&
      !is.na(probability$time_unit) && nzchar(probability$time_unit) &&
      identical(probability$source, "xgboost_aft_normal_native") &&
      identical(probability$ties, "not_applicable")
    if (!valid) {
      sb_bad_learner_state(
        "The saved XGBoost-AFT survival-probability state is invalid."
      )
    }
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- sb_xgboost_aft_survival_probability(
      prediction,
      probability$prediction_horizon,
      model$params$aft_loss_distribution_scale
    )
  }
  result
}

sb_xgboost_learner_adapter <- function(
  learner,
  params = sb_xgboost_params(learner)
) {
  params <- sb_validate_xgboost_params(learner, params)
  fit <- function(x, endpoint, event = NULL) {
    if (learner %in% sb_xgboost_binary_learner_ids()) {
      return(sb_fit_xgboost(x, endpoint, learner, params))
    }
    sb_fit_xgboost(
      x,
      list(time = endpoint, event = event),
      learner,
      params
    )
  }
  adapter <- list(
    id = learner,
    params = params,
    fit = fit,
    predict = sb_xgboost_prediction,
    extract_state = function(model, context = list()) {
      sb_extract_xgboost_learner_state(model, context)
    },
    predict_state = sb_predict_xgboost_learner_state
  )
  if (learner %in% sb_xgboost_cox_learner_ids()) {
    adapter$predict_survival <- sb_predict_xgboost_cox_survival
  } else if (learner %in% sb_xgboost_aft_learner_ids()) {
    adapter$predict_survival <- sb_predict_xgboost_aft_survival
  }
  adapter
}
