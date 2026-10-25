sb_catboost_aft_survival_params <- function(
    iterations = 100L,
    learning_rate = 0.1,
    depth = 6L,
    l2_leaf_reg = 3,
    seed = 1L) {
  list(
    iterations = sb_check_positive_integer(
      iterations,
      "catboost_iterations",
      minimum = 1L
    ),
    learning_rate = sb_validate_catboost_learning_rate(learning_rate),
    depth = sb_validate_catboost_depth(depth),
    l2_leaf_reg = sb_validate_catboost_l2_leaf_reg(l2_leaf_reg),
    random_seed = sb_check_positive_integer(
      seed,
      "catboost_seed",
      minimum = 1L
    ),
    loss_function = "SurvivalAft:dist=Normal;scale=1",
    eval_metric = "SurvivalAft:dist=Normal;scale=1",
    task_type = "CPU",
    thread_count = 1L,
    bootstrap_type = "Bayesian",
    bagging_temperature = 1,
    random_strength = 1,
    rsm = 1,
    grow_policy = "SymmetricTree",
    boosting_type = "Plain",
    use_best_model = FALSE,
    allow_writing_files = FALSE,
    save_snapshot = FALSE,
    logging_level = "Silent"
  )
}

sb_validate_catboost_aft_survival_inputs <- function(x, time, event) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time <= 0)) {
    sb_abort(
      paste0(
        "The fixed CatBoost AFT learner requires strictly positive ",
        "survival times."
      ),
      "stormbreaker_catboost_aft_nonpositive_time"
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
      "A fixed CatBoost AFT analysis fold must contain at least one event.",
      "stormbreaker_catboost_aft_no_event"
    )
  }
  list(x = x, time = time, event = as.integer(event))
}

sb_fit_catboost_aft_survival <- function(x, time, event, params) {
  inputs <- sb_validate_catboost_aft_survival_inputs(x, time, event)
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for the fixed CatBoost AFT learner.",
      "stormbreaker_missing_catboost"
    )
  }

  training_pool <- catboost::catboost.load_pool(
    data = as.matrix(inputs$x),
    label = cbind(
      inputs$time,
      ifelse(inputs$event == 1L, inputs$time, -1)
    )
  )
  booster <- catboost::catboost.train(
    learn_pool = training_pool,
    test_pool = NULL,
    params = params
  )
  tree_count <- as.integer(booster$tree_count)
  if (!identical(tree_count, params$iterations)) {
    sb_abort(
      "Fixed CatBoost AFT did not fit the declared number of iterations.",
      "stormbreaker_incomplete_catboost_fit"
    )
  }

  structure(
    list(
      features = names(inputs$x),
      params = params,
      iterations = params$iterations,
      distribution = "Normal",
      scale = 1,
      booster = booster
    ),
    class = "stormbreaker_catboost_aft_survival"
  )
}

sb_catboost_aft_raw_margin <- function(model, x) {
  if (!inherits(model, "stormbreaker_catboost_aft_survival")) {
    sb_abort(
      "The supplied model is not a fitted fixed CatBoost AFT learner.",
      "stormbreaker_bad_catboost_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fixed CatBoost AFT schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for fixed CatBoost AFT prediction.",
      "stormbreaker_missing_catboost"
    )
  }

  margin <- as.numeric(catboost::catboost.predict(
    model$booster,
    catboost::catboost.load_pool(as.matrix(x)),
    prediction_type = "RawFormulaVal",
    ntree_start = 0L,
    ntree_end = model$iterations,
    thread_count = 1L
  ))
  if (length(margin) != nrow(x) || anyNA(margin) || any(!is.finite(margin))) {
    sb_abort(
      "Fixed CatBoost AFT produced invalid log-time predictions.",
      "stormbreaker_invalid_catboost_prediction"
    )
  }
  unname(margin)
}

sb_predict_catboost_aft_survival <- function(model, x) {
  -sb_catboost_aft_raw_margin(model, x)
}

sb_catboost_aft_survival_probability <- function(
    risk, prediction_horizon, scale = 1) {
  valid <- is.numeric(risk) && length(risk) > 0L && !anyNA(risk) &&
    all(is.finite(risk)) && is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L && !is.na(prediction_horizon) &&
    is.finite(prediction_horizon) && prediction_horizon > 0 &&
    is.numeric(scale) && length(scale) == 1L && !is.na(scale) &&
    is.finite(scale) && scale > 0
  if (!valid) {
    sb_abort(
      "Fixed CatBoost AFT received invalid prediction inputs.",
      "stormbreaker_invalid_catboost_aft_prediction"
    )
  }
  probability <- stats::pnorm(
    (-risk - log(prediction_horizon)) / scale
  )
  if (anyNA(probability) || any(!is.finite(probability)) ||
    any(probability < 0 | probability > 1)) {
    sb_abort(
      "Fixed CatBoost AFT produced invalid survival probabilities.",
      "stormbreaker_invalid_catboost_aft_prediction"
    )
  }
  unname(probability)
}

sb_predict_catboost_aft_survival_survival <- function(
    model, x_analysis, time, event, x_assessment, prediction_horizon) {
  inputs <- sb_validate_catboost_aft_survival_inputs(
    x_analysis,
    time,
    event
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fixed CatBoost AFT schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  risk <- sb_predict_catboost_aft_survival(model, x_assessment)
  list(
    prediction_horizon = prediction_horizon,
    source = "catboost_aft_normal_native",
    ties = "not_applicable",
    probability = sb_catboost_aft_survival_probability(
      risk,
      prediction_horizon,
      model$scale
    )
  )
}

sb_extract_catboost_aft_survival_learner_state <- function(
    model, context = list()) {
  extracted <- sb_extract_catboost_tree_model(model)
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_catboost_aft_survival_survival(
      model,
      context[["x"]],
      context[["time"]],
      context[["event"]],
      context[["x"]][1L, , drop = FALSE],
      context[["prediction_horizon"]]
    )
    survival_probability <- list(
      prediction_horizon = probability$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context[["time_unit"]]
      } else {
        NA_character_
      },
      source = probability$source,
      ties = probability$ties
    )
  }
  structure(
    c(
      list(
        learner = "catboost_aft_survival",
        distribution = model$distribution,
        scale = model$scale,
        risk_reference = "catboost_aft_normal_native",
        risk_formula = "negative_catboost_aft_log_time_margin"
      ),
      extracted,
      list(survival_probability = survival_probability)
    ),
    class = c(
      "stormbreaker_catboost_aft_survival_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_catboost_aft_survival_learner_state <- function(state, x) {
  x <- sb_validate_learner_state_x(
    state,
    x,
    "catboost_aft_survival"
  )
  valid_state <- is.raw(state[["model_raw"]]) &&
    length(state[["model_raw"]]) > 0L &&
    is.numeric(state[["iterations"]]) &&
    length(state[["iterations"]]) == 1L &&
    !is.na(state[["iterations"]]) && state[["iterations"]] >= 1 &&
    identical(state[["distribution"]], "Normal") &&
    is.numeric(state[["scale"]]) && length(state[["scale"]]) == 1L &&
    !is.na(state[["scale"]]) && is.finite(state[["scale"]]) &&
    state[["scale"]] > 0
  if (!valid_state) {
    sb_bad_learner_state("The saved CatBoost AFT learner state is invalid.")
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for fixed CatBoost AFT prediction.",
      "stormbreaker_missing_catboost"
    )
  }
  booster <- tryCatch(
    unserialize(state[["model_raw"]]),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved CatBoost AFT state could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_catboost_state",
        parent = error
      )
    }
  )
  model <- structure(
    list(
      features = state[["features"]],
      params = state[["params"]],
      iterations = as.integer(state[["iterations"]]),
      distribution = state[["distribution"]],
      scale = state[["scale"]],
      booster = booster
    ),
    class = "stormbreaker_catboost_aft_survival"
  )
  result <- list(risk = sb_predict_catboost_aft_survival(model, x))
  probability <- state[["survival_probability"]]
  if (!is.null(probability)) {
    valid_probability <- is.list(probability) &&
      is.numeric(probability[["prediction_horizon"]]) &&
      length(probability[["prediction_horizon"]]) == 1L &&
      !is.na(probability[["prediction_horizon"]]) &&
      is.finite(probability[["prediction_horizon"]]) &&
      probability[["prediction_horizon"]] > 0 &&
      is.character(probability[["time_unit"]]) &&
      length(probability[["time_unit"]]) == 1L &&
      !is.na(probability[["time_unit"]]) &&
      nzchar(probability[["time_unit"]]) &&
      identical(probability[["source"]], "catboost_aft_normal_native") &&
      identical(probability[["ties"]], "not_applicable")
    if (!valid_probability) {
      sb_bad_learner_state(
        "The saved CatBoost AFT survival-probability state is invalid."
      )
    }
    result$prediction_horizon <- probability[["prediction_horizon"]]
    result$time_unit <- probability[["time_unit"]]
    result$survival_probability <- sb_catboost_aft_survival_probability(
      result$risk,
      probability[["prediction_horizon"]],
      model$scale
    )
  }
  result
}

sb_catboost_aft_survival_learner_adapter <- function(
    iterations = 100L,
    learning_rate = 0.1,
    depth = 6L,
    l2_leaf_reg = 3,
    seed = 1L) {
  params <- sb_catboost_aft_survival_params(
    iterations = iterations,
    learning_rate = learning_rate,
    depth = depth,
    l2_leaf_reg = l2_leaf_reg,
    seed = seed
  )
  list(
    id = "catboost_aft_survival",
    params = params,
    fit = function(x, time, event) {
      sb_fit_catboost_aft_survival(x, time, event, params)
    },
    predict = sb_predict_catboost_aft_survival,
    predict_survival = sb_predict_catboost_aft_survival_survival,
    extract_tree_model = sb_extract_catboost_tree_model,
    extract_state = function(model, context = list()) {
      sb_extract_catboost_aft_survival_learner_state(model, context)
    },
    predict_state = sb_predict_catboost_aft_survival_learner_state
  )
}
