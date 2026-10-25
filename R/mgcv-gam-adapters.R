
sb_mgcv_gam_binary_learner_ids <- function() character()

sb_mgcv_gam_survival_learner_ids <- function() character()

sb_mgcv_gam_learner_ids <- function() {
  c(
    sb_mgcv_gam_binary_learner_ids(),
    sb_mgcv_gam_survival_learner_ids()
  )
}

sb_mgcv_gam_params <- function(learner) {
  if (!learner %in% sb_mgcv_gam_learner_ids()) {
    sb_abort(
      paste0("Unknown mgcv GAM learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  common <- list(
    basis = "cs",
    basis_dimension = 5L,
    method = "REML",
    threads = 1L
  )
  if (learner %in% sb_mgcv_gam_binary_learner_ids()) {
    return(c(common, list(family = "binomial", link = "logit")))
  }
  c(common, list(family = "Cox PH", link = "identity", ties = "Peto"))
}

sb_mgcv_available <- function() {
  requireNamespace("mgcv", quietly = TRUE)
}

sb_require_mgcv <- function() {
  if (!sb_mgcv_available()) {
    sb_abort(
      "Package `mgcv` is required for GAM learners.",
      "stormbreaker_missing_mgcv"
    )
  }
  invisible(TRUE)
}

sb_mgcv_internal_feature_names <- function(n_features) {
  paste0(".stormbreaker_mgcv_x", seq_len(n_features))
}

sb_mgcv_formula_data <- function(x, response, params) {
  engine_features <- sb_mgcv_internal_feature_names(ncol(x))
  response_name <- ".stormbreaker_mgcv_response"
  analysis <- x
  names(analysis) <- engine_features
  analysis[[response_name]] <- response
  smooth_terms <- paste0(
    "s(", engine_features,
    ", bs = '", params$basis,
    "', k = ", params$basis_dimension, ")"
  )
  formula <- stats::as.formula(
    paste(response_name, "~", paste(smooth_terms, collapse = " + ")),
    env = asNamespace("mgcv")
  )
  list(
    data = analysis,
    formula = formula,
    engine_features = engine_features,
    response_name = response_name
  )
}

sb_mgcv_fit_call <- function(formula_data, endpoint, learner, params) {
  sb_reject_retired_learner(learner)
}


sb_validate_mgcv_engine <- function(engine, learner, engine_features,
                                    require_baseline = TRUE) {
  params <- sb_mgcv_gam_params(learner)
  coefficients <- tryCatch(
    stats::coef(engine),
    error = function(error) NULL
  )
  smooth_terms <- if (inherits(engine, "gam") && is.list(engine$smooth)) {
    vapply(engine$smooth, function(smooth) {
      as.character(smooth$term)[[1L]]
    }, "")
  } else {
    character()
  }
  smooth_dimensions <- if (inherits(engine, "gam") &&
      is.list(engine$smooth)) {
    vapply(engine$smooth, function(smooth) {
      as.integer(smooth$bs.dim)[[1L]]
    }, integer(1))
  } else {
    integer()
  }
  smooth_ok <- inherits(engine, "gam") &&
    identical(unname(smooth_terms), engine_features) &&
    identical(
      unname(smooth_dimensions),
      rep(params$basis_dimension, length(engine_features))
    )
  convergence_ok <- if (learner %in% sb_mgcv_gam_binary_learner_ids()) {
    isTRUE(engine$converged)
  } else {
    is.null(engine$converged) || isTRUE(engine$converged)
  }
  common_ok <- smooth_ok && convergence_ok &&
    identical(engine$method, params$method) &&
    is.list(engine$outer.info) &&
    identical(engine$outer.info$conv, "full convergence") &&
    is.numeric(coefficients) && length(coefficients) > 0L &&
    !anyNA(coefficients) && all(is.finite(coefficients)) &&
    is.numeric(engine$sp) && length(engine$sp) == length(engine_features) &&
    !anyNA(engine$sp) && all(is.finite(engine$sp)) && all(engine$sp >= 0)
  family_ok <- identical(engine$family$family, params$family) &&
    identical(engine$family$link, params$link)
  baseline_ok <- TRUE
  if (learner %in% sb_mgcv_gam_survival_learner_ids() && require_baseline) {
    baseline <- engine$family$data
    baseline_ok <- is.list(baseline) &&
      is.numeric(baseline$tr) && is.numeric(baseline$h) &&
      length(baseline$tr) > 0L &&
      length(baseline$tr) == length(baseline$h) &&
      !anyNA(baseline$tr) && all(is.finite(baseline$tr)) &&
      all(baseline$tr > 0) &&
      !anyNA(baseline$h) && all(is.finite(baseline$h)) &&
      all(baseline$h >= 0)
  }
  if (!isTRUE(common_ok && family_ok && baseline_ok)) {
    sb_abort(
      paste0("Learner `", learner, "` returned an incomplete mgcv GAM."),
      "stormbreaker_incomplete_mgcv_fit"
    )
  }
  engine
}

sb_mgcv_prediction_engine <- function(engine, learner,
                                      retain_cox_baseline = TRUE) {
  patient_fields <- c(
    "model", "y", "z", "weights", "prior.weights", "working.weights",
    "fitted.values", "residuals", "linear.predictors", "offset", "hat",
    "call"
  )
  engine[intersect(patient_fields, names(engine))] <- NULL
  if (learner %in% sb_mgcv_gam_survival_learner_ids()) {
    if (isTRUE(retain_cox_baseline)) {
      engine$family$data <- engine$family$data[c("tr", "h")]
    } else {
      engine$family$data <- NULL
    }
  }
  engine
}

sb_new_mgcv_gam_model <- function(learner, features, engine_features,
                                  response_name, params, engine) {
  structure(
    list(
      learner = learner,
      features = features,
      engine_features = engine_features,
      response_name = response_name,
      params = params,
      engine = engine
    ),
    class = "stormbreaker_mgcv_gam"
  )
}

sb_validate_mgcv_binary_inputs <- function(x, outcome, learner) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      paste0("A `", learner, "` analysis fold must contain both classes."),
      "stormbreaker_mgcv_binary_single_class"
    )
  }
  list(x = x, outcome = outcome)
}

sb_validate_mgcv_survival_inputs <- function(x, time, event, learner) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time <= 0)) {
    sb_abort(
      paste0("Learner `", learner, "` requires strictly positive times."),
      "stormbreaker_mgcv_nonpositive_time"
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
      paste0("Learner `", learner, "` requires at least one event."),
      "stormbreaker_mgcv_no_event"
    )
  }
  list(x = x, time = time, event = as.integer(event))
}

sb_fit_mgcv_gam <- function(x, endpoint, learner, params) {
  sb_reject_retired_learner(learner)
}


sb_validate_mgcv_model_x <- function(model, x) {
  valid_model <- inherits(model, "stormbreaker_mgcv_gam") &&
    is.character(model$learner) && length(model$learner) == 1L &&
    model$learner %in% sb_mgcv_gam_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.character(model$engine_features) &&
    length(model$engine_features) == length(model$features) &&
    is.character(model$response_name) && length(model$response_name) == 1L &&
    is.list(model$params)
  if (!valid_model) {
    sb_abort(
      "The fitted mgcv GAM is incomplete or invalid.",
      "stormbreaker_bad_mgcv_model"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted mgcv GAM.",
      "stormbreaker_mgcv_schema_mismatch"
    )
  }
  x
}

sb_mgcv_predict_call <- function(engine, newdata, type, learner) {
  tryCatch(
    stats::predict(engine, newdata = newdata, type = type),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) {
        stop(error)
      }
      sb_abort(
        paste0(
          "Learner `", learner, "` could not predict from its mgcv GAM: ",
          conditionMessage(error)
        ),
        "stormbreaker_mgcv_prediction_error",
        parent = error
      )
    }
  )
}

sb_predict_mgcv_gam <- function(model, x) {
  sb_require_mgcv()
  x <- sb_validate_mgcv_model_x(model, x)
  sb_validate_mgcv_engine(
    model$engine,
    model$learner,
    model$engine_features,
    require_baseline = FALSE
  )
  newdata <- x
  names(newdata) <- model$engine_features
  type <- if (model$learner %in% sb_mgcv_gam_binary_learner_ids()) {
    "response"
  } else {
    "link"
  }
  values <- unname(sb_mgcv_predict_call(
    model$engine,
    newdata,
    type,
    model$learner
  ))
  valid <- length(values) == nrow(x) && !anyNA(values) &&
    all(is.finite(values))
  if (model$learner %in% sb_mgcv_gam_binary_learner_ids()) {
    valid <- valid && all(values >= 0 & values <= 1)
  }
  if (!valid) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid predictions."),
      if (model$learner %in% sb_mgcv_gam_binary_learner_ids()) {
        "stormbreaker_invalid_binary_prediction"
      } else {
        "stormbreaker_invalid_survival_prediction"
      }
    )
  }
  values
}

sb_mgcv_cox_baseline_hazard <- function(model, prediction_horizon) {
  baseline <- model$engine$family$data
  valid <- is.list(baseline) && is.numeric(baseline$tr) &&
    is.numeric(baseline$h) && length(baseline$tr) > 0L &&
    length(baseline$tr) == length(baseline$h) &&
    !anyNA(baseline$tr) && all(is.finite(baseline$tr)) &&
    !anyNA(baseline$h) && all(is.finite(baseline$h)) &&
    all(baseline$h >= 0)
  if (!valid) {
    sb_abort(
      "The fitted mgcv Cox baseline hazard is incomplete or invalid.",
      "stormbreaker_invalid_baseline_hazard"
    )
  }
  ordering <- order(baseline$tr, seq_along(baseline$tr))
  times <- baseline$tr[ordering]
  hazards <- baseline$h[ordering]
  index <- findInterval(prediction_horizon, times)
  cumulative_hazard <- if (index == 0L) 0 else hazards[[index]]
  if (!is.finite(cumulative_hazard) || cumulative_hazard < 0) {
    sb_abort(
      "The fitted mgcv Cox baseline hazard is incomplete or invalid.",
      "stormbreaker_invalid_baseline_hazard"
    )
  }
  unname(cumulative_hazard)
}

sb_predict_mgcv_gam_survival <- function(model, x_analysis, time, event,
                                         x_assessment,
                                         prediction_horizon) {
  if (!inherits(model, "stormbreaker_mgcv_gam") ||
      !model$learner %in% sb_mgcv_gam_survival_learner_ids()) {
    sb_abort(
      "The fitted mgcv Cox GAM is incomplete or invalid.",
      "stormbreaker_bad_mgcv_model"
    )
  }
  inputs <- sb_validate_mgcv_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted mgcv Cox GAM.",
      "stormbreaker_mgcv_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  risk <- sb_predict_mgcv_gam(model, x_assessment)
  baseline_hazard <- sb_mgcv_cox_baseline_hazard(model, prediction_horizon)
  sb_survival_prediction_result(
    prediction_horizon,
    baseline_hazard,
    source = "mgcv_cox_ph_native_baseline",
    ties = "Peto",
    risk = risk
  )
}

sb_extract_mgcv_gam_learner_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_mgcv_gam")) {
    sb_bad_learner_state("The fitted mgcv GAM learner state is invalid.")
  }
  baseline_hazard <- NULL
  if (model$learner %in% sb_mgcv_gam_survival_learner_ids() &&
      "prediction_horizon" %in% names(context)) {
    prediction <- sb_predict_mgcv_gam_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    baseline_hazard <- list(
      prediction_horizon = prediction$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      cumulative_baseline_hazard = prediction$cumulative_baseline_hazard,
      baseline_survival = prediction$baseline_survival,
      source = prediction$source,
      ties = prediction$ties
    )
  }
  engine <- sb_mgcv_prediction_engine(
    model$engine,
    model$learner,
    retain_cox_baseline = FALSE
  )
  structure(
    list(
      learner = model$learner,
      features = model$features,
      engine_features = model$engine_features,
      response_name = model$response_name,
      params = model$params,
      engine = engine,
      baseline_hazard = baseline_hazard
    ),
    class = c(
      "stormbreaker_mgcv_gam_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_mgcv_gam_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
      !learner %in% sb_mgcv_gam_learner_ids()) {
    sb_bad_learner_state("The saved mgcv GAM learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_new_mgcv_gam_model(
    learner,
    state$features,
    state$engine_features,
    state$response_name,
    state$params,
    state$engine
  )
  values <- sb_predict_mgcv_gam(model, x)
  if (learner %in% sb_mgcv_gam_binary_learner_ids()) {
    return(list(probability = values))
  }
  result <- list(risk = values)
  baseline <- state$baseline_hazard
  if (!is.null(baseline)) {
    valid_baseline <- is.list(baseline) &&
      is.numeric(baseline$prediction_horizon) &&
      length(baseline$prediction_horizon) == 1L &&
      !is.na(baseline$prediction_horizon) &&
      is.finite(baseline$prediction_horizon) &&
      baseline$prediction_horizon > 0 &&
      is.character(baseline$time_unit) && length(baseline$time_unit) == 1L &&
      !is.na(baseline$time_unit) && nzchar(baseline$time_unit)
    if (!valid_baseline) {
      sb_bad_learner_state(
        "The saved mgcv Cox baseline-hazard state is invalid."
      )
    }
    result$prediction_horizon <- baseline$prediction_horizon
    result$time_unit <- baseline$time_unit
    result$survival_probability <- sb_survival_probability_from_hazard(
      values,
      baseline$cumulative_baseline_hazard
    )
  }
  result
}

sb_mgcv_gam_learner_adapter <- function(learner,
                                        params = sb_mgcv_gam_params(learner)) {
  if (!learner %in% sb_mgcv_gam_learner_ids()) {
    sb_abort(
      paste0("Unknown mgcv GAM learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  adapter <- list(
    id = learner,
    params = params,
    fit = if (learner %in% sb_mgcv_gam_binary_learner_ids()) {
      function(x, outcome) {
        sb_fit_mgcv_gam(x, list(outcome = outcome), learner, params)
      }
    } else {
      function(x, time, event) {
        sb_fit_mgcv_gam(
          x,
          list(time = time, event = event),
          learner,
          params
        )
      }
    },
    predict = sb_predict_mgcv_gam,
    extract_state = sb_extract_mgcv_gam_learner_state,
    predict_state = sb_predict_mgcv_gam_learner_state
  )
  if (learner %in% sb_mgcv_gam_survival_learner_ids()) {
    adapter$predict_survival <- sb_predict_mgcv_gam_survival
  }
  adapter
}
