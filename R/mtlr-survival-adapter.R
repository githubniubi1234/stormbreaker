sb_mtlr_survival_learner_ids <- function() {
  "mtlr_survival"
}

sb_mtlr_survival_params <- function(learner = sb_mtlr_survival_learner_ids()) {
  if (!identical(learner, sb_mtlr_survival_learner_ids())) {
    sb_abort(
      paste0("Unknown MTLR survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    interval_rule = "ceiling_sqrt_analysis_n",
    normalize = TRUE,
    C1 = 1,
    train_biases = TRUE,
    train_uncensored = TRUE,
    threshold = 1e-5,
    maximum_iterations = 5000L,
    lower_weight_bound = -15,
    upper_weight_bound = 15
  )
}

sb_mtlr_survival_parameter_semantics <- function() {
  paste0(
    "multi-task logistic regression with ceiling-sqrt analysis-fold time ",
    "intervals, normalization, C1 1, and bounded optimization fixed before ",
    "resampling"
  )
}

sb_mtlr_survival_available <- function() {
  requireNamespace("MTLR", quietly = TRUE)
}

sb_require_mtlr_survival_dependency <- function() {
  if (!sb_mtlr_survival_available()) {
    sb_abort(
      "Package `MTLR` is required for learner `mtlr_survival`.",
      "stormbreaker_missing_mtlr_survival_dependency"
    )
  }
  invisible(TRUE)
}

sb_mtlr_internal_features <- function(p) {
  paste0(".stormbreaker_mtlr_x", seq_len(p))
}

sb_mtlr_formula_data <- function(x, time, event) {
  internal_features <- sb_mtlr_internal_features(ncol(x))
  analysis <- x
  names(analysis) <- internal_features
  time_name <- ".stormbreaker_mtlr_time"
  event_name <- ".stormbreaker_mtlr_event"
  analysis[[time_name]] <- time
  analysis[[event_name]] <- event
  response <- paste0(
    "survival::Surv(", time_name, ", ", event_name, ")"
  )
  list(
    data = analysis,
    formula = stats::reformulate(
      internal_features,
      response = response,
      env = baseenv()
    ),
    internal_features = internal_features
  )
}

sb_mtlr_prediction_engine <- function(engine) {
  engine[c("x", "y", "response", "Call")] <- NULL
  attr(engine$Terms, ".Environment") <- baseenv()
  engine
}

sb_new_mtlr_survival_model <- function(
  features, internal_features, params, requested_intervals,
  effective_time_points, engine
) {
  structure(
    list(
      learner = sb_mtlr_survival_learner_ids(),
      features = features,
      internal_features = internal_features,
      params = params,
      requested_intervals = requested_intervals,
      effective_time_points = effective_time_points,
      engine = engine
    ),
    class = "stormbreaker_mtlr_survival_model"
  )
}

sb_fit_mtlr_survival <- function(x, time, event, learner, params) {
  inputs <- sb_validate_classical_survival_inputs(x, time, event, learner)
  sb_require_mtlr_survival_dependency()
  formula_data <- sb_mtlr_formula_data(
    inputs$x,
    inputs$time,
    inputs$event
  )
  intervals <- as.integer(ceiling(sqrt(nrow(inputs$x))))
  engine <- MTLR::mtlr(
    formula = formula_data$formula,
    data = formula_data$data,
    nintervals = intervals,
    normalize = params$normalize,
    C1 = params$C1,
    train_biases = params$train_biases,
    train_uncensored = params$train_uncensored,
    threshold = params$threshold,
    maxit = params$maximum_iterations,
    lower = params$lower_weight_bound,
    upper = params$upper_weight_bound
  )
  engine <- sb_mtlr_prediction_engine(engine)
  valid <- inherits(engine, "mtlr") &&
    is.matrix(engine$weight_matrix) &&
    nrow(engine$weight_matrix) == length(engine$time_points) &&
    ncol(engine$weight_matrix) == length(formula_data$internal_features) + 1L &&
    !anyNA(engine$weight_matrix) && all(is.finite(engine$weight_matrix)) &&
    is.numeric(engine$time_points) && length(engine$time_points) >= 1L &&
    !anyNA(engine$time_points) && all(is.finite(engine$time_points))
  if (!valid) {
    sb_abort(
      "MTLR did not return finite weights and a survival time grid.",
      "stormbreaker_incomplete_mtlr_survival_fit"
    )
  }
  sb_new_mtlr_survival_model(
    names(inputs$x),
    formula_data$internal_features,
    params,
    intervals,
    length(engine$time_points),
    engine
  )
}

sb_mtlr_prediction_data <- function(model, x) {
  if (!inherits(model, "stormbreaker_mtlr_survival_model")) {
    sb_abort(
      "The supplied model is not a fitted MTLR survival learner.",
      "stormbreaker_bad_mtlr_survival_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted MTLR learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  names(x) <- model$internal_features
  x
}

sb_predict_mtlr_survival_risk <- function(model, x) {
  sb_require_mtlr_survival_dependency()
  prediction_data <- sb_mtlr_prediction_data(model, x)
  mean_time <- as.numeric(stats::predict(
    model$engine,
    newdata = prediction_data,
    type = "mean_time"
  ))
  risk <- -mean_time
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "MTLR produced invalid native mean-survival risk scores.",
      "stormbreaker_invalid_mtlr_survival_prediction"
    )
  }
  unname(risk)
}

sb_mtlr_survival_probability <- function(model, x, prediction_horizon) {
  sb_require_mtlr_survival_dependency()
  prediction_data <- sb_mtlr_prediction_data(model, x)
  curves <- stats::predict(
    model$engine,
    newdata = prediction_data,
    type = "prob_times",
    times = prediction_horizon
  )
  probability <- if (is.data.frame(curves) && nrow(curves) >= 1L) {
    unname(as.numeric(curves[nrow(curves), -1L, drop = TRUE]))
  } else {
    numeric()
  }
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      "MTLR produced invalid native survival probabilities.",
      "stormbreaker_invalid_mtlr_survival_probability"
    )
  }
  probability
}

sb_mtlr_direct_survival_source <- function() {
  "MTLR_native_multitask_survival_curve"
}

sb_predict_mtlr_survival_probability <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  inputs <- sb_validate_classical_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted MTLR learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  list(
    prediction_horizon = prediction_horizon,
    source = sb_mtlr_direct_survival_source(),
    ties = "not_applicable",
    probability = sb_mtlr_survival_probability(
      model,
      x_assessment,
      prediction_horizon
    )
  )
}

sb_extract_mtlr_survival_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_mtlr_survival_model")) {
    sb_bad_learner_state("The fitted MTLR survival learner state is invalid.")
  }
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    horizon <- sb_validate_prediction_horizon(
      context$prediction_horizon,
      context$time
    )
    survival_probability <- list(
      prediction_horizon = horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      source = sb_mtlr_direct_survival_source(),
      ties = "not_applicable"
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      internal_features = model$internal_features,
      params = model$params,
      requested_intervals = model$requested_intervals,
      effective_time_points = model$effective_time_points,
      engine = model$engine,
      training_data_retention = "none",
      risk_reference = "MTLR_native_time_grid",
      risk_formula = "negative_native_mean_predicted_survival_time",
      baseline_hazard = NULL,
      survival_probability = survival_probability
    ),
    class = c(
      "stormbreaker_mtlr_survival_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_mtlr_survival_model <- function(state) {
  valid <- inherits(state, "stormbreaker_mtlr_survival_state") &&
    identical(state$learner, sb_mtlr_survival_learner_ids()) &&
    is.character(state$internal_features) &&
    length(state$internal_features) == length(state$features) &&
    is.list(state$params) && inherits(state$engine, "mtlr") &&
    identical(state$training_data_retention, "none")
  if (!isTRUE(valid)) {
    sb_bad_learner_state("The saved MTLR survival learner state is invalid.")
  }
  sb_new_mtlr_survival_model(
    state$features,
    state$internal_features,
    state$params,
    state$requested_intervals,
    state$effective_time_points,
    state$engine
  )
}

sb_predict_mtlr_survival_state <- function(state, x) {
  x <- sb_validate_learner_state_x(
    state,
    x,
    sb_mtlr_survival_learner_ids()
  )
  model <- sb_restore_mtlr_survival_model(state)
  result <- list(risk = sb_predict_mtlr_survival_risk(model, x))
  if (!is.null(state$survival_probability)) {
    probability <- state$survival_probability
    valid <- is.numeric(probability$prediction_horizon) &&
      length(probability$prediction_horizon) == 1L &&
      !is.na(probability$prediction_horizon) &&
      is.finite(probability$prediction_horizon) &&
      probability$prediction_horizon > 0 &&
      is.character(probability$time_unit) &&
      length(probability$time_unit) == 1L &&
      !is.na(probability$time_unit) && nzchar(probability$time_unit) &&
      identical(probability$source, sb_mtlr_direct_survival_source()) &&
      identical(probability$ties, "not_applicable")
    if (!isTRUE(valid)) {
      sb_bad_learner_state(
        "The saved MTLR direct survival-probability state is invalid."
      )
    }
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- sb_mtlr_survival_probability(
      model,
      x,
      probability$prediction_horizon
    )
  }
  result
}

sb_mtlr_survival_learner_adapter <- function(
  learner = sb_mtlr_survival_learner_ids(),
  params = sb_mtlr_survival_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_mtlr_survival(x, time, event, learner, params)
    },
    predict = sb_predict_mtlr_survival_risk,
    predict_survival = sb_predict_mtlr_survival_probability,
    extract_state = sb_extract_mtlr_survival_state,
    predict_state = sb_predict_mtlr_survival_state
  )
}
