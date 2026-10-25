sb_pamr_survival_learner_ids <- function() {
  "pamr_survival"
}

sb_pamr_survival_params <- function(
  learner = sb_pamr_survival_learner_ids()
) {
  if (!identical(learner, sb_pamr_survival_learner_ids())) {
    sb_abort(
      paste0("Unknown PAMR survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    threshold = 1,
    scale_sd = TRUE,
    offset_percent = 50,
    remove_zero_thresholds = FALSE,
    prior = "native_equal_survival_groups",
    survival_groups = 2L,
    ties = "breslow"
  )
}

sb_pamr_survival_parameter_semantics <- function() {
  paste0(
    "two native censoring-aware survival groups, nearest shrunken-centroid ",
    "threshold 1, variance scaling, equal group priors, and offset 50 fixed ",
    "before resampling"
  )
}

sb_pamr_survival_available <- function() {
  requireNamespace("pamr", quietly = TRUE)
}

sb_require_pamr_survival_dependency <- function() {
  if (!sb_pamr_survival_available()) {
    sb_abort(
      "Package `pamr` is required for learner `pamr_survival`.",
      "stormbreaker_missing_pamr_survival_dependency"
    )
  }
  invisible(TRUE)
}

sb_new_pamr_survival_model <- function(features, params, audit, engine) {
  structure(
    list(
      learner = sb_pamr_survival_learner_ids(),
      features = features,
      params = params,
      audit = audit,
      engine = engine
    ),
    class = "stormbreaker_pamr_survival_model"
  )
}

sb_fit_pamr_survival <- function(x, time, event, learner, params) {
  inputs <- sb_validate_classical_survival_inputs(x, time, event, learner)
  sb_require_pamr_survival_dependency()
  if (ncol(inputs$x) < 2L) {
    sb_abort(
      paste0(
        "PAMR survival training requires at least 2 predictors because ",
        "the native subset path does not preserve a one-row feature matrix."
      ),
      "stormbreaker_pamr_survival_too_few_features"
    )
  }
  data <- list(
    x = t(data.matrix(inputs$x)),
    survival.time = inputs$time,
    censoring.status = inputs$event,
    geneid = names(inputs$x),
    genenames = names(inputs$x)
  )
  engine <- NULL
  invisible(utils::capture.output(
    engine <- pamr::pamr.train(
      data,
      threshold = params$threshold,
      scale.sd = params$scale_sd,
      offset.percent = params$offset_percent,
      prior = NULL,
      remove.zeros = params$remove_zero_thresholds,
      sign.contrast = "both",
      ngroup.survival = params$survival_groups
    ),
    type = "output"
  ))
  cutoff <- unname(as.numeric(engine$cutoffs.survival))
  valid <- inherits(engine, "pamrtrained") &&
    identical(engine$problem.type, "surv.km") &&
    length(cutoff) == 1L && is.finite(cutoff) && cutoff > 0 &&
    identical(as.integer(engine$ngroup.survival), params$survival_groups)
  if (!valid) {
    sb_abort(
      "PAMR did not return a two-group survival centroid fit.",
      "stormbreaker_incomplete_pamr_survival_fit"
    )
  }
  model <- sb_new_pamr_survival_model(
    names(inputs$x),
    params,
    list(
      survival_cutoff = cutoff,
      early_survival_class = "1",
      late_survival_class = "2"
    ),
    sb_pamr_binary_prediction_engine(engine)
  )
  sb_predict_pamr_survival_risk(model, inputs$x)
  model
}

sb_validate_pamr_survival_model_x <- function(model, x) {
  if (!inherits(model, "stormbreaker_pamr_survival_model")) {
    sb_abort(
      "The supplied model is not a fitted PAMR survival learner.",
      "stormbreaker_bad_pamr_survival_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted PAMR survival schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  list(model = model, x = x)
}

sb_predict_pamr_survival_risk <- function(model, x) {
  sb_require_pamr_survival_dependency()
  validated <- sb_validate_pamr_survival_model_x(model, x)
  posterior <- pamr::pamr.predict(
    validated$model$engine,
    newx = t(data.matrix(validated$x)),
    threshold = validated$model$params$threshold,
    type = "posterior"
  )
  risk <- if (is.matrix(posterior) && "1" %in% colnames(posterior)) {
    unname(posterior[, "1"])
  } else {
    numeric()
  }
  if (length(risk) != nrow(validated$x) || anyNA(risk) ||
    any(!is.finite(risk)) || any(risk < 0 | risk > 1)) {
    sb_abort(
      "PAMR produced invalid early-survival posterior risk scores.",
      "stormbreaker_invalid_pamr_survival_prediction"
    )
  }
  risk
}

sb_pamr_survival_source <- function() {
  "pamr_survival_analysis_cox_breslow_calibration"
}

sb_fit_pamr_survival_calibration <- function(
  time, event, risk, prediction_horizon
) {
  calibration <- sb_classical_survival_stop_on_warning(
    survival::coxph(
      survival::Surv(time, event) ~ risk,
      data = data.frame(time = time, event = event, risk = risk),
      ties = "breslow",
      singular.ok = FALSE,
      model = FALSE,
      x = FALSE,
      y = FALSE
    ),
    "PAMR survival Cox calibration",
    "stormbreaker_pamr_survival_calibration_warning"
  )
  coefficient <- unname(as.numeric(stats::coef(calibration)))
  if (length(coefficient) != 1L || is.na(coefficient) ||
    !is.finite(coefficient) || coefficient < 0) {
    sb_abort(
      "PAMR survival risk could not be calibrated in its native direction.",
      "stormbreaker_invalid_pamr_survival_calibration"
    )
  }
  list(
    coefficient = coefficient,
    cumulative_baseline_hazard = sb_breslow_baseline_hazard(
      time,
      event,
      coefficient * risk,
      prediction_horizon
    )
  )
}

sb_predict_pamr_survival_probability <- function(
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
      "Analysis columns do not match the fitted PAMR survival schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  analysis_risk <- sb_predict_pamr_survival_risk(model, inputs$x)
  assessment_risk <- sb_predict_pamr_survival_risk(model, x_assessment)
  calibration <- sb_fit_pamr_survival_calibration(
    inputs$time,
    inputs$event,
    analysis_risk,
    prediction_horizon
  )
  result <- sb_survival_prediction_result(
    prediction_horizon,
    calibration$cumulative_baseline_hazard,
    source = sb_pamr_survival_source(),
    ties = model$params$ties,
    risk = calibration$coefficient * assessment_risk
  )
  result$calibration_coefficient <- calibration$coefficient
  result
}

sb_extract_pamr_survival_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_pamr_survival_model")) {
    sb_bad_learner_state("The fitted PAMR survival learner state is invalid.")
  }
  baseline_hazard <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_pamr_survival_probability(
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
      ties = probability$ties,
      calibration_coefficient = probability$calibration_coefficient
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      audit = model$audit,
      engine = model$engine,
      training_data_retention = "none",
      risk_reference = "pamr_native_survival_groups",
      risk_formula = "native_posterior_for_early_survival_group",
      baseline_hazard = baseline_hazard
    ),
    class = c(
      "stormbreaker_pamr_survival_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_pamr_survival_state <- function(state, x) {
  if (!inherits(state, "stormbreaker_pamr_survival_state") ||
    !identical(state$learner, sb_pamr_survival_learner_ids()) ||
    !inherits(state$engine, "pamrtrained") ||
    !identical(state$training_data_retention, "none")) {
    sb_bad_learner_state("The saved PAMR survival learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, state$learner)
  model <- sb_new_pamr_survival_model(
    state$features,
    state$params,
    state$audit,
    state$engine
  )
  risk <- sb_predict_pamr_survival_risk(model, x)
  result <- list(risk = risk)
  baseline <- state$baseline_hazard
  if (!is.null(baseline)) {
    coefficient <- baseline$calibration_coefficient
    valid <- is.numeric(coefficient) && length(coefficient) == 1L &&
      is.finite(coefficient) && coefficient >= 0 &&
      identical(baseline$source, sb_pamr_survival_source()) &&
      identical(baseline$ties, state$params$ties)
    if (!valid) {
      sb_bad_learner_state("The saved PAMR survival horizon state is invalid.")
    }
    result$prediction_horizon <- baseline$prediction_horizon
    result$time_unit <- baseline$time_unit
    result$survival_probability <- sb_survival_probability_from_hazard(
      coefficient * risk,
      baseline$cumulative_baseline_hazard
    )
  }
  result
}

sb_pamr_survival_learner_adapter <- function(
  learner = sb_pamr_survival_learner_ids(),
  params = sb_pamr_survival_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_pamr_survival(x, time, event, learner, params)
    },
    predict = sb_predict_pamr_survival_risk,
    predict_survival = sb_predict_pamr_survival_probability,
    extract_state = sb_extract_pamr_survival_state,
    predict_state = sb_predict_pamr_survival_state
  )
}
