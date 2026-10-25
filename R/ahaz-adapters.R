sb_ahaz_survival_learner_ids <- function() {
  c(
    "ahaz_lin_ying",
    "ahaz_lasso_bic",
    "ahaz_sscad_bic"
  )
}

sb_ahaz_available <- function() {
  requireNamespace("ahaz", quietly = TRUE)
}

sb_require_ahaz_dependency <- function() {
  if (!sb_ahaz_available()) {
    sb_abort(
      "Package `ahaz` is required for additive-hazards learners.",
      "stormbreaker_missing_ahaz"
    )
  }
  invisible(TRUE)
}

sb_ahaz_learner_params <- function(learner) {
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_ahaz_survival_learner_ids()) {
    sb_abort(
      "The requested additive-hazards learner is not supported.",
      "stormbreaker_unsupported_learner"
    )
  }
  if (identical(learner, "ahaz_lin_ying")) {
    return(list(
      estimator = "lin_ying",
      robust = FALSE
    ))
  }
  if (identical(learner, "ahaz_lasso_bic")) {
    return(list(
      estimator = "penalized",
      penalty = "lasso",
      alpha = 1,
      standardize = TRUE,
      tuning = "bic",
      nlambda = 100L,
      dfmax_rule = "min_30_p"
    ))
  }
  list(
    estimator = "penalized",
    penalty = "sscad",
    a = 3.7,
    nsteps = 1L,
    standardize = TRUE,
    tuning = "bic",
    nlambda = 100L,
    dfmax_rule = "min_30_p"
  )
}

sb_ahaz_parameter_semantics <- function(learner) {
  params <- sb_ahaz_learner_params(learner)
  if (identical(params$estimator, "lin_ying")) {
    return(paste0(
      "unpenalized Lin-Ying constant-coefficient additive-hazards fit ",
      "with model-based variance"
    ))
  }
  penalty <- if (identical(params$penalty, "lasso")) {
    "Lasso alpha 1"
  } else {
    "one-step SCAD with a=3.7"
  }
  paste0(
    penalty, ", standardization, 100-lambda path, and support cap ",
    "min(30,p) fixed; lambda selected by BIC inside each analysis fold"
  )
}

sb_ahaz_stop_on_warning <- function(code, stage) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      sb_abort(
        paste0(
          "The ahaz ", stage, " raised a warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_ahaz_warning",
        parent = warning
      )
    }
  )
}

sb_ahaz_wrap_error <- function(code, stage, subclass) {
  tryCatch(
    force(code),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) {
        stop(error)
      }
      sb_abort(
        paste0(
          "The ahaz ", stage, " failed: ", conditionMessage(error)
        ),
        subclass,
        parent = error
      )
    }
  )
}

sb_ahaz_engine_x <- function(x) {
  engine_x <- as.matrix(x)
  rownames(engine_x) <- NULL
  colnames(engine_x) <- paste0("x", seq_len(ncol(engine_x)))
  engine_x
}

sb_fit_ahaz_engine <- function(surv, engine_x, learner, params) {
  if (identical(params$estimator, "lin_ying")) {
    return(ahaz::ahaz(
      surv,
      engine_x,
      robust = params$robust
    ))
  }
  penalty <- if (identical(params$penalty, "lasso")) {
    ahaz::lasso.control(alpha = params$alpha)
  } else {
    ahaz::sscad.control(a = params$a, nsteps = params$nsteps)
  }
  ahaz::tune.ahazpen(
    surv,
    engine_x,
    standardize = params$standardize,
    penalty = penalty,
    tune = ahaz::bic.control(),
    dfmax = min(30L, ncol(engine_x)),
    nlambda = params$nlambda
  )
}

sb_ahaz_fit_audit <- function(engine, learner, coefficient) {
  if (identical(learner, "ahaz_lin_ying")) {
    return(list(
      selected_lambda = NA_real_,
      path_length = NA_integer_,
      selected_df = as.integer(sum(coefficient != 0))
    ))
  }
  selected_lambda <- as.numeric(engine$lambda.min)
  path <- as.numeric(engine$lambda)
  valid <- length(selected_lambda) == 1L && !is.na(selected_lambda) &&
    is.finite(selected_lambda) && selected_lambda > 0 &&
    length(path) > 0L && !anyNA(path) && all(is.finite(path)) &&
    all(path > 0) && is.numeric(engine$df) &&
    length(engine$df) == length(path)
  if (!valid) {
    sb_abort(
      paste0("Learner `", learner, "` returned an invalid BIC path."),
      "stormbreaker_invalid_ahaz_fit"
    )
  }
  selected_index <- which.min(abs(path - selected_lambda))
  list(
    selected_lambda = selected_lambda,
    path_length = as.integer(length(path)),
    selected_df = as.integer(engine$df[[selected_index]])
  )
}

sb_validate_ahaz_model <- function(model) {
  valid <- inherits(model, "stormbreaker_ahaz_survival_model") &&
    is.list(model) &&
    is.character(model$learner) && length(model$learner) == 1L &&
    !is.na(model$learner) &&
    model$learner %in% sb_ahaz_survival_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.numeric(model$coefficient) &&
    length(model$coefficient) == length(model$features) &&
    !anyNA(model$coefficient) && all(is.finite(model$coefficient)) &&
    identical(model$params, sb_ahaz_learner_params(model$learner)) &&
    is.list(model$audit)
  if (!valid) {
    sb_abort(
      "The supplied model is not a fitted ahaz survival learner.",
      "stormbreaker_bad_ahaz_fit"
    )
  }
  model
}

sb_ahaz_fit_times <- function(time, event, seed = 10101L) {
  event_times <- time[event == 1L]
  tied_times <- sort(unique(event_times[duplicated(event_times)]))
  fitted_time <- time
  audit <- list(
    method = "seeded_tied_event_offsets_fit_only",
    seed = seed,
    tied_event_groups = length(tied_times),
    tied_events = 0L,
    changed_events = 0L,
    max_abs_shift = 0
  )
  if (!length(tied_times)) return(list(time = time, audit = audit))
  distinct <- sort(unique(c(0, time)))
  withr::with_seed(seed, {
    for (value in tied_times) {
      ii <- which(time == value & event == 1L)
      previous <- distinct[match(value, distinct) - 1L]
      width <- min(value * 1e-8, (value - previous) / 4)
      # Keep deaths at/before their recorded time, so censoring at that time
      # remains after them. Random ranks avoid collisions from random draws.
      offsets <- sample.int(length(ii)) - 1L
      fitted_time[ii] <- value - width * offsets / length(ii)
      if (anyDuplicated(fitted_time[ii]) || any(fitted_time[ii] <= previous)) {
        sb_abort(
          "Tied event times cannot be separated at this numeric precision.",
          "stormbreaker_ahaz_time_precision"
        )
      }
      audit$tied_events <- audit$tied_events + length(ii)
    }
  })
  audit$changed_events <- sum(fitted_time != time)
  audit$max_abs_shift <- max(abs(fitted_time - time))
  list(time = fitted_time, audit = audit)
}

sb_fit_ahaz_survival <- function(x, time, event, learner, params) {
  inputs <- sb_validate_classical_survival_inputs(x, time, event, learner)
  sb_require_ahaz_dependency()
  engine_x <- sb_ahaz_engine_x(inputs$x)
  fit_times <- sb_ahaz_fit_times(inputs$time, inputs$event)
  surv <- survival::Surv(fit_times$time, inputs$event)
  engine <- sb_ahaz_wrap_error(
    sb_ahaz_stop_on_warning(
      sb_fit_ahaz_engine(surv, engine_x, learner, params),
      paste0("fit for `", learner, "`")
    ),
    paste0("fit for `", learner, "`"),
    "stormbreaker_ahaz_fit_failed"
  )
  coefficient <- unname(as.numeric(stats::coef(engine)))
  if (length(coefficient) != ncol(inputs$x) || anyNA(coefficient) ||
    any(!is.finite(coefficient))) {
    sb_abort(
      paste0(
        "Learner `", learner,
        "` did not return one finite coefficient per predictor."
      ),
      "stormbreaker_invalid_ahaz_fit"
    )
  }
  model <- structure(
    list(
      learner = learner,
      features = names(inputs$x),
      params = params,
      coefficient = coefficient,
      audit = c(
        sb_ahaz_fit_audit(engine, learner, coefficient),
        list(time_ties = fit_times$audit)
      )
    ),
    class = "stormbreaker_ahaz_survival_model"
  )
  sb_predict_ahaz_survival(model, inputs$x)
  if (!identical(learner, "ahaz_lin_ying")) model$diagnostic_path <- list(
    beta = engine$fit$beta, lambda = engine$lambda, df = engine$df,
    criterion = engine$tunem, selected_lambda = engine$lambda.min,
    criterion_type = engine$tune$type)
  model
}

sb_validate_ahaz_model_x <- function(model, x) {
  model <- sb_validate_ahaz_model(model)
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted ahaz learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  list(model = model, x = x)
}

sb_predict_ahaz_survival <- function(model, x) {
  validated <- sb_validate_ahaz_model_x(model, x)
  risk <- unname(as.vector(
    as.matrix(validated$x) %*% validated$model$coefficient
  ))
  if (length(risk) != nrow(validated$x) || anyNA(risk) ||
    any(!is.finite(risk))) {
    sb_abort(
      "The ahaz learner produced incomplete or non-finite risk scores.",
      "stormbreaker_invalid_ahaz_prediction"
    )
  }
  risk
}

sb_extract_ahaz_coefficients <- function(model) {
  model <- sb_validate_ahaz_model(model)
  data.frame(
    feature = model$features,
    coefficient = model$coefficient,
    nonzero = model$coefficient != 0,
    direction = ifelse(
      model$coefficient > 0,
      "higher_additive_hazard",
      ifelse(
        model$coefficient < 0,
        "lower_additive_hazard",
        "zero"
      )
    ),
    stringsAsFactors = FALSE
  )
}

sb_ahaz_survival_source <- function(learner) {
  paste0(learner, "_analysis_cox_breslow_calibration")
}

sb_fit_ahaz_calibration <- function(time, event, risk, prediction_horizon) {
  calibration <- sb_ahaz_wrap_error(
    sb_ahaz_stop_on_warning(
      survival::coxph(
        survival::Surv(time, event) ~ risk,
        data = data.frame(time = time, event = event, risk = risk),
        ties = "breslow",
        singular.ok = FALSE,
        model = FALSE,
        x = FALSE,
        y = FALSE
      ),
      "analysis-fold Cox/Breslow calibration"
    ),
    "analysis-fold Cox/Breslow calibration",
    "stormbreaker_ahaz_calibration_failed"
  )
  coefficient <- unname(as.numeric(stats::coef(calibration)))
  if (length(coefficient) != 1L || is.na(coefficient) ||
    !is.finite(coefficient)) {
    sb_abort(
      "The ahaz analysis-fold calibration coefficient is invalid.",
      "stormbreaker_ahaz_calibration_failed"
    )
  }
  if (coefficient < 0) {
    sb_abort(
      paste0(
        "The ahaz analysis-fold calibration reversed the native ",
        "higher-risk direction."
      ),
      "stormbreaker_ahaz_calibration_direction"
    )
  }
  calibrated_risk <- coefficient * risk
  if (anyNA(calibrated_risk) || any(!is.finite(calibrated_risk))) {
    sb_abort(
      "The ahaz calibrated analysis-fold risk is invalid.",
      "stormbreaker_ahaz_calibration_failed"
    )
  }
  list(
    coefficient = coefficient,
    cumulative_baseline_hazard = sb_breslow_baseline_hazard(
      time,
      event,
      calibrated_risk,
      prediction_horizon
    )
  )
}

sb_predict_ahaz_survival_probability <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  model <- sb_validate_ahaz_model(model)
  inputs <- sb_validate_classical_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted ahaz learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  analysis_risk <- sb_predict_ahaz_survival(model, inputs$x)
  assessment_risk <- sb_predict_ahaz_survival(model, x_assessment)
  calibration <- sb_fit_ahaz_calibration(
    inputs$time,
    inputs$event,
    analysis_risk,
    prediction_horizon
  )
  result <- sb_survival_prediction_result(
    prediction_horizon,
    calibration$cumulative_baseline_hazard,
    source = sb_ahaz_survival_source(model$learner),
    ties = "breslow",
    risk = calibration$coefficient * assessment_risk
  )
  result$calibration_coefficient <- calibration$coefficient
  result
}

sb_extract_ahaz_survival_state <- function(model, context = list()) {
  model <- sb_validate_ahaz_model(model)
  baseline_hazard <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_ahaz_survival_probability(
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
      cumulative_baseline_hazard =
        probability$cumulative_baseline_hazard,
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
      coefficient = model$coefficient,
      params = model$params,
      audit = model$audit,
      risk_reference = "ahaz_native",
      risk_formula = "x_beta",
      baseline_hazard = baseline_hazard
    ),
    class = c(
      "stormbreaker_ahaz_survival_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_ahaz_survival_state <- function(state, x) {
  learner <- if (is.list(state)) state$learner else NULL
  if (!inherits(state, "stormbreaker_ahaz_survival_state") ||
    !is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_ahaz_survival_learner_ids()) {
    sb_bad_learner_state("The saved ahaz learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  coefficient <- state$coefficient
  valid_coefficient <- is.numeric(coefficient) &&
    length(coefficient) == ncol(x) && !anyNA(coefficient) &&
    all(is.finite(coefficient))
  if (!valid_coefficient) {
    sb_bad_learner_state("The saved ahaz coefficients are invalid.")
  }
  risk <- unname(as.vector(as.matrix(x) %*% coefficient))
  if (anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "The saved ahaz state produced invalid risk predictions.",
      "stormbreaker_invalid_ahaz_prediction"
    )
  }
  result <- list(risk = risk)
  baseline <- state$baseline_hazard
  if (!is.null(baseline)) {
    calibration <- baseline$calibration_coefficient
    valid <- is.list(baseline) &&
      is.numeric(baseline$prediction_horizon) &&
      length(baseline$prediction_horizon) == 1L &&
      !is.na(baseline$prediction_horizon) &&
      is.finite(baseline$prediction_horizon) &&
      baseline$prediction_horizon > 0 &&
      is.numeric(calibration) && length(calibration) == 1L &&
      !is.na(calibration) && is.finite(calibration) && calibration >= 0 &&
      identical(baseline$source, sb_ahaz_survival_source(learner)) &&
      identical(baseline$ties, "breslow")
    if (!valid) {
      sb_bad_learner_state("The saved ahaz horizon state is invalid.")
    }
    result$prediction_horizon <- baseline$prediction_horizon
    result$time_unit <- baseline$time_unit
    result$survival_probability <- sb_survival_probability_from_hazard(
      calibration * risk,
      baseline$cumulative_baseline_hazard
    )
  }
  result
}

sb_ahaz_survival_learner_adapter <- function(
  learner, params = sb_ahaz_learner_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_ahaz_survival(x, time, event, learner, params)
    },
    predict = sb_predict_ahaz_survival,
    predict_survival = sb_predict_ahaz_survival_probability,
    extract_coefficients = sb_extract_ahaz_coefficients,
    extract_state = sb_extract_ahaz_survival_state,
    predict_state = sb_predict_ahaz_survival_state
  )
}
