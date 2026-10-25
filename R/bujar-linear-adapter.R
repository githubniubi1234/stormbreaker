sb_bujar_linear_learner_ids <- function() {
  "bujar_linear_boost"
}

sb_bujar_available <- function() {
  requireNamespace("bujar", quietly = TRUE)
}

sb_require_bujar_dependency <- function() {
  if (!sb_bujar_available()) {
    sb_abort(
      "Package `bujar` is required for Buckley-James linear boosting.",
      "stormbreaker_missing_bujar"
    )
  }
  invisible(TRUE)
}

sb_bujar_linear_params <- function(
  learner = sb_bujar_linear_learner_ids()
) {
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !identical(learner, sb_bujar_linear_learner_ids())) {
    sb_abort(
      "The requested Buckley-James learner is not supported.",
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    response_scale = "log_time",
    learner = "linear.regression",
    degree = 1,
    center = TRUE,
    mimpu = NULL,
    iter_bj = 20L,
    max_cycle = 10L,
    nu = 0.1,
    mstop = 50L,
    twin = FALSE,
    mstop2 = 100L,
    tuning = FALSE,
    cv = FALSE,
    vimpint = FALSE,
    rel_inf = FALSE,
    tolerance = 1e-4,
    cores = 1L,
    rng = 1L,
    seed = 1L,
    trace = FALSE
  )
}

sb_bujar_linear_parameter_semantics <- function(
  learner = sb_bujar_linear_learner_ids()
) {
  params <- sb_bujar_linear_params(learner)
  paste0(
    "Buckley-James log-time iterative imputation with componentwise ",
    "linear L2 boosting; ", params$iter_bj, " BJ iterations, ",
    params$mstop, " boosting steps, shrinkage ", params$nu,
    ", tolerance 1e-4, seed 1, and one-core execution fixed before ",
    "resampling"
  )
}

sb_bujar_stop_on_warning <- function(code, stage) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      sb_abort(
        paste0(
          "The Bujar ", stage, " raised a warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_bujar_warning",
        parent = warning
      )
    }
  )
}

sb_bujar_wrap_error <- function(code, stage, subclass) {
  tryCatch(
    force(code),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) {
        stop(error)
      }
      sb_abort(
        paste0(
          "The Bujar ", stage, " failed: ", conditionMessage(error)
        ),
        subclass,
        parent = error
      )
    }
  )
}

sb_bujar_engine_x <- function(x) {
  engine_x <- as.matrix(x)
  rownames(engine_x) <- NULL
  colnames(engine_x) <- paste0("x", seq_len(ncol(engine_x)))
  engine_x
}

sb_validate_bujar_linear_model <- function(model) {
  valid <- inherits(model, "stormbreaker_bujar_linear_model") &&
    is.list(model) &&
    identical(model$learner, sb_bujar_linear_learner_ids()) &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.numeric(model$intercept) && length(model$intercept) == 1L &&
    !is.na(model$intercept) && is.finite(model$intercept) &&
    is.numeric(model$coefficient) &&
    length(model$coefficient) == length(model$features) &&
    !anyNA(model$coefficient) && all(is.finite(model$coefficient)) &&
    identical(model$params, sb_bujar_linear_params(model$learner))
  if (!valid) {
    sb_abort(
      "The supplied model is not a fitted Bujar linear learner.",
      "stormbreaker_bad_bujar_fit"
    )
  }
  model
}

sb_fit_bujar_linear <- function(x, time, event, learner, params) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_bujar_dependency()
  engine_x <- sb_bujar_engine_x(inputs$x)
  engine <- sb_bujar_wrap_error(
    withr::with_seed(
      params$seed,
      bujar::bujar(
        y = log(inputs$time),
        cens = inputs$event,
        x = engine_x,
        degree = params$degree,
        learner = params$learner,
        center = params$center,
        mimpu = params$mimpu,
        iter.bj = params$iter_bj,
        max.cycle = params$max_cycle,
        nu = params$nu,
        mstop = params$mstop,
        twin = params$twin,
        mstop2 = params$mstop2,
        tuning = params$tuning,
        cv = params$cv,
        vimpint = params$vimpint,
        rel.inf = params$rel_inf,
        tol = params$tolerance,
        n.cores = params$cores,
        rng = params$rng,
        trace = params$trace
      )
    ),
    "linear-boosting fit",
    "stormbreaker_bujar_fit_failed"
  )
  converged <- is.list(engine) &&
    is.numeric(engine$contype) && length(engine$contype) == 1L &&
    !is.na(engine$contype) && engine$contype == 0 &&
    is.logical(engine$nonconv) && length(engine$nonconv) == 1L &&
    !is.na(engine$nonconv) && !engine$nonconv
  if (!converged) {
    sb_abort(
      paste0(
        "The Bujar linear learner did not converge (contype=",
        paste(engine$contype, collapse = ","), "; nonconv=",
        paste(engine$nonconv, collapse = ","), ")."
      ),
      "stormbreaker_bujar_nonconvergence"
    )
  }
  coefficient <- unname(as.numeric(engine$coef.bj))
  if (length(coefficient) != ncol(engine_x) + 1L ||
    anyNA(coefficient) || any(!is.finite(coefficient))) {
    sb_abort(
      "The Bujar linear learner did not return portable coefficients.",
      "stormbreaker_invalid_bujar_fit"
    )
  }
  intercept <- coefficient[[1L]]
  coefficient <- coefficient[-1L]
  portable <- unname(intercept + as.vector(engine_x %*% coefficient))
  native <- sb_bujar_wrap_error(
    as.numeric(stats::predict(engine, newx = engine_x)),
    "linear-boosting prediction",
    "stormbreaker_bujar_prediction_failed"
  )
  if (!isTRUE(all.equal(
    portable,
    unname(native),
    tolerance = 0,
    check.attributes = FALSE
  ))) {
    sb_abort(
      paste0(
        "Portable Bujar coefficients did not reproduce the native ",
        "training prediction."
      ),
      "stormbreaker_invalid_bujar_fit"
    )
  }
  structure(
    list(
      learner = learner,
      features = names(inputs$x),
      params = params,
      intercept = intercept,
      coefficient = coefficient
    ),
    class = "stormbreaker_bujar_linear_model"
  )
}

sb_validate_bujar_model_x <- function(model, x) {
  model <- sb_validate_bujar_linear_model(model)
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted Bujar learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  list(model = model, x = x)
}

sb_bujar_predicted_log_time <- function(model, x) {
  validated <- sb_validate_bujar_model_x(model, x)
  link <- unname(
    validated$model$intercept +
      as.vector(as.matrix(validated$x) %*% validated$model$coefficient)
  )
  if (length(link) != nrow(validated$x) || anyNA(link) ||
    any(!is.finite(link))) {
    sb_abort(
      "The Bujar linear learner produced invalid predicted log times.",
      "stormbreaker_invalid_bujar_prediction"
    )
  }
  link
}

sb_predict_bujar_linear <- function(model, x) {
  -sb_bujar_predicted_log_time(model, x)
}

sb_bujar_residual_curve <- function(time, event, predicted_log_time) {
  residual <- log(time) - predicted_log_time
  curve <- sb_bujar_wrap_error(
    sb_bujar_stop_on_warning(
      survival::survfit(survival::Surv(residual, event) ~ 1),
      "analysis-fold residual Kaplan-Meier calibration"
    ),
    "analysis-fold residual Kaplan-Meier calibration",
    "stormbreaker_bujar_residual_calibration_failed"
  )
  event_steps <- curve$n.event > 0
  curve_time <- unname(as.numeric(curve$time[event_steps]))
  curve_survival <- unname(as.numeric(curve$surv[event_steps]))
  valid <- length(curve_time) > 0L &&
    length(curve_time) == length(curve_survival) &&
    !anyNA(curve_time) && all(is.finite(curve_time)) &&
    !is.unsorted(curve_time) && !anyNA(curve_survival) &&
    all(is.finite(curve_survival)) &&
    all(curve_survival >= 0 & curve_survival <= 1)
  if (!valid) {
    sb_abort(
      "The Bujar residual Kaplan-Meier curve is invalid.",
      "stormbreaker_bujar_residual_calibration_failed"
    )
  }
  list(time = curve_time, survival = curve_survival)
}

sb_bujar_residual_probability <- function(
  curve_time, curve_survival, predicted_log_time, prediction_horizon
) {
  threshold <- log(prediction_horizon) - predicted_log_time
  index <- findInterval(threshold, curve_time)
  probability <- rep(1, length(index))
  use <- index > 0L
  probability[use] <- curve_survival[index[use]]
  if (anyNA(probability) || any(!is.finite(probability)) ||
    any(probability < 0 | probability > 1)) {
    sb_abort(
      "The Bujar residual calibration produced invalid survival values.",
      "stormbreaker_invalid_bujar_prediction"
    )
  }
  unname(probability)
}

sb_bujar_direct_survival_source <- function() {
  "bujar_linear_boost_analysis_residual_kaplan_meier"
}

sb_predict_bujar_linear_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  model <- sb_validate_bujar_linear_model(model)
  inputs <- sb_validate_additional_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted Bujar learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  analysis_link <- sb_bujar_predicted_log_time(model, inputs$x)
  assessment_link <- sb_bujar_predicted_log_time(model, x_assessment)
  curve <- sb_bujar_residual_curve(
    inputs$time,
    inputs$event,
    analysis_link
  )
  list(
    prediction_horizon = prediction_horizon,
    source = sb_bujar_direct_survival_source(),
    ties = "not_applicable",
    probability = sb_bujar_residual_probability(
      curve$time,
      curve$survival,
      assessment_link,
      prediction_horizon
    ),
    residual_time = curve$time,
    residual_survival = curve$survival
  )
}

sb_extract_bujar_linear_state <- function(model, context = list()) {
  model <- sb_validate_bujar_linear_model(model)
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_bujar_linear_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
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
      ties = probability$ties,
      residual_time = probability$residual_time,
      residual_survival = probability$residual_survival
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      intercept = model$intercept,
      coefficient = model$coefficient,
      params = model$params,
      risk_reference = "bujar_native",
      risk_formula = "negative_predicted_log_time",
      survival_probability = survival_probability
    ),
    class = c(
      "stormbreaker_bujar_linear_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_validate_bujar_probability_state <- function(state) {
  probability <- state$survival_probability
  valid <- is.list(probability) &&
    is.numeric(probability$prediction_horizon) &&
    length(probability$prediction_horizon) == 1L &&
    !is.na(probability$prediction_horizon) &&
    is.finite(probability$prediction_horizon) &&
    probability$prediction_horizon > 0 &&
    identical(probability$source, sb_bujar_direct_survival_source()) &&
    identical(probability$ties, "not_applicable") &&
    is.numeric(probability$residual_time) &&
    length(probability$residual_time) > 0L &&
    length(probability$residual_time) ==
      length(probability$residual_survival) &&
    !anyNA(probability$residual_time) &&
    all(is.finite(probability$residual_time)) &&
    !is.unsorted(probability$residual_time) &&
    is.numeric(probability$residual_survival) &&
    !anyNA(probability$residual_survival) &&
    all(is.finite(probability$residual_survival)) &&
    all(probability$residual_survival >= 0 &
      probability$residual_survival <= 1)
  if (!valid) {
    sb_bad_learner_state("The saved Bujar horizon state is invalid.")
  }
  probability
}

sb_predict_bujar_linear_state <- function(state, x) {
  learner <- if (is.list(state)) state$learner else NULL
  if (!inherits(state, "stormbreaker_bujar_linear_state") ||
    !identical(learner, sb_bujar_linear_learner_ids())) {
    sb_bad_learner_state("The saved Bujar learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  valid <- is.numeric(state$intercept) && length(state$intercept) == 1L &&
    !is.na(state$intercept) && is.finite(state$intercept) &&
    is.numeric(state$coefficient) &&
    length(state$coefficient) == ncol(x) &&
    !anyNA(state$coefficient) && all(is.finite(state$coefficient))
  if (!valid) {
    sb_bad_learner_state("The saved Bujar coefficients are invalid.")
  }
  link <- unname(
    state$intercept + as.vector(as.matrix(x) %*% state$coefficient)
  )
  if (anyNA(link) || any(!is.finite(link))) {
    sb_abort(
      "The saved Bujar state produced invalid predictions.",
      "stormbreaker_invalid_bujar_prediction"
    )
  }
  result <- list(risk = -link)
  if (!is.null(state$survival_probability)) {
    probability <- sb_validate_bujar_probability_state(state)
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- sb_bujar_residual_probability(
      probability$residual_time,
      probability$residual_survival,
      link,
      probability$prediction_horizon
    )
  }
  result
}

sb_bujar_linear_learner_adapter <- function(
  learner = sb_bujar_linear_learner_ids(),
  params = sb_bujar_linear_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_bujar_linear(x, time, event, learner, params)
    },
    predict = sb_predict_bujar_linear,
    predict_survival = sb_predict_bujar_linear_survival,
    extract_state = sb_extract_bujar_linear_state,
    predict_state = sb_predict_bujar_linear_state
  )
}
