sb_empty_survival_weak_calibration_patients <- function() {
  data.frame(
    row_instance = integer(),
    patient_id = character(),
    fold = integer(),
    time = numeric(),
    event = integer(),
    survival_probability = numeric(),
    predicted_event_probability = numeric(),
    transformed_prediction = numeric(),
    administrative_time = numeric(),
    event_by_horizon = integer(),
    stringsAsFactors = FALSE
  )
}

sb_empty_survival_weak_calibration_folds <- function() {
  sb_empty_survival_calibration_folds()
}

sb_survival_weak_calibration_empty_mean <- function() {
  list(
    status = "not_evaluated",
    method = "kaplan_meier_observed_over_mean_prediction",
    ratio_direction = "observed_over_expected",
    horizon_support = "not_checked",
    observed_event_probability = NA_real_,
    mean_predicted_event_probability = NA_real_,
    mean_calibration_ratio = NA_real_,
    n_followed_to_horizon = 0L,
    maximum_follow_up = NA_real_,
    km_survival_probability = NA_real_,
    event_curve = sb_empty_event_curve(),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_survival_weak_calibration_empty_slope <- function() {
  list(
    status = "not_fitted",
    method = "secondary_cox_proportional_hazards",
    formula = paste0(
      "Surv(administrative_time, event_by_horizon) ~ ",
      "transformed_prediction"
    ),
    ties = "efron",
    administrative_censoring_at_horizon = TRUE,
    proportional_hazards_assumption = TRUE,
    censoring_assumption = "independent_noninformative_censoring",
    coefficient = NA_real_,
    n_events_by_horizon = 0L,
    n_unique_transformed_prediction = 0L,
    iterations = NA_integer_,
    log_likelihood = c(initial = NA_real_, final = NA_real_),
    warnings = character(),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_survival_weak_calibration_empty_result <- function(
    prediction_horizon,
    time_unit,
    scope,
    n = 0L) {
  safe_horizon <- if (is.numeric(prediction_horizon) &&
      length(prediction_horizon) >= 1L) {
    as.numeric(prediction_horizon[[1L]])
  } else {
    NA_real_
  }
  safe_unit <- if (is.character(time_unit) && length(time_unit) >= 1L) {
    as.character(time_unit[[1L]])
  } else {
    NA_character_
  }
  safe_scope <- if (is.character(scope) && length(scope) >= 1L) {
    as.character(scope[[1L]])
  } else {
    NA_character_
  }
  list(
    metric = "survival_weak_calibration",
    calibration_level = "weak",
    estimand = "fixed_horizon_mean_calibration_and_slope",
    scope = safe_scope,
    purpose = "descriptive_evaluation_only",
    prediction_horizon = safe_horizon,
    time_unit = safe_unit,
    prediction_scale = "survival_probability",
    event_probability_transform = "complementary_log_log",
    transform_formula = "log(-log(S(h)))",
    probability_clipping = "none",
    status = "incomplete",
    input_status = "not_checked",
    mean_calibration_ratio = NA_real_,
    calibration_slope = NA_real_,
    mean_calibration_target = 1,
    slope_target = 1,
    n = as.integer(n),
    n_scored = 0L,
    n_events = 0L,
    n_censored = 0L,
    n_events_by_horizon = 0L,
    n_boundary_probability = 0L,
    n_followed_to_horizon = 0L,
    n_folds = 0L,
    n_evaluable_folds = 0L,
    mean_calibration = sb_survival_weak_calibration_empty_mean(),
    slope_fit = sb_survival_weak_calibration_empty_slope(),
    patient_predictions = sb_empty_survival_weak_calibration_patients(),
    per_fold = sb_empty_survival_weak_calibration_folds(),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_survival_weak_calibration_error_metric <- function(
    error,
    prediction_horizon,
    time_unit,
    scope,
    n,
    fold = NULL) {
  result <- sb_survival_weak_calibration_empty_result(
    prediction_horizon,
    time_unit,
    scope,
    n
  )
  if (!is.null(fold) && is.numeric(fold)) {
    fold_ids <- sort(unique(as.integer(fold[is.finite(fold)])))
    result$n_folds <- as.integer(length(fold_ids))
  }
  result$status <- "metric_failed"
  result$input_status <- "unexpected_error"
  result$error_class <- class(error)[[1L]]
  result$error_message <- conditionMessage(error)
  result
}

sb_survival_weak_calibration_patient_table <- function(
    patient_id,
    time,
    event,
    survival_probability,
    fold,
    prediction_horizon,
    log_cumulative_hazard = NULL) {
  patients <- sb_survival_calibration_patient_table(
    patient_id,
    time,
    event,
    survival_probability,
    fold,
    log_cumulative_hazard
  )
  n <- nrow(patients)
  valid_horizon <- is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) && is.finite(prediction_horizon)
  valid_time <- is.numeric(time) && length(time) == n
  valid_event <- (is.numeric(event) || is.logical(event)) &&
    length(event) == n
  patients$administrative_time <- if (valid_horizon && valid_time) {
    pmin(as.numeric(time), as.numeric(prediction_horizon))
  } else {
    rep(NA_real_, n)
  }
  patients$event_by_horizon <- if (
      valid_horizon && valid_time && valid_event) {
    as.integer(as.integer(event) == 1L & time <= prediction_horizon)
  } else {
    rep(NA_integer_, n)
  }
  patients[c(
    "row_instance",
    "patient_id",
    "fold",
    "time",
    "event",
    "survival_probability",
    "predicted_event_probability",
    "transformed_prediction",
    "administrative_time",
    "event_by_horizon"
  )]
}

sb_survival_weak_calibration_fit_slope <- function(fit_data) {
  survival::coxph(
    survival::Surv(administrative_time, event_by_horizon) ~
      transformed_prediction,
    data = fit_data,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
}

sb_survival_weak_calibration_component_error <- function(
    component,
    status,
    error_class,
    error_message) {
  component$status <- status
  component$error_class <- error_class
  component$error_message <- error_message
  component
}

sb_survival_weak_calibration_finish <- function(result) {
  if (identical(result$mean_calibration$status, "success") &&
      identical(result$slope_fit$status, "success")) {
    result$status <- "success"
    return(result)
  }
  failed_component <- if (!identical(
      result$mean_calibration$status,
      "success"
    )) {
    result$mean_calibration
  } else {
    result$slope_fit
  }
  result$status <- "incomplete"
  result$error_class <- failed_component$error_class
  result$error_message <- failed_component$error_message
  result
}

sb_survival_weak_calibration_metric <- function(
    patient_id,
    time,
    event,
    survival_probability,
    prediction_status,
    prediction_horizon,
    time_unit,
    scope,
    fold = NULL,
    log_cumulative_hazard = NULL) {
  n <- length(time)
  result <- sb_survival_weak_calibration_empty_result(
    prediction_horizon,
    time_unit,
    scope,
    n
  )
  patients <- sb_survival_weak_calibration_patient_table(
    patient_id,
    time,
    event,
    survival_probability,
    fold,
    prediction_horizon,
    log_cumulative_hazard
  )
  result$patient_predictions <- patients
  if (!is.null(log_cumulative_hazard)) result$transform_formula <- "saved_log_cumulative_hazard"

  valid_horizon <- is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) && is.finite(prediction_horizon) &&
    prediction_horizon > 0
  valid_unit <- is.character(time_unit) && length(time_unit) == 1L &&
    !is.na(time_unit) && time_unit %in% c("day", "month", "year")
  valid_scope <- is.character(scope) && length(scope) == 1L &&
    !is.na(scope) && nzchar(scope)
  valid_patient_id <- n > 0L && is.character(patient_id) &&
    length(patient_id) == n && !anyNA(patient_id) &&
    all(nzchar(trimws(patient_id)))
  valid_endpoint <- is.numeric(time) && length(time) == n &&
    (is.numeric(event) || is.logical(event)) && length(event) == n &&
    !anyNA(time) && !anyNA(event) &&
    all(is.finite(time)) && all(time >= 0) &&
    (!is.numeric(event) || all(is.finite(event))) &&
    all(event %in% c(0, 1))
  valid_probability <- is.numeric(survival_probability) &&
    length(survival_probability) == n &&
    !anyNA(survival_probability) &&
    all(is.finite(survival_probability)) &&
    all(survival_probability >= 0 & survival_probability <= 1)
  valid_status <- is.character(prediction_status) &&
    length(prediction_status) == n && !anyNA(prediction_status)
  fold_supplied <- !is.null(fold)
  valid_fold <- !fold_supplied || (
    is.numeric(fold) && length(fold) == n && !anyNA(fold) &&
      all(is.finite(fold)) && all(fold == trunc(fold))
  )

  if (valid_horizon && valid_endpoint) {
    result$n_followed_to_horizon <- as.integer(sum(time >= prediction_horizon))
    result$n_events_by_horizon <- as.integer(sum(
      event == 1L & time <= prediction_horizon
    ))
  }
  if (valid_endpoint) {
    result$n_events <- as.integer(sum(event == 1L))
    result$n_censored <- as.integer(sum(event == 0L))
  }
  if (valid_probability) {
    result$n_boundary_probability <- as.integer(sum(
      survival_probability %in% c(0, 1)
    ))
  }
  result$per_fold <- if (valid_fold && valid_horizon) {
    sb_survival_calibration_fold_audit(
      patients,
      prediction_status,
      prediction_horizon,
      fold_supplied
    )
  } else {
    sb_empty_survival_weak_calibration_folds()
  }
  result$n_folds <- as.integer(nrow(result$per_fold))
  result$n_evaluable_folds <- as.integer(sum(
    result$per_fold$status == "success"
  ))

  if (!valid_horizon || !valid_unit || !valid_scope ||
      !valid_patient_id || !valid_endpoint || !valid_probability ||
      !valid_status || !valid_fold) {
    result$input_status <- "input_failed"
    result$error_class <- "stormbreaker_bad_survival_weak_calibration_input"
    result$error_message <- paste0(
      "Survival weak calibration requires valid aligned endpoint, ",
      "probability, horizon, unit, scope, status, and optional fold vectors."
    )
    return(result)
  }
  if (any(prediction_status != "success")) {
    result$input_status <- "prediction_failed"
    result$error_class <-
      "stormbreaker_invalid_survival_weak_calibration_prediction"
    result$error_message <- paste0(
      "Complete successful survival probabilities are required for weak ",
      "calibration."
    )
    return(result)
  }

  result$input_status <- "success"
  result$n_scored <- as.integer(n)

  event_curve <- sb_fit_event_km(time, event)
  maximum_follow_up <- max(time)
  terminal_survival <- utils::tail(event_curve$survival_after, 1L)
  horizon_support <- if (maximum_follow_up >= prediction_horizon) {
    "observed_to_horizon"
  } else if (terminal_survival == 0) {
    "survival_zero_before_horizon"
  } else {
    "unsupported"
  }
  result$mean_calibration$event_curve <- event_curve
  result$mean_calibration$maximum_follow_up <- maximum_follow_up
  result$mean_calibration$n_followed_to_horizon <-
    result$n_followed_to_horizon
  result$mean_calibration$horizon_support <- horizon_support
  if (horizon_support == "unsupported") {
    result$mean_calibration <- sb_survival_weak_calibration_component_error(
      result$mean_calibration,
      "not_evaluable",
      "stormbreaker_survival_weak_calibration_horizon_unsupported",
      "The prediction horizon exceeds all observed endpoint follow-up."
    )
  } else {
    km_survival <- sb_event_survival_at(event_curve, prediction_horizon)
    observed_event_probability <- -expm1(log(km_survival))
    mean_predicted_event_probability <- mean(
      patients$predicted_event_probability
    )
    result$mean_calibration$km_survival_probability <- km_survival
    result$mean_calibration$observed_event_probability <-
      observed_event_probability
    result$mean_calibration$mean_predicted_event_probability <-
      mean_predicted_event_probability
    if (mean_predicted_event_probability == 0) {
      result$mean_calibration <- sb_survival_weak_calibration_component_error(
        result$mean_calibration,
        "not_evaluable",
        "stormbreaker_survival_weak_calibration_zero_mean_prediction",
        "The mean predicted event probability is zero."
      )
    } else {
      ratio <- observed_event_probability / mean_predicted_event_probability
      if (is.finite(ratio)) {
        result$mean_calibration$status <- "success"
        result$mean_calibration$mean_calibration_ratio <- unname(ratio)
        result$mean_calibration_ratio <- unname(ratio)
      } else {
        result$mean_calibration <- sb_survival_weak_calibration_component_error(
          result$mean_calibration,
          "not_evaluable",
          "stormbreaker_survival_weak_calibration_nonfinite_ratio",
          "The observed-to-expected calibration ratio is not finite."
        )
      }
    }
  }

  transformed_prediction <- patients$transformed_prediction
  result$slope_fit$n_events_by_horizon <- result$n_events_by_horizon
  result$slope_fit$n_unique_transformed_prediction <- as.integer(length(
    unique(transformed_prediction[is.finite(transformed_prediction)])
  ))
  if (result$n_boundary_probability > 0L && is.null(log_cumulative_hazard)) {
    result$slope_fit <- sb_survival_weak_calibration_component_error(
      result$slope_fit,
      "not_evaluable",
      "stormbreaker_survival_weak_calibration_boundary_probability",
      paste0(
        "Exact survival probabilities zero or one do not have a finite ",
        "complementary-log-log value; probabilities are not clipped."
      )
    )
    return(sb_survival_weak_calibration_finish(result))
  }
  if (result$n_events_by_horizon == 0L) {
    result$slope_fit <- sb_survival_weak_calibration_component_error(
      result$slope_fit,
      "not_evaluable",
      "stormbreaker_survival_weak_calibration_no_event",
      "At least one observed event by the prediction horizon is required."
    )
    return(sb_survival_weak_calibration_finish(result))
  }
  if (result$slope_fit$n_unique_transformed_prediction < 2L) {
    result$slope_fit <- sb_survival_weak_calibration_component_error(
      result$slope_fit,
      "not_evaluable",
      "stormbreaker_survival_weak_calibration_constant_prediction",
      "The transformed survival prediction is constant."
    )
    return(sb_survival_weak_calibration_finish(result))
  }

  fit_data <- data.frame(
    administrative_time = patients$administrative_time,
    event_by_horizon = patients$event_by_horizon,
    transformed_prediction = transformed_prediction
  )
  captured_warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      sb_survival_weak_calibration_fit_slope(fit_data),
      warning = function(warning) {
        captured_warnings <<- c(
          captured_warnings,
          conditionMessage(warning)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) error
  )
  result$slope_fit$warnings <- unique(captured_warnings)
  if (inherits(fit, "condition")) {
    result$slope_fit <- sb_survival_weak_calibration_component_error(
      result$slope_fit,
      "failed",
      class(fit)[[1L]],
      conditionMessage(fit)
    )
    return(sb_survival_weak_calibration_finish(result))
  }
  if (length(captured_warnings) > 0L) {
    result$slope_fit <- sb_survival_weak_calibration_component_error(
      result$slope_fit,
      "failed",
      "stormbreaker_survival_weak_calibration_cox_warning",
      paste(unique(captured_warnings), collapse = " ")
    )
    return(sb_survival_weak_calibration_finish(result))
  }

  coefficients <- stats::coef(fit)
  coefficient <- if (
      "transformed_prediction" %in% names(coefficients)) {
    unname(coefficients[["transformed_prediction"]])
  } else {
    NA_real_
  }
  result$slope_fit$iterations <- if (!is.null(fit$iter) &&
      length(fit$iter) >= 1L) {
    as.integer(fit$iter[[1L]])
  } else {
    NA_integer_
  }
  if (!is.null(fit$loglik) && length(fit$loglik) >= 2L) {
    result$slope_fit$log_likelihood <- c(
      initial = unname(fit$loglik[[1L]]),
      final = unname(fit$loglik[[2L]])
    )
  }
  if (length(coefficient) != 1L || is.na(coefficient) ||
      !is.finite(coefficient)) {
    result$slope_fit <- sb_survival_weak_calibration_component_error(
      result$slope_fit,
      "failed",
      "stormbreaker_survival_weak_calibration_singular_fit",
      "The secondary Cox fit did not return a finite named coefficient."
    )
    return(sb_survival_weak_calibration_finish(result))
  }

  result$slope_fit$status <- "success"
  result$slope_fit$coefficient <- coefficient
  result$calibration_slope <- coefficient
  sb_survival_weak_calibration_finish(result)
}
