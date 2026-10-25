sb_empty_survival_calibration_knots <- function() {
  data.frame(
    role = character(),
    quantile_probability = numeric(),
    transformed_value = numeric(),
    survival_probability = numeric(),
    event_probability = numeric(),
    stringsAsFactors = FALSE
  )
}

sb_empty_survival_calibration_patients <- function() {
  data.frame(
    row_instance = integer(),
    patient_id = character(),
    fold = integer(),
    time = numeric(),
    event = integer(),
    survival_probability = numeric(),
    predicted_event_probability = numeric(),
    transformed_prediction = numeric(),
    calibrated_survival_probability = numeric(),
    calibrated_event_probability = numeric(),
    absolute_error = numeric(),
    stringsAsFactors = FALSE
  )
}

sb_empty_survival_calibration_curve <- function() {
  data.frame(
    survival_probability = numeric(),
    predicted_event_probability = numeric(),
    transformed_prediction = numeric(),
    calibrated_survival_probability = numeric(),
    calibrated_event_probability = numeric(),
    n_row_instances = integer(),
    stringsAsFactors = FALSE
  )
}

sb_empty_survival_calibration_folds <- function() {
  data.frame(
    fold = integer(),
    status = character(),
    n = integer(),
    n_scored = integer(),
    n_events = integer(),
    n_censored = integer(),
    n_boundary_probability = integer(),
    n_followed_to_horizon = integer(),
    error_class = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
}

sb_survival_calibration_empty_fit <- function() {
  list(
    status = "not_fitted",
    method = "cox_ph_three_knot_restricted_natural_cubic_spline",
    formula = paste0(
      "Surv(time, event) ~ calibration_ns1 + calibration_ns2"
    ),
    ties = "efron",
    baseline_reference = "coxph_mean_basis_centered",
    coefficients = data.frame(
      term = c("calibration_ns1", "calibration_ns2"),
      estimate = c(NA_real_, NA_real_),
      stringsAsFactors = FALSE
    ),
    means = c(calibration_ns1 = NA_real_, calibration_ns2 = NA_real_),
    baseline_cumulative_hazard_horizon = NA_real_,
    iterations = NA_integer_,
    log_likelihood = c(initial = NA_real_, final = NA_real_),
    warnings = character(),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_survival_calibration_empty_result <- function(
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
    metric = "survival_calibration_curve",
    calibration_level = "moderate",
    estimand = "fixed_horizon_smoothed_calibration",
    method = "cox_ph_three_knot_restricted_natural_cubic_spline",
    scope = safe_scope,
    purpose = "descriptive_evaluation_only",
    prediction_horizon = safe_horizon,
    time_unit = safe_unit,
    prediction_scale = "survival_probability",
    event_probability_transform = "complementary_log_log",
    transform_formula = "log(-log(S(h)))",
    probability_clipping = "none",
    spline = list(
      type = "restricted_natural_cubic_spline",
      n_knots = 3L,
      knot_probabilities = c(0.10, 0.50, 0.90),
      quantile_type = 7L,
      outcome_tuning = "none"
    ),
    censoring_handling = "cox_partial_likelihood_complete_follow_up",
    censoring_assumption = "independent_noninformative_censoring",
    proportional_hazards_assumption = TRUE,
    administrative_censoring_at_horizon = FALSE,
    ties = "efron",
    baseline_reference = "coxph_mean_basis_centered",
    status = "incomplete",
    input_status = "not_checked",
    ici = NA_real_,
    e50 = NA_real_,
    e90 = NA_real_,
    n = as.integer(n),
    n_scored = 0L,
    n_events = 0L,
    n_censored = 0L,
    n_boundary_probability = 0L,
    n_followed_to_horizon = 0L,
    n_unique_survival_probability = 0L,
    n_folds = 0L,
    n_evaluable_folds = 0L,
    knots = sb_empty_survival_calibration_knots(),
    fit = sb_survival_calibration_empty_fit(),
    curve = sb_empty_survival_calibration_curve(),
    patient_calibration = sb_empty_survival_calibration_patients(),
    per_fold = sb_empty_survival_calibration_folds(),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_survival_calibration_fail <- function(
    result,
    status,
    input_status,
    error_class,
    error_message) {
  result$status <- status
  result$input_status <- input_status
  result$error_class <- error_class
  result$error_message <- error_message
  result
}

sb_survival_calibration_error_metric <- function(
    error,
    prediction_horizon,
    time_unit,
    scope,
    n,
    fold = NULL) {
  result <- sb_survival_calibration_empty_result(
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

sb_survival_calibration_patient_table <- function(
    patient_id,
    time,
    event,
    survival_probability,
    fold,
    log_cumulative_hazard = NULL) {
  n <- length(time)
  safe_patient_id <- if (is.character(patient_id) &&
      length(patient_id) == n) {
    patient_id
  } else {
    rep(NA_character_, n)
  }
  safe_event <- if ((is.numeric(event) || is.logical(event)) &&
      length(event) == n) {
    as.integer(event)
  } else {
    rep(NA_integer_, n)
  }
  safe_probability <- if (is.numeric(survival_probability) &&
      length(survival_probability) == n) {
    as.numeric(survival_probability)
  } else {
    rep(NA_real_, n)
  }
  safe_fold <- if (is.null(fold)) {
    rep(NA_integer_, n)
  } else if (is.numeric(fold) && length(fold) == n) {
    as.integer(fold)
  } else {
    rep(NA_integer_, n)
  }
  valid_probability <- !is.na(safe_probability) &
    is.finite(safe_probability) &
    safe_probability >= 0 & safe_probability <= 1
  event_probability <- rep(NA_real_, n)
  transformed <- rep(NA_real_, n)
  event_probability[valid_probability] <- -expm1(
    log(safe_probability[valid_probability])
  )
  transformed[valid_probability] <- log(
    -log(safe_probability[valid_probability])
  )
  if (!is.null(log_cumulative_hazard)) {
    transformed <- sb_validate_stable_survival_scale(
      log_cumulative_hazard, safe_probability)
  }
  data.frame(
    row_instance = seq_len(n),
    patient_id = safe_patient_id,
    fold = safe_fold,
    time = if (is.numeric(time) && length(time) == n) {
      as.numeric(time)
    } else {
      rep(NA_real_, n)
    },
    event = safe_event,
    survival_probability = safe_probability,
    predicted_event_probability = event_probability,
    transformed_prediction = transformed,
    calibrated_survival_probability = rep(NA_real_, n),
    calibrated_event_probability = rep(NA_real_, n),
    absolute_error = rep(NA_real_, n),
    stringsAsFactors = FALSE
  )
}

sb_survival_calibration_fold_audit <- function(
    patients,
    prediction_status,
    prediction_horizon,
    fold_supplied) {
  if (!fold_supplied) {
    return(sb_empty_survival_calibration_folds())
  }
  fold_ids <- sort(unique(patients$fold))
  rows <- lapply(fold_ids, function(fold_id) {
    index <- patients$fold == fold_id
    valid_status <- is.character(prediction_status) &&
      length(prediction_status) == nrow(patients) &&
      !anyNA(prediction_status[index]) &&
      all(prediction_status[index] == "success")
    valid_rows <- valid_status &&
      !anyNA(patients$time[index]) &&
      all(is.finite(patients$time[index])) &&
      all(patients$time[index] >= 0) &&
      !anyNA(patients$event[index]) &&
      all(patients$event[index] %in% c(0L, 1L)) &&
      !anyNA(patients$survival_probability[index]) &&
      all(is.finite(patients$survival_probability[index])) &&
      all(
        patients$survival_probability[index] >= 0 &
          patients$survival_probability[index] <= 1
      )
    data.frame(
      fold = as.integer(fold_id),
      status = if (valid_rows) "success" else "prediction_failed",
      n = as.integer(sum(index)),
      n_scored = if (valid_rows) as.integer(sum(index)) else 0L,
      n_events = as.integer(sum(patients$event[index] == 1L, na.rm = TRUE)),
      n_censored = as.integer(sum(
        patients$event[index] == 0L,
        na.rm = TRUE
      )),
      n_boundary_probability = as.integer(sum(
        patients$survival_probability[index] %in% c(0, 1),
        na.rm = TRUE
      )),
      n_followed_to_horizon = as.integer(sum(
        patients$time[index] >= prediction_horizon,
        na.rm = TRUE
      )),
      error_class = if (valid_rows) {
        NA_character_
      } else {
        "stormbreaker_invalid_survival_calibration_prediction"
      },
      error_message = if (valid_rows) {
        NA_character_
      } else {
        "OOF survival endpoint or probability rows are failed or invalid."
      },
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

sb_survival_calibration_knots <- function(transformed_prediction) {
  values <- unname(stats::quantile(
    transformed_prediction,
    probs = c(0.10, 0.50, 0.90),
    names = FALSE,
    type = 7
  ))
  cumulative_hazard <- exp(values)
  data.frame(
    role = c("lower_boundary", "interior", "upper_boundary"),
    quantile_probability = c(0.10, 0.50, 0.90),
    transformed_value = values,
    survival_probability = exp(-cumulative_hazard),
    event_probability = -expm1(-cumulative_hazard),
    stringsAsFactors = FALSE
  )
}

sb_survival_calibration_basis <- function(transformed_prediction, knots) {
  basis <- splines::ns(
    transformed_prediction,
    knots = knots[[2L]],
    Boundary.knots = knots[c(1L, 3L)],
    intercept = FALSE
  )
  colnames(basis) <- c("calibration_ns1", "calibration_ns2")
  basis
}

sb_survival_calibration_fit_cox <- function(fit_data) {
  survival::coxph(
    survival::Surv(time, event) ~ calibration_ns1 + calibration_ns2,
    data = fit_data,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
}

sb_survival_probability_from_log_cumulative_hazard <- function(value) {
  exp(-exp(value))
}

sb_survival_calibration_curve_table <- function(patients) {
  survival_values <- sort(
    unique(patients$survival_probability),
    decreasing = TRUE
  )
  group <- match(patients$survival_probability, survival_values)
  first <- match(survival_values, patients$survival_probability)
  data.frame(
    survival_probability = survival_values,
    predicted_event_probability =
      patients$predicted_event_probability[first],
    transformed_prediction = patients$transformed_prediction[first],
    calibrated_survival_probability =
      patients$calibrated_survival_probability[first],
    calibrated_event_probability =
      patients$calibrated_event_probability[first],
    n_row_instances = as.integer(tabulate(
      group,
      nbins = length(survival_values)
    )),
    stringsAsFactors = FALSE
  )
}

sb_survival_calibration_metric <- function(
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
  result <- sb_survival_calibration_empty_result(
    prediction_horizon,
    time_unit,
    scope,
    n
  )
  patients <- sb_survival_calibration_patient_table(
    patient_id,
    time,
    event,
    survival_probability,
    fold,
    log_cumulative_hazard
  )
  result$patient_calibration <- patients
  if (!is.null(log_cumulative_hazard)) result$transform_formula <- "saved_log_cumulative_hazard"

  valid_horizon <- is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) && is.finite(prediction_horizon) &&
    prediction_horizon > 0
  valid_unit <- is.character(time_unit) && length(time_unit) == 1L &&
    !is.na(time_unit) && time_unit %in% c("day", "month", "year")
  valid_scope <- is.character(scope) && length(scope) == 1L &&
    !is.na(scope) && nzchar(scope)
  valid_patient_id <- is.character(patient_id) && length(patient_id) == n &&
    !anyNA(patient_id) && all(nzchar(trimws(patient_id)))
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

  if (valid_horizon) {
    result$n_followed_to_horizon <- as.integer(sum(
      patients$time >= prediction_horizon,
      na.rm = TRUE
    ))
  }
  if (valid_endpoint) {
    result$n_events <- as.integer(sum(as.integer(event) == 1L))
    result$n_censored <- as.integer(sum(as.integer(event) == 0L))
  }
  if (valid_probability) {
    result$n_boundary_probability <- as.integer(sum(
      survival_probability %in% c(0, 1)
    ))
    result$n_unique_survival_probability <- as.integer(length(unique(
      survival_probability
    )))
  }
  result$per_fold <- if (valid_fold && valid_horizon) {
    sb_survival_calibration_fold_audit(
      patients,
      prediction_status,
      prediction_horizon,
      fold_supplied
    )
  } else {
    sb_empty_survival_calibration_folds()
  }
  result$n_folds <- as.integer(nrow(result$per_fold))
  result$n_evaluable_folds <- as.integer(sum(
    result$per_fold$status == "success"
  ))

  if (!valid_horizon || !valid_unit || !valid_scope ||
      !valid_patient_id || !valid_endpoint || !valid_probability ||
      !valid_status || !valid_fold) {
    return(sb_survival_calibration_fail(
      result,
      "incomplete",
      "input_failed",
      "stormbreaker_bad_survival_calibration_input",
      paste0(
        "Survival calibration requires valid aligned endpoint, probability, ",
        "horizon, unit, scope, status, and optional fold vectors."
      )
    ))
  }
  if (any(prediction_status != "success")) {
    return(sb_survival_calibration_fail(
      result,
      "incomplete",
      "prediction_failed",
      "stormbreaker_invalid_survival_calibration_prediction",
      "Complete successful survival probabilities are required for calibration."
    ))
  }
  result$n_scored <- as.integer(n)
  if (result$n_events == 0L) {
    return(sb_survival_calibration_fail(
      result,
      "incomplete",
      "not_evaluable",
      "stormbreaker_survival_calibration_no_event",
      "At least one observed event is required for survival calibration."
    ))
  }
  if (result$n_followed_to_horizon == 0L) {
    return(sb_survival_calibration_fail(
      result,
      "incomplete",
      "not_evaluable",
      "stormbreaker_survival_calibration_horizon_unsupported",
      paste0(
        "The prediction horizon exceeds all observed endpoint follow-up."
      )
    ))
  }
  if (result$n_boundary_probability > 0L && is.null(log_cumulative_hazard)) {
    return(sb_survival_calibration_fail(
      result,
      "incomplete",
      "not_evaluable",
      "stormbreaker_survival_calibration_boundary_probability",
      paste0(
        "Exact survival probabilities zero or one cannot enter a finite ",
        "complementary-log-log calibration fit without clipping, which is ",
        "not performed."
      )
    ))
  }

  transformed_prediction <- patients$transformed_prediction
  knots <- sb_survival_calibration_knots(transformed_prediction)
  result$knots <- knots
  knot_values <- knots$transformed_value
  if (anyNA(knot_values) || any(!is.finite(knot_values)) ||
      any(diff(knot_values) <= 0)) {
    return(sb_survival_calibration_fail(
      result,
      "incomplete",
      "not_evaluable",
      "stormbreaker_survival_calibration_duplicate_knots",
      paste0(
        "The fixed 10th, 50th, and 90th percentile spline knots are not ",
        "finite and strictly increasing."
      )
    ))
  }

  basis_attempt <- tryCatch(
    sb_survival_calibration_basis(transformed_prediction, knot_values),
    error = function(error) error
  )
  if (inherits(basis_attempt, "condition")) {
    return(sb_survival_calibration_fail(
      result,
      "fit_failed",
      "success",
      class(basis_attempt)[[1L]],
      conditionMessage(basis_attempt)
    ))
  }
  basis <- basis_attempt
  fit_data <- data.frame(
    time = as.numeric(time),
    event = as.integer(event),
    calibration_ns1 = basis[, 1L],
    calibration_ns2 = basis[, 2L]
  )
  captured_warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      sb_survival_calibration_fit_cox(fit_data),
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
  result$fit$warnings <- unique(captured_warnings)
  if (inherits(fit, "condition")) {
    result$fit$status <- "failed"
    result$fit$error_class <- class(fit)[[1L]]
    result$fit$error_message <- conditionMessage(fit)
    result$status <- "fit_failed"
    result$input_status <- "success"
    result$error_class <- result$fit$error_class
    result$error_message <- result$fit$error_message
    return(result)
  }

  coefficients <- stats::coef(fit)
  means <- fit$means
  result$fit$coefficients$estimate <- unname(coefficients[
    c("calibration_ns1", "calibration_ns2")
  ])
  result$fit$means <- unname(means[
    c("calibration_ns1", "calibration_ns2")
  ])
  names(result$fit$means) <- c("calibration_ns1", "calibration_ns2")
  result$fit$iterations <- as.integer(fit$iter[[1L]])
  if (length(fit$loglik) >= 2L) {
    result$fit$log_likelihood <- c(
      initial = unname(fit$loglik[[1L]]),
      final = unname(fit$loglik[[2L]])
    )
  }
  if (length(captured_warnings) > 0L) {
    result$fit$status <- "failed"
    result$fit$error_class <- "stormbreaker_survival_calibration_cox_warning"
    result$fit$error_message <- paste(unique(captured_warnings), collapse = " ")
    result$status <- "fit_failed"
    result$input_status <- "success"
    result$error_class <- result$fit$error_class
    result$error_message <- result$fit$error_message
    return(result)
  }
  valid_fit_parameters <- length(coefficients) == 2L &&
    all(c("calibration_ns1", "calibration_ns2") %in% names(coefficients)) &&
    !anyNA(coefficients) && all(is.finite(coefficients)) &&
    length(means) == 2L &&
    all(c("calibration_ns1", "calibration_ns2") %in% names(means)) &&
    !anyNA(means) && all(is.finite(means))
  if (!valid_fit_parameters) {
    result$fit$status <- "failed"
    result$fit$error_class <- "stormbreaker_survival_calibration_singular_fit"
    result$fit$error_message <- paste0(
      "The calibration Cox fit returned missing, singular, or non-finite ",
      "basis parameters."
    )
    result$status <- "fit_failed"
    result$input_status <- "success"
    result$error_class <- result$fit$error_class
    result$error_message <- result$fit$error_message
    return(result)
  }

  baseline <- tryCatch(
    withCallingHandlers(
      survival::basehaz(fit, centered = TRUE),
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
  result$fit$warnings <- unique(captured_warnings)
  if (inherits(baseline, "condition") || length(captured_warnings) > 0L) {
    result$fit$status <- "failed"
    result$fit$error_class <- if (inherits(baseline, "condition")) {
      class(baseline)[[1L]]
    } else {
      "stormbreaker_survival_calibration_baseline_warning"
    }
    result$fit$error_message <- if (inherits(baseline, "condition")) {
      conditionMessage(baseline)
    } else {
      paste(unique(captured_warnings), collapse = " ")
    }
    result$status <- "fit_failed"
    result$input_status <- "success"
    result$error_class <- result$fit$error_class
    result$error_message <- result$fit$error_message
    return(result)
  }
  valid_baseline <- is.data.frame(baseline) &&
    all(c("hazard", "time") %in% names(baseline)) &&
    !anyNA(baseline$hazard) && !anyNA(baseline$time) &&
    all(is.finite(baseline$hazard)) && all(baseline$hazard >= 0) &&
    all(is.finite(baseline$time))
  if (!valid_baseline) {
    result$fit$status <- "failed"
    result$fit$error_class <- "stormbreaker_survival_calibration_bad_baseline"
    result$fit$error_message <- paste0(
      "The centered calibration baseline cumulative hazard is invalid."
    )
    result$status <- "fit_failed"
    result$input_status <- "success"
    result$error_class <- result$fit$error_class
    result$error_message <- result$fit$error_message
    return(result)
  }
  eligible <- baseline$time <= prediction_horizon
  baseline_horizon <- if (any(eligible)) {
    unname(utils::tail(baseline$hazard[eligible], 1L))
  } else {
    0
  }
  result$fit$baseline_cumulative_hazard_horizon <- baseline_horizon

  centered_basis <- sweep(
    basis,
    2L,
    means[c("calibration_ns1", "calibration_ns2")],
    FUN = "-"
  )
  centered_lp <- drop(centered_basis %*% coefficients[
    c("calibration_ns1", "calibration_ns2")
  ])
  if (baseline_horizon == 0) {
    calibrated_survival <- rep(1, n)
    calibrated_event <- rep(0, n)
  } else {
    log_cumulative_hazard <- log(baseline_horizon) + centered_lp
    calibrated_survival <-
      sb_survival_probability_from_log_cumulative_hazard(
        log_cumulative_hazard
      )
    calibrated_event <- -expm1(-exp(log_cumulative_hazard))
  }
  valid_calibrated <- !anyNA(calibrated_survival) &&
    !anyNA(calibrated_event) &&
    all(is.finite(calibrated_survival)) &&
    all(is.finite(calibrated_event)) &&
    all(calibrated_survival >= 0 & calibrated_survival <= 1) &&
    all(calibrated_event >= 0 & calibrated_event <= 1)
  if (!valid_calibrated) {
    result$fit$status <- "failed"
    result$fit$error_class <-
      "stormbreaker_survival_calibration_invalid_probability"
    result$fit$error_message <- paste0(
      "The calibration fit produced invalid horizon probabilities."
    )
    result$status <- "fit_failed"
    result$input_status <- "success"
    result$error_class <- result$fit$error_class
    result$error_message <- result$fit$error_message
    return(result)
  }

  absolute_error <- abs(calibrated_survival - survival_probability)
  result$patient_calibration$calibrated_survival_probability <-
    calibrated_survival
  result$patient_calibration$calibrated_event_probability <- calibrated_event
  result$patient_calibration$absolute_error <- absolute_error
  result$curve <- sb_survival_calibration_curve_table(
    result$patient_calibration
  )
  result$ici <- unname(mean(absolute_error))
  result$e50 <- unname(stats::quantile(
    absolute_error,
    probs = 0.50,
    names = FALSE,
    type = 7
  ))
  result$e90 <- unname(stats::quantile(
    absolute_error,
    probs = 0.90,
    names = FALSE,
    type = 7
  ))
  result$fit$status <- "success"
  result$status <- "success"
  result$input_status <- "success"
  result
}
