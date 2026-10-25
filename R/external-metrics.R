sb_external_binary_auc_metric <- function(scored) {
  stratum <- sb_binary_auc_stratum(
    outcome = scored$outcome,
    probability = scored$probability
  )
  result <- list(
    metric = "roc_auc",
    scope = "pooled_external_cohort",
    aggregation = "all_external_positive_negative_pairs",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    status = "not_evaluable",
    estimate = NA_real_,
    n = stratum$n,
    n_positive = stratum$n_positive,
    n_negative = stratum$n_negative,
    correctly_ordered_pairs = 0,
    tied_pairs = 0,
    comparable_pairs = stratum$comparable_pairs,
    error_class = NA_character_,
    error_message = "Both outcome classes are required."
  )
  if (!identical(stratum$status, "success")) {
    return(result)
  }

  result$status <- "success"
  result$estimate <- stratum$estimate
  result$correctly_ordered_pairs <- stratum$correctly_ordered_pairs
  result$tied_pairs <- stratum$tied_pairs
  result$error_message <- NA_character_
  result
}

sb_external_binary_average_precision_metric <- function(scored) {
  outcome <- as.integer(scored$outcome)
  probability <- scored$probability
  n_positive <- as.integer(sum(outcome == 1L))
  n_negative <- as.integer(sum(outcome == 0L))
  empty_curve <- sb_empty_binary_precision_recall_curve()
  empty_curve$fold <- NULL
  result <- list(
    metric = "average_precision",
    scope = "pooled_external_cohort",
    aggregation = "positive_weighted_precision",
    definition = "sum((recall_i - recall_previous) * precision_i)",
    interpolation = "none",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    status = "not_evaluable",
    estimate = NA_real_,
    n = as.integer(length(outcome)),
    n_positive = n_positive,
    n_negative = n_negative,
    n_thresholds = 0L,
    positive_weighted_sum = 0,
    curve = empty_curve,
    error_class = NA_character_,
    error_message = "Both outcome classes are required."
  )
  if (n_positive == 0L || n_negative == 0L) {
    return(result)
  }
  curve <- sb_binary_precision_recall_curve(outcome, probability, 1L)
  curve$fold <- NULL
  positive_weighted_sum <- sum(
    curve$positive_at_threshold * curve$precision
  )
  result$status <- "success"
  result$estimate <- unname(positive_weighted_sum / n_positive)
  result$n_thresholds <- as.integer(nrow(curve))
  result$positive_weighted_sum <- unname(positive_weighted_sum)
  result$curve <- curve
  result$error_message <- NA_character_
  result
}

sb_external_binary_brier_metric <- function(scored) {
  squared_error <- (as.integer(scored$outcome) - scored$probability)^2
  list(
    metric = "brier_score",
    scope = "patient_level_external_cohort",
    aggregation = "mean_squared_error_over_external_patients",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    metric_direction = "lower_is_better",
    status = "success",
    estimate = unname(mean(squared_error)),
    n = as.integer(nrow(scored)),
    n_scored = as.integer(nrow(scored)),
    squared_error_sum = unname(sum(squared_error))
  )
}

sb_external_binary_log_loss_metric <- function(scored) {
  case_loss <- sb_binary_case_log_loss(
    scored$outcome,
    scored$probability
  )
  list(
    metric = "log_loss",
    scope = "patient_level_external_cohort",
    aggregation = "mean_case_loss_over_external_patients",
    definition = "negative_binary_log_likelihood",
    probability_clipping = "none",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    metric_direction = "lower_is_better",
    status = "success",
    estimate = unname(mean(case_loss)),
    n = as.integer(nrow(scored)),
    n_scored = as.integer(nrow(scored)),
    n_finite_loss = as.integer(sum(is.finite(case_loss))),
    n_infinite_loss = as.integer(sum(is.infinite(case_loss))),
    loss_sum = unname(sum(case_loss))
  )
}

sb_external_binary_mean_calibration_metric <- function(scored) {
  outcome <- as.integer(scored$outcome)
  probability <- scored$probability
  observed_events <- as.integer(sum(outcome == 1L))
  expected_events <- sum(probability)
  result <- list(
    metric = "mean_calibration_ratio",
    calibration_level = "mean",
    estimand = "binary_observed_over_expected_event_ratio",
    scope = "pooled_external_cohort",
    purpose = "descriptive_evaluation_only",
    ratio_direction = "observed_over_expected",
    probability_clipping = "none",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    status = "not_evaluable",
    estimate = NA_real_,
    mean_calibration_ratio = NA_real_,
    target = 1,
    n = as.integer(nrow(scored)),
    observed_events = observed_events,
    expected_events = unname(expected_events),
    observed_event_probability = unname(mean(outcome)),
    mean_predicted_event_probability = unname(mean(probability)),
    error_class = NA_character_,
    error_message = NA_character_
  )
  if (expected_events == 0) {
    result$error_class <-
      "stormbreaker_binary_mean_calibration_zero_expected"
    result$error_message <- paste0(
      "The total expected number of events is zero, so the observed-to-",
      "expected ratio is undefined."
    )
    return(result)
  }

  ratio <- observed_events / expected_events
  if (!is.finite(ratio)) {
    result$error_class <-
      "stormbreaker_binary_mean_calibration_nonfinite_ratio"
    result$error_message <- paste0(
      "The observed-to-expected ratio is not finite on the available ",
      "numeric scale."
    )
    return(result)
  }
  result$status <- "success"
  result$estimate <- unname(ratio)
  result$mean_calibration_ratio <- unname(ratio)
  result
}

sb_external_binary_calibration_metric <- function(scored) {
  outcome <- as.integer(scored$outcome)
  probability <- scored$probability
  n_positive <- as.integer(sum(outcome == 1L))
  n_negative <- as.integer(sum(outcome == 0L))
  n_boundary <- as.integer(sum(probability %in% c(0, 1)))
  intercept_fit <- sb_binary_empty_calibration_fit(
    "calibration_in_the_large",
    "(Intercept)"
  )
  slope_fit <- sb_binary_empty_calibration_fit(
    "logistic_recalibration",
    c("(Intercept)", "logit_probability")
  )
  result <- list(
    metric = "calibration_intercept_slope",
    scope = "pooled_external_cohort",
    purpose = "descriptive_evaluation_only",
    probability_clipping = "none",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    status = "incomplete",
    input_status = "success",
    calibration_intercept = NA_real_,
    calibration_slope = NA_real_,
    intercept_target = 0,
    slope_target = 1,
    n = as.integer(nrow(scored)),
    n_positive = n_positive,
    n_negative = n_negative,
    n_boundary_probability = n_boundary,
    intercept_fit = intercept_fit,
    slope_fit = slope_fit,
    error_class = NA_character_,
    error_message = NA_character_
  )
  if (n_positive == 0L || n_negative == 0L) {
    result$input_status <- "not_evaluable"
    result$error_message <- "Both outcome classes are required for calibration."
    return(result)
  }
  if (n_boundary > 0L) {
    result$input_status <- "not_evaluable"
    result$error_class <-
      "stormbreaker_binary_calibration_boundary_probability"
    result$error_message <- paste0(
      "Exact boundary probabilities cannot enter finite-logit calibration ",
      "without clipping, which is not performed."
    )
    return(result)
  }

  logit_probability <- stats::qlogis(probability)
  intercept_fit <- sb_binary_fit_calibration(
    outcome,
    logit_probability,
    "calibration_in_the_large"
  )
  slope_fit <- sb_binary_fit_calibration(
    outcome,
    logit_probability,
    "logistic_recalibration"
  )
  result$intercept_fit <- intercept_fit
  result$slope_fit <- slope_fit
  if (identical(intercept_fit$status, "success")) {
    result$calibration_intercept <- intercept_fit$coefficients$estimate[[1L]]
  }
  if (identical(slope_fit$status, "success")) {
    result$calibration_slope <- slope_fit$coefficients$estimate[[2L]]
  }
  if (identical(intercept_fit$status, "success") &&
      identical(slope_fit$status, "success")) {
    result$status <- "success"
  } else {
    failed_fit <- if (!identical(intercept_fit$status, "success")) {
      intercept_fit
    } else {
      slope_fit
    }
    result$error_class <- failed_fit$error_class
    result$error_message <- paste(
      stats::na.omit(c(
        intercept_fit$error_message,
        slope_fit$error_message
      )),
      collapse = " "
    )
  }
  result
}

sb_external_harrell_c_metric <- function(scored) {
  result <- list(
    metric = "harrell_c",
    scope = "external_cohort",
    aggregation = "all_external_comparable_pairs",
    time_weight = "n",
    risk_direction = "higher_is_higher_hazard",
    status = "not_evaluable",
    estimate = NA_real_,
    n = as.integer(nrow(scored)),
    n_events = as.integer(sum(scored$event == 1L)),
    concordant = 0,
    discordant = 0,
    tied_risk = 0,
    tied_time = 0,
    tied_both = 0,
    comparable_pairs = 0,
    error_class = NA_character_,
    error_message = "No comparable patient pairs in the external endpoint."
  )
  if (nrow(scored) < 2L || sum(scored$event == 1L) == 0L) {
    return(result)
  }
  stratum <- sb_harrell_c_stratum(
    time = scored$time,
    event = scored$event,
    risk = scored$risk
  )
  result$concordant <- stratum$concordant
  result$discordant <- stratum$discordant
  result$tied_risk <- stratum$tied_risk
  result$tied_time <- stratum$tied_time
  result$tied_both <- stratum$tied_both
  result$comparable_pairs <- stratum$comparable_pairs
  if (identical(stratum$status, "metric_failed")) {
    result$status <- "metric_failed"
    result$error_class <- stratum$error_class
    result$error_message <- stratum$error_message
    return(result)
  }
  if (identical(stratum$status, "not_evaluable")) {
    return(result)
  }
  result$status <- "success"
  result$estimate <- unname(
    (result$concordant + 0.5 * result$tied_risk) /
      result$comparable_pairs
  )
  result$error_message <- NA_character_
  result
}

sb_external_ipcw_brier_metric <- function(scored, endpoint) {
  horizon <- unique(as.numeric(scored$prediction_horizon))
  time <- scored$time
  event <- as.integer(scored$event)
  probability <- scored$survival_probability
  n <- as.integer(nrow(scored))
  curve_result <- tryCatch(
    sb_fit_censoring_km(time, event),
    error = function(error) error
  )
  if (inherits(curve_result, "condition")) {
    return(list(
      metric = "ipcw_brier_score",
      scope = "patient_level_external_cohort_at_locked_horizon",
      censoring_model = "external_reverse_kaplan_meier",
      censoring_assumption = "marginal_noninformative_censoring",
      status = "incomplete",
      estimate = NA_real_,
      prediction_horizon = horizon,
      time_unit = endpoint$time_unit,
      n = n,
      n_scored = 0L,
      error_class = class(curve_result)[[1L]],
      error_message = conditionMessage(curve_result),
      censoring_curve = sb_empty_censoring_curve()
    ))
  }
  curve <- curve_result
  event_by_horizon <- time <= horizon & event == 1L
  surviving_at_horizon <- time > horizon |
    (time == horizon & event == 0L)
  n_scored <- as.integer(sum(event_by_horizon | surviving_at_horizon))
  horizon_survival <- sb_censoring_survival_at(
    curve,
    horizon,
    left_limit = TRUE
  )
  result <- list(
    metric = "ipcw_brier_score",
    scope = "patient_level_external_cohort_at_locked_horizon",
    censoring_model = "external_reverse_kaplan_meier",
    censoring_assumption = "marginal_noninformative_censoring",
    censoring_ties = "events_before_censoring",
    censoring_evaluation = "left_continuous_before_censoring_at_time",
    aggregation = "summed_weighted_loss_over_all_external_patients",
    denominator = "all_external_patients_n",
    prediction_horizon = horizon,
    time_unit = endpoint$time_unit,
    status = "incomplete",
    estimate = NA_real_,
    n = n,
    n_scored = n_scored,
    n_events_by_horizon = as.integer(sum(event_by_horizon)),
    n_surviving_at_horizon = as.integer(sum(surviving_at_horizon)),
    n_censored_before_horizon = as.integer(sum(
      time < horizon & event == 0L
    )),
    n_censored_at_horizon = as.integer(sum(
      time == horizon & event == 0L
    )),
    censoring_survival_horizon = unname(horizon_survival),
    weight_sum = 0,
    weighted_loss_sum = 0,
    censoring_curve = curve,
    error_class = NA_character_,
    error_message = NA_character_
  )
  if (n_scored == 0L) {
    result$error_class <- "stormbreaker_no_evaluable_horizon_status"
    result$error_message <- paste0(
      "No external patient has observed event status at the locked horizon."
    )
    return(result)
  }

  event_survival <- vapply(
    time[event_by_horizon],
    function(value) sb_censoring_survival_at(
      curve,
      value,
      left_limit = TRUE
    ),
    numeric(1)
  )
  valid_event_survival <- all(is.finite(event_survival)) &&
    all(event_survival > 0)
  valid_horizon_survival <- !any(surviving_at_horizon) ||
    (is.finite(horizon_survival) && horizon_survival > 0)
  if (!valid_event_survival || !valid_horizon_survival) {
    result$error_class <- "stormbreaker_zero_external_censoring_survival"
    result$error_message <- paste0(
      "A reverse-Kaplan-Meier denominator required for external IPCW is ",
      "non-positive or non-finite."
    )
    return(result)
  }

  weights <- numeric(n)
  weights[event_by_horizon] <- 1 / event_survival
  if (any(surviving_at_horizon)) {
    weights[surviving_at_horizon] <- 1 / horizon_survival
  }
  losses <- numeric(n)
  losses[event_by_horizon] <- probability[event_by_horizon]^2
  losses[surviving_at_horizon] <-
    (1 - probability[surviving_at_horizon])^2
  weighted_loss <- weights * losses
  if (anyNA(weighted_loss) || any(!is.finite(weighted_loss))) {
    result$error_class <- "stormbreaker_invalid_external_ipcw_loss"
    result$error_message <- "External IPCW losses are missing or non-finite."
    return(result)
  }

  result$status <- "success"
  result$weight_sum <- unname(sum(weights))
  result$weighted_loss_sum <- unname(sum(weighted_loss))
  result$estimate <- unname(result$weighted_loss_sum / n)
  result
}

sb_external_time_auc_metric <- function(scored, endpoint, censoring_curve) {
  audit <- sb_time_auc_stratum(
    time = scored$time,
    event = scored$event,
    survival_probability = scored$survival_probability,
    log_cumulative_hazard = scored$log_cumulative_hazard,
    prediction_ok = !is.na(scored$status) & scored$status == "success",
    censoring_curve = censoring_curve,
    prediction_horizon = unique(as.numeric(scored$prediction_horizon)),
    fold = NA_integer_
  )
  list(
    metric = "time_auc",
    estimand = "cumulative_dynamic_auc",
    scope = "external_cohort_at_locked_horizon",
    marker = if (is.null(scored$log_cumulative_hazard)) {
      "one_minus_survival_probability"
    } else "log_cumulative_hazard",
    marker_direction = "higher_is_higher_event_risk",
    case_definition = "observed_event_time_less_than_or_equal_to_horizon",
    control_definition = "observed_follow_up_strictly_greater_than_horizon",
    excluded = "censored_at_or_before_horizon",
    censoring_model = "external_reverse_kaplan_meier",
    censoring_assumption = "marginal_noninformative_censoring",
    censoring_ties = "events_before_censoring",
    case_weight = "one_over_censoring_survival_at_event_left_limit",
    control_weight = "constant_at_horizon_cancelled_not_evaluated",
    marker_ties = "exact_half_credit",
    aggregation = "all_external_case_control_pairs",
    prediction_horizon = audit$prediction_horizon,
    time_unit = endpoint$time_unit,
    status = audit$status,
    estimate = audit$estimate,
    n = audit$n,
    n_cases = audit$n_cases,
    n_controls = audit$n_controls,
    n_excluded_censored = audit$n_excluded_censored,
    case_weight_sum = audit$case_weight_sum,
    weighted_concordant_pairs = audit$weighted_concordant_pairs,
    weighted_tied_pairs = audit$weighted_tied_pairs,
    weighted_pair_credit = audit$weighted_pair_credit,
    pair_denominator = audit$pair_denominator,
    error_class = audit$error_class,
    error_message = audit$error_message
  )
}

sb_external_uno_c_metric <- function(scored, endpoint, censoring_curve) {
  horizon <- unique(as.numeric(scored$prediction_horizon))
  audit <- tryCatch(
    sb_uno_c_stratum(
      time = scored$time,
      event = scored$event,
      risk = scored$risk,
      prediction_ok = !is.na(scored$status) & scored$status == "success",
      censoring_curve = censoring_curve,
      prediction_horizon = horizon,
      fold = NA_integer_
    ),
    error = function(error) error
  )
  if (inherits(audit, "condition")) {
    audit <- sb_uno_c_fail(
      sb_uno_c_audit(NA_integer_, horizon, scored$time, scored$event),
      "metric_failed",
      class(audit)[[1L]],
      conditionMessage(audit)
    )
  }

  list(
    metric = "uno_c",
    estimand = "truncated_concordance_probability",
    scope = "external_cohort_at_locked_horizon",
    marker = "risk",
    marker_direction = "higher_is_higher_event_risk",
    case_definition = "observed_event_time_strictly_less_than_horizon",
    comparator_definition = paste0(
      "later_observed_time_or_censoring_at_same_observed_time"
    ),
    censoring_model = "external_reverse_kaplan_meier",
    censoring_assumption = "marginal_random_censoring",
    censoring_ties = "events_before_censoring",
    case_weight = "one_over_squared_censoring_survival_at_event_left_limit",
    time_ties = "event_censor_comparable_event_event_not_comparable",
    marker_ties = "exact_half_credit",
    aggregation = "all_external_comparable_pairs",
    prediction_horizon = audit$prediction_horizon,
    time_unit = endpoint$time_unit,
    status = audit$status,
    estimate = audit$estimate,
    n = audit$n,
    n_eligible_events = audit$n_eligible_events,
    n_events_at_horizon_excluded = audit$n_events_at_horizon_excluded,
    n_comparable_pairs = audit$n_comparable_pairs,
    n_tied_time_comparable_pairs = audit$n_tied_time_comparable_pairs,
    censoring_survival_horizon = audit$censoring_survival_horizon,
    weighted_concordant_pairs = audit$weighted_concordant_pairs,
    weighted_tied_risk_pairs = audit$weighted_tied_risk_pairs,
    weighted_pair_credit = audit$weighted_pair_credit,
    weighted_pair_denominator = audit$weighted_pair_denominator,
    error_class = audit$error_class,
    error_message = audit$error_message
  )
}
