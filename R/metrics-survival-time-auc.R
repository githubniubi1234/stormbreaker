sb_time_auc_audit <- function(fold, prediction_horizon, time, event) {
  valid_endpoint <- is.numeric(time) &&
    (is.numeric(event) || is.logical(event)) &&
    length(event) == length(time) &&
    !anyNA(time) && !anyNA(event) &&
    all(is.finite(time)) && all(time >= 0) &&
    (!is.numeric(event) || all(is.finite(event))) &&
    all(event %in% c(0, 1))
  valid_horizon <- is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) && is.finite(prediction_horizon) &&
    prediction_horizon > 0
  if (valid_endpoint && valid_horizon) {
    event <- as.integer(event)
    cases <- event == 1L & time <= prediction_horizon
    controls <- time > prediction_horizon
    excluded_censored <- event == 0L & time <= prediction_horizon
  } else {
    cases <- controls <- excluded_censored <- logical(length(time))
  }
  data.frame(
    fold = as.integer(fold),
    status = "input_failed",
    prediction_horizon = as.numeric(prediction_horizon)[[1L]],
    n = as.integer(length(time)),
    n_cases = as.integer(sum(cases)),
    n_controls = as.integer(sum(controls)),
    n_excluded_censored = as.integer(sum(excluded_censored)),
    case_weight_sum = 0,
    weighted_concordant_pairs = 0,
    weighted_tied_pairs = 0,
    weighted_pair_credit = 0,
    pair_denominator = 0,
    estimate = NA_real_,
    error_class = NA_character_,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}

sb_time_auc_fail <- function(audit, status, error_class, error_message) {
  audit$status <- status
  audit$error_class <- error_class
  audit$error_message <- error_message
  audit
}

sb_time_auc_stratum <- function(time, event, survival_probability,
                                 prediction_ok, censoring_curve,
                                 prediction_horizon, fold = NA_integer_,
                                 log_cumulative_hazard = NULL) {
  audit <- sb_time_auc_audit(fold, prediction_horizon, time, event)
  valid_endpoint <- is.numeric(time) &&
    (is.numeric(event) || is.logical(event)) &&
    length(event) == length(time) &&
    !anyNA(time) && !anyNA(event) &&
    all(is.finite(time)) && all(time >= 0) &&
    (!is.numeric(event) || all(is.finite(event))) &&
    all(event %in% c(0, 1))
  valid_horizon <- is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) && is.finite(prediction_horizon) &&
    prediction_horizon > 0
  if (!valid_endpoint || !valid_horizon) {
    return(sb_time_auc_fail(
      audit,
      "input_failed",
      "stormbreaker_bad_time_auc_input",
      "Time AUC requires a valid survival endpoint and prediction horizon."
    ))
  }

  valid_prediction <- is.numeric(survival_probability) &&
    length(survival_probability) == length(time) &&
    !anyNA(survival_probability) &&
    all(is.finite(survival_probability)) &&
    all(survival_probability >= 0 & survival_probability <= 1) &&
    is.logical(prediction_ok) && length(prediction_ok) == length(time) &&
    !anyNA(prediction_ok) && all(prediction_ok)
  if (!valid_prediction) {
    return(sb_time_auc_fail(
      audit,
      "prediction_failed",
      "stormbreaker_invalid_survival_prediction",
      paste0(
        "Locked-horizon survival probabilities are failed, missing, or ",
        "invalid."
      )
    ))
  }
  if (audit$n_cases == 0L) {
    return(sb_time_auc_fail(
      audit,
      "not_evaluable",
      "stormbreaker_no_time_auc_case",
      "No observed event occurred at or before the prediction horizon."
    ))
  }
  if (audit$n_controls == 0L) {
    return(sb_time_auc_fail(
      audit,
      "not_evaluable",
      "stormbreaker_no_time_auc_control",
      "No patient has observed follow-up beyond the prediction horizon."
    ))
  }

  valid_curve <- is.data.frame(censoring_curve) &&
    nrow(censoring_curve) > 0L &&
    all(c("time", "survival_after") %in% names(censoring_curve))
  if (!valid_curve) {
    return(sb_time_auc_fail(
      audit,
      "censoring_failed",
      "stormbreaker_missing_censoring_curve",
      "The required marginal censoring curve is unavailable."
    ))
  }

  event <- as.integer(event)
  cases <- event == 1L & time <= prediction_horizon
  controls <- time > prediction_horizon
  case_survival <- vapply(
    time[cases],
    function(value) {
      sb_censoring_survival_at(
        censoring_curve,
        value,
        left_limit = TRUE
      )
    },
    numeric(1)
  )
  if (anyNA(case_survival) || any(!is.finite(case_survival)) ||
      any(case_survival <= 0)) {
    return(sb_time_auc_fail(
      audit,
      "censoring_failed",
      "stormbreaker_zero_censoring_survival",
      paste0(
        "The marginal censoring survival is non-positive or non-finite ",
        "before an observed case event."
      )
    ))
  }
  case_weights <- 1 / case_survival
  if (anyNA(case_weights) || any(!is.finite(case_weights))) {
    return(sb_time_auc_fail(
      audit,
      "censoring_failed",
      "stormbreaker_invalid_time_auc_weight",
      "A required time-AUC case weight is missing or non-finite."
    ))
  }

  marker <- if (is.null(log_cumulative_hazard)) {
    1 - survival_probability
  } else sb_validate_stable_survival_scale(log_cumulative_hazard, survival_probability)
  case_marker <- marker[cases]
  control_marker <- marker[controls]
  marker_values <- sort(unique(c(case_marker, control_marker)))
  case_groups <- match(case_marker, marker_values)
  control_groups <- match(control_marker, marker_values)
  control_counts <- tabulate(
    control_groups,
    nbins = length(marker_values)
  )
  case_weight_by_group <- vapply(
    seq_along(marker_values),
    function(group) sum(case_weights[case_groups == group]),
    numeric(1)
  )
  controls_below <- cumsum(control_counts) - control_counts
  weighted_concordant <- unname(sum(
    case_weight_by_group * controls_below
  ))
  weighted_tied <- unname(sum(case_weight_by_group * control_counts))
  weighted_pair_credit <- unname(
    weighted_concordant + 0.5 * weighted_tied
  )
  case_weight_sum <- unname(sum(case_weights))
  pair_denominator <- unname(case_weight_sum * audit$n_controls)
  valid_totals <- all(is.finite(c(
    weighted_concordant,
    weighted_tied,
    weighted_pair_credit,
    case_weight_sum,
    pair_denominator
  ))) && pair_denominator > 0
  if (!valid_totals) {
    return(sb_time_auc_fail(
      audit,
      "censoring_failed",
      "stormbreaker_invalid_time_auc_weight",
      "The weighted time-AUC pair totals are non-finite or non-positive."
    ))
  }

  audit$status <- "success"
  audit$case_weight_sum <- case_weight_sum
  audit$weighted_concordant_pairs <- weighted_concordant
  audit$weighted_tied_pairs <- weighted_tied
  audit$weighted_pair_credit <- weighted_pair_credit
  audit$pair_denominator <- pair_denominator
  audit$estimate <- unname(weighted_pair_credit / pair_denominator)
  audit
}

sb_survival_time_auc_metric <- function(predictions, analysis_inputs,
                                         prediction_horizon,
                                         censoring_curves) {
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  input_folds <- vapply(analysis_inputs, function(input) {
    fold <- input[["fold"]]
    if (is.numeric(fold) && length(fold) == 1L && !is.na(fold) &&
        is.finite(fold)) {
      as.integer(fold)
    } else {
      NA_integer_
    }
  }, integer(1))
  curve_names <- names(censoring_curves)
  if (is.null(curve_names)) {
    curve_names <- rep("", length(censoring_curves))
  }

  fold_results <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      ,
      drop = FALSE
    ]
    prediction_ok <- !is.na(fold_predictions$status) &
      fold_predictions$status == "success" &
      !is.na(fold_predictions$survival_probability_status) &
      fold_predictions$survival_probability_status == "success"
    empty_audit <- function() {
      sb_time_auc_stratum(
        fold_predictions$time,
        fold_predictions$event,
        fold_predictions$survival_probability,
        prediction_ok,
        sb_empty_censoring_curve(),
        prediction_horizon,
        fold_id
      )
    }

    matched_input <- which(input_folds == fold_id)
    if (length(matched_input) != 1L) {
      return(sb_time_auc_fail(
        empty_audit(),
        "analysis_input_failed",
        "stormbreaker_bad_time_auc_analysis_input",
        paste0(
          "Exactly one outer-analysis censoring input is required for this ",
          "fold."
        )
      ))
    }
    matched_curve <- which(curve_names == as.character(fold_id))
    if (length(matched_curve) == 0L) {
      return(sb_time_auc_fail(
        empty_audit(),
        "censoring_failed",
        "stormbreaker_missing_censoring_curve",
        "The matching outer-analysis censoring curve is unavailable."
      ))
    }
    if (length(matched_curve) != 1L) {
      return(sb_time_auc_fail(
        empty_audit(),
        "censoring_failed",
        "stormbreaker_bad_time_auc_censoring_curve",
        "Exactly one matching outer-analysis censoring curve is required."
      ))
    }

    sb_time_auc_stratum(
      fold_predictions$time,
      fold_predictions$event,
      fold_predictions$survival_probability,
      prediction_ok,
      censoring_curves[[matched_curve]],
      prediction_horizon,
      fold_id
    )
  })
  per_fold <- do.call(rbind, fold_results)
  rownames(per_fold) <- NULL
  evaluable <- per_fold$status == "success"
  complete <- length(fold_ids) > 0L && all(evaluable)
  weighted_concordant <- unname(sum(
    per_fold$weighted_concordant_pairs[evaluable]
  ))
  weighted_tied <- unname(sum(per_fold$weighted_tied_pairs[evaluable]))
  weighted_pair_credit <- unname(sum(
    per_fold$weighted_pair_credit[evaluable]
  ))
  pair_denominator <- unname(sum(per_fold$pair_denominator[evaluable]))

  list(
    metric = "time_auc",
    estimand = "cumulative_dynamic_auc",
    scope = "fold_stratified_out_of_fold_at_declared_horizon",
    marker = "one_minus_survival_probability",
    marker_direction = "higher_is_higher_event_risk",
    case_definition = "observed_event_time_less_than_or_equal_to_horizon",
    control_definition = "observed_follow_up_strictly_greater_than_horizon",
    excluded = "censored_at_or_before_horizon",
    censoring_model = "outer_analysis_reverse_kaplan_meier",
    censoring_ties = "events_before_censoring",
    case_weight = "one_over_censoring_survival_at_event_left_limit",
    control_weight = "constant_at_horizon_cancelled_not_evaluated",
    marker_ties = "exact_half_credit",
    aggregation = "summed_within_fold_weighted_pair_credit",
    prediction_horizon = as.numeric(prediction_horizon),
    status = if (complete) "success" else "incomplete",
    estimate = if (complete) {
      unname(weighted_pair_credit / pair_denominator)
    } else {
      NA_real_
    },
    n = as.integer(sum(per_fold$n)),
    n_cases = as.integer(sum(per_fold$n_cases)),
    n_controls = as.integer(sum(per_fold$n_controls)),
    n_excluded_censored = as.integer(sum(
      per_fold$n_excluded_censored
    )),
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    case_weight_sum = unname(sum(per_fold$case_weight_sum[evaluable])),
    weighted_concordant_pairs = weighted_concordant,
    weighted_tied_pairs = weighted_tied,
    weighted_pair_credit = weighted_pair_credit,
    pair_denominator = pair_denominator,
    per_fold = per_fold
  )
}
