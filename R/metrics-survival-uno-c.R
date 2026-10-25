sb_uno_c_audit <- function(fold, prediction_horizon, time, event) {
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
    eligible <- event == 1L & time < prediction_horizon
    excluded_at_horizon <- event == 1L & time == prediction_horizon
  } else {
    eligible <- excluded_at_horizon <- logical(length(time))
  }

  data.frame(
    fold = as.integer(fold),
    status = "input_failed",
    prediction_horizon = if (valid_horizon) {
      as.numeric(prediction_horizon)
    } else {
      NA_real_
    },
    n = as.integer(length(time)),
    n_eligible_events = as.integer(sum(eligible)),
    n_events_at_horizon_excluded = as.integer(sum(excluded_at_horizon)),
    n_comparable_pairs = 0L,
    n_tied_time_comparable_pairs = 0L,
    censoring_survival_horizon = NA_real_,
    weighted_concordant_pairs = 0,
    weighted_tied_risk_pairs = 0,
    weighted_pair_credit = 0,
    weighted_pair_denominator = 0,
    estimate = NA_real_,
    error_class = NA_character_,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}

sb_uno_c_fail <- function(audit, status, error_class, error_message) {
  audit$status <- status
  audit$error_class <- error_class
  audit$error_message <- error_message
  audit
}

sb_uno_c_stratum <- function(time, event, risk, prediction_ok,
                              censoring_curve, prediction_horizon,
                              fold = NA_integer_) {
  audit <- sb_uno_c_audit(fold, prediction_horizon, time, event)
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
    return(sb_uno_c_fail(
      audit,
      "input_failed",
      "stormbreaker_bad_uno_c_input",
      "Uno C requires a valid survival endpoint and prediction horizon."
    ))
  }

  valid_prediction <- is.numeric(risk) && length(risk) == length(time) &&
    !anyNA(risk) && all(is.finite(risk)) &&
    is.logical(prediction_ok) && length(prediction_ok) == length(time) &&
    !anyNA(prediction_ok) && all(prediction_ok)
  if (!valid_prediction) {
    return(sb_uno_c_fail(
      audit,
      "prediction_failed",
      "stormbreaker_invalid_uno_c_risk",
      "Uno C requires one finite successful risk score per patient."
    ))
  }
  if (audit$n_eligible_events == 0L) {
    return(sb_uno_c_fail(
      audit,
      "not_evaluable",
      "stormbreaker_no_uno_c_case",
      "No observed event occurred strictly before the prediction horizon."
    ))
  }

  valid_curve <- is.data.frame(censoring_curve) &&
    nrow(censoring_curve) > 0L &&
    all(c("time", "survival_after") %in% names(censoring_curve))
  if (!valid_curve) {
    return(sb_uno_c_fail(
      audit,
      "censoring_failed",
      "stormbreaker_missing_uno_c_censoring_curve",
      "The required marginal censoring curve is unavailable."
    ))
  }

  horizon_survival <- sb_censoring_survival_at(
    censoring_curve,
    prediction_horizon,
    left_limit = TRUE
  )
  audit$censoring_survival_horizon <- horizon_survival
  if (length(horizon_survival) != 1L || is.na(horizon_survival) ||
      !is.finite(horizon_survival) || horizon_survival <= 0) {
    return(sb_uno_c_fail(
      audit,
      "censoring_failed",
      "stormbreaker_zero_uno_c_horizon_censoring_survival",
      paste0(
        "The marginal censoring survival immediately before the prediction ",
        "horizon is non-positive or non-finite."
      )
    ))
  }

  event <- as.integer(event)
  eligible <- which(event == 1L & time < prediction_horizon)
  event_survival <- vapply(
    time[eligible],
    function(value) {
      sb_censoring_survival_at(
        censoring_curve,
        value,
        left_limit = TRUE
      )
    },
    numeric(1)
  )
  if (anyNA(event_survival) || any(!is.finite(event_survival)) ||
      any(event_survival <= 0)) {
    return(sb_uno_c_fail(
      audit,
      "censoring_failed",
      "stormbreaker_zero_uno_c_event_censoring_survival",
      paste0(
        "The marginal censoring survival immediately before an eligible ",
        "event is non-positive or non-finite."
      )
    ))
  }
  case_weights <- 1 / event_survival^2
  if (anyNA(case_weights) || any(!is.finite(case_weights))) {
    return(sb_uno_c_fail(
      audit,
      "censoring_failed",
      "stormbreaker_invalid_uno_c_weight",
      "A required squared Uno IPCW weight is missing or non-finite."
    ))
  }

  weighted_concordant <- 0
  weighted_tied_risk <- 0
  weighted_denominator <- 0
  n_comparable <- 0L
  n_tied_time <- 0L
  for (position in seq_along(eligible)) {
    index <- eligible[[position]]
    comparators <- which(
      time > time[[index]] |
        (time == time[[index]] & event == 0L)
    )
    if (length(comparators) == 0L) {
      next
    }
    weight <- case_weights[[position]]
    n_comparable <- n_comparable + length(comparators)
    n_tied_time <- n_tied_time + sum(
      time[comparators] == time[[index]]
    )
    weighted_concordant <- weighted_concordant + weight * sum(
      risk[[index]] > risk[comparators]
    )
    weighted_tied_risk <- weighted_tied_risk + weight * sum(
      risk[[index]] == risk[comparators]
    )
    weighted_denominator <- weighted_denominator +
      weight * length(comparators)
  }

  audit$n_comparable_pairs <- as.integer(n_comparable)
  audit$n_tied_time_comparable_pairs <- as.integer(n_tied_time)
  if (n_comparable == 0L) {
    return(sb_uno_c_fail(
      audit,
      "not_evaluable",
      "stormbreaker_no_uno_c_comparable_pair",
      "No patient pair is comparable for Uno C before the horizon."
    ))
  }

  weighted_credit <- weighted_concordant + 0.5 * weighted_tied_risk
  valid_totals <- all(is.finite(c(
    weighted_concordant,
    weighted_tied_risk,
    weighted_credit,
    weighted_denominator
  ))) && weighted_denominator > 0
  if (!valid_totals) {
    return(sb_uno_c_fail(
      audit,
      "censoring_failed",
      "stormbreaker_invalid_uno_c_weight",
      "The weighted Uno C pair totals are non-finite or non-positive."
    ))
  }

  audit$status <- "success"
  audit$weighted_concordant_pairs <- unname(weighted_concordant)
  audit$weighted_tied_risk_pairs <- unname(weighted_tied_risk)
  audit$weighted_pair_credit <- unname(weighted_credit)
  audit$weighted_pair_denominator <- unname(weighted_denominator)
  audit$estimate <- unname(weighted_credit / weighted_denominator)
  audit
}

sb_survival_uno_c_metric <- function(predictions, analysis_inputs,
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
      fold_predictions$status == "success"
    empty_audit <- function() {
      sb_uno_c_audit(
        fold_id,
        prediction_horizon,
        fold_predictions$time,
        fold_predictions$event
      )
    }

    matched_input <- which(input_folds == fold_id)
    if (length(matched_input) != 1L) {
      return(sb_uno_c_fail(
        empty_audit(),
        "analysis_input_failed",
        "stormbreaker_bad_uno_c_analysis_input",
        paste0(
          "Exactly one outer-analysis censoring input is required for this ",
          "fold."
        )
      ))
    }
    matched_curve <- which(curve_names == as.character(fold_id))
    if (length(matched_curve) == 0L) {
      return(sb_uno_c_fail(
        empty_audit(),
        "censoring_failed",
        "stormbreaker_missing_uno_c_censoring_curve",
        "The matching outer-analysis censoring curve is unavailable."
      ))
    }
    if (length(matched_curve) != 1L) {
      return(sb_uno_c_fail(
        empty_audit(),
        "censoring_failed",
        "stormbreaker_bad_uno_c_censoring_curve",
        "Exactly one matching outer-analysis censoring curve is required."
      ))
    }

    result <- tryCatch(
      sb_uno_c_stratum(
        fold_predictions$time,
        fold_predictions$event,
        fold_predictions$risk,
        prediction_ok,
        censoring_curves[[matched_curve]],
        prediction_horizon,
        fold_id
      ),
      error = function(error) error
    )
    if (inherits(result, "condition")) {
      return(sb_uno_c_fail(
        empty_audit(),
        "metric_failed",
        class(result)[[1L]],
        conditionMessage(result)
      ))
    }
    result
  })

  per_fold <- do.call(rbind, fold_results)
  rownames(per_fold) <- NULL
  evaluable <- per_fold$status == "success"
  complete <- length(fold_ids) > 0L && all(evaluable)
  weighted_concordant <- unname(sum(
    per_fold$weighted_concordant_pairs[evaluable]
  ))
  weighted_tied_risk <- unname(sum(
    per_fold$weighted_tied_risk_pairs[evaluable]
  ))
  weighted_credit <- unname(sum(
    per_fold$weighted_pair_credit[evaluable]
  ))
  weighted_denominator <- unname(sum(
    per_fold$weighted_pair_denominator[evaluable]
  ))

  list(
    metric = "uno_c",
    estimand = "truncated_concordance_probability",
    scope = "fold_stratified_out_of_fold_at_declared_horizon",
    marker = "risk",
    marker_direction = "higher_is_higher_event_risk",
    case_definition = "observed_event_time_strictly_less_than_horizon",
    comparator_definition = paste0(
      "later_observed_time_or_censoring_at_same_observed_time"
    ),
    censoring_model = "outer_analysis_reverse_kaplan_meier",
    censoring_assumption = "marginal_random_censoring",
    censoring_ties = "events_before_censoring",
    case_weight = "one_over_squared_censoring_survival_at_event_left_limit",
    time_ties = "event_censor_comparable_event_event_not_comparable",
    marker_ties = "exact_half_credit",
    aggregation = "summed_within_fold_weighted_pair_credit",
    prediction_horizon = as.numeric(prediction_horizon),
    status = if (complete) "success" else "incomplete",
    estimate = if (complete && weighted_denominator > 0) {
      unname(weighted_credit / weighted_denominator)
    } else {
      NA_real_
    },
    n = as.integer(sum(per_fold$n)),
    n_eligible_events = as.integer(sum(per_fold$n_eligible_events)),
    n_events_at_horizon_excluded = as.integer(sum(
      per_fold$n_events_at_horizon_excluded
    )),
    n_comparable_pairs = as.integer(sum(per_fold$n_comparable_pairs)),
    n_tied_time_comparable_pairs = as.integer(sum(
      per_fold$n_tied_time_comparable_pairs
    )),
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    weighted_concordant_pairs = weighted_concordant,
    weighted_tied_risk_pairs = weighted_tied_risk,
    weighted_pair_credit = weighted_credit,
    weighted_pair_denominator = weighted_denominator,
    per_fold = per_fold
  )
}
