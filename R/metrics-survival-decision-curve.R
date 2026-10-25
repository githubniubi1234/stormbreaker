sb_survival_dca_horizon_km <- function(time, event, prediction_horizon) {
  if (length(time) == 0L) {
    return(list(
      status = "success",
      horizon_support = "empty_selection",
      survival_probability = 1,
      event_probability = 0,
      max_follow_up = NA_real_,
      n_followed_to_horizon = 0L
    ))
  }

  curve <- sb_fit_event_km(time, event)
  survival_probability <- sb_event_survival_at(
    curve,
    prediction_horizon
  )
  max_follow_up <- max(time)
  observed_to_horizon <- max_follow_up >= prediction_horizon
  zero_before_horizon <- !observed_to_horizon &&
    isTRUE(survival_probability == 0)

  list(
    status = if (observed_to_horizon || zero_before_horizon) {
      "success"
    } else {
      "not_evaluable"
    },
    horizon_support = if (observed_to_horizon) {
      "observed_to_horizon"
    } else if (zero_before_horizon) {
      "survival_zero_before_horizon"
    } else {
      "insufficient_follow_up"
    },
    survival_probability = survival_probability,
    event_probability = 1 - survival_probability,
    max_follow_up = unname(max_follow_up),
    n_followed_to_horizon = as.integer(sum(time >= prediction_horizon))
  )
}

sb_survival_dca_empty_curve <- function(decision_thresholds, n) {
  threshold_odds <- decision_thresholds / (1 - decision_thresholds)
  data.frame(
    threshold_probability = decision_thresholds,
    threshold_odds = threshold_odds,
    n = rep.int(as.integer(n), length(decision_thresholds)),
    n_selected = rep.int(NA_integer_, length(decision_thresholds)),
    selected_fraction = rep.int(NA_real_, length(decision_thresholds)),
    n_selected_events_by_horizon = rep.int(
      NA_integer_,
      length(decision_thresholds)
    ),
    n_selected_followed_to_horizon = rep.int(
      NA_integer_,
      length(decision_thresholds)
    ),
    max_selected_follow_up = rep.int(
      NA_real_,
      length(decision_thresholds)
    ),
    horizon_support = rep.int(NA_character_, length(decision_thresholds)),
    selected_survival_probability = rep.int(
      NA_real_,
      length(decision_thresholds)
    ),
    selected_event_probability = rep.int(
      NA_real_,
      length(decision_thresholds)
    ),
    estimated_true_positive_rate = rep.int(
      NA_real_,
      length(decision_thresholds)
    ),
    estimated_false_positive_rate = rep.int(
      NA_real_,
      length(decision_thresholds)
    ),
    model_net_benefit = rep.int(NA_real_, length(decision_thresholds)),
    treat_all_net_benefit = rep.int(NA_real_, length(decision_thresholds)),
    treat_none_net_benefit = rep.int(0, length(decision_thresholds)),
    status = rep.int("input_failed", length(decision_thresholds)),
    error_class = rep.int(NA_character_, length(decision_thresholds)),
    error_message = rep.int(NA_character_, length(decision_thresholds)),
    stringsAsFactors = FALSE
  )
}

sb_survival_dca_valid_input <- function(predictions, prediction_horizon,
                                        time_unit, scope) {
  required <- c("time", "event", "survival_probability", "status")
  valid_frame <- is.data.frame(predictions) && nrow(predictions) > 0L &&
    all(vapply(required, function(name) {
      sum(names(predictions) == name) == 1L
    }, logical(1)))
  if (!valid_frame) {
    return(FALSE)
  }

  n <- nrow(predictions)
  time <- predictions$time
  event <- predictions$event
  probability <- predictions$survival_probability
  status <- predictions$status
  valid_probability_status <- TRUE
  if ("survival_probability_status" %in% names(predictions)) {
    valid_probability_status <-
      sum(names(predictions) == "survival_probability_status") == 1L &&
      is.character(predictions$survival_probability_status) &&
      length(predictions$survival_probability_status) == n &&
      !anyNA(predictions$survival_probability_status) &&
      all(predictions$survival_probability_status == "success")
  }

  is.numeric(time) && length(time) == n && !anyNA(time) &&
    all(is.finite(time)) && all(time >= 0) &&
    (is.numeric(event) || is.logical(event)) && length(event) == n &&
    !anyNA(event) && (!is.numeric(event) || all(is.finite(event))) &&
    all(event %in% c(0, 1)) &&
    is.numeric(probability) && length(probability) == n &&
    !anyNA(probability) && all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1) &&
    is.character(status) && length(status) == n && !anyNA(status) &&
    all(status == "success") && valid_probability_status &&
    is.numeric(prediction_horizon) && length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) && is.finite(prediction_horizon) &&
    prediction_horizon > 0 &&
    is.character(time_unit) && length(time_unit) == 1L &&
    !is.na(time_unit) && time_unit %in% c("day", "month", "year") &&
    is.character(scope) && length(scope) == 1L && !is.na(scope) &&
    nzchar(scope)
}

sb_survival_decision_curve_metric <- function(
    predictions,
    decision_thresholds,
    prediction_horizon,
    time_unit,
    scope) {
  decision_thresholds <- sb_validate_decision_thresholds(
    decision_thresholds
  )
  n <- if (is.data.frame(predictions)) nrow(predictions) else 0L
  curve <- sb_survival_dca_empty_curve(decision_thresholds, n)
  result <- list(
    metric = "decision_curve",
    task_type = "survival",
    scope = scope,
    purpose = "descriptive_clinical_utility_evaluation",
    prediction_target = "event_by_declared_horizon",
    threshold_source = "study_declared_survival_task",
    decision_rule = paste0(
      "one_minus_survival_probability_greater_than_or_equal_to_threshold"
    ),
    event_probability_estimator =
      "kaplan_meier_within_threshold_selected_group",
    net_benefit_formula = paste0(
      "selected_fraction * event_probability - selected_fraction * ",
      "survival_probability * threshold / (1 - threshold)"
    ),
    reference_strategies = c("treat_all", "treat_none"),
    censoring_assumption = paste0(
      "independent_noninformative_censoring_within_relevant_risk_groups"
    ),
    censoring_predictor_assumption =
      "censoring_independent_of_model_predictors",
    survival_ties = "events_before_censoring",
    prediction_horizon = as.numeric(prediction_horizon)[[1L]],
    time_unit = as.character(time_unit)[[1L]],
    status = "incomplete",
    n = as.integer(n),
    n_scored = 0L,
    n_thresholds = as.integer(length(decision_thresholds)),
    n_evaluable_thresholds = 0L,
    overall_horizon_support = NA_character_,
    overall_max_follow_up = NA_real_,
    overall_survival_probability = NA_real_,
    overall_event_probability = NA_real_,
    curve = curve,
    error_class = "stormbreaker_invalid_survival_prediction",
    error_message = paste0(
      "Survival endpoints, locked-horizon probabilities, or prediction ",
      "statuses are incomplete or invalid."
    )
  )

  if (!sb_survival_dca_valid_input(
    predictions,
    prediction_horizon,
    time_unit,
    scope
  )) {
    result$curve$status <- "prediction_failed"
    result$curve$error_class <- result$error_class
    result$curve$error_message <- result$error_message
    return(result)
  }

  time <- predictions$time
  event <- as.integer(predictions$event)
  event_probability <- 1 - predictions$survival_probability
  overall <- sb_survival_dca_horizon_km(
    time,
    event,
    prediction_horizon
  )
  result$n_scored <- as.integer(n)
  result$overall_horizon_support <- overall$horizon_support
  result$overall_max_follow_up <- overall$max_follow_up
  result$overall_survival_probability <- overall$survival_probability
  result$overall_event_probability <- overall$event_probability

  rows <- lapply(seq_along(decision_thresholds), function(index) {
    threshold <- decision_thresholds[[index]]
    selected <- event_probability >= threshold
    selected_n <- as.integer(sum(selected))
    selected_fraction <- selected_n / n
    selected_km <- sb_survival_dca_horizon_km(
      time[selected],
      event[selected],
      prediction_horizon
    )
    true_positive_rate <- if (identical(selected_km$status, "success")) {
      selected_fraction * selected_km$event_probability
    } else {
      NA_real_
    }
    false_positive_rate <- if (identical(selected_km$status, "success")) {
      selected_fraction * selected_km$survival_probability
    } else {
      NA_real_
    }
    model_net_benefit <- if (identical(selected_km$status, "success")) {
      true_positive_rate - false_positive_rate *
        curve$threshold_odds[[index]]
    } else {
      NA_real_
    }
    treat_all_net_benefit <- if (identical(overall$status, "success")) {
      overall$event_probability - overall$survival_probability *
        curve$threshold_odds[[index]]
    } else {
      NA_real_
    }
    row_success <- identical(selected_km$status, "success") &&
      identical(overall$status, "success")
    failed_support <- if (!identical(selected_km$status, "success")) {
      "The threshold-selected group has insufficient follow-up at the horizon."
    } else {
      "The complete evaluation sample has insufficient follow-up at the horizon."
    }

    data.frame(
      threshold_probability = threshold,
      threshold_odds = curve$threshold_odds[[index]],
      n = as.integer(n),
      n_selected = selected_n,
      selected_fraction = selected_fraction,
      n_selected_events_by_horizon = as.integer(sum(
        selected & time <= prediction_horizon & event == 1L
      )),
      n_selected_followed_to_horizon =
        selected_km$n_followed_to_horizon,
      max_selected_follow_up = selected_km$max_follow_up,
      horizon_support = selected_km$horizon_support,
      selected_survival_probability =
        selected_km$survival_probability,
      selected_event_probability = selected_km$event_probability,
      estimated_true_positive_rate = true_positive_rate,
      estimated_false_positive_rate = false_positive_rate,
      model_net_benefit = model_net_benefit,
      treat_all_net_benefit = treat_all_net_benefit,
      treat_none_net_benefit = 0,
      status = if (row_success) "success" else "not_evaluable",
      error_class = if (row_success) {
        NA_character_
      } else {
        "stormbreaker_insufficient_survival_dca_follow_up"
      },
      error_message = if (row_success) NA_character_ else failed_support,
      stringsAsFactors = FALSE
    )
  })
  result$curve <- do.call(rbind, rows)
  rownames(result$curve) <- NULL
  result$n_evaluable_thresholds <- as.integer(sum(
    result$curve$status == "success"
  ))
  complete <- result$n_evaluable_thresholds == result$n_thresholds
  result$status <- if (complete) "success" else "incomplete"
  result$error_class <- if (complete) {
    NA_character_
  } else {
    "stormbreaker_insufficient_survival_dca_follow_up"
  }
  result$error_message <- if (complete) {
    NA_character_
  } else {
    "At least one declared survival decision threshold is not evaluable."
  }
  result
}
