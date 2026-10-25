sb_validate_prediction_horizon <- function(prediction_horizon, analysis_time) {
  valid <- is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) &&
    is.finite(prediction_horizon) &&
    prediction_horizon > 0
  if (!valid) {
    sb_abort(
      "`prediction_horizon` must be one finite positive number.",
      "stormbreaker_bad_prediction_horizon"
    )
  }
  if (prediction_horizon > max(analysis_time)) {
    sb_abort(
      paste0(
        "The prediction horizon exceeds this outer analysis fold's ",
        "observed follow-up."
      ),
      "stormbreaker_prediction_horizon_unsupported"
    )
  }
  as.numeric(prediction_horizon)
}

sb_survival_probability_from_hazard <- function(risk,
                                                 cumulative_baseline_hazard) {
  valid_hazard <- is.numeric(cumulative_baseline_hazard) &&
    length(cumulative_baseline_hazard) == 1L &&
    !is.na(cumulative_baseline_hazard) &&
    is.finite(cumulative_baseline_hazard) &&
    cumulative_baseline_hazard >= 0
  if (!valid_hazard) {
    sb_abort(
      "The fitted baseline cumulative hazard is invalid.",
      "stormbreaker_invalid_baseline_hazard"
    )
  }
  if (length(risk) < 1L || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "The fitted linear predictors are incomplete or non-finite.",
      "stormbreaker_invalid_survival_prediction"
    )
  }

  probability <- if (cumulative_baseline_hazard == 0) {
    rep(1, length(risk))
  } else {
    log_patient_hazard <- log(cumulative_baseline_hazard) + risk
    patient_hazard <- exp(log_patient_hazard)
    exp(-patient_hazard)
  }
  if (length(probability) != length(risk) || anyNA(probability) ||
      any(!is.finite(probability)) ||
      any(probability < 0 | probability > 1)) {
    sb_abort(
      "The fitted survival probabilities are incomplete or invalid.",
      "stormbreaker_invalid_survival_prediction"
    )
  }
  unname(probability)
}

sb_breslow_baseline_hazard <- function(time, event, risk,
                                       prediction_horizon) {
  event_times <- sort(unique(time[event == 1L & time <= prediction_horizon]))
  if (length(event_times) == 0L) {
    return(0)
  }

  increments <- vapply(event_times, function(event_time) {
    risk_set <- risk[time >= event_time]
    maximum <- max(risk_set)
    log_denominator <- maximum + log(sum(exp(risk_set - maximum)))
    n_events <- sum(event == 1L & time == event_time)
    exp(log(n_events) - log_denominator)
  }, numeric(1))
  if (anyNA(increments) || any(!is.finite(increments)) ||
      any(increments <= 0)) {
    sb_abort(
      "The fixed-coefficient Breslow baseline hazard is invalid.",
      "stormbreaker_invalid_baseline_hazard"
    )
  }
  cumulative_hazard <- unname(sum(increments))
  if (!is.finite(cumulative_hazard) || cumulative_hazard < 0) {
    sb_abort(
      "The fixed-coefficient Breslow baseline hazard is invalid.",
      "stormbreaker_invalid_baseline_hazard"
    )
  }
  cumulative_hazard
}

sb_survival_prediction_result <- function(prediction_horizon,
                                           cumulative_baseline_hazard,
                                           source, ties, risk) {
  probability <- sb_survival_probability_from_hazard(
    risk,
    cumulative_baseline_hazard
  )
  list(
    prediction_horizon = as.numeric(prediction_horizon),
    cumulative_baseline_hazard = cumulative_baseline_hazard,
    baseline_survival = unname(exp(-cumulative_baseline_hazard)),
    source = source,
    ties = ties,
    probability = probability
  )
}
