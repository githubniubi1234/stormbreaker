sb_fit_censoring_km <- function(time, event) {
  time <- sb_validate_endpoint_vector(time, "time", length(time))
  event <- sb_validate_endpoint_vector(event, "event", length(time))
  if (any(time < 0) || !all(event %in% c(0, 1))) {
    sb_abort(
      "The censoring estimate requires non-negative time and 0/1 event.",
      "stormbreaker_bad_censoring_input"
    )
  }

  times <- sort(unique(time))
  n_risk <- vapply(times, function(value) sum(time >= value), integer(1))
  n_event <- vapply(
    times,
    function(value) sum(time == value & event == 1L),
    integer(1)
  )
  n_censoring <- vapply(
    times,
    function(value) sum(time == value & event == 0L),
    integer(1)
  )
  n_risk_for_censoring <- n_risk - n_event
  survival_before <- numeric(length(times))
  survival_after <- numeric(length(times))
  current <- 1
  for (index in seq_along(times)) {
    survival_before[index] <- current
    if (n_censoring[index] > 0L) {
      current <- current * (
        1 - n_censoring[index] / n_risk_for_censoring[index]
      )
    }
    survival_after[index] <- current
  }

  data.frame(
    time = times,
    n_risk = n_risk,
    n_event = n_event,
    n_censoring = n_censoring,
    n_risk_for_censoring = n_risk_for_censoring,
    survival_before = survival_before,
    survival_after = survival_after,
    stringsAsFactors = FALSE
  )
}

sb_censoring_survival_at <- function(curve, time, left_limit = FALSE) {
  eligible <- if (left_limit) curve$time < time else curve$time <= time
  if (!any(eligible)) {
    return(1)
  }
  unname(utils::tail(curve$survival_after[eligible], 1L))
}

sb_empty_censoring_curve <- function() {
  data.frame(
    time = numeric(),
    n_risk = integer(),
    n_event = integer(),
    n_censoring = integer(),
    n_risk_for_censoring = integer(),
    survival_before = numeric(),
    survival_after = numeric(),
    stringsAsFactors = FALSE
  )
}

sb_fit_event_km <- function(time, event) {
  time <- sb_validate_endpoint_vector(time, "time", length(time))
  event <- sb_validate_endpoint_vector(event, "event", length(time))
  if (any(time < 0) || !all(event %in% c(0, 1))) {
    sb_abort(
      "The survival reference requires non-negative time and 0/1 event.",
      "stormbreaker_bad_survival_reference_input"
    )
  }

  times <- sort(unique(time))
  n_risk <- vapply(times, function(value) sum(time >= value), integer(1))
  n_event <- vapply(
    times,
    function(value) sum(time == value & event == 1L),
    integer(1)
  )
  n_censoring <- vapply(
    times,
    function(value) sum(time == value & event == 0L),
    integer(1)
  )
  survival_before <- numeric(length(times))
  survival_after <- numeric(length(times))
  current <- 1
  for (index in seq_along(times)) {
    survival_before[index] <- current
    if (n_event[index] > 0L) {
      current <- current * (1 - n_event[index] / n_risk[index])
    }
    survival_after[index] <- current
  }

  data.frame(
    time = times,
    n_risk = n_risk,
    n_event = n_event,
    n_censoring = n_censoring,
    survival_before = survival_before,
    survival_after = survival_after,
    stringsAsFactors = FALSE
  )
}

sb_event_survival_at <- function(curve, time) {
  eligible <- curve$time <= time
  if (!any(eligible)) {
    return(1)
  }
  unname(utils::tail(curve$survival_after[eligible], 1L))
}

sb_empty_event_curve <- function() {
  data.frame(
    time = numeric(),
    n_risk = integer(),
    n_event = integer(),
    n_censoring = integer(),
    survival_before = numeric(),
    survival_after = numeric(),
    stringsAsFactors = FALSE
  )
}

sb_survival_ipcw_brier_metric <- function(predictions, analysis_inputs,
                                          prediction_horizon) {
  prediction_horizon <- as.numeric(prediction_horizon)
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  fold_results <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      ,
      drop = FALSE
    ]
    audit <- data.frame(
      fold = fold_id,
      status = "analysis_input_failed",
      prediction_horizon = prediction_horizon,
      n = as.integer(nrow(fold_predictions)),
      n_scored = 0L,
      n_events_by_horizon = as.integer(sum(
        fold_predictions$time <= prediction_horizon &
          fold_predictions$event == 1L
      )),
      n_surviving_past_horizon = as.integer(sum(
        fold_predictions$time > prediction_horizon
      )),
      n_censored_at_horizon = as.integer(sum(
        fold_predictions$time == prediction_horizon &
          fold_predictions$event == 0L
      )),
      n_censored_by_horizon = as.integer(sum(
        fold_predictions$time <= prediction_horizon &
          fold_predictions$event == 0L
      )),
      censoring_survival_horizon = NA_real_,
      weight_sum = 0,
      weighted_loss_sum = 0,
      estimate = NA_real_,
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    reference_audit <- data.frame(
      fold = fold_id,
      status = "analysis_input_failed",
      prediction_horizon = prediction_horizon,
      reference_survival_horizon = NA_real_,
      n = as.integer(nrow(fold_predictions)),
      n_scored = 0L,
      censoring_survival_horizon = NA_real_,
      weight_sum = 0,
      weighted_loss_sum = 0,
      estimate = NA_real_,
      source = "outer_analysis_event_only",
      ties = "events_before_censoring",
      horizon_evaluation = "right_continuous_after_events",
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    curve <- sb_empty_censoring_curve()
    reference_curve <- sb_empty_event_curve()
    input_folds <- vapply(
      analysis_inputs,
      function(input) as.integer(input$fold),
      integer(1)
    )
    matched <- which(input_folds == fold_id)
    if (length(matched) != 1L) {
      audit$error_class <- "stormbreaker_bad_ipcw_analysis_input"
      audit$error_message <- paste0(
        "Exactly one outer-analysis censoring input is required for this fold."
      )
      reference_audit$error_class <- audit$error_class
      reference_audit$error_message <- audit$error_message
      return(list(
        audit = audit,
        curve = curve,
        reference_audit = reference_audit,
        reference_curve = reference_curve
      ))
    }

    analysis <- analysis_inputs[[matched]]
    curve_result <- tryCatch(
      sb_fit_censoring_km(analysis$time, analysis$event),
      error = function(error) error
    )
    if (inherits(curve_result, "error")) {
      audit$error_class <- class(curve_result)[1]
      audit$error_message <- conditionMessage(curve_result)
      reference_audit$error_class <- audit$error_class
      reference_audit$error_message <- audit$error_message
      return(list(
        audit = audit,
        curve = curve,
        reference_audit = reference_audit,
        reference_curve = reference_curve
      ))
    }
    curve <- curve_result
    reference_curve_result <- tryCatch(
      sb_fit_event_km(analysis$time, analysis$event),
      error = function(error) error
    )
    if (inherits(reference_curve_result, "error")) {
      reference_audit$status <- "reference_failed"
      reference_audit$error_class <- class(reference_curve_result)[1]
      reference_audit$error_message <- conditionMessage(
        reference_curve_result
      )
    } else {
      reference_curve <- reference_curve_result
    }
    if (prediction_horizon > max(analysis$time)) {
      audit$status <- "horizon_unsupported"
      audit$error_class <- "stormbreaker_prediction_horizon_unsupported"
      audit$error_message <- paste0(
        "The prediction horizon exceeds this outer analysis fold's ",
        "observed follow-up."
      )
      reference_audit$status <- audit$status
      reference_audit$error_class <- audit$error_class
      reference_audit$error_message <- audit$error_message
      return(list(
        audit = audit,
        curve = curve,
        reference_audit = reference_audit,
        reference_curve = reference_curve
      ))
    }

    censoring_survival_horizon <- sb_censoring_survival_at(
      curve,
      prediction_horizon,
      left_limit = TRUE
    )
    audit$censoring_survival_horizon <- censoring_survival_horizon
    reference_audit$censoring_survival_horizon <-
      censoring_survival_horizon
    if (!is.finite(censoring_survival_horizon) ||
        censoring_survival_horizon <= 0) {
      audit$status <- "censoring_failed"
      audit$error_class <- "stormbreaker_zero_censoring_survival"
      audit$error_message <- paste0(
        "The outer-analysis censoring survival is zero at the prediction ",
        "horizon."
      )
      reference_audit$status <- audit$status
      reference_audit$error_class <- audit$error_class
      reference_audit$error_message <- audit$error_message
      return(list(
        audit = audit,
        curve = curve,
        reference_audit = reference_audit,
        reference_curve = reference_curve
      ))
    }

    event_by_horizon <- fold_predictions$time <= prediction_horizon &
      fold_predictions$event == 1L
    surviving_at_horizon <- fold_predictions$time > prediction_horizon |
      (
        fold_predictions$time == prediction_horizon &
          fold_predictions$event == 0L
      )
    event_survival <- vapply(
      fold_predictions$time[event_by_horizon],
      function(value) sb_censoring_survival_at(
        curve,
        value,
        left_limit = TRUE
      ),
      numeric(1)
    )
    if (any(!is.finite(event_survival)) || any(event_survival <= 0)) {
      audit$status <- "censoring_failed"
      audit$error_class <- "stormbreaker_zero_censoring_survival"
      audit$error_message <- paste0(
        "The outer-analysis censoring survival is zero before an observed ",
        "assessment event."
      )
      reference_audit$status <- audit$status
      reference_audit$error_class <- audit$error_class
      reference_audit$error_message <- audit$error_message
      return(list(
        audit = audit,
        curve = curve,
        reference_audit = reference_audit,
        reference_curve = reference_curve
      ))
    }

    weights <- numeric(nrow(fold_predictions))
    weights[event_by_horizon] <- 1 / event_survival
    weights[surviving_at_horizon] <- 1 / censoring_survival_horizon
    if (anyNA(weights) || any(!is.finite(weights))) {
      audit$status <- "censoring_failed"
      audit$error_class <- "stormbreaker_invalid_ipcw_weight"
      audit$error_message <- "The IPCW weights are missing or non-finite."
      reference_audit$status <- audit$status
      reference_audit$error_class <- audit$error_class
      reference_audit$error_message <- audit$error_message
      return(list(
        audit = audit,
        curve = curve,
        reference_audit = reference_audit,
        reference_curve = reference_curve
      ))
    }

    n_scored <- as.integer(sum(event_by_horizon | surviving_at_horizon))
    weight_sum <- unname(sum(weights))
    if (!inherits(reference_curve_result, "error")) {
      reference_survival <- sb_event_survival_at(
        reference_curve,
        prediction_horizon
      )
      valid_reference <- length(reference_survival) == 1L &&
        !is.na(reference_survival) &&
        is.finite(reference_survival) &&
        reference_survival >= 0 && reference_survival <= 1
      if (valid_reference) {
        reference_losses <- numeric(nrow(fold_predictions))
        reference_losses[event_by_horizon] <- reference_survival^2
        reference_losses[surviving_at_horizon] <- (1 - reference_survival)^2
        reference_weighted_loss <- weights * reference_losses
        if (anyNA(reference_weighted_loss) ||
            any(!is.finite(reference_weighted_loss))) {
          reference_audit$status <- "reference_failed"
          reference_audit$error_class <-
            "stormbreaker_invalid_reference_weighted_loss"
          reference_audit$error_message <- paste0(
            "The reference IPCW weighted losses are missing or non-finite."
          )
        } else {
          reference_audit$status <- "success"
          reference_audit$reference_survival_horizon <- reference_survival
          reference_audit$n_scored <- n_scored
          reference_audit$weight_sum <- weight_sum
          reference_audit$weighted_loss_sum <- unname(sum(
            reference_weighted_loss
          ))
          reference_audit$estimate <- unname(
            reference_audit$weighted_loss_sum / reference_audit$n
          )
        }
      } else {
        reference_audit$status <- "reference_failed"
        reference_audit$error_class <-
          "stormbreaker_invalid_reference_survival"
        reference_audit$error_message <- paste0(
          "The outer-analysis reference survival is invalid."
        )
      }
    }

    prediction_failed <-
      anyNA(fold_predictions$status) ||
      any(fold_predictions$status != "success") ||
      anyNA(fold_predictions$survival_probability_status) ||
      any(fold_predictions$survival_probability_status != "success") ||
      anyNA(fold_predictions$survival_probability) ||
      any(!is.finite(fold_predictions$survival_probability)) ||
      any(
        fold_predictions$survival_probability < 0 |
          fold_predictions$survival_probability > 1
      )
    if (prediction_failed) {
      audit$status <- "prediction_failed"
      audit$error_class <- "stormbreaker_invalid_survival_prediction"
      audit$error_message <- paste0(
        "OOF survival probabilities for this fold are failed or invalid."
      )
      return(list(
        audit = audit,
        curve = curve,
        reference_audit = reference_audit,
        reference_curve = reference_curve
      ))
    }

    losses <- numeric(nrow(fold_predictions))
    losses[event_by_horizon] <-
      fold_predictions$survival_probability[event_by_horizon]^2
    losses[surviving_at_horizon] <- (
      1 - fold_predictions$survival_probability[surviving_at_horizon]
    )^2
    weighted_loss <- weights * losses
    if (anyNA(weighted_loss) || any(!is.finite(weighted_loss))) {
      audit$status <- "prediction_failed"
      audit$error_class <- "stormbreaker_invalid_ipcw_weighted_loss"
      audit$error_message <- paste0(
        "The model IPCW weighted losses are missing or non-finite."
      )
      return(list(
        audit = audit,
        curve = curve,
        reference_audit = reference_audit,
        reference_curve = reference_curve
      ))
    }

    audit$status <- "success"
    audit$n_scored <- n_scored
    audit$weight_sum <- weight_sum
    audit$weighted_loss_sum <- unname(sum(weighted_loss))
    audit$estimate <- unname(audit$weighted_loss_sum / audit$n)
    list(
      audit = audit,
      curve = curve,
      reference_audit = reference_audit,
      reference_curve = reference_curve
    )
  })

  per_fold <- do.call(rbind, lapply(fold_results, `[[`, "audit"))
  rownames(per_fold) <- NULL
  censoring_curves <- lapply(fold_results, `[[`, "curve")
  names(censoring_curves) <- as.character(fold_ids)
  reference_per_fold <- do.call(
    rbind,
    lapply(fold_results, `[[`, "reference_audit")
  )
  rownames(reference_per_fold) <- NULL
  reference_curves <- lapply(fold_results, `[[`, "reference_curve")
  names(reference_curves) <- as.character(fold_ids)
  evaluable <- per_fold$status == "success"
  complete <- length(fold_ids) > 0L && all(evaluable)
  n <- as.integer(sum(per_fold$n))
  weighted_loss_sum <- unname(sum(per_fold$weighted_loss_sum[evaluable]))
  reference_evaluable <- reference_per_fold$status == "success"
  reference_complete <- length(fold_ids) > 0L && all(reference_evaluable)
  reference_weighted_loss_sum <- unname(sum(
    reference_per_fold$weighted_loss_sum[reference_evaluable]
  ))
  reference <- list(
    metric = "ipcw_brier_reference",
    type = "outer_analysis_kaplan_meier",
    source = "outer_analysis_event_only",
    scope = "patient_level_out_of_fold_at_declared_horizon",
    censoring_model = "outer_analysis_reverse_kaplan_meier",
    aggregation = "summed_weighted_loss_over_all_assessment_patients",
    prediction_horizon = prediction_horizon,
    ties = "events_before_censoring",
    horizon_evaluation = "right_continuous_after_events",
    status = if (reference_complete) "success" else "incomplete",
    estimate = if (reference_complete) {
      reference_weighted_loss_sum / n
    } else {
      NA_real_
    },
    n = n,
    n_scored = as.integer(sum(
      reference_per_fold$n_scored[reference_evaluable]
    )),
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(reference_evaluable)),
    weight_sum = unname(sum(
      reference_per_fold$weight_sum[reference_evaluable]
    )),
    weighted_loss_sum = reference_weighted_loss_sum,
    per_fold = reference_per_fold,
    survival_curves = reference_curves
  )
  skill_complete <- complete && reference_complete &&
    is.finite(reference_weighted_loss_sum) &&
    reference_weighted_loss_sum > 0
  skill_score <- list(
    metric = "ipcw_brier_skill_score",
    reference = reference$type,
    formula = paste0(
      "1 - model_weighted_loss_sum / reference_weighted_loss_sum"
    ),
    status = if (skill_complete) "success" else "incomplete",
    estimate = if (skill_complete) {
      unname(1 - weighted_loss_sum / reference_weighted_loss_sum)
    } else {
      NA_real_
    }
  )

  list(
    metric = "ipcw_brier_score",
    scope = "patient_level_out_of_fold_at_declared_horizon",
    censoring_model = "outer_analysis_reverse_kaplan_meier",
    censoring_ties = "events_before_censoring",
    censoring_evaluation = "left_continuous_before_censoring_at_time",
    aggregation = "summed_weighted_loss_over_all_assessment_patients",
    prediction_horizon = prediction_horizon,
    status = if (complete) "success" else "incomplete",
    estimate = if (complete) weighted_loss_sum / n else NA_real_,
    n = n,
    n_scored = as.integer(sum(per_fold$n_scored[evaluable])),
    n_events_by_horizon = as.integer(sum(per_fold$n_events_by_horizon)),
    n_surviving_past_horizon = as.integer(sum(
      per_fold$n_surviving_past_horizon
    )),
    n_censored_at_horizon = as.integer(sum(
      per_fold$n_censored_at_horizon
    )),
    n_censored_by_horizon = as.integer(sum(
      per_fold$n_censored_by_horizon
    )),
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    weight_sum = unname(sum(per_fold$weight_sum[evaluable])),
    weighted_loss_sum = weighted_loss_sum,
    per_fold = per_fold,
    censoring_curves = censoring_curves,
    reference = reference,
    skill_score = skill_score
  )
}
