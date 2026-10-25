sb_external_bootstrap_survival_decision_curve_comparisons <- function(metric) {
  curve <- data.frame(
    threshold_index = seq_len(nrow(metric$curve)),
    metric$curve,
    stringsAsFactors = FALSE
  )
  curve$overall_horizon_support <- metric$overall_horizon_support
  curve$overall_max_follow_up <- metric$overall_max_follow_up
  curve$overall_survival_probability <-
    metric$overall_survival_probability
  curve$overall_event_probability <- metric$overall_event_probability

  all_selected <- !is.na(curve$n_selected) &
    !is.na(curve$n) & curve$n_selected == curve$n
  none_selected <- !is.na(curve$n_selected) & curve$n_selected == 0L
  model_finite <- is.finite(curve$model_net_benefit)
  all_finite <- is.finite(curve$treat_all_net_benefit)

  model_minus_all <- rep.int(NA_real_, nrow(curve))
  model_minus_all[all_selected] <- 0
  ordinary_all <- !all_selected & model_finite & all_finite
  model_minus_all[ordinary_all] <-
    (curve$estimated_true_positive_rate[ordinary_all] -
       curve$overall_event_probability[ordinary_all]) -
    (curve$estimated_false_positive_rate[ordinary_all] -
       curve$overall_survival_probability[ordinary_all]) *
    curve$threshold_odds[ordinary_all]

  model_minus_none <- rep.int(NA_real_, nrow(curve))
  model_minus_none[none_selected] <- 0
  ordinary_none <- !none_selected & model_finite
  model_minus_none[ordinary_none] <- curve$model_net_benefit[ordinary_none]

  all_success <- is.finite(model_minus_all)
  none_success <- is.finite(model_minus_none)
  curve$model_minus_treat_all <- model_minus_all
  curve$model_minus_treat_none <- model_minus_none
  curve$model_minus_treat_all_status <- ifelse(
    all_success,
    "success",
    "not_evaluable"
  )
  curve$model_minus_treat_none_status <- ifelse(
    none_success,
    "success",
    "not_evaluable"
  )
  curve$model_minus_treat_all_error_class <- ifelse(
    all_success,
    NA_character_,
    curve$error_class
  )
  curve$model_minus_treat_all_error_message <- ifelse(
    all_success,
    NA_character_,
    curve$error_message
  )
  curve$model_minus_treat_none_error_class <- ifelse(
    none_success,
    NA_character_,
    curve$error_class
  )
  curve$model_minus_treat_none_error_message <- ifelse(
    none_success,
    NA_character_,
    curve$error_message
  )
  curve
}

sb_external_bootstrap_survival_decision_curve_replicate <- function(
    scored,
    decision_thresholds,
    replicate) {
  metric <- sb_survival_decision_curve_metric(
    scored,
    decision_thresholds,
    prediction_horizon = scored$prediction_horizon[[1L]],
    time_unit = scored$time_unit[[1L]],
    scope = "patient_level_external_bootstrap_sample"
  )
  curve <- sb_external_bootstrap_survival_decision_curve_comparisons(
    metric
  )
  data.frame(
    replicate = rep.int(as.integer(replicate), nrow(curve)),
    curve,
    stringsAsFactors = FALSE
  )
}

sb_external_bootstrap_survival_decision_point_comparisons <- function(curve) {
  comparisons <- c(
    "model_minus_treat_all",
    "model_minus_treat_none"
  )
  rows <- lapply(seq_len(nrow(curve)), function(index) {
    do.call(rbind, lapply(comparisons, function(comparison) {
      data.frame(
        threshold_index = curve$threshold_index[[index]],
        threshold_probability = curve$threshold_probability[[index]],
        threshold_odds = curve$threshold_odds[[index]],
        comparison = comparison,
        estimate = curve[[comparison]][[index]],
        status = curve[[paste0(comparison, "_status")]][[index]],
        error_class = curve[[
          paste0(comparison, "_error_class")
        ]][[index]],
        error_message = curve[[
          paste0(comparison, "_error_message")
        ]][[index]],
        stringsAsFactors = FALSE
      )
    }))
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

sb_external_bootstrap_survival_decision_intervals <- function(
    point_comparisons,
    replicate_curves,
    confidence_level,
    reps) {
  metric_key <- paste0(
    "survival_dca_",
    point_comparisons$threshold_index,
    "_",
    point_comparisons$comparison
  )
  point <- data.frame(
    metric = metric_key,
    estimate = point_comparisons$estimate,
    status = point_comparisons$status,
    error_class = point_comparisons$error_class,
    error_message = point_comparisons$error_message,
    stringsAsFactors = FALSE
  )
  replicate_rows <- lapply(seq_len(nrow(point_comparisons)), function(index) {
    comparison <- point_comparisons$comparison[[index]]
    selected <- replicate_curves$threshold_index ==
      point_comparisons$threshold_index[[index]]
    rows <- replicate_curves[selected, , drop = FALSE]
    data.frame(
      replicate = rows$replicate,
      metric = rep.int(metric_key[[index]], nrow(rows)),
      estimate = rows[[comparison]],
      status = rows[[paste0(comparison, "_status")]],
      error_class = rows[[paste0(comparison, "_error_class")]],
      error_message = rows[[paste0(comparison, "_error_message")]],
      stringsAsFactors = FALSE
    )
  })
  replicates <- do.call(rbind, replicate_rows)
  rownames(replicates) <- NULL
  intervals <- sb_external_percentile_intervals(
    point,
    replicates,
    confidence_level,
    reps
  )
  result <- data.frame(
    threshold_index = point_comparisons$threshold_index,
    threshold_probability = point_comparisons$threshold_probability,
    threshold_odds = point_comparisons$threshold_odds,
    comparison = point_comparisons$comparison,
    intervals[, setdiff(names(intervals), "metric"), drop = FALSE],
    stringsAsFactors = FALSE
  )
  target <- paste0(
    "`", result$comparison, "` at threshold ",
    format(result$threshold_probability, trim = TRUE)
  )
  too_few <- !is.na(result$error_class) &
    result$error_class == "stormbreaker_insufficient_bootstrap_replicates"
  result$error_message[too_few] <- paste0(
    "At least two finite successful bootstrap replicates are required for ",
    target[too_few], "."
  )
  partial <- !is.na(result$error_class) &
    result$error_class == "stormbreaker_incomplete_bootstrap_replicates"
  result$error_message[partial] <- paste0(
    "The displayed interval for ", target[partial],
    " uses only its finite successful replicates."
  )
  rownames(result) <- NULL
  result
}

sb_external_bootstrap_survival_decision_curve_result <- function(
    scored,
    decision_thresholds,
    replicate_curves,
    confidence_level,
    reps) {
  prediction_horizon <- scored$prediction_horizon[[1L]]
  time_unit <- scored$time_unit[[1L]]
  point_metric <- sb_survival_decision_curve_metric(
    scored,
    decision_thresholds,
    prediction_horizon,
    time_unit,
    scope = "patient_level_external_bootstrap_point"
  )
  point_curve <-
    sb_external_bootstrap_survival_decision_curve_comparisons(point_metric)
  point_comparisons <-
    sb_external_bootstrap_survival_decision_point_comparisons(point_curve)
  intervals <- sb_external_bootstrap_survival_decision_intervals(
    point_comparisons,
    replicate_curves,
    confidence_level,
    reps
  )
  complete <- nrow(intervals) >= 1L &&
    all(intervals$status == "success")
  list(
    status = if (complete) "success" else "incomplete",
    task_type = "survival",
    scope = "paired_patient_bootstrap_external_survival_decision_curve",
    interval_scope = "pointwise_not_simultaneous",
    event_probability_estimator = paste0(
      "kaplan_meier_within_threshold_selected_group_each_resample"
    ),
    prediction = "fixed",
    threshold_source = "study_declared_survival_task",
    comparison_targets = c(
      "model_minus_treat_all",
      "model_minus_treat_none"
    ),
    prediction_horizon = prediction_horizon,
    time_unit = time_unit,
    n_thresholds = as.integer(length(decision_thresholds)),
    point_curve = point_curve,
    point_comparisons = point_comparisons,
    intervals = intervals,
    replicates = replicate_curves
  )
}
