sb_horizon_view_metrics <- function() c("time_auc", "uno_c", "ipcw_brier")

sb_scalar_view_metrics <- function() {
  c("roc_auc", "average_precision", "brier_score", "log_loss", "harrell_c",
    sb_horizon_view_metrics())
}

sb_horizon_metric_row <- function(result, metric, task_type, external = FALSE,
                                 cohort_label = NULL) {
  id <- if (external) result$pipeline_id else result$pipeline$id
  if (!identical(task_type, "survival")) {
    sb_abort(paste0("Metric `", metric, "` is not saved in pipeline `", id,
      "`. No metric is recomputed by this function."),
      "stormbreaker_result_metric_not_saved")
  }
  value_or <- function(value, default) if (is.null(value)) default else value
  descriptor <- if (external) result$endpoint else result$task
  value <- result$metrics[[metric]]
  horizon <- value_or(value$prediction_horizon,
    value_or(descriptor$prediction_horizon, NA_real_))
  unit <- value_or(value$time_unit, value_or(descriptor$time_unit, NA_character_))
  status <- if (!is.null(value)) value$status else if (
    identical(result$status, "success") && is.na(horizon)) "not_requested" else
      "not_saved"
  output <- data.frame(
    pipeline_id = id, task_type = task_type, metric = metric,
    estimate = value_or(value$estimate, NA_real_), status = status,
    run_status = result$status, scope = value_or(value$scope, NA_character_),
    aggregation = value_or(value$aggregation, NA_character_),
    prediction_horizon = as.numeric(horizon), time_unit = unit,
    metric_error_class = value_or(value$error_class, NA_character_),
    metric_error_message = value_or(value$error_message, NA_character_),
    stringsAsFactors = FALSE
  )
  if (external) {
    output$cohort_label <- sb_result_cohort_label(cohort_label)
    output$patient_match_status <- result$patient_match$status
    output$error_stage <- result$error_stage
    output$error_class <- result$error_class
    output$error_message <- result$error_message
  }
  output
}

sb_horizon_view_labels <- function(tab) {
  horizon <- vapply(tab$prediction_horizon, function(value) {
    if (is.na(value)) "horizon not saved" else format(value, digits = 17, trim = TRUE)
  }, character(1))
  unit <- tab$time_unit
  unit[is.na(unit)] <- "unit not saved"
  paste0(as.character(tab$pipeline_id), " [", horizon, " ", unit, "]")
}
