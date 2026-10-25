sb_result_cohort_label <- function(cohort_label) {
  if (!is.null(cohort_label) &&
      (!is.character(cohort_label) || length(cohort_label) != 1L ||
       is.na(cohort_label) || !nzchar(trimws(cohort_label)))) {
    sb_abort("`cohort_label` must be NULL or one non-empty string.",
             "stormbreaker_bad_cohort_label")
  }
  if (is.null(cohort_label)) NA_character_ else cohort_label
}

sb_external_result_table <- function(x, component, metric, cohort_label) {
  is_external <- function(value) inherits(value, "stormbreaker_external_evaluation")
  results <- if (is_external(x)) list(x) else x$evaluations
  if (!is.list(results) || length(results) == 0L ||
      !all(vapply(results, is_external, logical(1)))) {
    sb_abort("The comparison must contain supported external evaluations.",
             "stormbreaker_bad_result_view")
  }
  label <- sb_result_cohort_label(cohort_label)
  supported <- sb_scalar_view_metrics()
  if (component == "metrics") {
    if (!is.character(metric) || length(metric) != 1L ||
        is.na(metric) || !metric %in% supported) {
      sb_abort(paste0("Unsupported `metric`; explicitly choose one of: ",
                      paste(supported, collapse = ", "), "."),
               "stormbreaker_bad_result_metric")
    }
  } else if (!is.null(metric)) {
    sb_abort('`metric` is only used for `component = "metrics"`; otherwise leave it NULL.',
             "stormbreaker_unused_result_metric")
  }
  value_or <- function(value, default) if (is.null(value)) default else value
  rows <- lapply(results, function(result) {
    if (component == "predictions") {
      tab <- result$scored_predictions
      tab$task_type <- rep(result$task_type, nrow(tab))
      tab$cohort_label <- rep(label, nrow(tab))
      return(tab)
    }
    if (component == "runs") {
      return(data.frame(
        pipeline_id = result$pipeline_id, task_type = result$task_type,
        status = result$status, cohort_label = label, n = result$n,
        patient_match_status = result$patient_match$status,
        error_stage = result$error_stage, error_class = result$error_class,
        error_message = result$error_message, stringsAsFactors = FALSE
      ))
    }
    if (metric %in% sb_horizon_view_metrics()) {
      return(sb_horizon_metric_row(result, metric, result$task_type,
        external = TRUE, cohort_label = cohort_label))
    }
    allowed <- if (identical(result$task_type, "binary")) supported[1:4] else "harrell_c"
    value <- result$metrics[[metric]]
    if (!metric %in% allowed ||
        (is.null(value) && !identical(result$status, "failed"))) {
      sb_abort(paste0("Metric `", metric, "` is not saved in pipeline `",
                      result$pipeline_id, "`. No metric is recomputed by this function."),
               "stormbreaker_result_metric_not_saved")
    }
    data.frame(
      pipeline_id = result$pipeline_id, task_type = result$task_type,
      metric = metric, estimate = value_or(value$estimate, NA_real_),
      status = value_or(value$status, "not_saved"), run_status = result$status,
      scope = value_or(value$scope, NA_character_),
      aggregation = value_or(value$aggregation, NA_character_),
      cohort_label = label,
      patient_match_status = result$patient_match$status,
      error_stage = result$error_stage, error_class = result$error_class,
      error_message = result$error_message,
      metric_error_class = value_or(value$error_class, NA_character_),
      metric_error_message = value_or(value$error_message, NA_character_),
      stringsAsFactors = FALSE
    )
  })
  if (component == "predictions" && any(vapply(rows, nrow, integer(1)) > 0L)) {
    rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
  }
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}
