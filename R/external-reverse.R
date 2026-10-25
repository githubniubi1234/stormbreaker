sb_external_reverse_evaluation <- function(prediction, endpoint) {
  original <- sb_evaluate_external(prediction, endpoint, reverse = FALSE)
  if (identical(original$status, "failed")) {
    return(original)
  }
  if (!identical(original$task_type, "binary")) {
    sb_abort("`reverse = TRUE` supports binary evaluation only.",
             "stormbreaker_bad_reverse")
  }
  auc <- original$metrics$roc_auc$estimate
  evaluable <- identical(original$metrics$roc_auc$status, "success") &&
    is.finite(auc)
  reversed <- evaluable && auc < 0.5
  result <- original
  if (reversed) {
    adjusted <- prediction
    adjusted$predictions$probability <- 1 - adjusted$predictions$probability
    result <- sb_evaluate_external(adjusted, endpoint, reverse = FALSE)
  }
  result$original_metrics <- original$metrics
  result$original_scored_predictions <- original$scored_predictions
  result$scored_predictions$probability_original <-
    original$scored_predictions$probability
  result$scored_predictions$reversed <- rep(reversed, nrow(result$scored_predictions))
  result$direction_correction <- data.frame(
    reverse = TRUE,
    reversed = reversed,
    status = if (evaluable) "success" else "not_evaluable",
    auc_original = auc,
    auc_processed = result$metrics$roc_auc$estimate,
    label_use = "post_hoc_evaluation_label_direction_selection",
    stringsAsFactors = FALSE
  )
  result$evaluation_boundary$direction_selection_uses_external_outcomes <- TRUE
  result$evaluation_boundary$outcome_use <- "metrics_and_post_hoc_direction_selection"
  result$evaluation_boundary$prediction_values_changed <- reversed
  result$provenance$direction_correction <- result$direction_correction
  result
}
