sb_binary_decision_curve_metric <- function(
    predictions,
    decision_thresholds,
    scope) {
  decision_thresholds <- sb_validate_decision_thresholds(
    decision_thresholds
  )
  n <- as.integer(nrow(predictions))
  outcome <- predictions$outcome
  probability <- predictions$probability
  prediction_status <- predictions$status
  n_positive <- as.integer(sum(outcome == 1L, na.rm = TRUE))
  n_negative <- as.integer(sum(outcome == 0L, na.rm = TRUE))
  threshold_odds <- decision_thresholds / (1 - decision_thresholds)

  curve <- data.frame(
    threshold_probability = decision_thresholds,
    threshold_odds = threshold_odds,
    n = rep.int(n, length(decision_thresholds)),
    n_positive = rep.int(n_positive, length(decision_thresholds)),
    n_negative = rep.int(n_negative, length(decision_thresholds)),
    n_selected = rep.int(NA_integer_, length(decision_thresholds)),
    true_positive = rep.int(NA_integer_, length(decision_thresholds)),
    false_positive = rep.int(NA_integer_, length(decision_thresholds)),
    true_negative = rep.int(NA_integer_, length(decision_thresholds)),
    false_negative = rep.int(NA_integer_, length(decision_thresholds)),
    sensitivity = rep.int(NA_real_, length(decision_thresholds)),
    specificity = rep.int(NA_real_, length(decision_thresholds)),
    positive_predictive_value = rep.int(
      NA_real_,
      length(decision_thresholds)
    ),
    negative_predictive_value = rep.int(
      NA_real_,
      length(decision_thresholds)
    ),
    model_net_benefit = rep.int(NA_real_, length(decision_thresholds)),
    treat_all_net_benefit = rep.int(
      NA_real_,
      length(decision_thresholds)
    ),
    treat_none_net_benefit = rep.int(NA_real_, length(decision_thresholds)),
    stringsAsFactors = FALSE
  )
  result <- list(
    metric = "decision_curve",
    task_type = "binary",
    scope = scope,
    purpose = "descriptive_clinical_utility_evaluation",
    threshold_source = "study_declared_binary_task",
    decision_rule = "probability_greater_than_or_equal_to_threshold",
    net_benefit_formula = paste0(
      "true_positive / n - false_positive / n * threshold / ",
      "(1 - threshold)"
    ),
    reference_strategies = c("treat_all", "treat_none"),
    status = "incomplete",
    n = n,
    n_scored = 0L,
    n_positive = n_positive,
    n_negative = n_negative,
    n_thresholds = as.integer(length(decision_thresholds)),
    curve = curve,
    error_class = "stormbreaker_incomplete_decision_curve_input",
    error_message = paste0(
      "Binary outcomes, probabilities, or prediction statuses are ",
      "incomplete or invalid."
    )
  )

  complete <- n >= 1L &&
    (is.numeric(outcome) || is.logical(outcome)) &&
    length(outcome) == n && !anyNA(outcome) &&
    (!is.numeric(outcome) || all(is.finite(outcome))) &&
    all(outcome %in% c(0, 1)) &&
    is.numeric(probability) && length(probability) == n &&
    !anyNA(probability) && all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1) &&
    is.character(prediction_status) && length(prediction_status) == n &&
    !anyNA(prediction_status) && all(prediction_status == "success")
  if (!complete) {
    return(result)
  }

  outcome <- as.integer(outcome)
  true_positive <- vapply(decision_thresholds, function(threshold) {
    as.integer(sum(probability >= threshold & outcome == 1L))
  }, integer(1))
  false_positive <- vapply(decision_thresholds, function(threshold) {
    as.integer(sum(probability >= threshold & outcome == 0L))
  }, integer(1))
  true_negative <- vapply(decision_thresholds, function(threshold) {
    as.integer(sum(probability < threshold & outcome == 0L))
  }, integer(1))
  false_negative <- vapply(decision_thresholds, function(threshold) {
    as.integer(sum(probability < threshold & outcome == 1L))
  }, integer(1))
  n_selected <- true_positive + false_positive
  n_not_selected <- true_negative + false_negative

  result$status <- "success"
  result$n_scored <- n
  result$curve$n_selected <- n_selected
  result$curve$true_positive <- true_positive
  result$curve$false_positive <- false_positive
  result$curve$true_negative <- true_negative
  result$curve$false_negative <- false_negative
  result$curve$sensitivity <- if (n_positive > 0L) {
    true_positive / n_positive
  } else {
    rep.int(NA_real_, length(decision_thresholds))
  }
  result$curve$specificity <- if (n_negative > 0L) {
    true_negative / n_negative
  } else {
    rep.int(NA_real_, length(decision_thresholds))
  }
  has_selected <- n_selected > 0L
  has_not_selected <- n_not_selected > 0L
  result$curve$positive_predictive_value[has_selected] <-
    true_positive[has_selected] / n_selected[has_selected]
  result$curve$negative_predictive_value[has_not_selected] <-
    true_negative[has_not_selected] / n_not_selected[has_not_selected]
  result$curve$model_net_benefit <- true_positive / n -
    false_positive / n * threshold_odds
  result$curve$treat_all_net_benefit <- n_positive / n -
    n_negative / n * threshold_odds
  result$curve$treat_none_net_benefit <- rep.int(
    0,
    length(decision_thresholds)
  )
  result$error_class <- NA_character_
  result$error_message <- NA_character_
  result
}
