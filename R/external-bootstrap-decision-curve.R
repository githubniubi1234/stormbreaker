sb_external_bootstrap_decision_thresholds <- function(evaluation) {
  metrics <- evaluation[["metrics"]]
  if (!is.list(metrics) || !"decision_curve" %in% names(metrics)) {
    return(NULL)
  }
  decision_curve <- metrics[["decision_curve"]]
  if (!is.list(decision_curve)) {
    sb_abort(
      "The external decision-curve evidence is malformed.",
      "stormbreaker_bad_external_evaluation"
    )
  }
  curve <- decision_curve[["curve"]]
  if (!is.data.frame(curve)) {
    sb_abort(
      "The external decision-curve evidence is malformed.",
      "stormbreaker_bad_external_evaluation"
    )
  }
  positions <- which(names(curve) == "threshold_probability")
  if (length(positions) != 1L) {
    sb_abort(
      "The external decision-curve threshold column is missing or duplicated.",
      "stormbreaker_bad_external_evaluation"
    )
  }
  sb_validate_decision_thresholds(curve[[positions]])
}

sb_external_bootstrap_decision_curve_comparisons <- function(curve) {
  curve$model_minus_treat_all <-
    (curve$true_negative * curve$threshold_odds -
       curve$false_negative) / curve$n
  curve$model_minus_treat_none <- curve$model_net_benefit
  curve
}

sb_external_bootstrap_decision_curve_replicate <- function(
    scored,
    decision_thresholds,
    replicate) {
  metric <- sb_binary_decision_curve_metric(
    scored,
    decision_thresholds,
    scope = "patient_level_external_bootstrap_sample"
  )
  curve <- sb_external_bootstrap_decision_curve_comparisons(metric$curve)
  data.frame(
    replicate = rep.int(as.integer(replicate), nrow(curve)),
    curve,
    stringsAsFactors = FALSE
  )
}

sb_external_bootstrap_decision_point_comparisons <- function(curve) {
  data.frame(
    threshold_probability = curve$threshold_probability,
    threshold_odds = curve$threshold_odds,
    model_minus_treat_all = curve$model_minus_treat_all,
    model_minus_treat_none = curve$model_minus_treat_none,
    stringsAsFactors = FALSE
  )
}

sb_external_bootstrap_decision_intervals <- function(
    point_comparisons,
    replicate_curves,
    confidence_level,
    reps) {
  probabilities <- sb_external_bootstrap_interval_method(
    confidence_level
  )$probabilities
  comparison_targets <- c(
    "model_minus_treat_all",
    "model_minus_treat_none"
  )
  rows <- lapply(seq_len(nrow(point_comparisons)), function(index) {
    threshold <- point_comparisons$threshold_probability[[index]]
    replicate_rows <-
      replicate_curves$threshold_probability == threshold
    lapply(comparison_targets, function(comparison) {
      bounds <- stats::quantile(
        replicate_curves[[comparison]][replicate_rows],
        probs = probabilities,
        names = FALSE,
        type = 7
      )
      data.frame(
        threshold_probability = threshold,
        threshold_odds = point_comparisons$threshold_odds[[index]],
        comparison = comparison,
        point_estimate = point_comparisons[[comparison]][[index]],
        lower = unname(bounds[[1L]]),
        upper = unname(bounds[[2L]]),
        confidence_level = confidence_level,
        n_requested = as.integer(reps),
        stringsAsFactors = FALSE
      )
    })
  })
  result <- do.call(rbind, unlist(rows, recursive = FALSE))
  rownames(result) <- NULL
  result
}

sb_external_bootstrap_decision_curve_result <- function(
    scored,
    decision_thresholds,
    replicate_curves,
    confidence_level,
    reps) {
  point_metric <- sb_binary_decision_curve_metric(
    scored,
    decision_thresholds,
    scope = "patient_level_external_bootstrap_point"
  )
  point_curve <- sb_external_bootstrap_decision_curve_comparisons(
    point_metric$curve
  )
  point_comparisons <- sb_external_bootstrap_decision_point_comparisons(
    point_curve
  )
  list(
    status = "success",
    scope = "paired_patient_bootstrap_external_decision_curve",
    interval_scope = "pointwise_not_simultaneous",
    comparison_targets = c(
      "model_minus_treat_all",
      "model_minus_treat_none"
    ),
    n_thresholds = as.integer(length(decision_thresholds)),
    point_curve = point_curve,
    point_comparisons = point_comparisons,
    intervals = sb_external_bootstrap_decision_intervals(
      point_comparisons,
      replicate_curves,
      confidence_level,
      reps
    ),
    replicates = replicate_curves
  )
}
