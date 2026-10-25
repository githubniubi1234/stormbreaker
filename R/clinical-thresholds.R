#' Evaluate explicit binary decision thresholds without selecting one
#'
#' @param x Saved binary CV or external results, or a named collection.
#' @param thresholds Explicit numeric probabilities in `[0, 1]`. Equal scores
#'   are classified together: probability greater than or equal to threshold
#'   means positive. No threshold is selected or optimized.
#' @param cohort_label Optional cohort display label.
#' @return A list with aggregate `table`, `diagnostics`, and saved metrics.
#'   Folds and cohorts remain separate. Undefined ratios are NA, not zero.
#'   Patient identifiers are not included in the returned tables.
#'   `kappa` is Cohen's chance-corrected agreement at each threshold; it can
#'   be negative and is NA when expected agreement is one.
#' @export
sb_threshold_metrics <- function(x, thresholds, cohort_label = NULL) {
  if (!is.numeric(thresholds) || !length(thresholds) ||
      any(!is.finite(thresholds)) || any(thresholds < 0 | thresholds > 1)) {
    sb_abort("Supply explicit finite thresholds in [0, 1].",
      "stormbreaker_bad_visualization_argument")
  }
  thresholds <- sort(unique(thresholds))
  leaves <- sb_viz_leaves(x, cohort_label)
  tables <- diagnostics <- saved <- list()
  ratio <- function(a, b) if (b == 0) NA_real_ else a / b
  for (leaf in leaves) {
    value <- leaf$value
    external <- inherits(value, "stormbreaker_external_evaluation")
    binary <- if (external) identical(value$task_type, "binary") else
      inherits(value, "stormbreaker_binary_cv_result")
    if (!binary) sb_abort("Decision thresholds require binary probabilities.",
      "stormbreaker_bad_visualization_task")
    pred <- if (external) value$scored_predictions else value$predictions
    saved[[leaf$context$view_id]] <- value$metrics
    folds <- if (!external && is.data.frame(pred) && nrow(pred) &&
      "fold" %in% names(pred)) unique(pred$fold) else NA_integer_
    for (fold in folds) {
      group <- if (!is.data.frame(pred)) data.frame() else if (external || is.na(fold)) pred else
        pred[!is.na(pred$fold) & pred$fold == fold, , drop = FALSE]
      status <- if (!identical(value$status, "success")) "failed_run" else if (!nrow(group)) "no_predictions" else "success"
      if (status == "success") {
        valid <- all(c("outcome", "probability", "status", "patient_id") %in% names(group)) &&
          is.numeric(group$outcome) && is.numeric(group$probability) &&
          all(is.finite(group$outcome)) && all(group$outcome %in% c(0, 1)) &&
          all(is.finite(group$probability)) && all(group$probability >= 0 & group$probability <= 1) &&
          !anyNA(group$status) && all(group$status == "success") &&
          !anyNA(group$patient_id) && !anyDuplicated(group$patient_id) &&
          (external || ("fold" %in% names(group) && !anyNA(group$fold)))
        if (!valid) status <- "invalid_predictions"
      }
      d <- data.frame(pipeline_id = sb_viz_pipeline_id(value), dataset = leaf$dataset,
        fold = as.integer(fold), panel = paste(sb_viz_panel(leaf),
          if (external) "External patients" else paste("Assessment fold", fold), sep = "\n"),
        status = status, n = nrow(group), stringsAsFactors = FALSE)
      d <- sb_view_annotate(d, leaf$context)
      d$panel <- paste(c(d$pipeline_id,
        if (leaf$dataset != "Unnamed cohort") leaf$dataset,
        if (!is.na(leaf$context$repeat_id)) paste("Repeat", leaf$context$repeat_id),
        if (!is.na(leaf$context$heldout_cohort_id)) leaf$context$heldout_cohort_id,
        if (external) "External" else paste("CV fold", fold)), collapse = "\n")
      diagnostics[[length(diagnostics) + 1L]] <- d
      if (status != "success") next
      rows <- lapply(thresholds, function(threshold) {
        positive <- group$probability >= threshold
        tp <- sum(positive & group$outcome == 1)
        fp <- sum(positive & group$outcome == 0)
        tn <- sum(!positive & group$outcome == 0)
        fn <- sum(!positive & group$outcome == 1)
        n <- nrow(group)
        observed <- (tp + tn) / n
        expected <- ((tp + fp) / n) * ((tp + fn) / n) +
          ((tn + fn) / n) * ((tn + fp) / n)
        data.frame(threshold = threshold, tp = tp, fp = fp, tn = tn, fn = fn,
          sensitivity = ratio(tp, tp + fn), specificity = ratio(tn, tn + fp),
          ppv = ratio(tp, tp + fp), npv = ratio(tn, tn + fn),
          accuracy = ratio(tp + tn, nrow(group)), f1 = ratio(2 * tp, 2 * tp + fp + fn),
          kappa = ratio(observed - expected, 1 - expected),
          positive_fraction = ratio(tp + fp, nrow(group)))
      })
      tab <- do.call(rbind, rows)
      for (name in names(d)) tab[[name]] <- rep(d[[name]], nrow(tab))
      tables[[length(tables) + 1L]] <- tab
    }
  }
  structure(list(table = if (length(tables)) sb_view_bind(tables) else data.frame(),
    diagnostics = if (length(diagnostics)) sb_view_bind(diagnostics) else data.frame(),
    saved_metrics = saved, threshold_rule = "probability >= threshold"),
    class = "stormbreaker_threshold_metrics")
}

#' Plot explicit threshold performance
#' @param x A result from [sb_threshold_metrics()].
#' @param metrics Metrics to display. Defaults to sensitivity, specificity,
#'   positive predictive value and negative predictive value.
#' @param ncol Number of panel columns.
#' @return A ggplot object with group availability attached. Undefined values
#'   remain missing. No model fitting or threshold optimization occurs.
#' @export
sb_plot_thresholds <- function(x,
                               metrics = c("sensitivity", "specificity", "ppv", "npv"),
                               ncol = 2L) {
  sb_require_curve_plotting()
  if (!inherits(x, "stormbreaker_threshold_metrics")) sb_abort(
    "Supply an sb_threshold_metrics result.", "stormbreaker_bad_visualization_argument")
  allowed <- c("sensitivity", "specificity", "ppv", "npv", "accuracy", "f1", "kappa", "positive_fraction")
  if (!is.character(metrics) || !length(metrics) || anyNA(metrics) ||
      anyDuplicated(metrics) || any(!metrics %in% allowed)) sb_abort(
    "Supply distinct supported threshold metrics.", "stormbreaker_bad_visualization_argument")
  ncol <- sb_viz_positive_integer(ncol, "ncol")
  if (!nrow(x$table)) {
    p <- sb_empty_display_plot("Threshold performance", "No usable binary predictions.",
      "Unavailable groups remain in diagnostics.")
  } else {
    data <- do.call(rbind, lapply(metrics, function(m) {
      d <- x$table
      d$metric <- m
      d$value <- d[[m]]
      d
    }))
    threshold <- value <- metric <- panel <- NULL
    p <- ggplot2::ggplot(data, ggplot2::aes(threshold, value, colour = metric)) +
      ggplot2::geom_line(na.rm = TRUE) + ggplot2::geom_point(size = 1, na.rm = TRUE) +
      ggplot2::facet_wrap(~panel, ncol = ncol, labeller = ggplot2::as_labeller(sb_view_wrap_label)) +
      ggplot2::scale_y_continuous(limits = c(if ("kappa" %in% metrics) -1 else 0, 1)) +
      ggplot2::labs(title = "How does the decision threshold change performance?",
        x = "Explicit probability threshold", y = "Metric value", colour = "Metric",
        caption = "Class 1 is positive; probability >= threshold. Separate folds/cohorts. No threshold selected. Undefined ratios are missing.") +
      sb_viz_theme()
  }
  attr(p, "stormbreaker_threshold_diagnostics") <- x$diagnostics
  p
}
