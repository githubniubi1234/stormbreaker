#' Plot a binary confusion matrix at an explicit threshold
#' @param x A result from [sb_threshold_metrics()].
#' @param threshold One threshold already evaluated in x.
#' @param normalize Show counts or proportions within each observed class.
#' @param ncol Panel columns.
#' @return A ggplot object. All four cells are retained. A missing observed
#'   class has undefined normalized proportions, not artificial zeros.
#' @export
sb_plot_confusion <- function(x, threshold, normalize = c("none", "observed"), ncol = 2L) {
  sb_require_curve_plotting()
  normalize <- match.arg(normalize)
  ncol <- sb_viz_positive_integer(ncol, "ncol")
  if (!inherits(x, "stormbreaker_threshold_metrics") || !is.numeric(threshold) ||
      length(threshold) != 1L || !is.finite(threshold) || threshold < 0 || threshold > 1) {
    sb_abort("Supply threshold results and one explicit threshold in [0, 1].",
      "stormbreaker_bad_visualization_argument")
  }
  if (!nrow(x$table)) {
    p <- sb_empty_display_plot("Confusion matrix", "No usable binary predictions.",
      "Unavailable groups remain in diagnostics.")
  } else {
    if (!threshold %in% x$table$threshold) sb_abort(
      "This threshold was not evaluated; call sb_threshold_metrics first.",
      "stormbreaker_bad_visualization_argument")
    tab <- x$table[x$table$threshold == threshold, , drop = FALSE]
    tab$panel_key <- interaction(tab$view_id, tab$pipeline_id,
      addNA(factor(tab$fold)), drop = TRUE)
    data <- do.call(rbind, lapply(seq_len(nrow(tab)), function(i) {
      r <- tab[i, , drop = FALSE]
      d <- r[rep(1L, 4L), , drop = FALSE]
      d$observed <- factor(c(0, 0, 1, 1), levels = c(0, 1))
      d$predicted <- factor(c(0, 1, 0, 1), levels = c(0, 1))
      d$count <- c(r$tn, r$fp, r$fn, r$tp)
      denominator <- c(rep(r$tn + r$fp, 2), rep(r$fn + r$tp, 2))
      d$value <- if (normalize == "none") d$count else
        ifelse(denominator == 0, NA_real_, d$count / denominator)
      d$label <- if (normalize == "none") as.character(d$count) else
        ifelse(is.na(d$value), paste0(d$count, "\nUndefined"),
          paste0(d$count, "\n", sprintf("%.1f%%", 100 * d$value)))
      d
    }))
    predicted <- observed <- value <- label <- panel_key <- NULL
    p <- ggplot2::ggplot(data, ggplot2::aes(predicted, observed, fill = value)) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::geom_text(ggplot2::aes(label = label), size = 3.5) +
      ggplot2::scale_fill_gradient(low = "#F2F6F9", high = "#80B1CF", na.value = "grey90") +
      ggplot2::facet_wrap(~panel_key, ncol = ncol,
        labeller = ggplot2::as_labeller(stats::setNames(tab$panel, tab$panel_key))) +
      ggplot2::labs(title = "Which patients are classified correctly?",
        subtitle = paste("Explicit threshold:", threshold), x = "Predicted class", y = "Observed class",
        fill = if (normalize == "none") "Count" else "Within observed class",
        caption = "Class 1 is positive; probability >= threshold. Separate folds/cohorts. No threshold selected.") +
      sb_viz_theme()
    attr(p, "stormbreaker_suggested_size") <- c(width = 11,
      height = max(4, 3 * ceiling(nlevels(tab$panel_key) / ncol)))
  }
  attr(p, "stormbreaker_threshold_diagnostics") <- x$diagnostics
  p
}
