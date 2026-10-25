#' Summarize saved binary probabilities by observed class
#'
#' @param x Saved binary CV or external results, or a named collection.
#' @param breaks Increasing bin boundaries, starting at zero and ending at one.
#' @param cohort_label Optional cohort display label.
#' @return Aggregate bin counts and within-class fractions, group diagnostics,
#'   and bin rules. Bins are left-closed, right-open except the last bin, which
#'   includes one. No patient identifiers, model refitting, or smoothing.
#' @export
sb_probability_distribution <- function(x, breaks = seq(0, 1, .1),
                                        cohort_label = NULL) {
  if (!is.numeric(breaks) || length(breaks) < 2L || any(!is.finite(breaks)) ||
      breaks[1] != 0 || utils::tail(breaks, 1) != 1 || any(diff(breaks) <= 0)) sb_abort(
    "Supply strictly increasing bin boundaries from zero to one.",
    "stormbreaker_bad_visualization_argument")
  thresholds <- sb_threshold_metrics(x, breaks, cohort_label)
  tab <- thresholds$table
  rows <- list()
  if (nrow(tab)) {
    for (start in seq.int(1L, nrow(tab), by = length(breaks))) {
      block <- tab[start + seq_along(breaks) - 1L, , drop = FALSE]
      for (observed in c(0L, 1L)) {
        cumulative <- if (observed == 1L) block$tp else block$fp
        # Last bin includes probability == 1 instead of subtracting it.
        count <- if (length(breaks) == 2L) cumulative[1L] else
          utils::head(cumulative, -1L) - c(cumulative[seq.int(2L, length(breaks) - 1L)], 0)
        metadata <- block[rep(1L, length(breaks) - 1L),
          names(thresholds$diagnostics), drop = FALSE]
        metadata$lower <- utils::head(breaks, -1L)
        metadata$upper <- utils::tail(breaks, -1L)
        metadata$observed_class <- observed
        metadata$count <- count
        metadata$class_n <- cumulative[1L]
        metadata$fraction <- if (cumulative[1L] == 0) NA_real_ else count / cumulative[1L]
        rows[[length(rows) + 1L]] <- metadata
      }
    }
  }
  structure(list(table = if (length(rows)) sb_view_bind(rows) else data.frame(),
    diagnostics = thresholds$diagnostics, breaks = breaks,
    bin_rule = "[lower, upper), except final bin [lower, 1]",
    saved_metrics = thresholds$saved_metrics),
    class = "stormbreaker_probability_distribution")
}

#' Plot saved probability distributions without patient identifiers
#' @param x An [sb_probability_distribution()] result.
#' @param scale Display count or fraction within each observed class.
#' @param ncol Number of panel columns.
#' @return A ggplot with adjacent class histograms for each fold or cohort.
#'   Missing classes have zero counts and undefined fractions, not fabricated
#'   density estimates. All group diagnostics are attached to the plot.
#' @export
sb_plot_probability_distribution <- function(x, scale = c("count", "fraction"),
                                             ncol = 2L) {
  sb_require_curve_plotting()
  if (!inherits(x, "stormbreaker_probability_distribution")) sb_abort(
    "Supply an sb_probability_distribution result.", "stormbreaker_bad_visualization_argument")
  scale <- match.arg(scale)
  ncol <- sb_viz_positive_integer(ncol, "ncol")
  d <- x$table
  if (!nrow(d)) {
    p <- sb_empty_display_plot("Prediction probability distribution",
      "No usable binary predictions.", "Unavailable groups remain in diagnostics.")
  } else {
    d$value <- d[[scale]]
    midpoint <- (d$lower + d$upper) / 2
    gap <- (d$upper - d$lower) * .03
    d$xmin <- ifelse(d$observed_class == 0, d$lower + gap, midpoint + gap)
    d$xmax <- ifelse(d$observed_class == 0, midpoint - gap, d$upper - gap)
    d$class <- factor(d$observed_class, levels = 0:1, labels = c("Observed 0", "Observed 1"))
    xmin <- xmax <- value <- class <- panel <- NULL
    p <- ggplot2::ggplot(d, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = value, fill = class)) +
      ggplot2::geom_rect(na.rm = TRUE) +
      ggplot2::facet_wrap(~panel, ncol = ncol, labeller = ggplot2::as_labeller(sb_view_wrap_label)) +
      ggplot2::scale_fill_manual(values = c("#2878A0", "#D77935"), drop = FALSE) +
      ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
      ggplot2::labs(title = "Prediction probability by observed outcome",
        x = "Predicted probability of class 1",
        y = if (scale == "count") "Patients per bin" else "Fraction within observed class",
        fill = NULL,
        caption = "Adjacent bars share each probability bin. Last bin includes 1.\nSeparate folds/cohorts; no smoothing. Absent-class fractions are undefined.") +
      sb_viz_theme()
  }
  attr(p, "stormbreaker_probability_diagnostics") <- x$diagnostics
  p
}
