#' Combine selected clinical views into one figure
#' @param x A result from [sb_clinician_views()].
#' @param views Explicit available plot names in the desired panel order.
#' @param ncol Number of columns in the combined figure.
#' @return A patchwork plot, suitable for ggplot2::ggsave(). The
#'   `stormbreaker_captions` attribute is a data frame of view names and separate
#'   figure legends in panel order. Requires the optional patchwork package.
#' @details Only layout changes: original plot data, layers, axes, legends and
#'   evaluation contexts are retained. Titles use the short selected view names;
#'   original titles and explanations stay in the separate caption table.
#'   Choose related views and export at a readable physical size. For many
#'   cohorts or complex faceted diagnostics, use fewer panels or separate figures
#'   rather than shrinking labels. No files are written and no models are fitted.
#' @export
sb_plot_clinician_panel <- function(x, views, ncol = 2L) {
  if (!inherits(x, "stormbreaker_clinician_views") || !is.character(views) ||
      !length(views) || anyNA(views) || anyDuplicated(views) ||
      !all(views %in% names(x$plots))) sb_abort(
    "Choose unique available plot names from names(x$plots); inspect x$status for unavailable views.",
    "stormbreaker_bad_report_input")
  if (!is.numeric(ncol) || length(ncol) != 1L || !is.finite(ncol) ||
      ncol < 1 || ncol != floor(ncol)) sb_abort(
    "ncol must be a positive integer.", "stormbreaker_bad_report_input")
  if (!requireNamespace("patchwork", quietly = TRUE)) sb_abort(
    "Install patchwork to combine clinical figures.", "stormbreaker_missing_dependency")
  indices <- match(views, names(x$plots))
  captions <- vapply(seq_along(views), function(i) {
    p <- x$plots[[indices[i]]]
    notes <- unique(c(p$labels$title, p$labels$subtitle, p$labels$caption, x$captions[indices[i]]))
    paste(notes[!is.na(notes) & nzchar(notes) & notes != views[i]], collapse = "\n")
  }, character(1))
  plots <- lapply(seq_along(views), function(i) {
    x$plots[[indices[i]]] + ggplot2::labs(title = views[i], subtitle = NULL, caption = NULL)
  })
  result <- patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(tag_levels = "A")
  attr(result, "stormbreaker_captions") <- data.frame(view = views, caption = captions)
  result
}
