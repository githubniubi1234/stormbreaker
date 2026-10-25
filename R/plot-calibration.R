#' Plot saved calibration curves
#'
#' Display existing descriptive calibration curves without fitting a new
#' calibration equation or changing model predictions.
#'
#' @param x A result from [sb_run_cv()], [sb_compare_cv()],
#'   [sb_evaluate_external()], or [sb_compare_external()]. Each comparison is
#'   displayed in supplied order with one panel per usable pipeline. Whole
#'   repeated CV, IECV and stacking CV are supported with separate context
#'   panels. IECV remains development validation, not final external validation.
#' @param cohort_label Optional non-empty display label for the cohort. The
#'   development/external evaluation boundary is retained in returned metadata.
#'   Default figures omit process/status text and explanatory footnotes. A label does
#'   not establish independence or external validity.
#' @param show_unavailable Logical; default FALSE hides unavailable panels.
#'   TRUE retains labeled audit panels for failed or unrequested curves.
#'
#' @return A `ggplot` object, requiring optional `ggplot2`. `plot$data` contains
#'   copied `x`/`y` coordinates, source row numbers, pipeline/cohort/scope and
#'   status information. The `stormbreaker_curve_status` attribute is a data
#'   frame with one diagnostic row per pipeline. No file is written.
#'   `stormbreaker_curve_data` retains all source coordinates, including hidden
#'   pipelines; `plot$data` contains only displayed pipelines by default.
#'
#' @details Binary coordinates come from the saved moderate-calibration curve.
#'   Survival coordinates are the saved predicted and calibrated *event*
#'   probabilities at the declared horizon, not survival probabilities. The
#'   axis target and panel horizon/unit identify this distinction. The
#'   calibrated coordinate is the existing smoothed descriptive estimate of
#'   observed probability, not a new model prediction to apply to patients.
#'
#'   Both metric and complete-run status must be successful for a curve to be
#'   drawn, and at least one saved line segment must be drawable. Otherwise its
#'   panel is hidden by default; source coordinates remain in diagnostic
#'   attributes and the original result. If none are usable, a no-results
#'   message is shown. Use [sb_export_failures()] for saved failure details.
#'   Missing/non-finite coordinates are not filled, and lines do not bridge
#'   excluded coordinates. An absent survival horizon is shown as not requested
#'   for an otherwise successful run. Missing information after a failure is
#'   not relabeled as an unrequested analysis.
#'
#'   The dashed identity line is an ideal-calibration reference. No confidence
#'   band, refitting, smoothing, ranking or threshold selection is performed.
#'   Curves remain development or external evaluation as recorded by the source.
#'   Comparisons can require a taller output image for readable panel labels.
#'
#' @seealso [sb_plot_dca()], [sb_plot_performance()]
#' @export
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   d <- sb_example_binary(120, seed = 1)
#'   task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
#'   result <- sb_run_cv(
#'     task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1
#'   )
#'   sb_plot_calibration(result)
#' }
sb_plot_calibration <- function(x, cohort_label = NULL, show_unavailable = FALSE) {
  sb_require_curve_plotting()
  original <- sb_curve_view_data(x, "calibration", cohort_label)
  view <- sb_curve_display_view(original, "calibration", show_unavailable)
  if (!nrow(view$status)) {
    p <- sb_empty_display_plot("Saved calibration curve",
      "No successful saved calibration curves are available.",
      "Use sb_export_failures() for details; unrequested curves are not failures.")
    p$data <- view$data
    attr(p, "stormbreaker_curve_status") <- original$status
    attr(p, "stormbreaker_curve_data") <- original$data
    return(p)
  }
  lines <- sb_curve_view_lines(view$data)
  annotations <- sb_curve_view_annotations(view)
  survival <- all(view$status$task_type == "survival")
  target <- if (survival) "event probability" else "probability"
  panel_label <- function(labels) {
    if (nrow(view$status) == 1L) return(rep(if (survival)
      paste(view$status$prediction_horizon, view$status$time_unit) else "", length(labels)))
    sub("\ncurve:.*", "", labels)
  }
  x <- y <- panel <- label <- line_group <- NULL
  p <- ggplot2::ggplot(view$data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey55") +
    ggplot2::geom_path(data = lines, ggplot2::aes(group = line_group),
                       colour = "#24547A") +
    ggplot2::geom_text(data = annotations,
                       ggplot2::aes(x = 0.5, y = 0.5, label = label),
                       inherit.aes = FALSE, colour = "grey35", size = 3.5) +
    ggplot2::facet_wrap(~panel, ncol = 1, drop = FALSE,
      labeller = ggplot2::as_labeller(panel_label)) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = "Calibration",
      x = paste("Predicted", target), y = paste("Observed", target)
    ) +
    sb_viz_theme()
  attr(p, "stormbreaker_curve_status") <- original$status
  attr(p, "stormbreaker_curve_data") <- original$data
  p
}
