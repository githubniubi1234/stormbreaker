#' Plot saved decision curves
#'
#' Display saved model, treat-all and treat-none net benefits at the original
#' study-declared threshold probabilities, without recalculating utility.
#'
#' @param x A result from [sb_run_cv()], [sb_compare_cv()],
#'   [sb_evaluate_external()], or [sb_compare_external()]. Comparisons keep one
#'   panel per usable pipeline in supplied order, with references from that same
#'   pipeline's saved curve. Whole repeated CV, IECV and stacking CV retain
#'   separate context panels and reference curves; they are never pooled.
#' @param cohort_label Optional non-empty cohort display label. Development
#'   versus external evaluation remains explicit; a name does not prove
#'   cohort independence.
#' @param show_unavailable Logical; default FALSE hides unavailable panels.
#' @param standardized Divide binary net benefit by the cohort event prevalence.
#'   Not implemented for survival outcomes. No model is refitted.
#' @param cost_benefit_axis Add a secondary threshold-odds cost:benefit scale.
#'   TRUE retains the previous labeled audit panels.
#'
#' @return A `ggplot` object using optional `ggplot2`. `plot$data` preserves
#'   original threshold (`x`) and net-benefit (`y`) coordinates for each
#'   `strategy`, with source row numbers and provenance. The attribute
#'   `stormbreaker_curve_status` contains one diagnostic row per pipeline.
#'   `stormbreaker_curve_data` retains all original coordinates, including
#'   hidden pipelines; `plot$data` contains displayed pipelines by default.
#'   No file is written.
#'
#' @details DCA must already have been requested using `decision_thresholds`
#'   when declaring the development task. This plotting function accepts no
#'   new thresholds and cannot compute missing curves. Survival thresholds
#'   apply to event probability by the saved horizon; its time unit and horizon
#'   are displayed. The saved censoring-aware estimator is not replaced.
#'
#'   Both curve and complete-run status must be successful for points to be
#'   drawn. A pipeline also needs at least one usable model net-benefit point.
#'   Unrequested and failed curves are hidden by default but retain diagnostic
#'   attributes. If no pipelines are usable, an explicit no-results message
#'   replaces the plot. Use [sb_export_failures()] for failure details.
#'   Non-finite or unsuccessful row-level coordinates are not
#'   drawn; separate line segments prevent bridging their gaps. Original values
#'   remain in `plot$data`, including negative net benefits. Empty panels use
#'   explicit missing-coordinate placeholders rather than invented values.
#'
#'   Points identify the actual saved thresholds; connecting segments are
#'   visual guides, not newly evaluated thresholds or a smoothed curve. No
#'   confidence interval, optimal threshold, ranking or clinical recommendation
#'   is inferred. Reference lines remain within their original pipeline panels.
#'   Use a taller output image for comparisons with many pipelines.
#'
#' @seealso [sb_plot_calibration()], [sb_plot_performance()]
#' @export
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   d <- sb_example_binary(120, seed = 1)
#'   task <- sb_binary_task(
#'     d, "patient_id", "outcome", "gene_a",
#'     decision_thresholds = c(0.1, 0.2, 0.3, 0.5)
#'   )
#'   result <- sb_run_cv(
#'     task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1
#'   )
#'   sb_plot_dca(result)
#' }
sb_plot_dca <- function(x, cohort_label = NULL, show_unavailable = FALSE,
                       standardized = FALSE, cost_benefit_axis = FALSE) {
  sb_require_curve_plotting()
  for (flag in list(standardized, cost_benefit_axis)) if (!is.logical(flag) ||
      length(flag) != 1L || is.na(flag)) stop("DCA display options must be TRUE or FALSE.")
  original <- sb_curve_view_data(x, "dca", cohort_label)
  models <- sb_clinician_models(sb_viz_leaves(x, cohort_label))
  if (standardized) {
    leaves <- sb_viz_leaves(x, cohort_label)
    if (!all(vapply(leaves, function(z) identical(z$value$task_type, "binary") ||
      inherits(z$value, "stormbreaker_binary_cv_result"), logical(1)))) stop("Standardized DCA currently supports binary outcomes only.")
    d <- original$data
    keys <- intersect(c("view_id", "pipeline_id", "panel"), names(d))
    key <- interaction(d[keys], drop = TRUE)
    for (k in levels(key)) {
      rows <- which(key == k)
      ref <- rows[d$strategy[rows] == "treat_all" & is.finite(d$y[rows]) & is.finite(d$x[rows])]
      prevalence <- d$y[ref] * (1 - d$x[ref]) + d$x[ref]
      if (!length(prevalence)) {
        if (any(is.finite(d$y[rows]))) stop("Standardization requires the saved treat-all reference.")
        next
      }
      if (any(prevalence <= 0) || diff(range(prevalence)) > 1e-8) stop("A positive common prevalence is required for standardized net benefit.")
      d$y[rows] <- d$y[rows] / mean(prevalence)
    }
    original$data <- d
  }
  view <- sb_curve_display_view(original, "dca", show_unavailable)
  if (!nrow(view$status)) {
    p <- sb_empty_display_plot("Saved decision curve",
      "No successful saved decision curves are available.",
      "Use sb_export_failures() for details; unrequested curves are not failures.")
    p$data <- view$data
    attr(p, "stormbreaker_curve_status") <- original$status
    attr(p, "stormbreaker_curve_data") <- original$data
    return(p)
  }
  lines <- sb_curve_view_lines(view$data)
  points <- view$data[view$data$draw, , drop = FALSE]
  annotations <- sb_curve_view_annotations(view)
  thresholds <- view$data$x[is.finite(view$data$x)]
  threshold_range <- if (length(unique(thresholds)) < 2L) c(0, 1) else range(thresholds)
  annotations$annotation_x <- rep(mean(threshold_range), nrow(annotations))
  x <- y <- strategy <- line_group <- panel <- label <- annotation_x <- NULL
  p <- ggplot2::ggplot(view$data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_path(data = lines,
                       ggplot2::aes(colour = strategy, linetype = strategy,
                                    group = line_group)) +
    ggplot2::geom_point(data = points, ggplot2::aes(colour = strategy), size = 1.5) +
    ggplot2::geom_text(data = annotations,
                       ggplot2::aes(x = annotation_x, y = 0, label = label),
                       inherit.aes = FALSE, colour = "grey35", size = 3.5) +
    ggplot2::facet_wrap(~panel, ncol = 1, drop = FALSE,
      labeller = if ("view_id" %in% names(view$status))
        ggplot2::as_labeller(sb_view_wrap_label) else "label_value") +
    ggplot2::scale_x_continuous(limits = threshold_range,
      sec.axis = if (cost_benefit_axis) ggplot2::dup_axis(name = "Cost:Benefit",
        breaks = c(.01, .1, .2, .5, .8), labels = c("1:99", "1:9", "1:4", "1:1", "4:1")) else ggplot2::waiver()) +
    ggplot2::scale_colour_manual(
      name = NULL, values = c(model = "#24547A", treat_all = "grey50", treat_none = "grey15"),
      limits = c("model", "treat_all", "treat_none"),
      breaks = c("model", "treat_all", "treat_none"),
      labels = c("Model", "Treat all", "Treat none"), drop = FALSE
    ) +
    ggplot2::scale_linetype_manual(
      name = NULL, values = c(model = 1, treat_all = 2, treat_none = 3),
      limits = c("model", "treat_all", "treat_none"),
      breaks = c("model", "treat_all", "treat_none"),
      labels = c("Model", "Treat all", "Treat none"), drop = FALSE
    ) +
    ggplot2::labs(
      title = "Saved decision curve",
      subtitle = sb_curve_view_scope_label(view),
      x = "Threshold probability", y = "Net benefit",
      caption = paste0(
        "Points: saved thresholds; lines: connecting guides. No threshold selection or confidence intervals.\n",
        sum(!view$data$draw & !is.na(view$data$source_row)),
        " saved coordinate(s) not drawn; inspect attr(plot, 'stormbreaker_curve_status')."
      )
    ) +
    sb_viz_theme() +
    ggplot2::theme(legend.position = "bottom")
  attr(p, "stormbreaker_curve_status") <- original$status
  attr(p, "stormbreaker_curve_data") <- original$data
  if (standardized) p <- p + ggplot2::labs(y = "Standardized net benefit")
  attr(p, "stormbreaker_figure_notes") <- p$labels$caption
  p <- sb_clinician_clean_plot(p, models)
  p <- p + ggplot2::labs(title = "Decision curve", subtitle = NULL, caption = NULL)
  p
}
