#' Plot saved CV or external performance without ranking
#'
#' Draw one selected scalar performance metric from saved results. This is
#' a descriptive display, not model selection or a new evaluation.
#'
#' @param x A CV result/comparison or external evaluation/single-cohort
#'   comparison, or whole repeated CV, IECV or stacking CV result. Contexts
#'   have distinct labels, without pooling. See [sb_result_table()].
#' @param metric One supported scalar metric, as documented in
#'   [sb_result_table()]. It must already exist in every supplied result.
#'   Failed external evaluations without a saved metric remain in diagnostics.
#' @param cohort_label Optional non-empty display label supplied by the caller,
#'   forwarded to [sb_result_table()]. No cohort identity is inferred.
#' @param show_unavailable Logical; default FALSE excludes unavailable pipeline
#'   labels. TRUE retains the previous audit view with labels but no points.
#'
#' @return A `ggplot` object. Its `data` retains every pipeline and the original
#'   metric values, with `pipeline_id` represented as an input-ordered factor.
#'   No file is written. The optional package `ggplot2` is required only when
#'   this function is called.
#'
#' @details Pipeline labels run from top to bottom in supplied order, not in
#'   performance order. Points require both the selected metric and the complete
#'   run to have status `"success"`, and the estimate to be finite. Otherwise
#'   the pipeline is omitted from the default display, and the caption counts
#'   omissions. If none are usable, an explicit no-results message is shown.
#'   Inspect [sb_result_table()] for metric and run status separately, or export
#'   [sb_export_failures()] for saved failure details.
#'   Estimates in the underlying table are never replaced with zero or clipped.
#'
#'   The subtitle reports the saved metric scope and aggregation. The plot does
#'   distinguish external evaluation from development CV and displays a supplied
#'   cohort label. An omitted external cohort label is explicitly unnamed. It does
#'   not fit models, recompute metrics, infer confidence intervals, or rank
#'   pipelines. The probability-scale metrics and concordance are displayed on
#'   a zero-to-one axis, except that log loss and IPCW Brier have a zero lower
#'   bound and no fixed upper bound. The IPCW weighting can yield estimates
#'   above one, which remain visible. No error bars are drawn.
#'   For `"time_auc"`, `"uno_c"`, and `"ipcw_brier"`, each pipeline label
#'   includes its saved horizon and time unit. Distinct targets remain distinct
#'   labels, not pooled or converted; they are not like-for-like comparisons.
#'   These plots also retain a `display_label` column in their data.
#'   Large comparisons may require a taller figure or fewer predeclared
#'   candidate pipelines for readable labels.
#'
#' @seealso [sb_result_table()]
#' @export
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   d <- sb_example_survival(n = 60, seed = 1)
#'   task <- sb_surv_task(
#'     d, "patient_id", "time", "event", "gene_a", "month"
#'   )
#'   result <- sb_run_cv(task, folds = 3, seed = 1)
#'   sb_plot_performance(result, "harrell_c")
#' }
sb_plot_performance <- function(x, metric, cohort_label = NULL,
                                show_unavailable = FALSE) {
  sb_validate_show_unavailable(show_unavailable)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    sb_abort(
      "Plotting requires ggplot2; install it with install.packages('ggplot2').",
      "stormbreaker_missing_plot_dependency"
    )
  }
  tab <- sb_result_table(x, metric = metric, cohort_label = cohort_label)
  external <- inherits(x, c("stormbreaker_external_evaluation",
                             "stormbreaker_external_comparison"))
  tab$pipeline_id <- factor(tab$pipeline_id, levels = unique(tab$pipeline_id))
  horizon_metric <- metric %in% sb_horizon_view_metrics()
  composite <- "view_id" %in% names(tab)
  separate_labels <- horizon_metric || composite
  if (separate_labels) {
    tab$display_label <- if (horizon_metric) sb_horizon_view_labels(tab) else as.character(tab$pipeline_id)
    if (composite) {
      tab$display_label <- vapply(seq_len(nrow(tab)), function(i) {
        paste(tab$display_label[i], sb_view_label(as.list(tab[i, ])),
              paste("Saved entry:", tab$view_id[i]), sep = "\n")
      }, character(1))
    }
    tab$display_label <- factor(tab$display_label, levels = unique(tab$display_label))
  }
  drawable <- sb_performance_display_rows(tab)
  points <- tab[drawable, , drop = FALSE]
  omitted <- sum(!drawable)
  title <- if (external) "External performance" else "CV performance"
  subtitle <- paste(
    c(if (!is.null(cohort_label)) cohort_label else if (external)
      "External cohort (not named)",
      paste(unique(tab$scope), unique(tab$aggregation), sep = "\n", collapse = "\n")),
    collapse = "\n"
  )
  if (composite) {
    title <- "Performance by context"
    subtitle <- paste(subtitle, paste(unique(vapply(seq_len(nrow(tab)), function(i) {
      sb_view_boundary(as.list(tab[i, ]))
    }, character(1))), collapse = "\n"), "No pooling across repetitions or cohorts.", sep = "\n")
  }
  caption <- paste0(
    omitted, " pipeline(s) omitted from points (run/metric status or nonfinite estimate).\n",
    "Use sb_export_failures() for details. No confidence intervals."
  )
  if (horizon_metric) {
    subtitle <- paste(subtitle,
      "Saved horizon and unit shown per pipeline; different targets are not equivalent.",
      sep = "\n")
  }
  if (!show_unavailable && !any(drawable)) {
    p <- sb_empty_display_plot(title, subtitle, caption)
    p$data <- tab
    return(p)
  }
  metric_label <- switch(metric,
    roc_auc = "ROC AUC",
    average_precision = "Average precision",
    brier_score = "Brier score",
    log_loss = "Log loss",
    harrell_c = "Harrell C",
    time_auc = "Cumulative/dynamic time-AUC",
    uno_c = "Truncated Uno C",
    ipcw_brier = "IPCW Brier score"
  )
  pipeline_id <- estimate <- display_label <- NULL
  mapping <- if (separate_labels) ggplot2::aes(x = display_label, y = estimate) else
    ggplot2::aes(x = pipeline_id, y = estimate)
  labels <- if (separate_labels) tab$display_label else tab$pipeline_id
  visible_labels <- if (separate_labels) points$display_label else points$pipeline_id
  p <- ggplot2::ggplot(tab, mapping) +
    ggplot2::geom_point(data = points, size = 2.5, colour = "#24547A") +
    ggplot2::scale_x_discrete(
      limits = if (show_unavailable) rev(levels(labels)) else
        rev(unique(as.character(visible_labels))), drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      limits = if (metric %in% c("log_loss", "ipcw_brier")) c(0, NA_real_) else c(0, 1)
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "Model",
      y = metric_label,
      title = title, subtitle = subtitle, caption = caption
    ) +
    sb_viz_theme()
  sb_clinician_clean_plot(p, sb_clinician_models(sb_viz_leaves(x, cohort_label)))
}
