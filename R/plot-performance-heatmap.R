#' Compare saved performance across models and cohorts in a heatmap
#'
#' Display one explicitly selected metric without averaging, ranking,
#' reevaluating patients or changing the saved estimates.
#'
#' @param x A saved result/comparison, supported composite CV object, or a
#'   uniquely named list of these; names identify cohorts or result sets.
#' @param metric A scalar metric supported by [sb_result_table()].
#' @param cohort_label Optional cohort label for one object, not a named list.
#' @param digits Number of decimal places printed in cells (zero to six).
#' @param annotate Logical; print estimates or unavailable status in cells.
#'
#' @return A ggplot object. `data` retains original estimates and status,
#'   a separate `fill_value` used only for successful cells, and explicit
#'   `not_supplied` rows for missing model/cohort combinations. No patient
#'   rows are included. Optional ggplot2 is required.
#'
#' @details Rows retain model input order. Columns retain cohort/context input
#'   order, including repetition and stacking role. Different tasks, horizons,
#'   time units, metric scopes or aggregations are separated into panels.
#'   Cohort names do not prove independent validation. No automatic winner or
#'   clustered ordering is inferred. A cell needs successful metric and run
#'   status and a finite estimate to be coloured. Unavailable cells are grey,
#'   never zero; original values remain in the data. Colour scales use the
#'   observed successful estimates and retain their numerical values; lighter
#'   colours indicate preferable performance. No confidence interval is
#'   inferred. Large grids require a suitably sized output image.
#'
#' @export
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   d <- sb_example_binary(60, seed = 1)
#'   task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
#'   cv <- sb_run_cv(task, folds = 3, seed = 1)
#'   sb_plot_performance_heatmap(list(Development = cv), "roc_auc")
#' }
sb_plot_performance_heatmap <- function(x, metric, cohort_label = NULL,
                                        digits = 3L, annotate = TRUE) {
  sb_require_curve_plotting()
  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits != floor(digits) || digits < 0 || digits > 6) {
    sb_abort("`digits` must be an integer from zero to six.", "stormbreaker_bad_visualization_argument")
  }
  if (!is.logical(annotate) || length(annotate) != 1L || is.na(annotate)) {
    sb_abort("`annotate` must be TRUE or FALSE.", "stormbreaker_bad_visualization_argument")
  }
  leaves <- sb_viz_leaves(x, cohort_label)
  rows <- lapply(leaves, function(leaf) {
    tab <- sb_result_table(leaf$value, metric = metric)
    tab <- sb_view_annotate(tab, leaf$context)
    tab$dataset <- leaf$dataset
    tab$column_label <- sb_viz_panel(leaf)
    horizon <- if (metric %in% sb_horizon_view_metrics()) {
      paste(tab$prediction_horizon, tab$time_unit)
    } else "horizon-free target"
    tab$target_label <- paste(tab$task_type, metric, horizon, tab$scope, tab$aggregation, sep = " | ")
    tab$cell_status <- ifelse(tab$run_status != "success" | is.na(tab$run_status), "failed_run",
      ifelse(tab$status != "success" | is.na(tab$status), "unavailable_metric",
        ifelse(is.finite(tab$estimate), "success", "nonfinite_estimate")))
    tab
  })
  data <- sb_view_bind(rows)
  model_order <- unique(data$pipeline_id)
  column_order <- unique(data$column_label)
  target_order <- unique(data$target_label)
  key <- paste(match(data$pipeline_id, model_order), match(data$column_label, column_order),
    match(data$target_label, target_order), sep = "/")
  if (anyDuplicated(key)) sb_abort("Multiple saved entries occupy the same model/cohort/target cell; give them distinct result-set labels.",
    "stormbreaker_ambiguous_heatmap_cell")
  additions <- list()
  for (target in target_order) {
    local <- data[data$target_label == target, , drop = FALSE]
    grid <- expand.grid(pipeline_id = unique(local$pipeline_id), column_label = unique(local$column_label),
      stringsAsFactors = FALSE)
    present <- paste(match(local$pipeline_id, model_order), match(local$column_label, column_order), sep = "/")
    missing <- !paste(match(grid$pipeline_id, model_order), match(grid$column_label, column_order), sep = "/") %in% present
    if (!any(missing)) next
    absent <- grid[missing, , drop = FALSE]
    absent$target_label <- target
    absent$cell_status <- "not_supplied"
    absent$estimate <- NA_real_
    absent$metric <- metric
    additions[[length(additions) + 1L]] <- absent
  }
  if (length(additions)) data <- sb_view_bind(c(list(data), additions))
  data$fill_value <- ifelse(data$cell_status == "success", data$estimate, NA_real_)
  data$annotation <- ifelse(data$cell_status == "success", formatC(data$estimate, format = "f", digits = digits),
    "")
  labels <- c(failed_run = "failed", unavailable_metric = "unavailable", nonfinite_estimate = "nonfinite", not_supplied = "not supplied")
  unavailable <- data$cell_status != "success"
  data$annotation[unavailable] <- unname(labels[data$cell_status[unavailable]])
  data$pipeline_id <- factor(data$pipeline_id, levels = rev(model_order))
  data$column_label <- factor(data$column_label, levels = column_order)
  data$target_label <- factor(data$target_label, levels = target_order)
  lower_better <- metric %in% c("brier_score", "log_loss", "ipcw_brier")
  pipeline_id <- column_label <- fill_value <- annotation <- target_label <- NULL
  p <- ggplot2::ggplot(data, ggplot2::aes(x = column_label, y = pipeline_id, fill = fill_value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::scale_fill_viridis_c(option = "C", direction = if (lower_better) -1 else 1,
      na.value = "grey90", name = metric) +
    ggplot2::facet_wrap(~target_label, scales = "free", ncol = 1,
      labeller = ggplot2::as_labeller(sb_view_wrap_label)) +
    ggplot2::scale_x_discrete(labels = sb_view_wrap_label) +
    ggplot2::scale_y_discrete(labels = sb_view_wrap_label) +
    ggplot2::labs(title = "Saved model-by-cohort performance",
      subtitle = if (all(unavailable)) "No usable successful estimates; availability is shown." else
        paste(if (lower_better) "Lower" else "Higher", "is better; input order, no model selection."),
      x = NULL, y = "Pipeline",
      caption = "Separate targets and validation contexts. Grey cells are unavailable, not zero. No pooling or inferred intervals.") +
    sb_viz_theme() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 9))
  if (annotate) p <- p + ggplot2::geom_label(ggplot2::aes(label = annotation),
    size = 2.8, fill = "white", alpha = 0.85)
  p
}
