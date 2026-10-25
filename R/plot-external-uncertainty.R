#' Plot saved external bootstrap confidence intervals
#'
#' Read existing uncertainty results without bootstrap resampling, metric
#' recomputation or model selection.
#'
#' @param x An external bootstrap object, or a uniquely named list of such
#'   objects. For `sb_plot_external_difference`, use paired external bootstrap
#'   comparison objects instead. List names identify cohorts/result sets.
#' @param metric One saved scalar metric, explicitly selected by name.
#' @param cohort_label Optional label when supplying a single object.
#'
#' @return A ggplot object whose data preserves the saved interval endpoints,
#'   availability, confidence level and successful-replicate denominator.
#'
#' @details Intervals describe external-patient sampling conditional on locked
#'   predictions, not development-fitting or model-selection uncertainty.
#'   Partial bootstrap intervals are shown only when the saved object marks
#'   them available, with their partial denominator explicitly labelled.
#'   Missing intervals are never reconstructed from a point estimate. Different
#'   task/metric targets are placed in separate panels. Numerical values and
#'   intervals are not clipped to theoretical metric bounds.
#'
#'   Differences retain the saved candidate-minus-reference direction. A
#'   positive difference is favourable only for a higher-is-better metric.
#'   Failed candidates remain visible without invented points. These scalar
#'   interfaces do not flatten nested decision-curve uncertainty results.
#'
#' @export
sb_plot_external_intervals <- function(x, metric, cohort_label = NULL) {
  sb_plot_bootstrap_intervals(x, metric, cohort_label, paired = FALSE)
}

#' @rdname sb_plot_external_intervals
#' @export
sb_plot_external_difference <- function(x, metric, cohort_label = NULL) {
  sb_plot_bootstrap_intervals(x, metric, cohort_label, paired = TRUE)
}

sb_viz_uncertainty_inputs <- function(x, classes, cohort_label) {
  if (inherits(x, classes)) {
    label <- sb_result_cohort_label(cohort_label)
    return(stats::setNames(list(x), if (is.na(label)) "Unnamed cohort" else label))
  }
  if (!is.null(cohort_label) || !is.list(x) || !length(x) || is.null(names(x)) ||
      anyNA(names(x)) || any(!nzchar(trimws(names(x)))) || anyDuplicated(names(x)) ||
      !all(vapply(x, inherits, logical(1), what = classes))) {
    sb_abort("Supply a supported uncertainty object or a uniquely named list of them.",
      "stormbreaker_bad_visualization_input")
  }
  x
}

sb_viz_uncertainty_target <- function(x, metric) {
  target <- x$metric_target
  if (is.null(target)) target <- x$prediction_target
  horizon_metric <- !metric %in% c("roc_auc", "average_precision", "brier_score", "log_loss", "harrell_c") &&
    identical(x$task_type, "survival")
  label <- if (!horizon_metric) "horizon-free target" else if (isTRUE(target$has_horizon)) {
    paste(target$prediction_horizon, target$time_unit)
  } else "horizon not saved"
  paste(x$task_type, metric, label, sep = " | ")
}

sb_viz_empty_interval <- function(metric, status = "interval_not_saved") {
  data.frame(metric = metric, point_estimate = NA_real_, lower = NA_real_, upper = NA_real_,
    confidence_level = NA_real_, status = status, interval_available = FALSE,
    replicate_completeness = "not_saved", n_requested = NA_integer_,
    n_finite_successful = NA_integer_, stringsAsFactors = FALSE)
}

sb_viz_interval_label <- function(tab) {
  level <- ifelse(is.finite(tab$confidence_level), paste0(format(100 * tab$confidence_level, trim = TRUE), "% CI"), "CI not saved")
  counts <- ifelse(is.finite(tab$n_requested) & is.finite(tab$n_finite_successful),
    paste0(tab$n_finite_successful, "/", tab$n_requested, " draws; ", tab$replicate_completeness),
    "replicate count not saved")
  paste(level, counts, tab$status, sep = " | ")
}

sb_plot_bootstrap_intervals <- function(x, metric, cohort_label, paired) {
  sb_require_curve_plotting()
  if (!is.character(metric) || length(metric) != 1L || is.na(metric) || !nzchar(metric)) {
    sb_abort("Explicitly select one saved metric.", "stormbreaker_bad_result_metric")
  }
  classes <- if (paired) "stormbreaker_external_bootstrap_comparison" else "stormbreaker_external_bootstrap"
  values <- sb_viz_uncertainty_inputs(x, classes, cohort_label)
  rows <- lapply(seq_along(values), function(i) {
    value <- values[[i]]
    all <- value$intervals
    tab <- if (is.data.frame(all) && "metric" %in% names(all)) all[!is.na(all$metric) & all$metric == metric, , drop = FALSE] else NULL
    if (!is.null(tab) && nrow(tab) > 0L) {
      required <- c("point_estimate", "lower", "upper", "confidence_level", "status", "interval_available",
        "replicate_completeness", "n_requested", "n_finite_successful")
      if (!all(required %in% names(tab))) sb_abort("Saved interval fields are incomplete.", "stormbreaker_bad_visualization_input")
    }
    if (paired) {
      candidates <- value$candidate_status
      ids <- candidates$candidate_pipeline_id
      if (!length(ids) || anyNA(ids) || anyDuplicated(ids)) {
        sb_abort("A paired comparison needs distinct saved candidate IDs.", "stormbreaker_bad_visualization_input")
      }
      pieces <- lapply(seq_along(ids), function(j) {
        row <- if (!is.null(tab) && nrow(tab)) tab[tab$candidate_pipeline_id == ids[[j]], , drop = FALSE] else NULL
        if (is.null(row) || !nrow(row)) {
          row <- sb_viz_empty_interval(metric, if (candidates$status[[j]] == "failed") "failed" else "interval_not_saved")
          row$candidate_pipeline_id <- ids[[j]]
        }
        if (nrow(row) != 1L) sb_abort("Duplicate candidate/metric interval rows.", "stormbreaker_bad_visualization_input")
        row$reference_pipeline_id <- value$reference_pipeline_id
        row$pipeline_id <- ids[[j]]
        row
      })
      tab <- sb_view_bind(pieces)
    } else {
      if (is.null(tab) || !nrow(tab)) tab <- sb_viz_empty_interval(metric,
        if (identical(value$status, "failed")) "failed" else "interval_not_saved")
      if (nrow(tab) != 1L) sb_abort("Duplicate model/metric interval rows.", "stormbreaker_bad_visualization_input")
      tab$pipeline_id <- value$pipeline_id
    }
    tab$dataset <- names(values)[[i]]
    tab$target_label <- sb_viz_uncertainty_target(value, metric)
    tab$interval_type <- "External patient-bootstrap CI"
    tab$row_label <- paste(tab$dataset, tab$pipeline_id,
      if (paired) paste("Reference:", value$reference_pipeline_id) else "Locked model",
      sb_viz_interval_label(tab), sep = "\n")
    tab
  })
  data <- sb_view_bind(rows)
  sb_viz_interval_plot(data,
    title = if (paired) "Paired external performance differences" else "External performance uncertainty",
    xlab = if (paired) paste0(metric, " (candidate - reference)") else metric,
    caption = paste("Conditional on locked predictions; no development-fitting uncertainty or model selection.",
      if (paired) "Same patient draws within each comparison; positive favours candidates only for higher-is-better metrics." else
        "Saved bootstrap intervals and replicate completeness; no new resampling."),
    reference = if (paired) 0 else NULL)
}

sb_viz_interval_plot <- function(data, title, xlab, caption, reference = NULL) {
  data$point_available <- is.finite(data$point_estimate) & !data$status %in% c("failed", "failed_run")
  data$draw_interval <- !is.na(data$interval_available) & data$interval_available &
    is.finite(data$lower) & is.finite(data$upper) & data$lower <= data$upper & data$point_available
  data$row_label <- factor(data$row_label, levels = rev(unique(data$row_label)))
  data$target_label <- factor(data$target_label, levels = unique(data$target_label))
  lower <- upper <- row_label <- point_estimate <- interval_type <- target_label <- NULL
  p <- ggplot2::ggplot(data, ggplot2::aes(y = row_label)) +
    ggplot2::geom_segment(data = data[data$draw_interval, , drop = FALSE],
      ggplot2::aes(x = lower, xend = upper, yend = row_label, linetype = interval_type),
      colour = "#24547A", linewidth = 0.7) +
    ggplot2::geom_point(data = data[data$point_available, , drop = FALSE],
      ggplot2::aes(x = point_estimate, shape = interval_type), colour = "#24547A", size = 2.4)
  if (!is.null(reference)) p <- p + ggplot2::geom_vline(xintercept = reference, linetype = 2, colour = "grey60")
  # Train row labels from all saved rows without adding estimates or marks.
  p + ggplot2::geom_blank() +
    ggplot2::facet_wrap(~target_label, ncol = 1, scales = "free",
      labeller = ggplot2::as_labeller(sb_view_wrap_label)) +
    ggplot2::scale_y_discrete(labels = sb_view_wrap_label) +
    ggplot2::labs(title = title, subtitle = if (!any(data$point_available)) "No usable point estimates; unavailable rows retained." else NULL,
      x = xlab, y = NULL, linetype = "Interval", shape = "Interval", caption = sb_viz_wrap(caption)) + sb_viz_theme()
}

#' Plot saved external Meta-analysis results as a forest plot
#'
#' @param x A saved external discrimination, mean-calibration,
#'   calibration-slope, discrimination-difference or prediction-error-difference
#'   Meta-analysis object.
#' @param show_prediction_interval Logical; show the saved between-cohort
#'   prediction interval on a separate row, distinct from the pooled CI.
#'
#' @return A ggplot object with cohort points, saved cohort bootstrap CIs when
#'   present, the saved pooled CI and optionally the saved prediction interval.
#'
#' @details No Meta model is refitted and no cohort CI is inferred from its
#'   standard error. Cohort interval levels are those in their source bootstrap
#'   objects and can differ from the pooled level; each row is labelled.
#'   Missing source intervals are explicit while cohort point estimates remain
#'   visible. Pooled and prediction endpoints are copied unchanged, including
#'   values beyond usual metric bounds. A prediction interval describes
#'   between-cohort variation, not an individual patient's prediction.
#'   The original small-cohort caution remains in the caption. Differences
#'   retain candidate-minus-reference direction; calibration target lines are
#'   one, and difference target lines are zero.
#'
#' @export
sb_plot_meta_forest <- function(x, show_prediction_interval = TRUE) {
  sb_require_curve_plotting()
  classes <- c("stormbreaker_external_discrimination_meta", "stormbreaker_external_mean_calibration_meta",
    "stormbreaker_external_calibration_slope_meta", "stormbreaker_external_discrimination_difference_meta",
    "stormbreaker_external_prediction_error_difference_meta")
  if (!inherits(x, classes) || !identical(x$status, "success")) {
    sb_abort("Supply a successful saved external Meta-analysis object.", "stormbreaker_bad_visualization_input")
  }
  if (!is.logical(show_prediction_interval) || length(show_prediction_interval) != 1L || is.na(show_prediction_interval)) {
    sb_abort("`show_prediction_interval` must be TRUE or FALSE.", "stormbreaker_bad_visualization_argument")
  }
  paired <- inherits(x, classes[4:5])
  cohorts <- x$cohorts
  sources <- if (paired) x$comparisons else x$bootstraps
  if (!is.data.frame(cohorts) || !nrow(cohorts) || anyDuplicated(cohorts$cohort_id)) {
    sb_abort("Saved cohort rows are missing or duplicated.", "stormbreaker_bad_visualization_input")
  }
  rows <- lapply(seq_len(nrow(cohorts)), function(i) {
    id <- cohorts$cohort_id[[i]]
    source <- sources[[id]]
    intervals <- source$intervals
    row <- NULL
    if (is.data.frame(intervals) && nrow(intervals)) {
      take <- !is.na(intervals$metric) & intervals$metric == x$metric
      if (paired) take <- take & intervals$candidate_pipeline_id == x$candidate_pipeline_id
      row <- intervals[take, , drop = FALSE]
    }
    if (is.null(row) || !nrow(row)) row <- sb_viz_empty_interval(x$metric)
    if (nrow(row) != 1L) sb_abort("Multiple saved intervals for one cohort and target.", "stormbreaker_bad_visualization_input")
    row$point_estimate <- if (paired) cohorts$point_difference[[i]] else cohorts$point_estimate[[i]]
    row$cohort_id <- id
    row$interval_type <- "Cohort bootstrap CI"
    row$row_label <- paste(id, sb_viz_interval_label(row), sep = "\n")
    row
  })
  summary <- x$summary
  if (!is.data.frame(summary) || nrow(summary) != 1L) sb_abort("One saved pooled summary is required.", "stormbreaker_bad_visualization_input")
  pooled <- sb_viz_empty_interval(x$metric, "success")
  pooled$point_estimate <- summary$estimate
  pooled$lower <- summary$lower
  pooled$upper <- summary$upper
  pooled$confidence_level <- x$confidence_level
  pooled$interval_available <- all(is.finite(c(pooled$lower, pooled$upper)))
  pooled$interval_type <- "Pooled CI"
  pooled$row_label <- paste0("Pooled estimate\n", 100 * x$confidence_level, "% CI")
  rows[[length(rows) + 1L]] <- pooled
  if (show_prediction_interval) {
    pi <- pooled
    pi$lower <- summary$prediction_lower
    pi$upper <- summary$prediction_upper
    pi$interval_available <- all(is.finite(c(pi$lower, pi$upper)))
    pi$interval_type <- "Prediction interval"
    pi$row_label <- paste0("Between-cohort prediction interval\n", 100 * x$confidence_level, "% PI")
    rows[[length(rows) + 1L]] <- pi
  }
  data <- sb_view_bind(rows)
  data$target_label <- sb_viz_uncertainty_target(x, x$metric)
  xlab <- if (paired) paste0(x$metric, " (candidate - reference)") else x$metric
  reference <- if (paired) 0 else if (x$metric %in% c("mean_calibration_ratio", "calibration_slope")) 1 else NULL
  sb_viz_interval_plot(data, "External validation forest plot", xlab,
    paste("Saved cohort bootstrap CIs, pooled CI and between-cohort PI; no refitting.",
      x$prediction_interval_caution), reference)
}
