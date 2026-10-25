#' Plot ROC and precision-recall curves from saved binary predictions
#'
#' Compute descriptive threshold coordinates without fitting models or
#' replacing saved performance metrics. ROC defaults to one pooled out-of-fold
#' curve per model/context; precision-recall retains separate assessment folds.
#'
#' @param x A binary CV or external evaluation/comparison, a supported
#'   repeated/IECV/stacking CV object, or a uniquely named list of these.
#' @param cohort_label Optional cohort display label for a single object.
#'   For a named list, labels come from the list names.
#' @param ncol Positive number of panel columns. Defaults to two.
#' @param by_fold For ROC, show separate assessment folds instead of a single
#'   pooled out-of-fold curve. Repetitions and external cohorts remain separate.
#' @param annotate Add displayed-curve AUC labels to ROC panels.
#' @param ci_reps Optional patient-bootstrap replicates for ROC AUC intervals;
#'   0 omits intervals. Use 1000 for automatic 95 percent interval labels.
#'   Requires independent patient rows, not repeated observations per patient.
#' @param seed Random seed for ROC intervals.
#'
#' @return A ggplot object. Its data contains aggregate threshold coordinates,
#'   not patient IDs. Attributes `stormbreaker_curve_status` and
#'   `stormbreaker_saved_metrics` retain availability and unchanged saved
#'   metrics for each entry. Optional ggplot2 is required.
#'
#' @details Both functions treat class 1 as positive and higher probabilities
#'   as more positive. Equal probabilities enter together. ROC joins the
#'   resulting points linearly; PR uses vertical-then-horizontal steps so
#'   areas correspond to noninterpolated average precision within a panel.
#'   The starting PR precision at recall zero is defined as one.
#'
#'   Repetitions, held-out IECV cohorts and stacking roles remain separate.
#'   External evaluation has no fictitious folds. Saved aggregate AUC/AP objects
#'   are available as an attribute and are not replaced by displayed-curve AUC.
#'   Pooled OOF is not apparent training performance. Bootstrap intervals
#'   condition on saved predictions; they do not include development uncertainty.
#'   There is no threshold selection or confidence band. A failed run or a
#'   group with invalid/incomplete predictions is not
#'   drawn; no invalid rows are silently removed. ROC needs both classes;
#'   PR needs at least one positive. Unavailable groups remain in diagnostics.
#'   Survival time-dependent ROC is not this binary interface.
#'
#' @export
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   d <- sb_example_binary(60, seed = 1)
#'   task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
#'   cv <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1)
#'   sb_plot_roc(cv)
#'   sb_plot_pr(cv)
#' }
sb_plot_roc <- function(x, cohort_label = NULL, ncol = 2L,
                        by_fold = FALSE, annotate = TRUE, ci_reps = 0L, seed = 1L) {
  sb_plot_binary_discrimination(x, "roc", cohort_label, ncol,
    by_fold, annotate, ci_reps, seed)
}

#' @rdname sb_plot_roc
#' @export
sb_plot_pr <- function(x, cohort_label = NULL, ncol = 2L) {
  sb_plot_binary_discrimination(x, "pr", cohort_label, ncol)
}

sb_binary_threshold_coordinates <- function(outcome, probability, kind) {
  ord <- order(probability, decreasing = TRUE)
  prob <- probability[ord]
  y <- outcome[ord]
  end <- c(which(diff(prob) != 0), length(prob))
  tp <- c(0, cumsum(y == 1)[end])
  fp <- c(0, cumsum(y == 0)[end])
  positive <- sum(y == 1)
  negative <- sum(y == 0)
  data.frame(threshold = c(Inf, prob[end]), true_positive = tp, false_positive = fp,
    x = if (kind == "roc") fp / negative else tp / positive,
    y = if (kind == "roc") tp / positive else c(1, utils::tail(tp, -1)/(utils::tail(tp, -1) + utils::tail(fp, -1))))
}

sb_plot_binary_discrimination <- function(x, kind, cohort_label, ncol,
                                         by_fold = TRUE, annotate = FALSE,
                                         ci_reps = 0L, seed = 1L) {
  sb_require_curve_plotting()
  ncol <- sb_viz_positive_integer(ncol, "ncol")
  leaves <- sb_viz_leaves(x, cohort_label)
  coordinates <- diagnostics <- saved <- annotations <- list()
  if (!is.logical(by_fold) || length(by_fold) != 1L || is.na(by_fold) ||
      !is.logical(annotate) || length(annotate) != 1L || is.na(annotate)) stop("by_fold and annotate must be TRUE or FALSE.")
  if (!is.numeric(ci_reps) || length(ci_reps) != 1L || !is.finite(ci_reps) ||
      ci_reps < 0 || ci_reps != as.integer(ci_reps) || ci_reps == 1) stop("ci_reps must be 0 or an integer of at least 2.")
  withr::local_seed(seed)
  for (leaf in leaves) {
    value <- leaf$value
    external <- inherits(value, "stormbreaker_external_evaluation")
    binary <- if (external) identical(value$task_type, "binary") else
      inherits(value, "stormbreaker_binary_cv_result")
    if (!binary) sb_abort("ROC/PR requires binary saved predictions, not survival scores.",
      "stormbreaker_bad_visualization_task")
    id <- sb_viz_pipeline_id(value)
    pred <- if (external) value$scored_predictions else value$predictions
    metric <- if (kind == "roc") "roc_auc" else "average_precision"
    saved[[leaf$context$view_id]] <- value$metrics[[metric]]
    groups <- if (by_fold && !external && is.data.frame(pred) && nrow(pred) && "fold" %in% names(pred)) {
      unique(pred$fold)
    } else NA_integer_
    for (fold in groups) {
      group <- if (!is.data.frame(pred)) data.frame() else if (external || is.na(fold)) pred else
        pred[!is.na(pred$fold) & pred$fold == fold, , drop = FALSE]
      panel <- paste(sb_viz_panel(leaf), if (external) "External patients" else
        if (by_fold) paste("Assessment fold", fold) else "Out-of-fold", sep = "\n")
      status <- if (!identical(value$status, "success")) "failed_run" else if (!nrow(group)) "no_predictions" else "success"
      if (status == "success") {
        columns <- c("outcome", "probability", "status", "patient_id")
        valid <- all(columns %in% names(group)) && is.numeric(group$outcome) &&
          is.numeric(group$probability) && all(is.finite(group$outcome)) &&
          all(group$outcome %in% c(0, 1)) && all(is.finite(group$probability)) &&
          all(group$probability >= 0 & group$probability <= 1) &&
          !anyNA(group$status) && all(group$status == "success") &&
          !anyNA(group$patient_id) && !anyDuplicated(group$patient_id) &&
          (external || ("fold" %in% names(group) && !anyNA(group$fold)))
        if (!valid) status <- "invalid_predictions" else if (!sum(group$outcome == 1) ||
            (kind == "roc" && !sum(group$outcome == 0))) status <- "single_class"
      }
      diagnostic <- data.frame(pipeline_id = id, dataset = leaf$dataset, fold = as.integer(fold),
        panel = panel, curve_status = status, run_status = value$status,
        n = nrow(group), stringsAsFactors = FALSE)
      diagnostic <- sb_view_annotate(diagnostic, leaf$context)
      diagnostics[[length(diagnostics) + 1L]] <- diagnostic
      if (status != "success") next
      tab <- sb_binary_threshold_coordinates(group$outcome, group$probability, kind)
      for (name in names(diagnostic)) tab[[name]] <- rep(diagnostic[[name]], nrow(tab))
      tab$series_id <- paste(leaf$context$view_id, fold, sep = "/")
      coordinates[[length(coordinates) + 1L]] <- tab
      if (annotate && kind == "roc") {
        a <- sb_roc_interval(group$outcome, group$probability, ci_reps)
        a$pipeline_id <- id
        a$panel <- panel
        annotations[[length(annotations) + 1L]] <- a
      }
    }
  }
  status <- sb_view_bind(diagnostics)
  caption <- paste0(if (by_fold) "Separate assessment folds/external cohorts; not pooled. " else
    "One pooled out-of-fold curve per model and resampling context; not apparent training performance. ", "Saved aggregate ",
    if (kind == "roc") "AUC" else "AP", " is unchanged. ",
    sum(status$curve_status != "success"), " unavailable curve group(s). ",
    "AUC labels describe the displayed curve. Patient-bootstrap intervals condition on saved predictions, not development uncertainty.")
  title <- if (kind == "roc") "ROC" else "Precision-recall"
  if (!length(coordinates)) {
    p <- sb_empty_display_plot(title, "No usable saved binary predictions for these curves.", caption)
    p$data <- data.frame()
  } else {
    data <- sb_view_bind(coordinates)
    data$panel <- factor(data$panel, levels = unique(status$panel[status$curve_status == "success"]))
    data$pipeline_id <- factor(data$pipeline_id, levels = unique(status$pipeline_id))
    x <- y <- pipeline_id <- series_id <- panel <- NULL
    p <- ggplot2::ggplot(data, ggplot2::aes(x = x, y = y, colour = pipeline_id, group = series_id))
    p <- p + if (kind == "roc") ggplot2::geom_path(linewidth = 0.65) else
      ggplot2::geom_step(direction = "vh", linewidth = 0.65)
    if (kind == "roc") p <- p + ggplot2::geom_abline(slope = 1, intercept = 0,
      colour = "grey65", linetype = 2, inherit.aes = FALSE)
    p <- p + ggplot2::facet_wrap(~panel, ncol = ncol, labeller = ggplot2::as_labeller(sb_view_wrap_label)) +
      ggplot2::scale_x_continuous(limits = c(0, 1)) + ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = title, x = if (kind == "roc") "False positive rate" else "Recall",
        y = if (kind == "roc") "Sensitivity" else "Precision", colour = "Pipeline", caption = caption) +
      sb_viz_theme()
    if (length(annotations)) {
      a <- sb_view_bind(annotations)
      a$label <- paste0("AUC ", sprintf("%.3f", a$auc),
        if (ci_reps > 0) paste0(" (95% CI ", sprintf("%.3f", a$lower), "\u2013", sprintf("%.3f", a$upper), ")") else "")
      a$label_y <- .04 + .055 * (stats::ave(seq_len(nrow(a)), a$panel, FUN = seq_along) - 1)
      label <- label_y <- NULL
      p <- p + ggplot2::geom_text(data = a,
        ggplot2::aes(x = .98, y = label_y, label = label, colour = pipeline_id),
        inherit.aes = FALSE, show.legend = FALSE, hjust = 1, size = 3)
      attr(p, "stormbreaker_roc_statistics") <- a
    }
  }
  attr(p, "stormbreaker_curve_status") <- status
  attr(p, "stormbreaker_saved_metrics") <- saved
  if (kind == "roc" && nrow(p$data)) {
    p <- sb_clinician_clean_plot(p, sb_clinician_models(leaves))
    attr(p, "stormbreaker_figure_notes") <- caption
    p <- p + ggplot2::labs(subtitle = NULL, caption = NULL)
  }
  p
}
