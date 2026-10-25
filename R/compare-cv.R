#' Compare admitted pipelines on the same outer folds
#'
#' `sb_compare_cv()` runs every supplied pipeline with the same survival or
#' binary task and outer-fold specification. A scalar count creates one shared
#' deterministic event- or outcome-stratified assignment. A multi-element
#' numeric vector named by every task patient ID supplies one exact assignment
#' to every pipeline; names, not vector position, bind patients to folds. The
#' package does not repair or rebalance a supplied plan, derive it from model
#' performance, or certify its design quality. The returned table copies each
#' run's fold-stratified Harrell C, or
#' its IPCW Brier score, fold-local Kaplan-Meier reference Brier, Brier skill
#' score, cumulative/dynamic time-dependent AUC, truncated Uno IPCW
#' concordance, pooled observed/expected mean calibration, calibration slope,
#' and moderate-calibration ICI/E50/E90 when the survival task declared a
#' prediction horizon, or binary
#' ROC AUC,
#' fold-stratified average precision, patient-level OOF Brier
#' score and log loss, their analysis-prevalence references, and their skill
#' scores, descriptive OOF calibration intercept and slope, a separate saved
#' binary moderate-calibration curve's ICI/E50/E90 summaries, an optional saved
#' task-appropriate decision curve, the saved fold-local selection-stability
#' status and mean pairwise Jaccard, and the separate glmnet nonzero active-set
#' stability audit. Pipelines without a learned selector are marked not
#' applicable for selector stability; Ridge and non-glmnet learners are marked
#' not applicable for active-set stability. The comparison does
#' not rank predictions from different folds against each other or recompute a
#' different metric. Weak or moderate calibration equations are not applied to
#' predictions. The binary moderate-calibration curve is not refitted here and
#' ICI/E50/E90 have no universal acceptable cutoff. Decision curves are copied
#' as complete per-pipeline objects; they are not flattened, integrated,
#' smoothed, ranked, or used to select a threshold.
#' Survival censoring distributions, reference survival, probabilities, losses,
#' within-fold time-AUC and Uno C pair quantities, both saved weak-calibration
#' summaries, the moderate-calibration curve, and its ICI/E50/E90 summaries are
#' likewise taken from the saved CV result and are never re-estimated by this
#' comparison.
#' The moderate curve is descriptive and is not applied to predictions. No
#' survival calibration intercept or universal calibration cutoff is reported.
#'
#' This is a descriptive comparison. It does not rank pipelines, select a
#' winner, estimate uncertainty, or turn development CV into external
#' validation. Both Jaccard summaries are descriptive: outer analysis folds
#' overlap, and no confidence interval or chance correction is inferred. A
#' not-evaluable glmnet active-set audit does not erase otherwise complete OOF
#' predictions or metrics. A
#' failed pipeline remains in the table with its complete audited CV result and
#' no replacement learner; later supplied pipelines are still run in their
#' original order. Task, fold, and pipeline-contract errors are validated
#' before this collection behavior applies.
#' Fixed LightGBM, CatBoost, and ranger pipelines are compared as their declared
#' configurations; this function does not reinterpret them as tuned algorithms.
#'
#' @param task A task returned by [sb_surv_task()] or [sb_binary_task()].
#' @param pipelines A non-empty list of distinct pipeline specifications
#'   returned by [sb_pipeline()]. Input order is preserved.
#' @param folds Either one positive whole-number shared outer-fold count or a
#'   multi-element numeric vector named by the complete unique task patient IDs.
#'   Supplied values must be consecutive positive whole-number labels
#'   `1, ..., K`, with `K >= 2`, and are preserved after exact ID alignment.
#' @param seed Positive integer used to generate scalar-count outer folds. It
#'   does not choose a supplied patient-ID assignment, whose effective outer
#'   seed is recorded as `NA`.
#'
#' @param reuse_selectors Reuse successful, warning-free fits of the seeded
#'   Boruta, varSelRF, survival Lasso, CoxBoost and RSF minimal-depth selectors,
#'   and deterministic StepCox, within this comparison when their complete
#'   analysis inputs and settings match exactly. Other selectors are unchanged.
#'   No fit is shared across different outer or inner analysis inputs.
#' @param progress Show model and outer-fold progress messages. Defaults to
#'   `interactive()`; set TRUE to include progress in script logs.
#' @return A `stormbreaker_cv_comparison` containing the task type, a
#'   descriptive `comparison` data frame, the complete named `results` from
#'   [sb_run_cv()], the shared fold count, effective seed, source, aligned
#'   patient-to-fold assignment, and, when declared, a named list of copied
#'   `decision_curves`.
#' @export
#'
#' @examples
#' d <- sb_example_survival(n = 50, seed = 1)
#' task <- sb_surv_task(
#'   d, "patient_id", "time", "event",
#'   c("gene_a", "gene_b"), "month"
#' )
#' compared <- sb_compare_cv(
#'   task,
#'   list(
#'     sb_pipeline(),
#'     sb_pipeline(
#'       selector = "univ_cox_top_k",
#'       selector_n_features = 1
#'     )
#'   ),
#'   folds = 5,
#'   seed = 2026
#' )
#' compared
#' outer_plan <- setNames(
#'   compared$fold_assignment$fold,
#'   compared$fold_assignment$patient_id
#' )
#' identical(unname(outer_plan), compared$fold_assignment$fold)
#'
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   binary_data <- sb_example_binary(n = 60, seed = 1)
#'   binary_task <- sb_binary_task(
#'     binary_data, "patient_id", "outcome",
#'     c("gene_a", "gene_b", "gene_c")
#'   )
#'   sb_compare_cv(
#'     binary_task,
#'     list(
#'       sb_pipeline(learner = "logistic"),
#'       sb_pipeline(learner = "gaussian_naive_bayes"),
#'       sb_pipeline(
#'         learner = "glmnet_logistic",
#'         inner_folds = 3,
#'         inner_seed = 1701
#'       )
#'     ),
#'     folds = 3,
#'     seed = 2026
#'   )
#' }
#'
#' if (requireNamespace("lightgbm", quietly = TRUE)) {
#'   binary_data <- sb_example_binary(n = 50, seed = 2)
#'   binary_task <- sb_binary_task(
#'     binary_data, "patient_id", "outcome",
#'     c("gene_a", "gene_b", "gene_c")
#'   )
#'   sb_compare_cv(
#'     binary_task,
#'     list(
#'       sb_pipeline(
#'         learner = "lightgbm_binary",
#'         lightgbm_nrounds = 10
#'       ),
#'       sb_pipeline(
#'         selector = "univ_binary_auc_top_k",
#'         learner = "lightgbm_binary",
#'         selector_n_features = 2,
#'         lightgbm_nrounds = 10
#'       )
#'     ),
#'     folds = 5,
#'     seed = 2026
#'   )
#' }
#'
#' if (requireNamespace("catboost", quietly = TRUE)) {
#'   binary_data <- sb_example_binary(n = 50, seed = 3)
#'   binary_task <- sb_binary_task(
#'     binary_data, "patient_id", "outcome",
#'     c("gene_a", "gene_b", "gene_c")
#'   )
#'   sb_compare_cv(
#'     binary_task,
#'     list(
#'       sb_pipeline(
#'         learner = "catboost_binary",
#'         catboost_iterations = 10
#'       ),
#'       sb_pipeline(
#'         selector = "univ_binary_auc_top_k",
#'         learner = "catboost_binary",
#'         selector_n_features = 2,
#'         catboost_iterations = 10
#'       )
#'     ),
#'     folds = 5,
#'     seed = 2026
#'   )
#' }
sb_compare_cv <- function(task, pipelines, folds = 5L, seed = 1L,
                          reuse_selectors = TRUE, progress = interactive()) {
  if (inherits(pipelines, "stormbreaker_stack")) {
    pipelines <- sb_upgrade_legacy_names(pipelines)
  } else if (is.list(pipelines)) pipelines <- lapply(pipelines, sb_upgrade_legacy_names)

  sb_check_progress(progress)
  task_type <- if (inherits(task, "stormbreaker_binary_task")) {
    "binary"
  } else if (inherits(task, "stormbreaker_surv_task")) {
    "survival"
  } else {
    sb_abort(
      "`task` must be created by `sb_surv_task()` or `sb_binary_task()`.",
      "stormbreaker_bad_task"
    )
  }
  valid_collection <- is.list(pipelines) && length(pipelines) > 0L &&
    all(vapply(
      pipelines,
      function(pipeline) inherits(pipeline, "stormbreaker_pipeline"),
      logical(1)
    ))
  if (!valid_collection) {
    sb_abort(
      "`pipelines` must be a non-empty list created with `sb_pipeline()`.",
      "stormbreaker_bad_pipeline_collection"
    )
  }

  pipeline_ids <- vapply(
    pipelines,
    function(pipeline) pipeline$id,
    character(1)
  )
  if (anyDuplicated(pipeline_ids) > 0L) {
    sb_abort(
      "Every compared pipeline must have a distinct pipeline ID.",
      "stormbreaker_duplicate_pipeline_id"
    )
  }
  matching_task <- vapply(
    pipelines,
    function(pipeline) identical(pipeline$task_type, task_type),
    logical(1)
  )
  if (!all(matching_task)) {
    sb_abort(
      "Every compared pipeline must match the supplied task type.",
      "stormbreaker_pipeline_task_mismatch"
    )
  }

  results <- sb_with_selector_reuse(lapply(
    seq_along(pipelines),
    function(i) {
      pipeline <- pipelines[[i]]
      sb_progress_message(progress, "CV model ", i, "/", length(pipelines), ": ",
        sb_display_selector(pipeline$selector), " + ",
        sb_display_learner(pipeline$learner, pipeline$learner_params))
      result <- sb_run_cv(task, pipeline, folds = folds, seed = seed, progress = progress)
      sb_progress_message(progress, "CV model ", i, "/", length(pipelines),
        " complete: ", result$status)
      result
    }
  ), enabled = reuse_selectors)
  names(results) <- pipeline_ids

  rows <- lapply(seq_along(results), function(index) {
    result <- results[[index]]
    pipeline <- pipelines[[index]]
    if (identical(task_type, "binary")) {
      roc_auc <- result$metrics$roc_auc
      average_precision <- result$metrics$average_precision
      brier <- result$metrics$brier_score
      log_loss <- result$metrics$log_loss
      calibration <- result$metrics$calibration
      moderate_calibration <- result$metrics$moderate_calibration
      selection_stability <- result$selection_stability
      active_set_stability <- result$glmnet_active_set_stability
      data.frame(
        pipeline_id = pipeline$id,
        selector = pipeline$selector,
        learner = pipeline$learner,
        status = result$status,
        complete_oof = result$complete_oof,
        roc_auc_status = roc_auc$status,
        roc_auc = roc_auc$estimate,
        comparable_pairs = roc_auc$comparable_pairs,
        average_precision_status = average_precision$status,
        average_precision = average_precision$estimate,
        brier_status = brier$status,
        brier_score = brier$estimate,
        n_scored = brier$n_scored,
        reference_brier_status = brier$reference$status,
        reference_brier = brier$reference$estimate,
        brier_skill_status = brier$skill_score$status,
        brier_skill_score = brier$skill_score$estimate,
        log_loss_status = log_loss$status,
        log_loss = log_loss$estimate,
        reference_log_loss_status = log_loss$reference$status,
        reference_log_loss = log_loss$reference$estimate,
        log_loss_skill_status = log_loss$skill_score$status,
        log_loss_skill_score = log_loss$skill_score$estimate,
        calibration_status = calibration$status,
        calibration_intercept = calibration$calibration_intercept,
        calibration_slope = calibration$calibration_slope,
        moderate_calibration_status = moderate_calibration$status,
        moderate_calibration_ici = moderate_calibration$ici,
        moderate_calibration_e50 = moderate_calibration$e50,
        moderate_calibration_e90 = moderate_calibration$e90,
        selection_stability_status = selection_stability$status,
        mean_pairwise_jaccard = selection_stability$mean_pairwise_jaccard,
        n_successful_selection_folds = selection_stability$n_successful_folds,
        n_folds = as.integer(nrow(result$runs)),
        n_successful_folds = as.integer(
          sum(result$runs$status == "success")
        ),
        glmnet_active_set_stability_status = active_set_stability$status,
        glmnet_active_set_mean_pairwise_jaccard =
          active_set_stability$mean_pairwise_jaccard,
        n_successful_glmnet_active_set_folds =
          active_set_stability$n_successful_folds,
        stringsAsFactors = FALSE
      )
    } else {
      metric <- result$metrics$harrell_c
      selection_stability <- result$selection_stability
      active_set_stability <- result$glmnet_active_set_stability
      row <- data.frame(
        pipeline_id = pipeline$id,
        selector = pipeline$selector,
        learner = pipeline$learner,
        status = result$status,
        complete_oof = result$complete_oof,
        harrell_c_status = metric$status,
        harrell_c = metric$estimate,
        comparable_pairs = metric$comparable_pairs,
        selection_stability_status = selection_stability$status,
        mean_pairwise_jaccard = selection_stability$mean_pairwise_jaccard,
        n_successful_selection_folds = selection_stability$n_successful_folds,
        n_folds = as.integer(nrow(result$runs)),
        n_successful_folds = as.integer(
          sum(result$runs$status == "success")
        ),
        glmnet_active_set_stability_status = active_set_stability$status,
        glmnet_active_set_mean_pairwise_jaccard =
          active_set_stability$mean_pairwise_jaccard,
        n_successful_glmnet_active_set_folds =
          active_set_stability$n_successful_folds,
        stringsAsFactors = FALSE
      )
      if (!is.null(result$metrics$ipcw_brier)) {
        ipcw_brier <- result$metrics$ipcw_brier
        weak_calibration <- result$metrics$weak_calibration
        row <- data.frame(
          pipeline_id = row$pipeline_id,
          selector = row$selector,
          learner = row$learner,
          status = row$status,
          complete_oof = row$complete_oof,
          complete_survival_probability_oof =
            result$complete_survival_probability_oof,
          prediction_horizon = ipcw_brier$prediction_horizon,
          harrell_c_status = row$harrell_c_status,
          harrell_c = row$harrell_c,
          comparable_pairs = row$comparable_pairs,
          ipcw_brier_status = ipcw_brier$status,
          ipcw_brier = ipcw_brier$estimate,
          reference_ipcw_brier_status = ipcw_brier$reference$status,
          reference_ipcw_brier = ipcw_brier$reference$estimate,
          ipcw_brier_skill_status = ipcw_brier$skill_score$status,
          ipcw_brier_skill_score = ipcw_brier$skill_score$estimate,
          time_auc_status = result$metrics$time_auc$status,
          time_auc = result$metrics$time_auc$estimate,
          time_auc_pair_denominator =
            result$metrics$time_auc$pair_denominator,
          uno_c_status = result$metrics$uno_c$status,
          uno_c = result$metrics$uno_c$estimate,
          uno_c_weighted_pair_denominator =
            result$metrics$uno_c$weighted_pair_denominator,
          calibration_status = result$metrics$calibration$status,
          calibration_ici = result$metrics$calibration$ici,
          calibration_e50 = result$metrics$calibration$e50,
          calibration_e90 = result$metrics$calibration$e90,
          mean_calibration_status = weak_calibration$mean_calibration$status,
          mean_calibration_ratio =
            weak_calibration$mean_calibration_ratio,
          calibration_slope_status = weak_calibration$slope_fit$status,
          calibration_slope = weak_calibration$calibration_slope,
          selection_stability_status = row$selection_stability_status,
          mean_pairwise_jaccard = row$mean_pairwise_jaccard,
          n_successful_selection_folds =
            row$n_successful_selection_folds,
          n_folds = row$n_folds,
          n_successful_folds = row$n_successful_folds,
          glmnet_active_set_stability_status =
            row$glmnet_active_set_stability_status,
          glmnet_active_set_mean_pairwise_jaccard =
            row$glmnet_active_set_mean_pairwise_jaccard,
          n_successful_glmnet_active_set_folds =
            row$n_successful_glmnet_active_set_folds,
          stringsAsFactors = FALSE
        )
      }
      row
    }
  })
  comparison <- do.call(rbind, rows)
  rownames(comparison) <- NULL
  comparison$prediction_variation <- vapply(results, function(result) {
    sb_prediction_variation(result$predictions, task_type)$status
  }, character(1))
  metric_status <- if (identical(task_type, "binary")) {
    comparison$roc_auc_status == "success" &
      comparison$average_precision_status == "success" &
      comparison$brier_status == "success" &
      comparison$reference_brier_status == "success" &
      comparison$brier_skill_status == "success" &
      comparison$log_loss_status == "success" &
      comparison$reference_log_loss_status == "success" &
      comparison$log_loss_skill_status == "success" &
      comparison$calibration_status == "success" &
      comparison$selection_stability_status %in% c(
        "success", "not_applicable"
      )
  } else {
    survival_metric_status <- comparison$harrell_c_status == "success" &
      comparison$selection_stability_status %in% c(
        "success", "not_applicable"
      )
    if ("ipcw_brier_status" %in% names(comparison)) {
      survival_metric_status <- survival_metric_status &
        comparison$ipcw_brier_status == "success" &
        comparison$reference_ipcw_brier_status == "success" &
        comparison$ipcw_brier_skill_status == "success" &
        comparison$time_auc_status == "success" &
        comparison$uno_c_status == "success"
    }
    survival_metric_status
  }
  complete <- all(comparison$status == "success") &&
    all(metric_status)

  compared <- list(
    status = if (complete) "success" else "incomplete",
    task_type = task_type,
    comparison = comparison,
    results = results,
    n_patients = results[[1L]]$task$n,
    folds = results[[1L]]$provenance$folds,
    seed = results[[1L]]$provenance$seed,
    fold_source = results[[1L]]$provenance$fold_source,
    fold_assignment = results[[1L]]$provenance$fold_assignment
  )
  if ("decision_thresholds" %in% names(task)) {
    decision_curves <- lapply(
      results,
      function(result) result$metrics$decision_curve
    )
    names(decision_curves) <- names(results)
    compared$decision_curves <- decision_curves
  }
  structure(compared, class = "stormbreaker_cv_comparison")
}

#' @export
print.stormbreaker_cv_comparison <- function(x, ...) {
  cat("<stormbreaker_cv_comparison:", x$status, ">\n")
  cat(
    "Pipelines:", nrow(x$comparison),
    "| Patients:", x$n_patients,
    "| Shared outer folds:", x$folds, "\n"
  )
  fold_label <- if (identical(x$fold_source, "preassigned_patient_id")) {
    "preassigned patient-ID folds"
  } else {
    "generated stratified folds"
  }
  cat("Outer resampling:", fold_label, "\n")
  metric_columns <- if (identical(x$task_type, "binary")) {
    c(
      "roc_auc", "average_precision", "brier_score", "brier_skill_score",
      "log_loss", "log_loss_skill_score",
      "calibration_intercept", "calibration_slope",
      "moderate_calibration_ici", "moderate_calibration_e50",
      "moderate_calibration_e90",
      "comparable_pairs", "mean_pairwise_jaccard",
      "glmnet_active_set_mean_pairwise_jaccard"
    )
  } else {
    survival_columns <- c(
      "harrell_c", "comparable_pairs", "mean_pairwise_jaccard",
      "glmnet_active_set_mean_pairwise_jaccard"
    )
    if ("ipcw_brier" %in% names(x$comparison)) {
      survival_columns <- c(
        "harrell_c", "ipcw_brier", "reference_ipcw_brier",
        "ipcw_brier_skill_score", "time_auc", "uno_c", "calibration_ici",
        "calibration_e50", "calibration_e90", "mean_calibration_ratio",
        "calibration_slope",
        "time_auc_pair_denominator", "uno_c_weighted_pair_denominator",
        "comparable_pairs",
        "mean_pairwise_jaccard",
        "glmnet_active_set_mean_pairwise_jaccard"
      )
    }
    survival_columns
  }
  print(
    x$comparison[c(
      "pipeline_id", "status", metric_columns
    )],
    row.names = FALSE,
    digits = 3
  )
  cat("No pipeline is selected by this comparison.\n")
  invisible(x)
}
