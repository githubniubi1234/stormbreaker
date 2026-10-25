# STORMBREAKER 0.3.0 visualization tutorial. Synthetic data, not clinical evidence.
# Requires stormbreaker and optional ggplot2 and metafor (for Meta-analysis).
# Source this file, then call stormbreaker_visualization_example().
# Pass a NEW output directory to export all figures and an offline gallery.
stormbreaker_visualization_example <- function(output_dir = NULL, reps = 40L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 to draw the gallery.")
  if (!requireNamespace("metafor", quietly = TRUE)) stop("Install metafor to run the five Meta-analysis examples.")
  library(stormbreaker)
  features <- c("gene_a", "gene_b", "gene_c")
  development <- bi_example(120, seed = 301)
  task <- bi_task(development, "patient_id", "outcome", features)
  pipelines <- list(Logistic = sb_pipeline(learner = "logistic"),
    Top2 = sb_pipeline(selector = "univ_binary_auc_top_k", learner = "logistic",
      selector_n_features = 2L))
  cv <- compare_cv(task, pipelines, folds = 3L, seed = 302)
  # No winning model is chosen from the comparison; both are predeclared.
  fits <- lapply(pipelines, function(p) fit_final(task, p))
  cohorts <- lapply(1:3, function(i) {
    d <- bi_example(100, seed = 310 + i)
    d$patient_id <- paste0("external_", i, "_", seq_len(nrow(d)))
    d
  })
  names(cohorts) <- c("Cohort A", "Cohort B", "Cohort C")
  evaluations <- lapply(cohorts, function(d) {
    predictions <- lapply(fits, function(fit)
      predict_final(fit, d[c("patient_id", features)], "patient_id"))
    compare_external(predictions, bi_endpoint(d, "patient_id", "outcome"))
  })
  bootstraps <- lapply(evaluations, function(e)
    sb_bootstrap_external(e$evaluations[[1]], reps = reps, seed = 320))
  paired <- lapply(evaluations, function(e)
    sb_bootstrap_external_comparison(e, pipelines[[1]]$id, reps = reps, seed = 321))
  meta <- list(
    Discrimination = sb_meta_external_discrimination(bootstraps),
    Calibration_slope = sb_meta_external_calibration_slope(bootstraps),
    Mean_calibration = sb_meta_external_mean_calibration(bootstraps),
    Discrimination_difference = sb_meta_external_discrimination_difference(paired, pipelines[[2]]$id),
    Prediction_error_difference = sb_meta_external_prediction_error_difference(paired, pipelines[[2]]$id)
  )
  explanation <- additive_shap(fits[[1]], cohorts[[1]][1:30, ],
    "patient_id", development[features])
  plots <- list(
    ROC_by_assessment_fold = plot_roc(cv),
    PR_by_assessment_fold = plot_pr(cv),
    External_AUC_heatmap = plot_heatmap(evaluations, "roc_auc"),
    External_AUC_intervals = plot_ext_ci(bootstraps, "roc_auc"),
    Paired_AUC_differences = plot_ext_diff(paired, "roc_auc"),
    Selection_frequency = plot_freq(cv$results[[2]]),
    Selection_Jaccard = plot_jaccard(cv$results[[2]]),
    Selection_intersections = plot_upset(cv$results[[2]]),
    SHAP_importance = plot_shap_imp(explanation),
    SHAP_beeswarm = plot_shap(explanation,
      feature_values = cohorts[[1]][1:30, ], id_col = "patient_id"),
    SHAP_waterfall = plot_shap_one(explanation, top_n = 2L)
  )
  for (name in names(meta)) plots[[paste0("Meta_", name)]] <- plot_forest(meta[[name]])
  # Fixed-horizon survival evidence stays separate from binary evidence.
  survival <- surv_example(120, seed = 330)
  survival_task <- surv_task(survival, "patient_id", "time", "event", features,
    time_unit = "month", prediction_horizon = 6)
  survival_fit <- fit_final(survival_task, sb_pipeline(learner = "coxph"))
  survival_external <- surv_example(100, seed = 331)
  survival_external$patient_id <- paste0("surv_external_", seq_len(nrow(survival_external)))
  survival_evaluation <- eval_external(
    predict_final(survival_fit, survival_external[c("patient_id", features)], "patient_id"),
    surv_endpoint(survival_external, "patient_id", "time", "event",
      time_unit = "month", prediction_horizon = 6))
  survival_bootstrap <- sb_bootstrap_external(survival_evaluation, reps = reps, seed = 332)
  plots$Survival_horizon_heatmap <- plot_heatmap(survival_evaluation, "ipcw_brier")
  # Bootstrap interval keys retain their saved name, ipcw_brier_score.
  plots$Survival_horizon_intervals <- plot_ext_ci(survival_bootstrap, "ipcw_brier_score")
  plots$No_selection_available <- plot_freq(cv$results[[1]])
  # The small number of bootstrap replicates is for a fast software tutorial,
  # not a recommended precision for research. Only three synthetic cohorts:
  # pooled intervals and heterogeneity estimates are highly uncertain.
  if (!is.null(output_dir)) {
    if (file.exists(output_dir)) stop("Choose a new output directory; existing output is never overwritten.")
    if (!dir.create(output_dir, recursive = TRUE)) stop("Cannot create output directory.")
    heights <- rep(7, length(plots))
    heights[1:2] <- 9
    for (i in seq_along(plots)) {
      ggplot2::ggsave(file.path(output_dir, paste0(names(plots)[i], ".png")),
        plots[[i]], width = 12, height = heights[i], dpi = 140, bg = "white")
    }
    sb_render_visual_report(plots, file.path(output_dir, "gallery.html"),
      title = "STORMBREAKER 0.3.0 | Synthetic visualization gallery",
      width = 12, height = heights)
  }
  invisible(list(plots = plots, cv = cv, evaluations = evaluations,
    bootstraps = bootstraps, paired = paired, meta = meta, explanation = explanation,
    survival_evaluation = survival_evaluation, survival_bootstrap = survival_bootstrap))
}
