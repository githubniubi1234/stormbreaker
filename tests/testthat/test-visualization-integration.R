test_that("all five real Meta outputs and the synthetic tutorial retain saved endpoints", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("metafor")
  path <- system.file("examples", "visualization-core.R", package = "stormbreaker")
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env)
  x <- env$stormbreaker_visualization_example(reps = 40L)
  expect_length(x$plots, 19L)
  for (p in x$plots) expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  for (m in x$meta) {
    before <- serialize(m, NULL)
    p <- sb_plot_meta_forest(m)
    expect_identical(tail(p$data$point_estimate, 2), rep(m$summary$estimate, 2))
    expect_identical(tail(p$data$lower, 2), c(m$summary$lower, m$summary$prediction_lower))
    expect_identical(tail(p$data$upper, 2), c(m$summary$upper, m$summary$prediction_upper))
    expect_identical(serialize(m, NULL), before)
  }
  expected <- subset(x$survival_bootstrap$intervals, metric == "ipcw_brier_score")
  p <- sb_plot_external_intervals(x$survival_bootstrap, "ipcw_brier_score")
  expect_identical(p$data$lower, expected$lower)
  expect_identical(p$data$upper, expected$upper)
  expect_match(as.character(p$data$target_label), "6 month", fixed = TRUE)
})

test_that("real repeated and IECV results remain separate plotting panels", {
  skip_if_not_installed("ggplot2")
  tasks <- lapply(1:2, function(i) {
    d <- sb_example_binary(60, seed = 500 + i)
    d$patient_id <- paste0(i, "_", d$patient_id)
    sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  })
  names(tasks) <- c("A", "B")
  pipelines <- list(sb_pipeline(learner = "logistic"))
  repeated <- sb_compare_repeated_cv(tasks[[1]], pipelines, folds = 3, seeds = c(501, 502))
  iecv <- sb_run_internal_external_cv(tasks, pipelines)
  for (x in list(repeated, iecv)) {
    before <- serialize(x, NULL)
    roc <- sb_plot_roc(x, by_fold = TRUE)
    pr <- sb_plot_pr(x)
    expect_equal(length(unique(roc$data$panel)), if (identical(x, repeated)) 6L else 2L)
    expect_equal(length(unique(pr$data$panel)), length(unique(roc$data$panel)))
    expect_s3_class(ggplot2::ggplot_build(sb_plot_performance_heatmap(x, "roc_auc")), "ggplot_built")
    expect_identical(serialize(x, NULL), before)
  }
})

test_that("real glmnet active sets and stack roles are usable without pooling", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("glmnet")
  d <- sb_example_binary(90, seed = 510)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  cv <- sb_run_cv(task, sb_pipeline(learner = "glmnet_logistic", glmnet_alpha = 1,
    inner_folds = 3, inner_seed = 511), folds = 3, seed = 511)
  p <- sb_plot_selection_frequency(cv, source = "glmnet_active_set")
  expect_identical(p$data$selection_frequency,
    cv$glmnet_active_set_stability$feature_frequency$selection_frequency)
  for (fun in list(sb_plot_selection_jaccard, sb_plot_selection_upset)) {
    expect_s3_class(ggplot2::ggplot_build(fun(cv, source = "glmnet_active_set")), "ggplot_built")
  }
  stack <- sb_stack(list(sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")), inner_folds = 3)
  result <- sb_run_stack_cv(task, stack, folds = 3, seed = 512)
  p <- sb_plot_roc(result, by_fold = TRUE)
  expect_setequal(unique(p$data$model_role), c("member", "stack"))
  expect_equal(length(unique(p$data$panel)), 6L)
  expect_s3_class(ggplot2::ggplot_build(sb_plot_performance_heatmap(result, "roc_auc")), "ggplot_built")
})

test_that("real Cox and naive Bayes SHAP contributions reconstruct plotted scores", {
  skip_if_not_installed("ggplot2")
  for (learner in c("coxph", "gaussian_naive_bayes")) {
    d <- if (learner == "coxph") sb_example_survival(80, seed = 520) else sb_example_binary(80, seed = 521)
    features <- c("gene_a", "gene_b", "gene_c")
    task <- if (learner == "coxph") sb_surv_task(d, "patient_id", "time", "event", features, time_unit = "month") else
      sb_binary_task(d, "patient_id", "outcome", features)
    fit <- sb_fit_final(task, sb_pipeline(learner = learner))
    explanation <- sb_additive_shap(fit, d[1:10, ], "patient_id", d[features])
    before <- serialize(explanation, NULL)
    for (fun in list(sb_plot_shap_importance, sb_plot_shap_beeswarm, sb_plot_shap_waterfall)) {
      expect_s3_class(ggplot2::ggplot_build(fun(explanation, top_n = 2)), "ggplot_built")
    }
    expect_identical(serialize(explanation, NULL), before)
  }
})
