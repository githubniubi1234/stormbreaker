test_that("SuperPC screening, loadings and CV are read from the original objects", {
  skip_if_not_installed("superpc")
  for (id in c("superpc", "superpc_cv")) {
    pipeline <- if (id == "superpc_cv") sb_pipeline(learner = id,
      learner_tuning_folds = 3L, learner_seed = 8903L) else sb_pipeline(learner = id)
    fit <- sb_fit_final(make_component_survival_task(), pipeline, capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    original <- serialize(fit, NULL)
    model <- fit$diagnostic_artifacts$learner
    engine <- model$engine
    screen <- sb_superpc_diagnostics(fit)
    loadings <- sb_superpc_diagnostics(fit, "loadings")
    expect_identical(screen$native_table, engine$trained$feature.scores)
    expect_equal(screen$table$selected, unname(engine$prediction$which.features))
    expect_identical(loadings$native_table, engine$prediction$u)
    expect_equal(loadings$table$loading, as.vector(engine$prediction$u))
    expect_equal(loadings$table$feature, model$features[engine$prediction$which.features])
    expect_equal(screen$table$selected, abs(screen$table$score) >= screen$threshold)
    if (requireNamespace("ggplot2", quietly = TRUE)) for (q in list(screen, loadings)) expect_s3_class(
      ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
    cv <- tryCatch(sb_superpc_diagnostics(fit, "cv_curve"), error = identity)
    if (id == "superpc_cv") {
      native <- engine$cross_validation
      expect_equal(cv$table$threshold, rep(native$thresholds, each = nrow(native$scor)))
      expect_equal(cv$table$score, as.vector(native$scor))
      expect_equal(cv$table$lower, as.vector(native$scor.lower))
      expect_equal(cv$table$upper, as.vector(native$scor.upper))
      expect_equal(cv$table$component[cv$table$selected], 1L)
      expect_equal(cv$table$threshold[cv$table$selected], native$thresholds[which.max(native$scor[1, ])])
      expect_equal(unname(cv$folds), c(model$audit$requested_folds, native$n.fold))
      if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(
        ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(cv)), "ggplot_built")
    } else expect_s3_class(cv, "stormbreaker_unavailable_diagnostic")
    expect_identical(serialize(fit, NULL), original)
    fit$diagnostic_artifacts <- NULL
    err <- tryCatch(sb_superpc_diagnostics(fit), error = identity)
    expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
  }
})

test_that("SuperPC preserves its one-selected-feature singular vector", {
  skip_if_not_installed("superpc")
  fit <- sb_fit_final(make_component_survival_task(),
    sb_pipeline(learner = "superpc", superpc_threshold = 2.5), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  q <- sb_superpc_diagnostics(fit, "loadings")
  expect_equal(nrow(q$table), 1L)
  expect_identical(q$native_table, fit$diagnostic_artifacts$learner$engine$prediction$u)
  expect_equal(q$table$feature, "gene_01")
  if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(
    ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
})
