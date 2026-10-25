test_that("five PLS learners retain exact native components without enlarging portable state", {
  skip_if_not_installed("mixOmics")
  skip_if_not_installed("plsRglm")
  skip_if_not_installed("plsRcox")
  d <- make_sb_binary_fixture(n = 100L, seed = 108L)
  binary <- sb_binary_task(d, "patient_id", "outcome", c("gene_c", "gene_a", "gene_b"))
  for (id in c("plsrglm_binary", "plsda_binary", "splsda_binary", "plsrcox", "plsrcox_cv")) {
    task <- if (grepl("binary", id)) binary else make_component_survival_task()
    pipeline <- if (id == "plsrcox_cv") sb_pipeline(learner = id, learner_tuning_folds = 3L,
      learner_seed = 8903L) else sb_pipeline(learner = id)
    fit <- sb_fit_final(task, pipeline, capture_diagnostics = TRUE)
    compact <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success")
    expect_identical(compact$status, "success")
    expect_identical(fit$learner_state, compact$learner_state)
    expect_null(fit$learner_state$diagnostic_components)
    # Version 2 compares the full saved state without ALTREP representation
    # changes when native deferred score row names are materialized by unname().
    original <- serialize(fit, NULL, version = 2L)
    loadings <- sb_pls_diagnostics(fit)
    scores <- sb_pls_diagnostics(fit, "score_distribution")
    model <- fit$diagnostic_artifacts$learner
    native <- if (startsWith(id, "plsrcox")) list(loadings = model$engine$pp,
      scores = model$engine$tt) else model$diagnostic_components
    expect_identical(loadings$native_table, native$loadings)
    expect_identical(scores$native_table, unname(native$scores))
    expect_equal(loadings$table$loading, as.vector(native$loadings))
    expect_equal(scores$table$score, as.vector(native$scores))
    expect_equal(loadings$table$feature[loadings$table$component == 1L], model$features)
    expect_identical(scores$sensitive_training_data, TRUE)
    expect_null(rownames(scores$native_table))
    if (requireNamespace("ggplot2", quietly = TRUE)) for (q in list(loadings, scores)) expect_s3_class(
      ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
    expect_true(identical(serialize(fit, NULL, version = 2L), original), info = id)
    err <- tryCatch(sb_pls_diagnostics(compact), error = identity)
    expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
  }
})

test_that("PLS diagnostics preserve sparse zeros and actual one-component fits", {
  skip_if_not_installed("mixOmics")
  set.seed(108L)
  x <- as.data.frame(matrix(rnorm(100L * 12L), 100L))
  names(x) <- paste0("gene_", seq_len(12L))
  d <- data.frame(id = seq_len(100L), outcome = rbinom(100L, 1L,
    plogis(x[[1]] - x[[2]])), x)
  task <- sb_binary_task(d, "id", "outcome", names(x))
  fit <- sb_fit_final(task, sb_pipeline(learner = "splsda_binary"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  q <- sb_pls_diagnostics(fit)
  expect_equal(nrow(q$table), 24L)
  expect_true(any(q$table$loading == 0))
  expect_identical(q$native_table, fit$diagnostic_artifacts$learner$diagnostic_components$loadings)
  task <- sb_binary_task(d, "id", "outcome", "gene_1")
  fit <- sb_fit_final(task, sb_pipeline(learner = "plsda_binary"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  q <- sb_pls_diagnostics(fit)
  expect_identical(unique(q$table$component), 1L)
  expect_equal(nrow(q$table), 1L)
  fit$diagnostic_artifacts$learner$diagnostic_components <- NULL
  expect_error(sb_pls_diagnostics(fit), class = "stormbreaker_missing_diagnostic_artifacts")
})
