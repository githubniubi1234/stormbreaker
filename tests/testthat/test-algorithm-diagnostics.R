test_that("optional diagnostic capture preserves original model predictions", {
  skip_if_not_installed("glmnet")
  d <- make_sb_binary_fixture(n = 90L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  pipeline <- sb_pipeline(learner = "glmnet_logistic", inner_folds = 3L, inner_seed = 5L)
  a <- sb_fit_final(task, pipeline)
  b <- sb_fit_final(task, pipeline, capture_diagnostics = TRUE)
  expect_identical(b$status, "success")
  expect_null(a$diagnostic_artifacts)
  expect_identical(a$learner_state, b$learner_state)
  expect_identical(sb_predict_final(a, d, "patient_id"), sb_predict_final(b, d, "patient_id"))
  q <- sb_algorithm_diagnostics(b)
  expect_equal(q$table$coefficient, as.vector(as.matrix(b$diagnostic_artifacts$learner$beta)))
  expect_equal(length(unique(q$table$feature)), 2)
  skip_if_not_installed("ggplot2")
  expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
})

test_that("ABESS diagnostics preserve the native GIC path and chosen support", {
  skip_if_not_installed("abess", minimum_version = "0.4.11")
  for (kind in c("binary", "survival")) {
    if (kind == "binary") {
      d <- make_sb_binary_fixture(n = 84L, seed = 8611L)
      task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
      learner <- "abess_binary_gic"
    } else {
      d <- make_sb_fixture(n = 90L, seed = 8612L)
      task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
      learner <- "abess_cox_gic"
    }
    fit <- sb_fit_final(task, sb_pipeline(learner = learner), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- fit$learner_state
    q <- sb_algorithm_diagnostics(fit, "support_criterion")
    native <- fit$diagnostic_artifacts$learner
    coefficient_path <- sb_algorithm_diagnostics(fit, "coefficient_path")
    expect_equal(coefficient_path$table$coefficient, as.vector(as.matrix(native$engine$beta)))
    expect_s3_class(ggplot2::ggplotGrob(sb_plot_algorithm_diagnostic(coefficient_path)), "gtable")
    expect_identical(q$table[names(native$criterion_path)], native$criterion_path)
    expect_identical(q$table$support_size[q$table$selected], native$selected_support_size)
    expect_identical(fit$learner_state, before)
    expect_equal(sum(q$table$selected), 1)
    damaged <- fit
    damaged$diagnostic_artifacts$learner$criterion_path$gic[1] <- NA_real_
    err <- tryCatch(sb_algorithm_diagnostics(damaged, "support_criterion"), error = identity)
    expect_s3_class(err, "stormbreaker_invalid_diagnostic")
    missing <- fit
    missing$diagnostic_artifacts <- NULL
    err <- tryCatch(sb_algorithm_diagnostics(missing, "support_criterion"), error = identity)
    expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
    }
  }
})
