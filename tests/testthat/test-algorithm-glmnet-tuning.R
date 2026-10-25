test_that("all six glmnet task/mixing variants expose unchanged saved tuning", {
  skip_if_not_installed("glmnet")
  for (kind in c("binary", "survival")) {
    d <- if (kind == "binary") make_sb_binary_fixture(n = 90L, seed = 108L) else
      make_sb_fixture(n = 90L, seed = 8612L)
    task <- if (kind == "binary") sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")) else
      sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
    for (alpha in c(0, .5, 1)) {
      learner <- if (kind == "binary") "glmnet_logistic" else "glmnet_cox"
      fit <- sb_fit_final(task, sb_pipeline(learner = learner,
        inner_folds = 3L, inner_seed = 5L, glmnet_alpha = alpha))
      expect_identical(fit$status, "success")
      expect_null(fit$diagnostic_artifacts)
      before <- fit
      for (type in c("cv_curve", "nonzero_path")) {
        q <- sb_algorithm_diagnostics(fit, type)
        expect_identical(q$table[names(fit$tuning$curve)], fit$tuning$curve)
        expect_equal(q$table$lambda[q$table$selected], fit$tuning$lambda)
        expect_identical(q$alpha, alpha)
        if (requireNamespace("ggplot2", quietly = TRUE)) {
          expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
        }
      }
      expect_identical(fit, before)
    }
  }
})

test_that("tuning diagnostic rejects missing or corrupt native uncertainty", {
  skip_if_not_installed("glmnet")
  d <- make_sb_binary_fixture(n = 90L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "glmnet_logistic", inner_folds = 3L, inner_seed = 5L))
  fit$tuning$curve$standard_error[1] <- -1
  err <- tryCatch(sb_algorithm_diagnostics(fit, "cv_curve"), error = identity)
  expect_s3_class(err, "stormbreaker_invalid_diagnostic")
  fit$tuning$curve <- NULL
  err <- tryCatch(sb_algorithm_diagnostics(fit, "nonzero_path"), error = identity)
  expect_s3_class(err, "stormbreaker_invalid_diagnostic")
})
