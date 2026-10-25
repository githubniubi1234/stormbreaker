test_that("pooled ROC has one curve and direct scores have reproducible intervals", {
  d <- sb_example_binary(80, seed = 42)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  cv <- sb_run_cv(task, folds = 3)
  p <- sb_plot_roc(cv, ci_reps = 20)
  expect_length(unique(p$data$panel), 1L)
  expect_equal(attr(p, "stormbreaker_roc_statistics")$auc,
    sb_roc_interval(cv$predictions$outcome, cv$predictions$probability)$auc)
  expect_length(unique(sb_plot_roc(cv, by_fold = TRUE)$data$panel), 3L)
  before <- cv
  cv$predictions$probability[1] <- NA_real_
  expect_identical(attr(sb_plot_roc(cv), "stormbreaker_curve_status")$curve_status, "invalid_predictions")
  y <- c(1, 0, 1, 0); score <- c(.9, .7, .7, .1)
  a <- sb_plot_prediction_roc(y, score, kind = "external", ci_reps = 40)
  b <- sb_plot_prediction_roc(y, score, kind = "external", ci_reps = 40)
  expect_equal(attr(a, "stormbreaker_roc_statistics")$auc, .875)
  expect_identical(attr(a, "stormbreaker_roc_statistics"), attr(b, "stormbreaker_roc_statistics"))
  expect_identical(attr(a, "prediction_origin"), "external")
  expect_equal(nrow(ggplot2::ggplot_build(a)$layout$layout), 1)
})

test_that("ranking uses only declared cohorts and leaves incomplete models unranked", {
  skip_if_not_installed("patchwork")
  d <- data.frame(model = rep(c("a", "b", "c"), each = 3),
    cohort = rep(c("Train", "Validation", "Test"), 3),
    auc = c(.9, .8, .5, .8, .9, .4, .99, NA, 1))
  z <- sb_plot_model_ranking(d, "Validation", top_n = 2)
  expect_identical(z$ranking$model, c("b", "a", "c"))
  expect_equal(z$ranking$mean_auc, c(.9, .8, NA))
  expect_identical(z$ranking$status[3], "incomplete")
  expect_s3_class(z$plot, "patchwork")
  expect_equal(nrow(z$data), 6)
})

test_that("bootstrap calibration refits each sample and computes the correction", {
  d <- sb_example_binary(100, seed = 42)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  original <- sb_fit_final
  calls <- 0L
  local_mocked_bindings(sb_fit_final = function(...) {
    calls <<- calls + 1L
    original(...)
  })
  z <- sb_bootstrap_calibration(task, sb_pipeline(learner = "logistic"), reps = 12, seed = 7)
  expect_identical(calls, 13L)
  expect_identical(z$status$status, rep("success", 12))
  expect_equal(z$curve$corrected, z$curve$apparent - z$curve$optimism)
  expect_equal(z$curve$effective_reps, colSums(is.finite(z$replicates)))
  expect_gt(sum(is.finite(z$curve$corrected)), 20)
  expect_gt(max(abs(z$curve$optimism), na.rm = TRUE), 0)
  again <- sb_bootstrap_calibration(task, sb_pipeline(learner = "logistic"), reps = 12, seed = 7)
  expect_identical(z$curve, again$curve)
  expect_s3_class(sb_plot_bootstrap_calibration(z), "ggplot")
})

test_that("nomogram density and patient marks preserve original probabilities", {
  d <- sb_example_binary(100, seed = 42)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  nom <- sb_nomogram(fit, d)
  p <- sb_plot_nomogram(nom, reference_data = d, patient = d[1, , drop = FALSE])
  expect_equal(attr(p, "stormbreaker_nomogram_patient")$probability,
    sb_predict_final(fit, d[1, c("patient_id", "gene_a", "gene_b")], "patient_id")$predictions$probability)
  expect_gt(nrow(attr(p, "stormbreaker_nomogram_density")), 100)
  expect_null(p$labels$caption)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

test_that("standardized binary DCA divides by prevalence without changing inputs", {
  d <- sb_example_binary(80, seed = 42)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"), decision_thresholds = c(.1, .2, .5, .8))
  cv <- sb_run_cv(task, folds = 3)
  before <- serialize(cv, NULL)
  ordinary <- sb_plot_dca(cv)
  standard <- sb_plot_dca(cv, standardized = TRUE, cost_benefit_axis = TRUE)
  expect_equal(standard$data$y, ordinary$data$y / mean(d$outcome))
  expect_identical(serialize(cv, NULL), before)
  expect_identical(standard$labels$y, "Standardized net benefit")
  expect_s3_class(ggplot2::ggplot_build(standard), "ggplot_built")
})
