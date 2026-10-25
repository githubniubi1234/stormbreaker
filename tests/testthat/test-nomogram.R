test_that("nomogram points exactly reconstruct the fitted logistic predictor", {
  d <- make_sb_binary_fixture(n = 60L, seed = 108L)
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b")),
    sb_pipeline(learner = "logistic"))
  original <- fit
  n <- sb_nomogram(fit, d)
  v <- sb_nomogram_values(n, d)
  expected <- drop(as.matrix(d[fit$learner_state$features]) %*% fit$learner_state$coefficient) + fit$learner_state$intercept
  expect_equal(v$score, expected, tolerance = 1e-12)
  expect_equal(v$probability, plogis(expected), tolerance = 1e-12)
  expect_identical(fit, original)
  expect_equal(sum(v$outside_display_range), 0)
  outside <- d[1, ]
  outside$gene_a <- max(d$gene_a) + 100
  vv <- sb_nomogram_values(n, outside)
  expect_identical(vv$outside_display_range, TRUE)
  expect_equal(vv$score, sum(unlist(outside[fit$learner_state$features]) * fit$learner_state$coefficient) + fit$learner_state$intercept)
  skip_if_not_installed("ggplot2")
  expect_s3_class(ggplot2::ggplot_build(sb_plot_nomogram(n)), "ggplot_built")
})

test_that("nomogram display preserves zero terms and spaces short axes", {
  d <- make_sb_binary_fixture(n = 60L, seed = 108L)
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b")),
    sb_pipeline(learner = "logistic"))
  # Synthetic coefficient configurations exercise exact point conversion.
  fit$learner_state$coefficient[] <- c(0, -.00001)
  n <- sb_nomogram(fit, d)
  p <- sb_plot_nomogram(n)
  expect_identical(attr(p, "stormbreaker_nomogram_omitted_zero_terms"), "gene_a")
  expect_equal(nrow(n$features), 2)
  expect_equal(nrow(attr(sb_plot_nomogram(n, show_zero = TRUE), "stormbreaker_nomogram_axes")), 5)
  expect_equal(sb_nomogram_values(n, d)$score,
    fit$learner_state$intercept + d$gene_b * -.00001, tolerance = 1e-12)
  fit$learner_state$coefficient[] <- 0
  z <- sb_nomogram(fit, d)
  expect_equal(z$total_points, 0)
  expect_equal(sb_nomogram_values(z, d)$probability,
    rep(plogis(fit$learner_state$intercept), nrow(d)))
  expect_s3_class(ggplot2::ggplot_build(sb_plot_nomogram(z)), "ggplot_built")
  constant <- d
  constant$gene_a <- 5
  constant$gene_b <- 10
  c <- sb_nomogram(fit, constant)
  expect_equal(sb_nomogram_values(c, constant)$total_points, rep(0, nrow(d)))
})

test_that("penalized nomogram reproduces the original portable predictor", {
  skip_if_not_installed("glmnet")
  d <- make_sb_binary_fixture(n = 90L, seed = 108L)
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b")),
    sb_pipeline(learner = "glmnet_logistic", inner_folds = 3L, inner_seed = 5L))
  expect_identical(fit$status, "success")
  n <- sb_nomogram(fit, d)
  v <- sb_nomogram_values(n, d)
  prediction <- sb_predict_final(fit, d, "patient_id")
  expect_identical(prediction$status, "success")
  expect_equal(v$probability, prediction$predictions$probability, tolerance = 1e-12)
})

test_that("nomogram preserves the original Cox horizon baseline", {
  d <- make_sb_fixture(n = 72L, seed = 9701L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month", prediction_horizon = 6)
  fit <- sb_fit_final(task, sb_pipeline())
  n <- sb_nomogram(fit, d)
  v <- sb_nomogram_values(n, d)
  expected <- drop(as.matrix(d[fit$learner_state$features]) %*% fit$learner_state$coefficient)
  expect_equal(v$score, expected, tolerance = 1e-12)
  expect_equal(v$probability, exp(-fit$learner_state$baseline_hazard$cumulative_baseline_hazard * exp(expected)), tolerance = 1e-12)
  expect_identical(n$baseline, fit$learner_state$baseline_hazard)
  expect_identical(v$probability_type, "survival_at_saved_horizon")
})
