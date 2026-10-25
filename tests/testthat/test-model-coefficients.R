test_that("coefficients retain prediction units and optional reference scaling", {
  d <- make_sb_binary_fixture(n = 60L, seed = 108L)
  fit <- sb_fit_final(sb_binary_task(
    d, "patient_id", "outcome",
    c("gene_a", "gene_b")
  ), sb_pipeline(learner = "logistic"))
  x <- sb_model_coefficients(fit, d)
  expect_equal(x$table$coefficient, fit$learner_state$coefficient)
  expect_equal(x$table$reference_sd, unname(vapply(d[fit$required_features], sd, numeric(1))))
  expect_equal(x$table$scaled_coefficient, x$table$coefficient * x$table$reference_sd)
  expect_equal(x$intercept, fit$learner_state$intercept)
  expect_identical(x$score_scale, "class_1_log_odds")
  expect_identical(x$reference_n, nrow(d))
  expect_equal(sb_model_coefficients(fit)$table$coefficient, x$table$coefficient)
  skip_if_not_installed("ggplot2")
  p <- sb_plot_coefficients(x, scale = "reference_sd")
  expect_s3_class(p, "ggplot")
  expect_equal(sort(p$data$value), sort(x$table$scaled_coefficient))
})
