test_that("Cancerclass views preserve centroid order and the saved probability map", {
  skip_if_not_installed("cancerclass")
  skip_if_not_installed("Biobase")
  d <- sb_example_binary(160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features), sb_pipeline(learner = "cancerclass_binary"))
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  centroids <- sb_cancerclass_diagnostics(fit)
  mapping <- sb_cancerclass_diagnostics(fit, "score_mapping", c(-2, 3))
  m <- matrix(centroids$table$mean, nrow = length(features))
  expect_identical(centroids$table$feature, rep(features, 2L))
  expect_equal(m, unname(fit$learner_state$predictor[features, c("0", "1")]))
  expect_identical(mapping$calibration_coefficients, fit$learner_state$calibration_coefficients)
  expect_equal(range(mapping$table$score), c(-2, 3))
  b <- mapping$calibration_coefficients
  expect_equal(mapping$table$probability, plogis(b[1L] + b[2L] * mapping$table$score))
  q <- d[c(1L, 8L, 13L), features]; q$gene_c <- q$gene_c + .71
  z <- (cancerclass::get.d(t(as.matrix(q)), m[, 1L], method = "cor") -
    cancerclass::get.d(t(as.matrix(q)), m[, 2L], method = "cor")) /
    cancerclass::get.d(m[, 1L], m[, 2L], method = "cor")
  expect_equal(as.numeric(plogis(b[1L] + b[2L] * z)),
    sb_predict_final(fit, data.frame(patient_id = 1:3, q), "patient_id")$predictions$probability,
    tolerance = 1e-12)
  for (view in list(centroids, mapping)) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(view)), "ggplot_built")
  expect_identical(sb_native_diagnostics(fit, complete = TRUE)$status$status, rep("available", 2L))
  expect_identical(serialize(fit, NULL, version = 2L), before)
})
