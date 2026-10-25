test_that("Gaussian-process kernel terms preserve native probabilities", {
  skip_if_not_installed("kernlab")
  d <- sb_example_binary(100L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features),
    sb_pipeline(learner = "gaussian_process_binary"))
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  e <- unserialize(fit$learner_state$engine_raw)
  weights <- sb_gaussian_process_diagnostics(fit)
  expect_equal(weights$table$weight, -as.numeric(e@alpha[[1L]]))
  expect_identical(weights$table$reference_row, e@alphaindex[[1L]])
  for (rows in list(8L, c(1L, 8L, 19L))) {
    q <- d[rows, ]; q$gene_c <- q$gene_c * 1.3 + .71
    z <- sb_gaussian_process_diagnostics(fit, "query_contributions", q)
    native <- kernlab::predict(e, as.matrix(q[features]), type = "probabilities")
    expect_equal(z$summary$probability, as.numeric(native[, "1"]), tolerance = 1e-12)
    expect_equal(plogis(z$summary$latent_score), z$summary$probability, tolerance = 1e-12)
    for (i in seq_along(rows)) {
      scaled <- (as.numeric(q[i, features]) - e@scaling$x.scale[["scaled:center"]]) /
        e@scaling$x.scale[["scaled:scale"]]
      terms <- vapply(seq_along(e@alphaindex[[1L]]), function(j) {
        as.numeric(e@kernelf(scaled, e@xmatrix[e@alphaindex[[1L]][j], ])) * -e@alpha[[1L]][j]
      }, numeric(1))
      expect_equal(z$table$contribution[z$table$query == i], terms)
    }
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    expect_identical(sb_native_diagnostics(fit, complete = TRUE, newdata = q)$status$status, rep("available", 2L))
  }
  expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(weights)), "ggplot_built")
  expect_identical(sb_native_diagnostics(fit, complete = TRUE)$status$status, c("available", "unavailable"))
  expect_identical(serialize(fit, NULL, version = 2L), before)
})
