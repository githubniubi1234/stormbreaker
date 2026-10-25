test_that("neural diagnostics preserve every native connection and forward prediction", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("neuralnet")
  d <- sb_example_binary(n = 160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  for (id in c("nnet_binary", "neuralnet_multilayer_binary")) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success", info = id)
    before <- serialize(fit, NULL, version = 2L)
    diagnostic <- sb_neural_diagnostics(fit)
    parameters <- diagnostic$parameters
    e <- unserialize(fit$learner_state$engine_raw)
    if (id == "nnet_binary") {
      expect_equal(unlist(parameters$weights, use.names = FALSE), e$model$wts, tolerance = 0)
    } else expect_equal(parameters$weights, e$model$weights[[1L]], tolerance = 0)
    expect_equal(nrow(diagnostic$table), sum(vapply(parameters$weights, length, integer(1))))
    expect_equal(sum(diagnostic$nodes$bias), length(parameters$weights))
    expect_identical(diagnostic$nodes$label[diagnostic$nodes$layer == 0 & !diagnostic$nodes$bias], features)
    value <- sweep(sweep(as.matrix(d[features]), 2L, parameters$center, "-"), 2L, parameters$scale, "/")
    for (weights in parameters$weights) {
      score <- cbind(1, value) %*% weights
      value <- if (id == "nnet_binary") ifelse(score < -15, 0,
        ifelse(score > 15, 1, 1 / (1 + exp(-score)))) else e$model$act.fct(score)
    }
    expect_equal(as.numeric(value), sb_predict_final(fit, d, "patient_id")$predictions$probability, tolerance = 1e-12)
    expect_error(sb_neural_diagnostics(fit, "learning_curve"), class = "stormbreaker_unsupported_diagnostic")
    expect_error(sb_nomogram(fit, d), class = "stormbreaker_unsupported_coefficients")
    bundle <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", 2L))
    if (requireNamespace("ggplot2", quietly = TRUE)) for (view in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(view)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_null(p$labels$subtitle)
      expect_null(p$labels$caption)
    }
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
