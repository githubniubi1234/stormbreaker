test_that("public final fits replay constrained predictions and original-unit views", {
  skip_if_not_installed("brglm2")
  withr::local_seed(410)
  x <- as.data.frame(matrix(stats::rnorm(320), 80, 4))
  names(x) <- paste0("gene", 1:4)
  x[] <- sweep(as.matrix(x), 1L, rowMeans(x), "-")
  d <- x
  d$id <- paste0("p", seq_len(nrow(d)))
  d$outcome <- rep(0:1, 40)
  task <- sb_binary_task(d, "id", "outcome", names(x))
  for (id in c("brglm2_logistic", "lda_binary")) {
    fit <- sb_fit_final(task, sb_model_pipelines("binary", id)[[1L]])
    expect_identical(fit$status, "success")
    prediction <- sb_predict_final(fit, d[c("id", names(x))], "id")
    expect_identical(prediction$status, "success")
    replay <- sb_predict_final(unserialize(serialize(fit, NULL)), d[c("id", names(x))], "id")
    expect_identical(replay$predictions, prediction$predictions)
    coefficients <- sb_model_coefficients(fit)
    expect_identical(coefficients$table$feature, names(x))
    expect_equal(stats::plogis(coefficients$intercept +
      as.vector(as.matrix(x) %*% coefficients$table$coefficient)),
      prediction$predictions$probability, tolerance = 1e-10)
    if (id == "lda_binary") {
      diagnostic <- sb_discriminant_diagnostics(fit)
      expect_equal(nrow(diagnostic$table), 2L * ncol(x))
      expect_equal(unname(diagnostic$parameters$means["0", ]),
        unname(colMeans(x[d$outcome == 0L, ])), tolerance = 1e-12)
    }
  }
})
