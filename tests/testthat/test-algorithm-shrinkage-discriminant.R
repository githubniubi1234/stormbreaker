test_that("SDA and SLDA diagnostics preserve native prediction geometry", {
  skip_if_not_installed("sda")
  skip_if_not_installed("ipred")
  d <- sb_example_binary(160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in c("sda_binary", "slda_binary")) {
    fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features), sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    e <- unserialize(fit$learner_state$engine_raw)
    first <- sb_shrinkage_discriminant_diagnostics(fit)
    second <- sb_shrinkage_discriminant_diagnostics(fit, if (id == "sda_binary") "regularization" else "class_means")
    if (id == "sda_binary") {
      expect_identical(first$table$feature, features)
      expect_equal(first$table$coefficient, as.numeric(e$beta[2L, ] - e$beta[1L, ]))
      expect_identical(second$table$value, as.numeric(e$regularization))
    } else {
      loadings <- matrix(first$table$loading, nrow = length(features))
      expect_equal(loadings, unname(e$scores))
      expect_identical(first$table$feature, rep(features, ncol(e$scores)))
      expect_equal(matrix(second$table$mean, nrow = 2L), unname(e$mylda$means))
    }
    for (rows in list(8L, c(1L, 8L, 13L))) {
      q <- d[rows, ]; q$gene_c <- q$gene_c + .71
      if (id == "sda_binary") {
        p <- plogis(as.numeric(as.matrix(q[features]) %*% first$table$coefficient) + first$intercept)
        reconstructed <- zapsmall(cbind(1 - p, p))[, 2L]
        native <- get("predict.sda", asNamespace("sda"))(e, as.matrix(q[features]), verbose = FALSE)$posterior[, "1"]
      } else {
        projected <- as.data.frame(as.matrix(q[features]) %*% loadings)
        names(projected) <- colnames(e$mylda$means)
        reconstructed <- predict(e$mylda, projected)$posterior[, "1"]
        native <- sb_predict_slda_binary(e, as.matrix(q[features]))
      }
      expect_equal(as.numeric(reconstructed), as.numeric(native), tolerance = 1e-12)
      expect_equal(as.numeric(reconstructed), sb_predict_final(fit, q, "patient_id")$predictions$probability, tolerance = 1e-12)
    }
    for (z in list(first, second)) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    expect_identical(sb_native_diagnostics(fit, complete = TRUE)$status$status, rep("available", 2L))
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
