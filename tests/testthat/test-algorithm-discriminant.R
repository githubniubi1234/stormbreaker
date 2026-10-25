test_that("discriminant diagnostics preserve native prediction equations and compact state", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("klaR")
  d <- sb_example_binary(n = 160L, seed = 108L)
  d$gene_a <- 10 + d$gene_a
  d$gene_b <- 3 * d$gene_b - 5
  for (features in list(c("gene_c", "gene_a", "gene_b"), "gene_a")) {
    task <- sb_binary_task(d, "patient_id", "outcome", features)
    for (id in c("lda_binary", "qda_binary", "rda_binary")) {
      fit <- sb_fit_final(task, sb_pipeline(learner = id))
      expect_identical(fit$status, "success", info = id)
      expect_null(fit$diagnostic_artifacts)
      before <- serialize(fit, NULL, version = 2L)
      diagnostic <- sb_discriminant_diagnostics(fit)
      p <- diagnostic$parameters
      x <- as.matrix(d[features])
      score <- p$intercept + drop(x %*% p$coefficient) + rowSums((x %*% p$quadratic) * x)
      native <- sb_predict_final(fit, d, "patient_id")$predictions$probability
      expect_equal(plogis(score), native, tolerance = 1e-12, info = id)
      e <- unserialize(fit$learner_state$engine_raw)
      expect_equal(p$means, if (id == "rda_binary") t(e$means) else e$means)
      expect_equal(diagnostic$table$value, as.vector(p$means))
      expect_equal(p$priors, e$prior)
      expect_identical(names(p$coefficient), features)
      expect_identical(dimnames(p$quadratic), list(features, features))
      if (id == "lda_binary") {
        expect_true(all(p$quadratic == 0))
        expect_error(sb_discriminant_diagnostics(fit, "class_correlation"),
          class = "stormbreaker_unavailable_diagnostic")
        coefficients <- sb_model_coefficients(fit, d)
        expect_equal(unname(coefficients$table$coefficient), unname(p$coefficient))
        nomogram <- sb_nomogram(fit, d)
        expect_equal(sb_nomogram_values(nomogram, d)$score, score, tolerance = 1e-12)
        shap <- sb_additive_shap(fit, d, "patient_id", d)
        expect_equal(unname(rowSums(shap$contributions)) + shap$reference_value, score, tolerance = 1e-12)
      } else {
        correlation <- sb_discriminant_diagnostics(fit, "class_correlation")
        expect_equal(nrow(correlation$table), 2 * length(features)^2)
        for (class in c("0", "1")) expect_equal(correlation$table$value[correlation$table$class == class],
          as.vector(stats::cov2cor(p$covariance[[class]])))
        expect_error(sb_nomogram(fit, d), class = "stormbreaker_unsupported_coefficients")
        expect_error(sb_additive_shap(fit, d, "patient_id", d), class = "stormbreaker_unsupported_additive_shap")
      }
      expect_error(sb_discriminant_diagnostics(fit, "unknown"), class = "stormbreaker_unsupported_diagnostic")
      bundle <- sb_native_diagnostics(fit, complete = TRUE)
      expect_identical(bundle$status$status, rep("available", 2L))
      if (requireNamespace("ggplot2", quietly = TRUE)) for (view in bundle$diagnostics) {
        plot <- sb_plot_algorithm_diagnostic(view)
        expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
        expect_null(plot$labels$subtitle)
        expect_null(plot$labels$caption)
      }
      expect_identical(serialize(fit, NULL, version = 2L), before)
    }
  }
})
