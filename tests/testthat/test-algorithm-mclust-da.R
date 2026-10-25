test_that("model-based discriminant views preserve full component distributions", {
  skip_if_not_installed("mclust")
  d <- sb_example_binary(180L, seed = 108L)
  d$gene_c <- d$gene_c + rep(c(-5, 5), 90L)
  d$gene_b <- .4 * d$gene_c + d$gene_b
  for (features in list("gene_c", c("gene_c", "gene_a", "gene_b"))) {
    for (id in sb_mclust_da_learner_ids()) {
      fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features), sb_pipeline(learner = id))
      expect_identical(fit$status, "success")
      before <- serialize(fit, NULL, version = 2L)
      e <- fit$learner_state$engine
      centers <- sb_mclust_da_diagnostics(fit)
      correlations <- sb_mclust_da_diagnostics(fit, "correlations")
      weights <- sb_mclust_da_diagnostics(fit, "mixture_weights")
      if (id == "mclustda_binary_bic") expect_gt(max(weights$table$component), 1L)
      q <- d[c(1L, 8L, 13L), ]; q$gene_c <- q$gene_c + .71
      mass <- matrix(0, nrow(q), 2L)
      for (cls in c("0", "1")) {
        par <- e$models[[cls]]$parameters
        w <- weights$table[weights$table$class == cls, ]
        expect_equal(w$weight, as.numeric(par$pro))
        for (j in seq_along(par$pro)) {
          rows <- centers$table$class == cls & centers$table$component == j
          mu <- centers$table$mean[rows]
          expect_identical(centers$table$feature[rows], features)
          expect_equal(mu, as.numeric(matrix(par$mean, nrow = length(features))[, j]))
          covariance <- if (length(features) == 1L) matrix(rep(par$variance$sigmasq, length.out = length(par$pro))[j]) else par$variance$sigma[, , j]
          cor_rows <- correlations$table$class == cls & correlations$table$component == j
          cor <- matrix(correlations$table$correlation[cor_rows], nrow = length(features))
          expect_equal(cor, unname(stats::cov2cor(covariance)))
          sd <- sqrt(diag(covariance))
          recovered_cov <- cor * outer(sd, sd)
          delta <- sweep(as.matrix(q[features]), 2L, mu)
          density <- exp(-rowSums((delta %*% solve(recovered_cov)) * delta) / 2) /
            sqrt((2 * pi)^length(features) * det(recovered_cov))
          k <- match(cls, c("0", "1"))
          mass[, k] <- mass[, k] + density * w$weight[j] * w$class_prior[j]
        }
      }
      expect_equal(mass[, 2L] / rowSums(mass), sb_predict_final(fit, q, "patient_id")$predictions$probability, ignore_attr = TRUE, tolerance = 1e-12)
      for (z in list(centers, correlations, weights)) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
      expect_identical(sb_native_diagnostics(fit, complete = TRUE)$status$status, rep("available", 3L))
      expect_identical(serialize(fit, NULL, version = 2L), before)
    }
  }
})
