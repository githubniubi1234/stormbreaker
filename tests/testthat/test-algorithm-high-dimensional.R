test_that("HDA and HDDA original geometry reproduces native probabilities", {
  skip_if_not_installed("hda")
  skip_if_not_installed("HDclassif")
  d <- sb_example_binary(160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in c("hda_binary", "hdda_binary_cattell")) {
    fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features), sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    e <- fit$learner_state$engine
    loadings <- sb_high_dimensional_diagnostics(fit)
    q <- d[c(1L, 8L, 13L), ]; q$gene_c <- q$gene_c + .71
    mass <- matrix(0, nrow(q), 2L)
    if (id == "hda_binary") {
      densities <- sb_high_dimensional_diagnostics(fit, "densities")
      m <- matrix(loadings$table$loading, nrow = length(features))
      expect_equal(m, unname(e$hda.loadings[, seq_len(e$reduced.dimension), drop = FALSE]))
      expect_identical(loadings$table$feature, rep(features, e$reduced.dimension))
      projected <- as.matrix(q[features]) %*% m
      native_projected <- get("predict.hda", asNamespace("hda"))(e, as.matrix(q[features]), task = "dr")
      expect_equal(unname(projected), unname(native_projected))
      for (i in 1:2) {
        par <- densities$parameters[densities$parameters$class == as.character(i - 1L), ]
        expect_equal(par$mean, unname(e$naivebayes$tables[[1L]][i, 1L]))
        mass[, i] <- dnorm(projected[, 1L], par$mean, par$sd) * as.numeric(e$naivebayes$apriori[i])
      }
      views <- list(loadings, densities)
    } else {
      spectrum <- sb_high_dimensional_diagnostics(fit, "spectrum")
      centers <- sb_high_dimensional_diagnostics(fit, "centers")
      for (i in 1:2) {
        cls <- e$kname[i]
        rows <- loadings$table$class == cls
        Q <- matrix(loadings$table$loading[rows], nrow = length(features))
        expect_equal(Q, e$Q[[i]])
        spec <- spectrum$table[spectrum$table$class == cls, ]
        expect_equal(spec$eigenvalue, as.numeric(e$ev[i, ]))
        expect_equal(sum(spec$retained), as.integer(e$d[i]))
        k <- ncol(Q)
        a <- spec$model_variance[seq_len(k)]
        b <- spec$model_variance[k + 1L]
        cov <- Q %*% diag(a - b, k) %*% t(Q) + diag(b, length(features))
        cov <- cov * outer(e$scaling$sd, e$scaling$sd)
        mu <- centers$table$mean[centers$table$class == cls]
        expect_equal(mu, as.numeric(e$mu[i, ]) * as.numeric(e$scaling$sd) + as.numeric(e$scaling$mu))
        delta <- sweep(as.matrix(q[features]), 2L, mu)
        mass[, i] <- exp(-rowSums((delta %*% solve(cov)) * delta) / 2) /
          sqrt((2 * pi)^length(features) * det(cov)) * as.numeric(e$prop[i])
      }
      views <- list(loadings, spectrum, centers)
    }
    expect_equal(as.numeric(mass[, 2L] / rowSums(mass)), sb_predict_final(fit, q, "patient_id")$predictions$probability, tolerance = 1e-12)
    for (z in views) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    expect_identical(sb_native_diagnostics(fit, complete = TRUE)$status$status, rep("available", length(views)))
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
