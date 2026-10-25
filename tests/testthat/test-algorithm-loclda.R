test_that("local-LDA neighborhoods reproduce native local posterior calculations", {
  skip_if_not_installed("klaR")
  d <- sb_example_binary(160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features), sb_pipeline(learner = "loclda_binary"))
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  e <- unserialize(fit$learner_state$engine_raw)
  for (rows in list(8L, c(1L, 8L, 13L))) {
    q <- d[rows, ]; q$gene_c <- q$gene_c + .71
    z <- sb_loclda_neighborhood(fit, q)
    for (i in seq_len(nrow(q))) {
      t <- z$table[z$table$query == i, ]
      distances <- as.numeric(dist(rbind(as.matrix(q[i, features]), e$learn)))[seq_len(nrow(e$learn))]
      ids <- which(distances <= sort(distances)[e$k])
      expect_identical(t$reference_row, ids)
      expect_equal(t$distance, distances[ids])
      expect_equal(t$weight, vapply(distances[ids] / sort(distances)[e$k], e$weight.func, numeric(1)))
      t <- t[t$used, ]
      pooled <- matrix(0, length(features), length(features))
      means <- matrix(0, 2L, length(features))
      prior <- numeric(2L)
      for (j in 1:2) {
        tt <- t[t$class == as.character(j - 1L), ]
        X <- e$learn[tt$reference_row, , drop = FALSE]
        means[j, ] <- colSums(X * tt$weight) / sum(tt$weight)
        pooled <- pooled + cov.wt(X, wt = tt$weight)$cov * nrow(X)
        prior[j] <- sum(tt$weight) / sum(t$weight)
      }
      pooled <- pooled / (nrow(t) - 2L)
      delta <- sweep(means, 2L, as.numeric(q[i, features]))
      score <- -.5 * rowSums((delta %*% solve(pooled)) * delta) + log(prior)
      expect_equal(plogis(score[2L] - score[1L]), z$predictions$probability[i], tolerance = 1e-12)
    }
    native <- get("predict.loclda", asNamespace("klaR"))(e, as.matrix(q[features]))$posterior[, 2L]
    expect_equal(z$predictions$probability, as.numeric(native))
    expect_equal(z$predictions$probability, sb_predict_final(fit, q, "patient_id")$predictions$probability)
    built <- ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z))
    expect_s3_class(built, "ggplot_built")
    expect_identical(as.character(built$plot$scales$get_scales("shape")$get_labels()), "Used")
    expect_identical(sb_native_diagnostics(fit, newdata = q)$status$status, "available")
  }
  expect_identical(sb_native_diagnostics(fit)$status$status, "unavailable")
  expect_identical(serialize(fit, NULL, version = 2L), before)
})
