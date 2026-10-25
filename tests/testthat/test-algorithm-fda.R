test_that("FDA/PDA/MDA native geometry reconstructs posterior probabilities", {
  skip_if_not_installed("mda")
  d <- sb_example_binary(160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in c("fda_binary", "pda_binary", "mda_binary")) {
    fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features), sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    e <- unserialize(fit$learner_state$engine_raw)
    centers <- sb_fda_diagnostics(fit)
    means <- matrix(centers$table$coordinate, nrow = nrow(e$means))
    expect_equal(means, unname(as.matrix(e$means)), ignore_attr = TRUE)
    mapping <- centers$table[centers$table$axis == 1L, ]
    weights <- if (id == "mda_binary") sb_fda_diagnostics(fit, "mixture_weights") else NULL
    if (!is.null(weights)) {
      expect_equal(weights$table$weight, as.numeric(unlist(e$sub.prior)))
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(weights)), "ggplot_built")
    }
    for (rows in list(8L, c(1L, 8L, 13L))) {
      q <- d[rows, ]; q$gene_c <- q$gene_c + .71
      z <- sb_fda_diagnostics(fit, "query_variates", q)
      coordinates <- matrix(z$table$coordinate, nrow = nrow(q))
      batch <- rbind(q[features], d[22L, features, drop = FALSE])
      native_coordinates <- predict(e, batch, type = "variates", dimension = e$dimension)[seq_len(nrow(q)), , drop = FALSE]
      expect_equal(coordinates, unname(as.matrix(native_coordinates)), ignore_attr = TRUE)
      mass <- matrix(0, nrow(q), 2L)
      for (i in seq_len(nrow(means))) {
        cls <- match(mapping$class[i], c("0", "1"))
        weight <- if (is.null(weights)) 1 else weights$table$weight[i]
        distance <- rowSums(sweep(coordinates, 2L, means[i, ], "-")^2)
        mass[, cls] <- mass[, cls] + exp(-distance / 2) * weight * as.numeric(e$prior[cls])
      }
      posterior <- mass[, 2L] / rowSums(mass)
      expect_equal(posterior, sb_predict_final(fit, q, "patient_id")$predictions$probability, tolerance = 1e-12)
      native_posterior <- predict(e, batch, type = "posterior")[seq_len(nrow(q)), "1"]
      expect_equal(posterior, unname(native_posterior), tolerance = 1e-12)
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    }
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(centers)), "ggplot_built")
    expect_identical(tail(sb_native_diagnostics(fit, complete = TRUE)$status$status, 1L), "unavailable")
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
