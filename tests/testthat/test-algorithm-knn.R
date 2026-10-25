test_that("explicit KNN neighborhoods reproduce original predictions", {
  skip_if_not_installed("kknn")
  d <- sb_example_binary(n = 100L, seed = 108L)
  d$gene_c <- 2 * d$gene_c + 4
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  for (id in c("knn_binary", "kernel_weighted_knn_binary")) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    e <- unserialize(fit$learner_state$engine_raw)
    for (rows in list(1L, c(1L, 7L, 11L))) {
      query <- d[rows, ]
      if (length(rows) > 1L) query$gene_c <- query$gene_c + 0.17
      z <- sb_knn_neighborhood(fit, query)
      expect_equal(nrow(z$table), length(rows) * e$k)
      expect_equal(z$predictions$probability,
        sb_predict_final(fit, query, "patient_id")$predictions$probability, tolerance = 1e-12)
      expect_equal(as.numeric(tapply(z$table$voting_weight * (z$table$class == "1"),
        z$table$query, sum)), z$predictions$probability, tolerance = 1e-12)
      if (id == "knn_binary") {
        for (i in seq_along(rows)) {
          standardized <- (as.numeric(query[i, features]) - e$center) / e$scale
          distances <- apply(e$training, 1L, function(row) sum((row - standardized)^2))
          chosen <- order(distances, seq_along(distances))[seq_len(e$k)]
          tab <- z$table[z$table$query == i, ]
          expect_identical(tab$training_row, chosen)
          expect_equal(tab$distance, sqrt(unname(distances[chosen])), tolerance = 1e-12)
          expect_equal(tab$voting_weight, rep(1 / e$k, e$k))
        }
      } else {
        training <- e$training
        training$y <- e$outcome
        native <- kknn::kknn(y ~ ., training, query[features], k = e$k,
          distance = e$distance, kernel = e$kernel, scale = e$scale)
        expect_equal(z$table$training_row, as.integer(t(native$C)))
        expect_equal(z$table$distance, as.numeric(t(native$D)))
        expect_equal(z$table$native_weight, as.numeric(t(native$W)))
        expect_identical(z$table$class, as.character(t(native$CL)))
      }
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
      expect_identical(sb_native_diagnostics(fit, newdata = query)$status$status, "available")
    }
    expect_identical(sb_native_diagnostics(fit)$status$status, "unavailable")
    expect_match(sb_native_diagnostics(fit)$status$reason, "explicit query")
    if (id == "knn_binary") {
      tied <- fit
      e$training[,] <- 0
      tied$learner_state$engine_raw <- serialize(e, NULL, version = 2L)
      tie_query <- as.data.frame(as.list(e$center))
      explained <- sb_knn_neighborhood(tied, tie_query)
      expect_identical(explained$table$training_row, seq_len(e$k))
      expect_equal(explained$table$distance, rep(0, e$k))
      expect_equal(explained$predictions$probability, mean(e$outcome[seq_len(e$k)]))
    }
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
