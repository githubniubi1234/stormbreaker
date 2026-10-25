test_that("fern split mapping reproduces native scores without changing fits", {
  skip_if_not_installed("rFerns")
  for (features in list("gene_c", c("gene_c", "gene_a", "gene_b"))) {
    d <- sb_example_binary(100L, seed = 108L)
    fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features),
      sb_pipeline(learner = "random_ferns_binary"))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    e <- unserialize(fit$learner_state$engine_raw)
    usage <- sb_ferns_diagnostics(fit)
    thresholds <- sb_ferns_diagnostics(fit, "thresholds")
    expect_equal(usage$table$splits, as.integer(table(factor(thresholds$table$feature, levels = features))))
    expect_identical(thresholds$table$threshold, e$model$threReal)
    query <- d[c(1L, 8L, 13L), features, drop = FALSE]
    query[[1L]] <- query[[1L]] + .137
    depth <- unname(e$parameters["depth"])
    nferns <- unname(e$parameters["ferns"])
    scores <- array(e$model$scores, c(length(e$classLabels), 2^depth, nferns))
    reconstructed <- matrix(0, nrow(query), length(e$classLabels))
    for (f in seq_len(nferns)) {
      z <- thresholds$table[thresholds$table$fern == f, ]
      leaves <- vapply(seq_len(nrow(query)), function(i) {
        sum((as.numeric(query[i, z$feature, drop = FALSE]) < z$threshold) * 2^(z$split - 1L)) + 1L
      }, numeric(1))
      reconstructed <- reconstructed + t(scores[, leaves, f, drop = FALSE][, , 1L])
    }
    native <- get("predict.rFerns", asNamespace("rFerns"))(e, query, scores = TRUE)
    expect_equal(reconstructed / nferns, unname(as.matrix(native)), tolerance = 1e-12)
    for (z in list(usage, thresholds)) {
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    }
    expect_identical(sb_native_diagnostics(fit, complete = TRUE)$status$status, rep("available", 2L))
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
