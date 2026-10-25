test_that("all LVQ diagnostic views retain original codebooks and prediction geometry", {
  skip_if_not_installed("class")
  d <- sb_example_binary(n = 160L, seed = 108L)
  d$gene_a <- d$gene_a + 10
  d$gene_b <- 4 * d$gene_b - 5
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  for (id in sb_lvq_binary_learner_ids()) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success", info = id)
    expect_null(fit$diagnostic_artifacts)
    before <- serialize(fit, NULL, version = 2L)
    prediction_before <- sb_predict_final(fit, d, "patient_id")
    e <- unserialize(fit$learner_state$engine_raw)
    diagnostic <- sb_lvq_diagnostics(fit)
    expect_identical(diagnostic$table$value, as.vector(e$prototypes))
    expect_equal(diagnostic$table$original_value,
      as.vector(sweep(sweep(e$prototypes, 2L, e$scale, "*"), 2L, e$center, "+")), tolerance = 0)
    expect_identical(diagnostic$table$feature, rep(features, each = nrow(e$prototypes)))
    expect_identical(diagnostic$table$class, rep(e$prototype_class, length(features)))
    expect_equal(diagnostic$prototype_counts$count,
      as.integer(table(factor(e$prototype_class, levels = c("0", "1")))))
    expect_identical(diagnostic$iterations, e$iterations)
    expect_identical(names(diagnostic$center), features)
    expect_equal(unname(diagnostic$scale), e$scale)
    distances <- sb_lvq_diagnostics(fit, "prototype_distances")
    expected <- vapply(seq_len(nrow(e$prototypes)), function(i)
      rowSums(sweep(e$prototypes, 2L, e$prototypes[i, ], "-")^2), numeric(nrow(e$prototypes)))
    expect_equal(distances$table$value, as.vector(expected), tolerance = 1e-12)
    expect_equal(sum(diag(matrix(distances$table$value, nrow(e$prototypes)))), 0)
    expect_error(sb_lvq_diagnostics(fit, "training_loss"), class = "stormbreaker_unsupported_diagnostic")
    bundle <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", 2L))
    if (requireNamespace("ggplot2", quietly = TRUE)) for (view in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(view)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_null(p$labels$subtitle)
      expect_null(p$labels$caption)
    }
    expect_identical(sb_predict_final(fit, d, "patient_id"), prediction_before)
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
