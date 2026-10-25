test_that("LIBLINEAR SVC and DWD weights preserve the original decision function", {
  skip_if_not_installed("LiblineaR")
  skip_if_not_installed("DWDLargeR")
  d <- sb_example_binary(n = 100L, seed = 108L)
  d$gene_a <- 20 + 5 * d$gene_a
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in sb_linear_margin_diagnostic_ids()) {
    task <- sb_binary_task(d, "patient_id", "outcome", features)
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    q <- sb_svm_diagnostics(fit, "coefficients")
    scores <- sb_svm_diagnostics(fit, "decision_distribution")
    expect_identical(q$table$feature, features)
    expect_equal(as.numeric(data.matrix(d[features]) %*% q$table$coefficient + q$intercept), scores$table$decision)
    e <- fit$diagnostic_artifacts$learner$engine
    if (id != "dwd_binary") {
      x <- sweep(sweep(data.matrix(d[features]), 2L, e$center, "-"), 2L, e$scale, "/")
      native <- stats::predict(e$model, x, decisionValues = TRUE)$decisionValues
      expect_equal(scores$table$decision, as.numeric(native[, "1"] - native[, "0"]), tolerance = 1e-10)
      p <- sb_predict_liblinear_binary(e, data.matrix(d[features]))
    } else {
      expect_identical(q$table$coefficient, as.numeric(e$coefficient))
      p <- plogis(2 * (e$intercept + data.matrix(d[features]) %*% e$coefficient))
    }
    factor <- if (id == "dwd_binary") 2 else 1
    expect_equal(plogis(factor * scores$table$decision), as.numeric(p), tolerance = 1e-10)
    for (plot in list(q, scores)) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(plot)), "ggplot_built")
    expect_true(scores$sensitive_training_data)
    compact <- fit
    compact$diagnostic_artifacts <- NULL
    expect_identical(sb_svm_diagnostics(compact, "coefficients"), q)
    expect_error(sb_svm_diagnostics(compact, "decision_distribution"), class = "stormbreaker_missing_diagnostic_artifacts")
    expect_error(sb_svm_diagnostics(fit, "support_counts"), class = "stormbreaker_unavailable_diagnostic")
    expect_true(identical(serialize(fit, NULL, version = 2L), before))
  }
})
