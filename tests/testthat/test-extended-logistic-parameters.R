test_that("five additional logistic predictors retain exact original-unit explanations", {
  d <- sb_example_binary(n = 160L, seed = 108L)
  d$gene_a <- d$gene_a + 10
  d$gene_b <- 3 * d$gene_b - 5
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  dependencies <- c(stepglm_binary = "stats", glmboost_binary = "mboost",
    plsrglm_binary = "plsRglm", brglm2_logistic = "brglm2", bayesglm_logistic_cauchy = "arm")
  for (id in names(dependencies)) {
    skip_if_not_installed(dependencies[[id]])
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success", info = id)
    before <- serialize(fit, NULL, version = 2L)
    p <- sb_model_coefficients(fit, d)
    score <- p$intercept + as.numeric(as.matrix(d[features]) %*% p$table$coefficient)
    prediction <- sb_predict_final(fit, d, "patient_id")
    expect_equal(plogis(score), prediction$predictions$probability, tolerance = 1e-12, info = id)
    nomogram <- sb_nomogram(fit, d)
    expect_equal(sb_nomogram_values(nomogram, d)$score, score, tolerance = 1e-12, info = id)
    shap <- sb_additive_shap(fit, d, "patient_id", d)
    expect_identical(shap$score_scale, "class_1_log_odds")
    expect_equal(unname(rowSums(shap$contributions)) + shap$reference_value, score, tolerance = 1e-12, info = id)
    expect_identical(serialize(fit, NULL, version = 2L), before)
    expect_equal(unname(p$table$scaled_coefficient), unname(p$table$coefficient * vapply(d[features], sd, numeric(1))))
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      expect_s3_class(ggplot2::ggplot_build(sb_plot_nomogram(nomogram)), "ggplot_built")
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(sb_linear_coefficients_diagnostic(fit))), "ggplot_built")
    }
  }
})
