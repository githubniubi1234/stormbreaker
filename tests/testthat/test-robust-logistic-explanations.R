test_that("robust logistic explanations preserve the native Mqle model", {
  skip_if_not_installed("robustbase")
  d <- sb_example_binary(n = 160L, seed = 108L)
  d$gene_c <- d$gene_c * 2 + 4
  d$gene_a <- d$gene_a * 0.3 - 7
  features <- c("gene_c", "gene_a", "gene_b")
  for (selected in list(features, "gene_c")) {
    task <- sb_binary_task(d, "patient_id", "outcome", selected)
    fit <- sb_fit_final(task, sb_pipeline(learner = "robust_logistic_mqle"))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    native_data <- d[selected]
    names(native_data) <- paste0(".stormbreaker_robust_x", seq_along(selected))
    native_data$.stormbreaker_robust_outcome <- factor(d$outcome, levels = c(0, 1))
    native <- robustbase::glmrob(
      reformulate(names(native_data)[seq_along(selected)], ".stormbreaker_robust_outcome"),
      family = binomial("logit"), data = native_data, method = "Mqle", weights.on.x = "none",
      control = robustbase::glmrobMqle.control(acc = 1e-4, test.acc = "coef", maxit = 50L, tcc = 3.5),
      model = FALSE, x = FALSE)
    expect_identical(native$converged, TRUE)
    parameters <- sb_model_coefficients(fit, d)
    expect_identical(parameters$table$feature, selected)
    expect_equal(unname(parameters$table$coefficient), unname(coef(native)[-1L]), tolerance = 1e-12)
    expect_equal(parameters$intercept, unname(coef(native)[1L]), tolerance = 1e-12)
    score <- as.numeric(as.matrix(d[selected]) %*% parameters$table$coefficient + parameters$intercept)
    expect_equal(plogis(score), unname(predict(native, native_data, type = "response")), tolerance = 1e-12)
    expect_equal(sb_predict_final(fit, d, "patient_id")$predictions$probability,
      plogis(score), tolerance = 1e-12)
    nomogram <- sb_nomogram(fit, d)
    expect_equal(sb_nomogram_values(nomogram, d)$score, score, tolerance = 1e-12)
    shap <- sb_additive_shap(fit, d, "patient_id", d)
    expected <- sweep(sweep(as.matrix(d[selected]), 2L, colMeans(d[selected]), "-"),
      2L, unname(coef(native)[-1L]), "*")
    expect_equal(unname(shap$contributions), unname(expected), tolerance = 1e-12)
    expect_equal(unname(rowSums(shap$contributions)) + shap$reference_value, score, tolerance = 1e-12)
    expect_identical(shap$score_scale, "class_1_log_odds")
    bundle <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(bundle$status$status, "available")
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(bundle$diagnostics[[1L]])), "ggplot_built")
    expect_s3_class(ggplot2::ggplot_build(sb_plot_nomogram(nomogram)), "ggplot_built")
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
