test_that("LIBLINEAR logistic explanations match native probabilities in input units", {
  skip_if_not_installed("LiblineaR")
  d <- sb_example_binary(n = 100L, seed = 108L)
  d$gene_c <- d$gene_c * 2 + 4
  d$gene_a <- d$gene_a * 0.3 - 7
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  for (id in c("liblinear_l2_logistic_binary", "liblinear_l1_logistic_binary")) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success", info = id)
    before <- serialize(fit, NULL, version = 2L)
    e <- unserialize(fit$learner_state$engine_raw)
    x <- as.matrix(d[features])
    native_x <- sweep(sweep(x, 2L, e$center, "-"), 2L, e$scale, "/")
    native <- predict(e$model, native_x, proba = TRUE)$probabilities[, "1"]
    parameters <- sb_model_coefficients(fit, d)
    score <- as.numeric(x %*% parameters$table$coefficient + parameters$intercept)
    expect_identical(parameters$table$feature, features)
    expect_equal(plogis(score), unname(native), tolerance = 1e-12)
    expect_equal(sb_predict_final(fit, d, "patient_id")$predictions$probability,
      unname(native), tolerance = 1e-12)
    nomogram <- sb_nomogram(fit, d)
    expect_equal(sb_nomogram_values(nomogram, d)$score, score, tolerance = 1e-12)
    shap <- sb_additive_shap(fit, d, "patient_id", d)
    expect_identical(shap$score_scale, "class_1_log_odds")
    expect_equal(unname(rowSums(shap$contributions)) + shap$reference_value,
      score, tolerance = 1e-12)
    expected_shap <- sweep(sweep(x, 2L, colMeans(x), "-"),
      2L, parameters$table$coefficient, "*")
    expect_equal(unname(shap$contributions), unname(expected_shap), tolerance = 1e-12)
    expect_s3_class(ggplot2::ggplot_build(sb_plot_nomogram(nomogram)), "ggplot_built")
    diagnostic <- sb_linear_coefficients_diagnostic(fit)
    expect_equal(diagnostic$table$coefficient, parameters$table$coefficient)
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(diagnostic)), "ggplot_built")
    expect_identical(serialize(fit, NULL, version = 2L), before)
    # Equivalent native encoding with the opposite first class tests sign handling.
    reversed <- fit
    e$model$ClassNames <- rev(e$model$ClassNames)
    e$model$W <- -e$model$W
    reversed$learner_state$engine_raw <- serialize(e, NULL, version = 2L)
    rp <- sb_model_coefficients(reversed)
    expect_equal(rp$table$coefficient, parameters$table$coefficient, tolerance = 1e-12)
    expect_equal(rp$intercept, parameters$intercept, tolerance = 1e-12)
    expect_equal(unname(predict(e$model, native_x, proba = TRUE)$probabilities[, "1"]),
      unname(native), tolerance = 1e-12)
  }
})
