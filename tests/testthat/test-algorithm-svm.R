test_that("five LIBSVM learners expose native counts and oriented decisions", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(n = 100L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in c("svm_linear_binary", "svm_radial_binary",
      "svm_radial_weighted_binary", "svm_polynomial_binary", "svm_sigmoid_binary")) {
    task <- sb_binary_task(d, "patient_id", "outcome", features)
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    original <- serialize(fit, NULL, version = 2L)
    e <- fit$diagnostic_artifacts$learner$engine
    counts <- sb_svm_diagnostics(fit)
    native_counts <- table(factor(d$outcome[e$index], levels = c(0, 1)))
    expect_equal(counts$table$support_vectors[match(c("0", "1"), counts$table$class)], as.integer(native_counts))
    q <- sb_svm_diagnostics(fit, "decision_distribution")
    native <- stats::predict(e, data.matrix(d[features]), decision.values = TRUE)
    expect_identical(q$native_table, unname(attr(native, "decision.values")))
    expect_equal(as.integer(q$table$decision > 0), as.integer(as.character(native)))
    expect_identical(q$table$class, as.character(d$outcome))
    expect_true(q$sensitive_training_data)
    expect_null(rownames(q$native_table))
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(counts)), "ggplot_built")
    compact <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(compact$status, "success")
    expect_identical(compact$learner_state, fit$learner_state)
    expect_identical(sb_svm_diagnostics(compact), counts)
    expect_error(sb_svm_diagnostics(compact, "decision_distribution"), class = "stormbreaker_missing_diagnostic_artifacts")
    bad <- fit
    bad$diagnostic_artifacts$predictors <- bad$diagnostic_artifacts$predictors[rev(features)]
    expect_error(sb_svm_diagnostics(bad, "decision_distribution"), class = "stormbreaker_invalid_diagnostic")
    if (id != "svm_linear_binary") expect_error(sb_svm_diagnostics(fit, "coefficients"),
      class = "stormbreaker_unavailable_diagnostic")
    expect_true(identical(serialize(fit, NULL, version = 2L), original))
  }
})

test_that("linear SVM weights reproduce native decisions for either native label order", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(n = 100L, seed = 108L)
  d$gene_a <- 30 + 15 * d$gene_a
  features <- c("gene_b", "gene_a", "gene_c")
  orientations <- numeric()
  for (descending in c(FALSE, TRUE)) {
    data <- d[order(d$outcome, decreasing = descending), ]
    task <- sb_binary_task(data, "patient_id", "outcome", features)
    fit <- sb_fit_final(task, sb_pipeline(learner = "svm_linear_binary"), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    q <- sb_svm_diagnostics(fit, "coefficients")
    scores <- sb_svm_diagnostics(fit, "decision_distribution")
    expect_equal(as.numeric(data.matrix(data[features]) %*% q$table$coefficient + q$intercept),
      scores$table$decision, tolerance = 1e-10)
    expect_identical(q$table$feature, features)
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
    orientations <- c(orientations, q$orientation)
  }
  expect_setequal(orientations, c(-1, 1))
})
