test_that("earth diagnostics retain logistic rather than least-squares coefficients", {
  skip_if_not_installed("earth")
  d <- sb_example_binary(n = 160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  fit <- sb_fit_final(task, sb_pipeline(learner = "earth_binary"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  native <- fit$diagnostic_artifacts$learner$engine
  diagnostic <- sb_earth_diagnostics(fit)
  expect_equal(diagnostic$table$coefficient, as.numeric(native$glm.coefficients), tolerance = 0)
  expect_false(isTRUE(all.equal(diagnostic$table$coefficient, as.numeric(native$coefficients))))
  expect_identical(diagnostic$table$term, rownames(native$glm.coefficients))
  expect_equal(diagnostic$table$basis_index, native$selected.terms)
  expect_equal(diagnostic$cuts, native$cuts, tolerance = 0)
  expect_equal(diagnostic$directions, native$dirs, tolerance = 0)
  basis <- get("model.matrix.earth", asNamespace("earth"))(native, d[features])
  expect_equal(plogis(drop(basis %*% diagnostic$table$coefficient)),
    sb_predict_final(fit, d, "patient_id")$predictions$probability, tolerance = 1e-12)
  shap <- sb_additive_shap(fit, d, "patient_id", d)
  expect_equal(unname(rowSums(shap$contributions)) + shap$reference_value,
    drop(basis %*% diagnostic$table$coefficient), tolerance = 1e-12)
  expect_identical(colnames(shap$contributions), features)
  expect_identical(shap$score_scale, "class_1_log_odds")
  pruning <- sb_earth_diagnostics(fit, "pruning_gcv")
  expect_equal(pruning$table$gcv, as.numeric(native$gcv.per.subset), tolerance = 0)
  expect_equal(pruning$table$terms[pruning$table$selected], length(native$selected.terms))
  compact <- fit
  compact$diagnostic_artifacts <- NULL
  expect_equal(sb_earth_diagnostics(compact)$table, diagnostic$table)
  expect_error(sb_earth_diagnostics(compact, "pruning_gcv"), class = "stormbreaker_missing_diagnostic_artifacts")
  expect_identical(sb_native_diagnostics(compact, complete = TRUE)$status$status, c("available", "missing_capture"))
  expect_error(sb_nomogram(fit, d), class = "stormbreaker_unsupported_coefficients")
  bundle <- sb_native_diagnostics(fit, complete = TRUE)
  expect_identical(bundle$status$status, rep("available", 2L))
  if (requireNamespace("ggplot2", quietly = TRUE)) for (view in bundle$diagnostics) {
    p <- sb_plot_algorithm_diagnostic(view)
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    expect_null(p$labels$subtitle)
    expect_null(p$labels$caption)
  }
  expect_identical(serialize(fit, NULL, version = 2L), before)
})
