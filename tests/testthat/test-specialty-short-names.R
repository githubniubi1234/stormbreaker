test_that("every specialty short name preserves the existing function and arguments", {
  pairs <- read.csv(system.file("guides", "specialty-name-map.csv",
    package = "stormbreaker", mustWork = TRUE), stringsAsFactors = FALSE)
  exports <- getNamespaceExports("stormbreaker")
  expect_gt(nrow(pairs), 0L)
  expect_contains(exports, pairs$short)
  expect_contains(exports, pairs$original)
  for (i in seq_len(nrow(pairs))) {
    short <- getExportedValue("stormbreaker", pairs$short[i])
    old <- getExportedValue("stormbreaker", pairs$original[i])
    expect_identical(short, old, info = pairs$short[i])
    expect_identical(formals(short), formals(old), info = pairs$short[i])
    expect_lt(nchar(pairs$short[i]), nchar(pairs$original[i]))
  }
})

test_that("all diagnostic and drawing exports have a preferred short name", {
  preferred <- read.csv(system.file("guides", "preferred-name-map.csv",
    package = "stormbreaker", mustWork = TRUE), stringsAsFactors = FALSE)
  exports <- getNamespaceExports("stormbreaker")
  scoped <- grep("^sb_plot_|_diagnostics$", exports, value = TRUE)
  expect_contains(preferred$original, scoped)
  expect_identical(anyDuplicated(preferred$short), 0L)
  expect_identical(plot_diag, plot_algorithm)
  expect_identical(plot_shap_imp, plot_shap_importance)
  expect_identical(plot_shap_dep, plot_shap_dependence)
  expect_identical(plot_shap_one, plot_shap_waterfall)
})

test_that("short Lasso interfaces preserve native tables and coordinates", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("ggplot2")
  data <- bi_example(80, seed = 42)
  task <- bi_task(data, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  fit <- fit_final(task, bi_pipeline("Lasso"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  cv <- path_diag(fit, "cv_curve")
  path <- path_diag(fit, "coefficient_path")
  expect_identical(cv, sb_algorithm_diagnostics(fit, "cv_curve"))
  expect_identical(path, sb_algorithm_diagnostics(fit, "coefficient_path"))
  expect_gt(nrow(path$table), 0L)
  expect_identical(plot_diag(cv)$data, sb_plot_algorithm_diagnostic(cv)$data)
  expect_identical(native_diag(fit, complete = TRUE)$status,
    sb_native_diagnostics(fit, complete = TRUE)$status)
})

test_that("short forest interface preserves original importance", {
  skip_if_not_installed("randomForest")
  skip_if_not_installed("ggplot2")
  data <- bi_example(80, seed = 42)
  task <- bi_task(data, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  fit <- fit_final(task, bi_pipeline("RF"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  diagnostic <- forest_diag(fit, "importance")
  expect_identical(diagnostic, sb_forest_diagnostics(fit, "importance"))
  expect_equal(nrow(diagnostic$table), length(fit$required_features))
  expect_identical(plot_diag(diagnostic)$data, sb_plot_algorithm_diagnostic(diagnostic)$data)
})

test_that("short PAMR interface preserves both centroid stages", {
  skip_if_not_installed("pamr")
  skip_if_not_installed("ggplot2")
  data <- bi_example(80, seed = 42)
  task <- bi_task(data, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  fit <- fit_final(task, bi_pipeline("pamr"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  diagnostic <- pamr_diag(fit)
  expect_identical(diagnostic, sb_pamr_diagnostics(fit))
  expect_setequal(diagnostic$table$stage, c("Before", "After"))
  expect_equal(nrow(diagnostic$table), length(fit$required_features) * 4L)
  expect_identical(plot_diag(diagnostic)$data, sb_plot_algorithm_diagnostic(diagnostic)$data)
})

test_that("short diagnostic interfaces preserve missing-capture failures", {
  data <- bi_example(60, seed = 42)
  task <- bi_task(data, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- fit_final(task, bi_pipeline("Logistic"))
  old <- tryCatch(sb_algorithm_diagnostics(fit), error = identity)
  new <- tryCatch(path_diag(fit), error = identity)
  expect_s3_class(new, "stormbreaker_missing_diagnostic_artifacts")
  expect_identical(class(new), class(old))
  expect_identical(conditionMessage(new), conditionMessage(old))
  expect_null(fit$diagnostic_artifacts)
})
