test_that("automatic extraction executes native routes and preserves unavailable reasons", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_binary(n = 100L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  for (id in c("gaussian_naive_bayes", "logistic", "decision_tree_binary", "knn_binary")) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    auto <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(serialize(fit, NULL, version = 2L), before)
    if (id == "gaussian_naive_bayes") {
      expect_identical(auto$status$status, c("available", "available"))
      expect_identical(auto$diagnostics[[1]], sb_naive_bayes_diagnostics(fit))
      expect_length(sb_native_diagnostics(fit)$diagnostics, 1L)
    } else if (id == "logistic") {
      expect_identical(auto$status$status, "available")
      expect_identical(auto$diagnostics[[1]]$table, sb_model_coefficients(fit)$table)
    } else {
      expected <- if (id == "knn_binary") "unavailable" else "missing_capture"
      expect_identical(auto$status$status, expected)
      expect_length(auto$diagnostics, 0L)
      result <- sb_run_cv(task, sb_pipeline(learner = id), folds = 3L, seed = 41L)
      views <- sb_clinician_views(result, native = auto)
      expect_true(expected %in% views$status$status)
      expect_identical(views$diagnostics[["Native extraction"]], auto$status)
    }
  }
  fit <- sb_fit_final(task, sb_pipeline(learner = "decision_tree_binary"), capture_diagnostics = TRUE)
  expect_identical(sb_native_diagnostics(fit)$diagnostics[[1]], sb_tree_structure(fit))
})

test_that("all catalog routes name real compatible extractors without claiming availability", {
  plans <- list()
  for (type in c("binary", "survival")) {
    catalog <- sb_model_catalog(type)
    expect_gt(nrow(catalog), 0L)
    plans[[type]] <- do.call(rbind, lapply(unique(catalog$learner), sb_native_diagnostic_plan))
  }
  plan <- do.call(rbind, plans)
  for (name in unique(plan$extractor)) {
    extractor <- get(name, envir = asNamespace("stormbreaker"))
    expect_true(is.function(extractor), info = name)
    expect_true("fit" %in% names(formals(extractor)), info = name)
    if (any(nzchar(plan$type[plan$extractor == name]))) {
      expect_true("type" %in% names(formals(extractor)), info = name)
    }
  }
})

test_that("automatic glmnet detail distinguishes compact tuning from missing native paths", {
  skip_if_not_installed("glmnet")
  d <- sb_example_binary(n = 100L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "glmnet_logistic", inner_folds = 3L, inner_seed = 5L))
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  auto <- sb_native_diagnostics(fit, complete = TRUE)
  expect_identical(auto$status$type, c("cv_curve", "nonzero_path", "coefficient_path"))
  expect_identical(auto$status$status, c("available", "available", "missing_capture"))
  expect_identical(auto$diagnostics[[1]], sb_algorithm_diagnostics(fit, "cv_curve"))
  expect_identical(serialize(fit, NULL, version = 2L), before)
})
