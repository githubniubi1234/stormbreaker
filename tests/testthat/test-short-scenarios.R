run_short_scenario <- function(file) {
  environment <- new.env(parent = globalenv())
  sys.source(system.file("examples", "scenarios", file,
    package = "stormbreaker", mustWork = TRUE), envir = environment)
  environment$example_result
}

test_that("binary scenario connects saved CV tables and three clinical plots", {
  skip_if_not_installed("ggplot2")
  x <- run_short_scenario("01-binary.R")
  expect_equal(nrow(x$input), 120L)
  expect_equal(nrow(x$metrics), 1L)
  expect_identical(x$metrics$status, "success")
  expect_setequal(names(x$figures), c("ROC", "Calibration", "DCA"))
  expect_gt(nrow(x$figures$ROC$data), 0L)
  expect_identical(x$metrics$estimate, x$comparison$results[[1]]$metrics$roc_auc$estimate)
})

test_that("comparison scenario returns both declared candidates", {
  skip_if_not_installed("ggplot2")
  x <- run_short_scenario("02-compare.R")
  expect_equal(nrow(x$metrics), 2L)
  expect_identical(x$metrics$status, rep("success", 2))
  expect_length(x$comparison$results, 2L)
  expect_s3_class(x$figures$Performance, "ggplot")
  expect_gt(nrow(x$figures$ROC$data), 0L)
})

test_that("external scenario keeps both evaluations and the original predictions", {
  skip_if_not_installed("ggplot2")
  x <- run_short_scenario("03-external.R")
  expect_identical(x$finished$fit$status, "success")
  expect_equal(nrow(x$finished$predictions$Validation_A), 80L)
  expect_equal(nrow(x$finished$predictions$Validation_B), 70L)
  expect_identical(x$metrics$status, rep("success", 2))
  expect_setequal(x$metrics$cohort_label, c("Validation_A", "Validation_B"))
  expect_setequal(names(x$finished$direction_correction), c("Validation_A", "Validation_B"))
  expect_identical(x$finished$direction_correction$Validation_A$reverse, TRUE)
  expect_identical(x$finished$direction_correction$Validation_B$reverse, TRUE)
  original_a <- predict_final(x$finished$fit, x$external$Validation_A,
    "patient_id")$predictions
  original_b <- predict_final(x$finished$fit, x$external$Validation_B,
    "patient_id")$predictions
  expect_identical(x$finished$predictions$Validation_A, original_a)
  expect_identical(x$finished$predictions$Validation_B, original_b)
  expect_identical(x$finished$predictions$Validation_A$probability,
    x$finished$evaluated_predictions$Validation_A$probability)
  expect_identical(x$finished$predictions$Validation_B$probability,
    x$finished$evaluated_predictions$Validation_B$probability)
  expect_s3_class(x$figures$ROC, "ggplot")
})

test_that("survival scenario keeps horizon units and complete external rows", {
  skip_if_not_installed("ggplot2")
  x <- run_short_scenario("04-survival.R")
  p <- x$prediction$predictions
  expect_identical(x$fit$status, "success")
  expect_identical(x$evaluation$status, "success")
  expect_equal(nrow(p), 100L)
  expect_equal(unique(p$prediction_horizon), 6)
  expect_identical(unique(p$time_unit), "month")
  expect_gte(min(p$survival_probability), 0)
  expect_lte(max(p$survival_probability), 1)
  expect_setequal(x$horizon_metrics$metric, c("time_auc", "ipcw_brier"))
  expect_equal(unique(x$horizon_metrics$prediction_horizon), 6)
  expect_null(x$evaluation$direction_correction)
  expect_gt(nrow(x$figures$Calibration$data), 0L)
})

test_that("selector scenario retains native CV and coefficient evidence", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("glmnet")
  x <- run_short_scenario("05-diagnostics.R")
  expect_identical(x$fit$status, "success")
  expect_setequal(x$genes$gene, x$fit$required_features)
  expect_contains(c("gene_a", "gene_b", "gene_c"), x$genes$gene)
  expect_gt(nrow(x$genes), 0L)
  expect_contains(names(x$figures), c("selector_1_cv", "selector_1_coefficient_path", "Coefficients"))
  expect_s3_class(x$figures$selector_1_cv, "ggplot")
  expect_length(x$diagnostics$unavailable, 0L)
})

test_that("SHAP scenario reconstructs probabilities from its declared background", {
  skip_if_not_installed("ggplot2")
  x <- run_short_scenario("06-shap.R")
  expect_identical(dim(x$shap$contributions), c(8L, 3L))
  expect_identical(rownames(x$shap$contributions), x$query$patient_id)
  expect_identical(x$shap$score_scale, "class_1_probability")
  expect_equal(unname(x$shap$reference_value + rowSums(x$shap$contributions)),
    x$shap$scores$score, tolerance = 1e-10)
  expect_setequal(names(x$figures), c("Importance", "Beeswarm", "Dependence", "Waterfall"))
  expect_s3_class(x$figures$Waterfall, "ggplot")
})

test_that("tutorial code can be extracted and replayed independently", {
  skip_if_not_installed("knitr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("glmnet")
  input <- test_path("..", "..", "vignettes", "stormbreaker-short-scenarios.Rmd")
  if (!file.exists(input)) {
    input <- system.file("doc", "stormbreaker-short-scenarios.Rmd",
      package = "stormbreaker", mustWork = TRUE)
  }
  directory <- withr::local_tempdir()
  output <- file.path(directory, "tutorial.R")
  knitr::purl(input, output = output, quiet = TRUE)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    paste0("root <- ", paste(deparse(find.package("stormbreaker")), collapse = "")),
    "if (file.exists(file.path(root, 'Meta', 'package.rds'))) {",
    "  library(stormbreaker, lib.loc = dirname(root))",
    "} else {",
    "  pkgload::load_all(root, helpers = FALSE)",
    "}",
    "set.seed(42)",
    "grDevices::pdf('tutorial.pdf')",
    "source('tutorial.R', echo = FALSE, print.eval = TRUE)",
    "grDevices::dev.off()",
    "saveRDS(list(metrics = metrics, contributions = shap$contributions,",
    "  reference_value = shap$reference_value, scores = shap$scores$score), 'result.rds')"
  ), file.path(directory, "replay.R"), useBytes = TRUE)
  withr::local_dir(directory)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "replay.R"),
    stdout = "replay.log", stderr = "replay.log")
  expect_identical(status, 0L,
    info = paste(tail(readLines("replay.log", warn = FALSE), 20), collapse = "\n"))
  if (status != 0L) return(invisible(NULL))
  result <- readRDS("result.rds")
  expect_s3_class(result$metrics, "data.frame")
  expect_identical(dim(result$contributions), c(8L, 3L))
  expect_equal(
    unname(result$reference_value + rowSums(result$contributions)),
    result$scores, tolerance = 1e-10)
})
