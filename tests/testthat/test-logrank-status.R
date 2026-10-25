test_that("log-rank status distinguishes a successful test from unavailable groups", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  d <- sb_example_survival(60, 42)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b"), "month")
  fit <- sb_fit_final(task, surv_pipeline("Cox"))
  reference <- surv_risk_reference(fit, task)

  result <- surv_risk_plots(reference, task)

  native <- survival::survdiff(survival::Surv(time, event) ~ group, result$patients)
  expect_identical(result$logrank_status, "success")
  expect_identical(result$logrank_error, NA_character_)
  expect_equal(result$logrank_p, stats::pchisq(native$chisq, 1, lower.tail = FALSE))
  reference$cutoff <- Inf
  single <- surv_risk_plots(reference, task)
  expect_identical(single$status, "single_group")
  expect_identical(single$logrank_status, "single_group")
  expect_identical(single$logrank_p, NA_real_)
})

test_that("log-rank reports no events without inventing a P value", {
  d <- data.frame(time = 1:6, event = 0L,
    group = factor(rep(c("Low", "High"), each = 3)))

  result <- sb_risk_logrank(d)

  expect_identical(result$status, "no_events")
  expect_identical(result$p, NA_real_)
  expect_identical(result$error, NA_character_)
})

test_that("a log-rank error remains visible without discarding usable panels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  d <- sb_example_survival(60, 42)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b"), "month")
  fit <- sb_fit_final(task, surv_pipeline("Cox"))
  reference <- surv_risk_reference(fit, task)
  local_mocked_bindings(survdiff = function(...) stop("SYNTHETIC_LOGRANK_FAILURE"),
    .package = "survival")
  warnings <- character()

  result <- withCallingHandlers(surv_risk_plots(reference, task), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })

  expect_identical(result$logrank_status, "failed")
  expect_identical(result$logrank_error, "SYNTHETIC_LOGRANK_FAILURE")
  expect_identical(result$logrank_p, NA_real_)
  expect_match(warnings, "SYNTHETIC_LOGRANK_FAILURE")
  expect_s3_class(result$plot, "patchwork")
  expect_equal(nrow(result$patients), nrow(d))
})

test_that("completion retains and exports a failed log-rank calculation", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  d <- sb_example_survival(60, 42)
  run <- surv_run(d, "time", "event", c("gene_a", "gene_b"),
    "patient_id", models = "Cox", folds = 3, plots = FALSE)
  external <- d
  external$patient_id <- paste0("EXT_", external$patient_id)
  local_mocked_bindings(survdiff = function(...) stop("SYNTHETIC_LOGRANK_FAILURE"),
    .package = "survival")
  warnings <- character()

  result <- withCallingHandlers(sb_complete_run(run, "Cox", list(Validation = external)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  expect_identical(result$evaluations$Validation$status, "success")
  expect_identical(result$risk_groups$Validation$logrank_status, "failed")
  expect_match(result$issues$risk$Validation, "SYNTHETIC_LOGRANK_FAILURE")
  expect_match(warnings, "SYNTHETIC_LOGRANK_FAILURE")
  expect_s3_class(result$plots$external$Validation, "patchwork")
  output <- file.path(withr::local_tempdir(), "export")
  sb_export_run(result, output)
  exported <- read.csv(file.path(output, "risk-groups.csv"))
  expect_identical(exported$logrank_status, "failed")
  expect_identical(exported$logrank_error, "SYNTHETIC_LOGRANK_FAILURE")
  expect_match(paste(read.csv(file.path(output, "issues.csv"))$message, collapse = " "),
    "SYNTHETIC_LOGRANK_FAILURE")
})

test_that("completion does not claim a log-rank test when plots were not requested", {
  d <- sb_example_survival(60, 42)
  run <- surv_run(d, "time", "event", c("gene_a", "gene_b"),
    "patient_id", models = "Cox", folds = 3, plots = FALSE)
  external <- d
  external$patient_id <- paste0("EXT_", external$patient_id)
  local_mocked_bindings(survdiff = function(...) stop("Unexpected log-rank test"),
    .package = "survival")

  result <- sb_complete_run(run, "Cox", list(Validation = external), plots = FALSE)

  expect_identical(result$evaluations$Validation$status, "success")
  expect_identical(result$risk_groups$Validation$logrank_status, "not_computed")
  expect_identical(result$risk_groups$Validation$logrank_p, NA_real_)
  expect_identical(result$risk_groups$Validation$logrank_error, NA_character_)
})
