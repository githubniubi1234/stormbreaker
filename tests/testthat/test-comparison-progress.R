test_that("CV progress describes actual folds without changing results", {
  d <- sb_example_binary(45, 10)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  pipelines <- list(sb_pipeline(learner = "logistic"))
  silent <- sb_compare_cv(task, pipelines, folds = 3, seed = 5, progress = FALSE)
  messages <- capture.output(
    verbose <- sb_compare_cv(task, pipelines, folds = 3, seed = 5, progress = TRUE),
    type = "message"
  )
  expect_identical(verbose, silent)
  expect_match(paste(messages, collapse = " "), "CV model 1/1")
  expect_match(paste(messages, collapse = " "), "Fold 3/3")
  expect_match(paste(messages, collapse = " "), "complete")
})

test_that("cohort progress includes final fit, prediction and evaluation", {
  d <- sb_example_binary(45, 10)
  v <- sb_example_binary(25, 11)
  v$patient_id <- paste0("v", v$patient_id)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  val <- sb_binary_task(v, "patient_id", "outcome", c("gene_a", "gene_b"))
  messages <- capture.output(
    result <- sb_compare_cohorts(task, list(Validation = val),
      roles = c(Validation = "validation"), models = "logistic", folds = 3,
      view = FALSE, progress = TRUE), type = "message"
  )
  expect_s3_class(result, "stormbreaker_cohort_comparison")
  text <- paste(messages, collapse = " ")
  expect_match(text, "Final fit 1/1")
  expect_match(text, "Predict Validation")
  expect_match(text, "Evaluate Validation")
})

test_that("survival progress preserves results and invalid flags are rejected", {
  d <- sb_example_survival(45, 3)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month")
  silent <- sb_run_cv(task, folds = 3, progress = FALSE)
  messages <- capture.output(verbose <- sb_run_cv(task, folds = 3, progress = TRUE), type = "message")
  expect_identical(verbose, silent)
  expect_match(paste(messages, collapse = " "), "Fold 3/3 complete")
  expect_error(sb_run_cv(task, progress = NA), "TRUE or FALSE", class = "stormbreaker_bad_progress")
})
