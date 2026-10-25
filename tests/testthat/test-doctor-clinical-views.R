clinical_binary_run <- function(thresholds = NULL) {
  bi_run(sb_example_binary(120, seed = 42), "outcome",
    c("gene_a", "gene_b"), "patient_id", models = "Logistic",
    folds = 3, plots = FALSE, decision_thresholds = thresholds
  )
}

test_that("short binary entrypoint retains declared decision thresholds", {
  thresholds <- c(.1, .3, .5)
  run <- clinical_binary_run(thresholds)

  expect_identical(run$development_task$decision_thresholds, thresholds)
  expect_identical(run$comparison$results[[1]]$metrics$decision_curve$status, "success")
})

test_that("clinical views accept a doctor run without fitting again", {
  run <- clinical_binary_run(c(.1, .3, .5))
  saved <- run$comparison$results[[1]]$metrics$decision_curve$curve
  original <- run$comparison$results[[1]]$predictions
  local_mocked_bindings(sb_fit_final = function(...) stop("Unexpected refit"))

  views <- sb_clinician_views(run)

  expect_setequal(names(views$plots), c("ROC", "Calibration", "Decision curve"))
  expect_equal(views$plots[["Decision curve"]]$data$x, rep(saved$threshold_probability, 3))
  expect_equal(sort(views$plots[["Decision curve"]]$data$y),
    sort(c(saved$model_net_benefit, saved$treat_all_net_benefit, saved$treat_none_net_benefit)))
  expect_identical(run$comparison$results[[1]]$predictions, original)
})

test_that("named external cohorts retain saved calibration and decision curves", {
  run <- clinical_binary_run(c(.1, .3, .5))
  external <- sb_example_binary(100, seed = 17)
  completed <- sb_complete_run(run, "Logistic", list(first = external),
    plots = FALSE, reverse = FALSE
  )
  evaluation <- completed$evaluations$first
  cohorts <- list(first = evaluation, second = evaluation)

  views <- sb_clinician_views(cohorts)

  expect_identical(views$status$status, rep("available", 3))
  expect_setequal(names(views$plots), c("ROC", "Calibration", "Decision curve"))
  curve <- sb_plot_calibration(evaluation)$data
  expect_equal(views$plots$Calibration$data$x, rep(curve$x, 2))
  expect_equal(views$plots$Calibration$data$y, rep(curve$y, 2))
  expect_match(paste(views$plots$Calibration$data$panel, collapse = " "), "first")
  expect_match(paste(views$plots$Calibration$data$panel, collapse = " "), "second")
})

test_that("completed doctor results display external rather than development evidence", {
  run <- clinical_binary_run(c(.1, .3, .5))
  completed <- sb_complete_run(run, "Logistic",
    list(validation = sb_example_binary(100, seed = 17)),
    plots = FALSE, reverse = FALSE
  )
  original <- completed$evaluations$validation$scored_predictions
  local_mocked_bindings(sb_fit_final = function(...) stop("Unexpected refit"))

  views <- sb_clinician_views(completed)

  expect_setequal(names(views$plots), c("ROC", "Calibration", "Decision curve"))
  expect_match(paste(views$plots$Calibration$data$panel, collapse = " "), "External evaluation")
  expect_identical(completed$evaluations$validation$scored_predictions, original)
})

test_that("short survival entrypoint retains the clinical horizon and thresholds", {
  thresholds <- c(.1, .3, .5)
  run <- surv_run(sb_example_survival(120, seed = 42), "time", "event",
    c("gene_a", "gene_b"), "patient_id", models = "Cox",
    folds = 3, plots = FALSE, prediction_horizon = 6,
    decision_thresholds = thresholds
  )

  expect_identical(run$development_task$prediction_horizon, 6)
  expect_identical(run$development_task$decision_thresholds, thresholds)
  expect_identical(run$comparison$results[[1]]$metrics$time_auc$status, "success")
  expect_identical(run$comparison$results[[1]]$metrics$decision_curve$status, "success")
})

test_that("missing external evaluations are not replaced with training evidence", {
  run <- clinical_binary_run()
  completed <- sb_complete_run(run, "Logistic", plots = FALSE)

  result <- tryCatch(sb_clinician_views(completed), error = identity)

  expect_s3_class(result, "stormbreaker_bad_visualization_input")
  expect_match(conditionMessage(result), "No saved external evaluations", fixed = TRUE)
})

test_that("doctor comparison plots show model names instead of parameter IDs", {
  run <- bi_run(sb_example_binary(120, seed = 42), "outcome",
    c("gene_a", "gene_b"), "patient_id", models = c("Logistic", "NaiveBayes"),
    folds = 3
  )
  plot <- run$plots$performance

  expect_identical(plot$labels$title, "Model comparison")
  expect_null(plot$labels$subtitle)
  expect_null(plot$labels$caption)
  expect_setequal(plot$scales$get_scales("x")$get_labels(),
    c("Logistic", "Gaussian NB"))
  expect_setequal(as.character(plot$data$pipeline_id),
    vapply(run$pipelines, function(pipeline) pipeline$id, character(1)))
})
