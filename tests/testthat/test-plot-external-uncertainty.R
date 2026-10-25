make_sb_plot_bootstrap <- function(id = "model", metric = "roc_auc") {
  structure(list(status = "success", pipeline_id = id, task_type = "binary",
    confidence_level = 0.95, prediction_target = list(has_horizon = FALSE),
    intervals = data.frame(metric = metric, point_estimate = 0.7, lower = 0.6, upper = 0.8,
      confidence_level = 0.95, status = "success", interval_available = TRUE,
      replicate_completeness = "complete", n_requested = 100L, n_finite_successful = 100L)),
    class = "stormbreaker_external_bootstrap")
}

test_that("external interval display copies endpoints and bootstrap completeness", {
  skip_if_not_installed("ggplot2")
  a <- make_sb_plot_bootstrap()
  b <- a
  b$status <- "incomplete"
  b$intervals$replicate_completeness <- "partial"
  b$intervals$n_finite_successful <- 80L
  b$intervals$status <- "incomplete"
  before <- b
  p <- sb_plot_external_intervals(list(A = a, B = b), "roc_auc")
  expect_s3_class(p, "ggplot")
  expect_identical(p$data$lower, c(0.6, 0.6))
  expect_identical(p$data$upper, c(0.8, 0.8))
  expect_match(as.character(p$data$row_label[2]), "80/100", fixed = TRUE)
  expect_identical(b, before)
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 2L)
})

test_that("paired differences retain reference direction and failed candidates", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_plot_bootstrap()
  class(x) <- "stormbreaker_external_bootstrap_comparison"
  x$reference_pipeline_id <- "reference"
  x$intervals$candidate_pipeline_id <- "candidate"
  x$intervals$reference_pipeline_id <- "reference"
  x$intervals$metric_direction <- "higher_is_better"
  x$candidate_status <- data.frame(candidate_pipeline_id = c("candidate", "broken"), status = c("success", "failed"))
  p <- sb_plot_external_difference(x, "roc_auc")
  expect_equal(nrow(p$data), 2L)
  expect_identical(p$data$reference_pipeline_id, c("reference", "reference"))
  expect_identical(p$data$status, c("success", "failed"))
  expect_match(p$labels$x, "candidate - reference", fixed = TRUE)
  expect_equal(ggplot2::ggplot_build(p)$data[[3]]$xintercept, 0)
})

test_that("forest plots distinguish saved cohort CI, pooled CI and prediction interval", {
  skip_if_not_installed("ggplot2")
  a <- make_sb_plot_bootstrap()
  b <- a
  b$intervals <- NULL
  x <- structure(list(status = "success", task_type = "binary", metric = "roc_auc",
    confidence_level = 0.95, pipeline_id = "model", metric_target = list(has_horizon = FALSE),
    cohorts = data.frame(cohort_id = c("A", "B"), point_estimate = c(0.7, 0.72)),
    summary = data.frame(estimate = 0.71, lower = 0.63, upper = 0.78,
      prediction_lower = 0.55, prediction_upper = 0.83),
    bootstraps = list(A = a, B = b), prediction_interval_caution = "Only two cohorts."),
    class = "stormbreaker_external_discrimination_meta")
  before <- x
  p <- sb_plot_meta_forest(x)
  expect_equal(p$data$point_estimate, c(0.7, 0.72, 0.71, 0.71))
  expect_equal(p$data$lower, c(0.6, NA, 0.63, 0.55))
  expect_identical(p$data$interval_type, c("Cohort bootstrap CI", "Cohort bootstrap CI", "Pooled CI", "Prediction interval"))
  expect_identical(p$data$status[2], "interval_not_saved")
  expect_match(p$labels$caption, "Only two cohorts", fixed = TRUE)
  expect_identical(x, before)
  expect_equal(nrow(sb_plot_meta_forest(x, show_prediction_interval = FALSE)$data), 3L)
})

test_that("real external bootstrap output can be plotted without re-evaluation", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_binary(60, seed = 17)
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  evaluation <- sb_evaluate_external(sb_predict_final(fit, d, "patient_id"),
    sb_binary_endpoint(d, "patient_id", "outcome"))
  b <- sb_bootstrap_external(evaluation, reps = 10, seed = 9)
  local_mocked_bindings(sb_bootstrap_external = function(...) stop("Do not bootstrap"),
    sb_evaluate_external = function(...) stop("Do not evaluate"))
  p <- sb_plot_external_intervals(b, "brier_score")
  expected <- b$intervals[b$intervals$metric == "brier_score", ]
  expect_identical(p$data$lower, expected$lower)
  expect_identical(p$data$point_estimate, expected$point_estimate)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

make_sb_plot_missing_bootstrap <- function(id = "missing", metric = "roc_auc") {
  x <- make_sb_plot_bootstrap(id, metric)
  x$status <- "incomplete"
  x$intervals[c("point_estimate", "lower", "upper")] <- NA_real_
  x$intervals$status <- "incomplete"
  x$intervals$interval_available <- FALSE
  x$intervals$replicate_completeness <- "none"
  x$intervals$n_finite_successful <- 0L
  x
}

test_that("all unavailable intervals retain visible cohort and status rows", {
  skip_if_not_installed("ggplot2")
  withr::local_pdf(NULL)
  inputs <- list(Alpha = make_sb_plot_missing_bootstrap("first"),
    Beta = make_sb_plot_missing_bootstrap("second"))
  before <- inputs
  p <- sb_plot_external_intervals(inputs, "roc_auc")
  built <- ggplot2::ggplot_build(p)
  labels <- as.character(built$layout$panel_params[[1]]$y$get_labels())
  expect_setequal(labels, sb_view_wrap_label(as.character(p$data$row_label)))
  expect_length(labels, 2L)
  expect_equal(sum(grepl("0/100", labels, fixed = TRUE)), 2L)
  expect_equal(sum(grepl("incomplete", labels, fixed = TRUE)), 2L)
  expect_equal(nrow(built$data[[1]]), 0L)
  expect_equal(nrow(built$data[[2]]), 0L)
  expect_identical(p$data$point_estimate, c(NA_real_, NA_real_))
  expect_identical(inputs, before)
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})

test_that("unavailable interval rows stay in their own target panels", {
  skip_if_not_installed("ggplot2")
  good <- make_sb_plot_bootstrap("available", "calibration_slope")
  missing <- make_sb_plot_missing_bootstrap("unavailable", "calibration_slope")
  survival <- missing
  survival$task_type <- "survival"
  survival$prediction_target <- list(has_horizon = TRUE,
    prediction_horizon = 24, time_unit = "month")
  p <- sb_plot_external_intervals(list(Binary = good,
    Binary_missing = missing, Survival_missing = survival), "calibration_slope")
  built <- ggplot2::ggplot_build(p)
  panels <- built$layout$layout
  expect_equal(nrow(panels), 2L)
  for (i in seq_len(nrow(panels))) {
    rows <- as.character(p$data$target_label) == as.character(panels$target_label[i])
    labels <- as.character(built$layout$panel_params[[i]]$y$get_labels())
    expect_setequal(labels, sb_view_wrap_label(as.character(p$data$row_label[rows])))
  }
  available <- ggplot2::ggplot_build(sb_plot_external_intervals(good, "calibration_slope"))
  expect_equal(built$layout$panel_params[[1]]$x.range,
    available$layout$panel_params[[1]]$x.range)
  expect_equal(built$data[[1]]$x, 0.6)
  expect_equal(built$data[[1]]$xend, 0.8)
  expect_equal(built$data[[2]]$x, 0.7)
})

test_that("all failed paired candidates keep labels and only the reference line", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_plot_missing_bootstrap()
  class(x) <- "stormbreaker_external_bootstrap_comparison"
  x$reference_pipeline_id <- "reference"
  x$intervals$candidate_pipeline_id <- "first"
  x$intervals$reference_pipeline_id <- "reference"
  x$intervals$status <- "failed"
  x$candidate_status <- data.frame(candidate_pipeline_id = c("first", "second"),
    status = c("failed", "failed"))
  p <- sb_plot_external_difference(x, "roc_auc")
  built <- ggplot2::ggplot_build(p)
  labels <- as.character(built$layout$panel_params[[1]]$y$get_labels())
  expect_setequal(labels, sb_view_wrap_label(as.character(p$data$row_label)))
  expect_length(labels, 2L)
  expect_equal(sum(grepl("Reference: reference", labels, fixed = TRUE)), 2L)
  expect_equal(nrow(built$data[[1]]), 0L)
  expect_equal(nrow(built$data[[2]]), 0L)
  expect_equal(built$data[[3]]$xintercept, 0)
})
