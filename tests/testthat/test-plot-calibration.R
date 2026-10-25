test_that("binary and survival calibration coordinates are copied exactly", {
  skip_if_not_installed("ggplot2")
  for (type in c("binary", "survival")) {
    fixture <- make_sb_curve_view_fixture(type)
    for (boundary in c("cv", "external")) {
      x <- fixture[[boundary]]
      before <- x
      p <- sb_plot_calibration(x, cohort_label = "Demonstration cohort")
      key <- if (type == "binary") "moderate_calibration" else "calibration"
      metric <- x$metrics[[key]]
      x_column <- if (type == "binary") "probability" else "predicted_event_probability"
      y_column <- if (type == "binary") "calibrated_probability" else "calibrated_event_probability"
      expect_s3_class(p, "ggplot")
      expect_null(p$labels$caption)
      expect_null(p$labels$subtitle)
      expect_identical(p$data$x, metric$curve[[x_column]])
      expect_identical(p$data$y, metric$curve[[y_column]])
      expect_identical(p$data$source_row, seq_len(nrow(metric$curve)))
      info <- attr(p, "stormbreaker_curve_status")
      expect_identical(info$scope, metric$scope)
      expect_identical(info$metric_status, metric$status)
      expect_identical(info$run_status, x$status)
      expect_identical(info$cohort, "Demonstration cohort")
      expect_identical(info$validation,
                       if (boundary == "cv") "development_cv" else "external_evaluation")
      if (type == "survival") {
        expect_identical(info$prediction_horizon, 6)
        expect_identical(info$time_unit, "month")
        expect_match(as.character(p$data$panel[[1]]), "6 month", fixed = TRUE)
        expect_match(p$labels$x, "event probability", fixed = TRUE)
      }
      built <- ggplot2::ggplot_build(p)
      expect_equal(built$data[[2]]$x, p$data$x)
      expect_equal(built$data[[2]]$y, p$data$y)
      expect_identical(x, before)
    }
  }
})

test_that("both comparison classes preserve panels and all original curves", {
  skip_if_not_installed("ggplot2")
  for (type in c("binary", "survival")) {
    f <- make_sb_curve_view_fixture(type, comparison = TRUE)
    for (boundary in c("cv", "external")) {
      x <- f[[boundary]]
      results <- if (boundary == "cv") x$results else x$evaluations
      ids <- vapply(results, function(z) {
        if (boundary == "cv") z$pipeline$id else z$pipeline_id
      }, character(1))
      p <- sb_plot_calibration(x)
      info <- attr(p, "stormbreaker_curve_status")
      expect_identical(info$pipeline_id, unname(ids))
      expect_identical(levels(p$data$panel), info$panel)
      expect_identical(nrow(ggplot2::ggplot_build(p)$layout$layout), 2L)
      key <- if (type == "binary") "moderate_calibration" else "calibration"
      column <- if (type == "binary") "probability" else "predicted_event_probability"
      for (i in seq_along(results)) {
        expect_identical(p$data$x[p$data$pipeline_id == ids[[i]]],
                         results[[i]]$metrics[[key]]$curve[[column]])
      }
    }
  }
})

test_that("explicit audit view keeps absent and failed calibration panels", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_survival(60, seed = 1)
  task <- sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month")
  x <- sb_run_cv(task, folds = 3, seed = 1)
  p <- sb_plot_calibration(x, show_unavailable = TRUE)
  expect_identical(attr(p, "stormbreaker_curve_status")$metric_status, "not_requested")
  expect_identical(nrow(ggplot2::ggplot_build(p)$layout$layout), 1L)
  expect_identical(nrow(ggplot2::ggplot_build(p)$data[[2]]), 0L)
  expect_identical(p$data$source_row, NA_integer_)

  x <- make_sb_curve_view_fixture()$cv
  original <- x$metrics$moderate_calibration$curve$calibrated_probability
  x$status <- "failed"
  p <- sb_plot_calibration(x, show_unavailable = TRUE)
  expect_identical(p$data$y, original)
  expect_identical(attr(p, "stormbreaker_curve_status")$run_status, "failed")
  expect_identical(nrow(ggplot2::ggplot_build(p)$data[[2]]), 0L)
  expect_match(p$layers[[3]]$data$label, "No curve", fixed = TRUE)
})

test_that("calibration drawing never refits a model or calibration equation", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture()$external
  before <- x
  local_mocked_bindings(
    sb_fit_final = function(...) stop("No fitting"),
    sb_evaluate_external = function(...) stop("No evaluation"),
    sb_binary_moderate_calibration_metric = function(...) stop("No calibration"),
    sb_survival_calibration_metric = function(...) stop("No calibration")
  )
  expect_s3_class(sb_plot_calibration(x), "ggplot")
  expect_identical(x, before)
})

test_that("curve plotting rejects unsupported containers and invalid cohort labels", {
  skip_if_not_installed("ggplot2")
  expect_error(sb_plot_calibration(list()), class = "stormbreaker_bad_curve_view")
  expect_error(sb_plot_calibration(structure(list(results = list()),
                                            class = "stormbreaker_cv_comparison")),
               class = "stormbreaker_bad_curve_view")
  x <- make_sb_curve_view_fixture()$cv
  for (label in list(NA_character_, character(), c("a", "b"), "", 1)) {
    expect_error(sb_plot_calibration(x, cohort_label = label),
                 class = "stormbreaker_bad_cohort_label")
  }
})
