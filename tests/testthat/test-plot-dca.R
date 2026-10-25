test_that("DCA preserves every threshold and all three saved strategies", {
  skip_if_not_installed("ggplot2")
  for (type in c("binary", "survival")) {
    f <- make_sb_curve_view_fixture(type)
    for (boundary in c("cv", "external")) {
      x <- f[[boundary]]
      before <- x
      p <- sb_plot_dca(x, cohort_label = "Example cohort")
      curve <- x$metrics$decision_curve$curve
      for (strategy in c("model", "treat_all", "treat_none")) {
        rows <- p$data$strategy == strategy
        expect_identical(p$data$x[rows], curve$threshold_probability)
        expect_identical(p$data$y[rows], curve[[paste0(strategy, "_net_benefit")]])
        expect_identical(p$data$source_row[rows], seq_len(nrow(curve)))
      }
      info <- attr(p, "stormbreaker_curve_status")
      expect_identical(info$scope, x$metrics$decision_curve$scope)
      expect_identical(info$cohort, "Example cohort")
      expect_identical(info$metric_status, "success")
      if (type == "survival") {
        expect_identical(info$prediction_horizon, 6)
        expect_identical(info$time_unit, "month")
      }
      built <- ggplot2::ggplot_build(p)
      expect_equal(sort(built$data[[1]]$y), sort(p$data$y))
      expect_equal(sort(built$data[[2]]$x), sort(p$data$x))
      expect_identical(x, before)
    }
  }
})

test_that("DCA retains negative benefits and single-threshold points", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture()$cv
  x$metrics$decision_curve$curve$model_net_benefit[[1]] <- -0.4
  p <- sb_plot_dca(x)
  expect_identical(p$data$y[[1]], -0.4)
  expect_equal(min(ggplot2::ggplot_build(p)$data[[1]]$y), min(p$data$y))
  x$metrics$decision_curve$curve <- x$metrics$decision_curve$curve[1, , drop = FALSE]
  p <- sb_plot_dca(x)
  expect_identical(nrow(p$data), 3L)
  expect_identical(p$scales$get_scales("x")$limits, c(0, 1))
  expect_no_warning(ggplot2::ggplot_build(p))
  expect_identical(nrow(ggplot2::ggplot_build(p)$data[[2]]), 3L)
})

test_that("DCA comparisons keep reference strategies within each pipeline panel", {
  skip_if_not_installed("ggplot2")
  for (type in c("binary", "survival")) {
    f <- make_sb_curve_view_fixture(type, comparison = TRUE)
    for (boundary in c("cv", "external")) {
      x <- f[[boundary]]
      results <- if (boundary == "cv") x$results else x$evaluations
      p <- sb_plot_dca(x)
      info <- attr(p, "stormbreaker_curve_status")
      expect_identical(nrow(info), 2L)
      expect_identical(levels(p$data$panel), info$panel)
      expect_identical(nrow(ggplot2::ggplot_build(p)$layout$layout), 2L)
      for (i in seq_along(results)) {
        rows <- p$data$pipeline_id == info$pipeline_id[[i]] & p$data$strategy == "treat_all"
        expect_identical(p$data$y[rows],
                         results[[i]]$metrics$decision_curve$curve$treat_all_net_benefit)
      }
    }
  }
})

test_that("audit view distinguishes not-requested DCA from failed curves", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_binary(60, seed = 1)
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  x <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1)
  p <- sb_plot_dca(x, show_unavailable = TRUE)
  expect_identical(attr(p, "stormbreaker_curve_status")$metric_status, "not_requested")
  expect_identical(p$scales$get_scales("x")$limits, c(0, 1))
  expect_no_warning(ggplot2::ggplot_build(p))
  expect_identical(nrow(ggplot2::ggplot_build(p)$data[[2]]), 0L)
  x$status <- "failed"
  p <- sb_plot_dca(x, show_unavailable = TRUE)
  expect_identical(attr(p, "stormbreaker_curve_status")$metric_status, "not_saved")

  x <- make_sb_curve_view_fixture("survival")$external
  x$metrics$decision_curve$status <- "incomplete"
  x$metrics$decision_curve$curve$model_net_benefit[[2]] <- NA_real_
  p <- sb_plot_dca(x, show_unavailable = TRUE)
  expect_identical(attr(p, "stormbreaker_curve_status")$metric_status, "incomplete")
  expect_no_warning(ggplot2::ggplot_build(p))
  expect_identical(p$data$y[[2]], NA_real_)
  expect_identical(nrow(ggplot2::ggplot_build(p)$data[[1]]), 0L)
})

test_that("unusable DCA rows do not create an interpolated line across a gap", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("survival")$cv
  x$metrics$decision_curve$curve$status[[3]] <- "not_evaluable"
  p <- sb_plot_dca(x)
  expect_identical(sum(p$data$draw), 12L)
  model_line <- p$layers[[1]]$data[p$layers[[1]]$data$strategy == "model", ]
  expect_identical(length(unique(model_line$segment)), 2L)
  expect_identical(p$data$y[3], x$metrics$decision_curve$curve$model_net_benefit[3])
})

test_that("DCA plotting cannot choose thresholds or recompute clinical utility", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture()$external
  before <- x
  local_mocked_bindings(
    sb_fit_final = function(...) stop("No fitting"),
    sb_evaluate_external = function(...) stop("No evaluation"),
    sb_binary_decision_curve_metric = function(...) stop("No recomputation"),
    sb_survival_decision_curve_metric = function(...) stop("No recomputation")
  )
  expect_s3_class(sb_plot_dca(x), "ggplot")
  expect_identical(x, before)
  expect_identical(names(formals(sb_plot_dca)), c("x", "cohort_label", "show_unavailable",
    "standardized", "cost_benefit_axis"))
})
