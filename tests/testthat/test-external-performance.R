test_that("external performance preserves metric values and cohort identity", {
  skip_if_not_installed("ggplot2")
  for (type in c("binary", "survival")) {
    x <- make_sb_curve_view_fixture(type, comparison = TRUE)$external
    metrics <- if (type == "binary") c("roc_auc", "average_precision", "brier_score", "log_loss") else "harrell_c"
    for (metric in metrics) {
      p <- sb_plot_performance(x, metric, cohort_label = "Independent cohort")
      expect_identical(p$data$estimate, sb_result_table(x, metric = metric)$estimate)
      expect_identical(levels(p$data$pipeline_id), names(x$evaluations))
      expect_identical(p$labels$title, "External performance")
      expect_match(p$labels$subtitle, "Independent cohort", fixed = TRUE)
      expect_match(p$labels$subtitle, x$evaluations[[1]]$metrics[[metric]]$aggregation, fixed = TRUE)
      expect_length(p$layers, 1L)
      expect_s3_class(p$layers[[1]]$geom, "GeomPoint")
      expect_equal(ggplot2::ggplot_build(p)$data[[1]]$y, p$data$estimate)
    }
  }
})

test_that("external failures and nonfinite values are retained but not drawn", {
  skip_if_not_installed("ggplot2")
  a <- make_sb_curve_view_fixture("binary")$external
  b <- a
  b$pipeline_id <- "failed"
  b$status <- "failed"
  b$metrics <- list()
  c <- a
  c$pipeline_id <- "infinite"
  c$metrics$log_loss$estimate <- Inf
  x <- structure(list(evaluations = list(b, a, c)), class = "stormbreaker_external_comparison")
  p <- sb_plot_performance(x, "log_loss")
  expect_identical(p$data$estimate, c(NA_real_, a$metrics$log_loss$estimate, Inf))
  expect_identical(nrow(p$layers[[1]]$data), 1L)
  expect_identical(levels(p$data$pipeline_id), c("failed", a$pipeline_id, "infinite"))
  expect_match(p$labels$caption, "2 pipeline", fixed = TRUE)
  x$evaluations <- list(b, c)
  empty <- sb_plot_performance(x, "log_loss", show_unavailable = TRUE)
  expect_identical(nrow(ggplot2::ggplot_build(empty)$data[[1]]), 0L)
  expect_identical(nrow(empty$data), 2L)
})

test_that("external performance never refits or evaluates and CV remains unchanged", {
  skip_if_not_installed("ggplot2")
  fixture <- make_sb_curve_view_fixture("survival")
  x <- fixture$external
  before <- x
  local_mocked_bindings(
    sb_fit_final = function(...) stop("must not fit"),
    sb_evaluate_external = function(...) stop("must not evaluate"),
    sb_external_harrell_c_metric = function(...) stop("must not compute")
  )
  p <- sb_plot_performance(x, "harrell_c")
  expect_identical(p$data$estimate, x$metrics$harrell_c$estimate)
  expect_identical(x, before)
  cv <- sb_plot_performance(fixture$cv, "harrell_c")
  expect_identical(cv$labels$title, "CV performance")
  expect_identical(names(sb_result_table(fixture$cv, metric = "harrell_c")), c("pipeline_id", "task_type", "metric", "estimate", "status", "run_status", "scope", "aggregation"))
})
