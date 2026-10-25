make_sb_heatmap_fixture <- function(id = "model A", estimate = 0.75, status = "success") {
  structure(list(status = status, task_type = "binary", pipeline_id = id,
    patient_match = list(status = "success"), error_stage = NA_character_,
    error_class = NA_character_, error_message = NA_character_,
    metrics = list(roc_auc = list(estimate = estimate, status = "success",
      scope = "external_cohort", aggregation = "patient_pairs"))),
    class = "stormbreaker_external_evaluation")
}

test_that("performance heatmaps preserve named cohorts, model order and missing cells", {
  skip_if_not_installed("ggplot2")
  a <- make_sb_heatmap_fixture()
  b <- make_sb_heatmap_fixture("model B", 0.6)
  comparison <- structure(list(evaluations = list(a, b)), class = "stormbreaker_external_comparison")
  x <- list(Cohort1 = comparison, Cohort2 = a)
  before <- x
  p <- sb_plot_performance_heatmap(x, "roc_auc")
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 4L)
  expect_identical(levels(p$data$pipeline_id), c("model B", "model A"))
  expect_equal(sum(p$data$cell_status == "not_supplied"), 1L)
  expect_equal(sort(p$data$estimate[is.finite(p$data$estimate)]), c(0.6, 0.75, 0.75))
  expect_identical(x, before)
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 4L)
})

test_that("failed estimates are retained but never painted as success", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_heatmap_fixture(status = "failed")
  p <- sb_plot_performance_heatmap(list(Cohort = x), "roc_auc")
  expect_equal(p$data$estimate, 0.75)
  expect_identical(p$data$fill_value, NA_real_)
  expect_identical(p$data$cell_status, "failed_run")
  expect_match(p$labels$subtitle, "No usable", fixed = TRUE)
  expect_error(sb_plot_performance_heatmap(list(a = x, a = x), "roc_auc"),
    class = "stormbreaker_bad_visualization_input")
})

test_that("time targets and validation scopes are not conflated", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_survival(45, seed = 14)
  results <- lapply(c(3, 6), function(horizon) {
    task <- sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month",
      prediction_horizon = horizon)
    sb_run_cv(task, folds = 3, seed = 9)
  })
  p <- sb_plot_performance_heatmap(list(ThreeMonths = results[[1]], SixMonths = results[[2]]), "time_auc")
  expect_equal(length(unique(p$data$target_label)), 2L)
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 2L)
  expect_identical(sort(unique(p$data$prediction_horizon)), c(3, 6))
})
