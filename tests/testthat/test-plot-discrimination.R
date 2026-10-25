make_sb_discrimination_fixture <- function(external = FALSE) {
  predictions <- data.frame(patient_id = paste0("p", 1:4), outcome = c(1L, 0L, 1L, 0L),
    probability = c(0.9, 0.7, 0.7, 0.1), status = "success", stringsAsFactors = FALSE)
  metrics <- list(
    roc_auc = list(estimate = 0.875, status = "success", scope = "fixture", aggregation = "fixture"),
    average_precision = list(estimate = 5/6, status = "success", scope = "fixture", aggregation = "fixture"))
  if (external) return(structure(list(status = "success", task_type = "binary",
    pipeline_id = "test", scored_predictions = predictions, metrics = metrics),
    class = "stormbreaker_external_evaluation"))
  predictions$fold <- 1L
  predictions$pipeline_id <- "test"
  structure(list(status = "success", task_type = "binary", pipeline = list(id = "test"),
    predictions = predictions, metrics = metrics), class = "stormbreaker_binary_cv_result")
}

test_that("ROC uses tied thresholds and preserves saved estimands", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_discrimination_fixture()
  before <- x
  p <- sb_plot_roc(x)
  expect_s3_class(p, "ggplot")
  expect_equal(p$data$x, c(0, 0, 0.5, 1))
  expect_equal(p$data$y, c(0, 0.5, 1, 1))
  expect_equal(sum(diff(p$data$x) * (head(p$data$y, -1) + tail(p$data$y, -1))/2), 0.875)
  expect_identical(attr(p, "stormbreaker_saved_metrics")[[1]]$estimate, 0.875)
  expect_identical(x, before)
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 1L)
})

test_that("PR uses noninterpolated AP step geometry", {
  skip_if_not_installed("ggplot2")
  p <- sb_plot_pr(make_sb_discrimination_fixture(external = TRUE))
  expect_equal(p$data$x, c(0, 0.5, 1, 1))
  expect_equal(p$data$y, c(1, 1, 2/3, 0.5))
  expect_equal(sum(diff(p$data$x) * tail(p$data$y, -1)), 5/6)
  expect_match(p$labels$caption, "not pooled", fixed = TRUE)
  expect_identical(p$layers[[1]]$geom_params$direction, "vh")
})

test_that("optional per-fold ROC preserves folds and invalid-row diagnostics", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_discrimination_fixture()
  second <- x$predictions
  second$patient_id <- paste0("q", 1:4)
  second$fold <- 2L
  second$probability <- rev(second$probability)
  x$predictions <- rbind(x$predictions, second)
  p <- sb_plot_roc(x, by_fold = TRUE)
  expect_equal(length(unique(p$data$fold)), 2L)
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 2L)
  x$predictions$probability[1] <- NA_real_
  p <- sb_plot_roc(x, by_fold = TRUE)
  status <- attr(p, "stormbreaker_curve_status")
  expect_identical(status$curve_status, c("invalid_predictions", "success"))
  expect_identical(unique(p$data$fold), 2L)
  expect_match(attr(p, "stormbreaker_figure_notes"), "1 unavailable", fixed = TRUE)
})

test_that("failed and single-class curves remain explicit", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_discrimination_fixture()
  x$status <- "failed"
  p <- sb_plot_roc(x)
  expect_identical(attr(p, "stormbreaker_curve_status")$curve_status, "failed_run")
  expect_match(p$labels$subtitle, "No usable", fixed = TRUE)
  x$status <- "success"
  x$predictions$outcome <- 1L
  expect_identical(attr(sb_plot_roc(x), "stormbreaker_curve_status")$curve_status, "single_class")
  expect_equal(tail(sb_plot_pr(x)$data$y, 1), 1)
})

test_that("named cohorts and composite repetitions retain their boundaries", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_discrimination_fixture()
  comparison <- structure(list(results = list(x)), class = "stormbreaker_cv_comparison")
  repeated <- structure(list(repetitions = list(comparison, comparison)),
    class = "stormbreaker_repeated_cv_comparison")
  p <- sb_plot_roc(list(Development = repeated, External = make_sb_discrimination_fixture(TRUE)))
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 3L)
  expect_identical(unique(p$data$repeat_id), c(1L, 2L, NA_integer_))
  expect_equal(length(unique(p$data$series_id)), 3L)
  expect_error(sb_plot_roc(list(x, x)), class = "stormbreaker_bad_visualization_input")
})
