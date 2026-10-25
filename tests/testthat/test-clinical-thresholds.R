clinical_threshold_fixture <- function() {
  structure(list(status = "success", task_type = "binary", pipeline_id = "clinical",
    scored_predictions = data.frame(patient_id = letters[1:4], outcome = c(1L, 0L, 1L, 0L),
      probability = c(.9, .7, .7, .1), status = "success"), metrics = list()),
    class = "stormbreaker_external_evaluation")
}

test_that("Kappa uses observed margins and preserves negative agreement", {
  x <- clinical_threshold_fixture()
  expect_equal(sb_threshold_metrics(x, c(0, .7, 1))$table$kappa, c(0, .5, 0))
  x$scored_predictions$probability <- x$scored_predictions$outcome
  expect_equal(sb_threshold_metrics(x, .5)$table$kappa, 1)
  x$scored_predictions$probability <- 1 - x$scored_predictions$outcome
  r <- sb_threshold_metrics(x, .5)
  expect_equal(r$table$kappa, -1)
  p <- sb_plot_thresholds(r, metrics = "kappa")
  expect_equal(ggplot2::ggplot_build(p)$data[[2]]$y, -1)
  x$scored_predictions$outcome[] <- 1L
  x$scored_predictions$probability[] <- 1
  expect_equal(sb_threshold_metrics(x, .5)$table$kappa, NA_real_)
})

test_that("clinical thresholds preserve ties and all confusion counts", {
  x <- clinical_threshold_fixture()
  before <- x
  result <- sb_threshold_metrics(x, c(.7, 1, 0, .7))
  expect_equal(result$table$threshold, c(0, .7, 1))
  expect_equal(result$table$tp, c(2, 2, 0))
  expect_equal(result$table$fp, c(2, 1, 0))
  expect_equal(result$table$tn, c(0, 1, 2))
  expect_equal(result$table$fn, c(0, 0, 2))
  expect_equal(result$table$ppv, c(.5, 2 / 3, NA_real_))
  expect_equal(result$table$npv, c(NA_real_, 1, .5))
  expect_identical(x, before)
  expect_length(intersect("patient_id", names(result$table)), 0)
})

test_that("single classes retain defined metrics and invalid rows are not dropped", {
  x <- clinical_threshold_fixture()
  x$scored_predictions$outcome[] <- 1L
  r <- sb_threshold_metrics(x, .5)
  expect_equal(r$table$sensitivity, .75)
  expect_equal(r$table$specificity, NA_real_)
  x$scored_predictions$probability[1] <- NA_real_
  r <- sb_threshold_metrics(x, .5)
  expect_equal(nrow(r$table), 0L)
  expect_identical(r$diagnostics$status, "invalid_predictions")
})

test_that("threshold plotting uses all requested points", {
  skip_if_not_installed("ggplot2")
  p <- sb_plot_thresholds(sb_threshold_metrics(clinical_threshold_fixture(), c(.2, .7, .95)))
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 12L)
  expect_length(unique(p$data$metric), 4L)
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 1L)
})

test_that("confusion cells preserve class orientation and counts", {
  skip_if_not_installed("ggplot2")
  x <- sb_threshold_metrics(clinical_threshold_fixture(), .7)
  p <- sb_plot_confusion(x, .7)
  expect_equal(p$data$count, c(1, 1, 0, 2))
  expect_equal(as.character(p$data$observed), c("0", "0", "1", "1"))
  q <- sb_plot_confusion(x, .7, normalize = "observed")
  expect_equal(q$data$value, c(.5, .5, 0, 1))
  expect_s3_class(ggplot2::ggplot_build(q), "ggplot_built")
})

test_that("different models never share confusion-cell text positions", {
  skip_if_not_installed("ggplot2")
  a <- clinical_threshold_fixture()
  b <- a
  b$pipeline_id <- "second_model"
  b$scored_predictions$probability <- 1 - b$scored_predictions$probability
  input <- structure(list(evaluations = list(a, b)), class = "stormbreaker_external_comparison")
  thresholds <- sb_threshold_metrics(input, .5)
  p <- sb_plot_confusion(thresholds, .5)
  built <- ggplot2::ggplot_build(p)
  labels <- built$data[[2L]]
  expect_equal(nrow(built$layout$layout), 2L)
  expect_equal(nrow(labels), 8L)
  expect_equal(anyDuplicated(labels[c("PANEL", "x", "y")]), 0L)
  expect_equal(sum(p$data$count), 8L)
})
