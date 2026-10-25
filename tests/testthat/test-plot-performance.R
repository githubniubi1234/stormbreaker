test_that("performance points use saved estimates in input order", {
  skip_if_not_installed("ggplot2")
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  x <- sb_compare_cv(task, list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), folds = 3, seed = 1)
  before <- x
  p <- sb_plot_performance(x, "roc_auc")
  expect_s3_class(p, "ggplot")
  expect_identical(p$data$estimate, sb_result_table(x, metric = "roc_auc")$estimate)
  expect_identical(levels(p$data$pipeline_id), names(x$results))
  expect_identical(p$scales$get_scales("x")$limits, rev(names(x$results)))
  expect_length(p$layers, 1L)
  expect_s3_class(p$layers[[1]]$geom, "GeomPoint")
  expect_match(p$labels$subtitle, x$results[[1]]$metrics$roc_auc$aggregation,
               fixed = TRUE)
  expect_match(p$labels$caption, "No confidence intervals", fixed = TRUE)
  expect_identical(x, before)
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$data[[1]]$y, p$data$estimate)
})

test_that("failed runs and unavailable metrics remain visible without zero filling", {
  skip_if_not_installed("ggplot2")
  d <- make_sb_fixture()
  task <- sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month")
  a <- sb_run_cv(task, folds = 3, seed = 1)
  b <- a
  b$pipeline$id <- "failed_run"
  b$status <- "failed"
  c <- a
  c$pipeline$id <- "unavailable_metric"
  c$metrics$harrell_c$status <- "not_evaluable"
  c$metrics$harrell_c$estimate <- NA_real_
  x <- structure(list(results = list(b, a, c)),
                 class = "stormbreaker_cv_comparison")
  p <- sb_plot_performance(x, "harrell_c")
  expect_identical(nrow(p$data), 3L)
  expect_identical(nrow(p$layers[[1]]$data), 1L)
  expect_identical(p$data$estimate,
                   c(a$metrics$harrell_c$estimate, a$metrics$harrell_c$estimate, NA_real_))
  expect_identical(levels(p$data$pipeline_id), c("failed_run", a$pipeline$id,
                                              "unavailable_metric"))
  expect_match(p$labels$caption, "2 pipeline", fixed = TRUE)
  x$results <- list(b, c)
  empty <- sb_plot_performance(x, "harrell_c", show_unavailable = TRUE)
  expect_identical(nrow(empty$layers[[1]]$data), 0L)
  expect_identical(nrow(empty$data), 2L)
  expect_identical(nrow(ggplot2::ggplot_build(empty)$data[[1]]), 0L)
})

test_that("plotting keeps infinite loss in its table but does not draw it", {
  skip_if_not_installed("ggplot2")
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  x <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1)
  x$metrics$log_loss$estimate <- Inf
  p <- sb_plot_performance(x, "log_loss", show_unavailable = TRUE)
  expect_identical(p$data$estimate, Inf)
  expect_identical(nrow(p$layers[[1]]$data), 0L)
  expect_match(p$labels$caption, "1 pipeline", fixed = TRUE)
})

test_that("plotting is read-only and works for each supported binary metric", {
  skip_if_not_installed("ggplot2")
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  x <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1)
  before <- x
  local_mocked_bindings(
    sb_run_cv = function(...) stop("must not run CV"),
    sb_fit_final = function(...) stop("must not fit"),
    sb_binary_auc_metric = function(...) stop("must not recompute")
  )
  for (metric in c("roc_auc", "average_precision", "brier_score", "log_loss")) {
    p <- sb_plot_performance(x, metric)
    expect_identical(p$data$estimate, x$metrics[[metric]]$estimate)
    expect_identical(nrow(ggplot2::ggplot_build(p)$data[[1]]), 1L)
  }
  expect_identical(x, before)
})
