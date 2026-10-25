test_that("PAMR displays native shrunken contrasts without dropping genes", {
  skip_if_not_installed("pamr")
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  fit <- sb_fit_multiclass(task, "pamr", 3, 1)
  result <- sb_multiclass_native_plots(fit)
  native <- pamr::pamr.predict(fit$engine, newx = NULL,
    threshold = fit$parameters$threshold, type = "centroid")
  expected <- (native - fit$engine$centroid.overall) / fit$engine$sd
  tab <- result$tables$shrunken_centroids
  expect_equal(tab$contrast, as.vector(expected))
  expect_equal(nrow(tab), 12L)
  expect_setequal(as.character(tab$feature), names(iris)[1:4])
  expect_setequal(tab$class, levels(iris$Species))
  expect_s3_class(ggplot2::ggplot_build(result$plots$shrunken_centroids), "ggplot_built")
  fit$parameters$threshold <- 1e6
  zero <- sb_multiclass_native_plots(fit)$tables$shrunken_centroids
  expect_equal(nrow(zero), 12L)
  expect_true(all(zero$contrast == 0))
})
