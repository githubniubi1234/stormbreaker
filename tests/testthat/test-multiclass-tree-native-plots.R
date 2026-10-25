test_that("C50 and regularized forest plots use saved native quantities", {
  skip_if_not_installed("C50")
  skip_if_not_installed("RRF")
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  for (method in c("c50", "c50_rules", "rrf")) {
    fit <- sb_fit_multiclass(task, method, 3, 1)
    p <- sb_multiclass_native_plots(fit)
    native <- if (method == "rrf") fit$engine$importance else C50::C5imp(fit$engine, metric = "usage", pct = TRUE)
    expect_equal(p$tables$importance$importance, as.numeric(native[, 1]))
    expect_setequal(as.character(p$tables$importance$feature), fit$features[match(rownames(native), fit$feature_codes)])
    expect_s3_class(ggplot2::ggplot_build(p$plots$importance), "ggplot_built")
    if (method == "rrf") {
      expect_equal(p$tables$oob_error$error, as.vector(fit$engine$err.rate))
      expect_setequal(p$tables$oob_error$series, c("OOB", fit$classes))
      expect_s3_class(ggplot2::ggplot_build(p$plots$oob_error), "ggplot_built")
    }
  }
})
