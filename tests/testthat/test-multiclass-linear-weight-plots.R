test_that("linear multiclass weight plots retain native parameters", {
  skip_if_not_installed("LiblineaR")
  skip_if_not_installed("xgboost")
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  for (method in c("liblinear_l2_logistic", "liblinear_l1_logistic", "xgboost_linear")) {
    f <- sb_fit_multiclass(task, method, 3, 1)
    p <- sb_multiclass_native_plots(f)
    tab <- p$tables$linear_weights
    expect_equal(nrow(tab), (length(f$features) + 1L) * length(f$classes))
    expect_setequal(tab$class, f$classes)
    if (method != "xgboost_linear") expect_equal(tab$weight, as.vector(f$engine$W)) else {
      b <- sb_load_xgboost_state_model(list(model_raw = f$engine, params = f$parameters),
        "xgboost_binary_linear")$booster
      x <- as.matrix(d[1:4])
      colnames(x) <- f$feature_codes
      native <- predict(b, x, outputmargin = TRUE, reshape = TRUE)
      beta <- matrix(tab$weight, nrow = 3)
      reconstructed <- cbind(x, 1) %*% t(beta)
      expect_equal(sweep(reconstructed, 1, reconstructed[, 1]),
        sweep(native, 1, native[, 1]), tolerance = 1e-5, ignore_attr = TRUE)
    }
    expect_s3_class(ggplot2::ggplot_build(p$plots$linear_weights), "ggplot_built")
  }
})
