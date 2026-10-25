test_that("PRE preserves native multiclass coefficients and coefficient-only CV", {
  skip_if_not_installed("pre")
  d <- iris
  d$id <- seq_len(nrow(d))
  f <- sb_fit_multiclass(sb_multiclass_task(d, "id", "Species", names(iris)[1:4]), "pre", 3, 1)
  out <- sb_multiclass_native_plots(f)
  native <- stats::coef(f$engine, penalty.par.val = "lambda.1se")
  expect_equal(nrow(out$tables$coefficients), nrow(native) * 3L)
  for (k in 1:3) expect_equal(out$tables$coefficients$coefficient[out$tables$coefficients$class == f$classes[k]],
    native[[paste0("coefficient.", f$codes[k])]])
  expect_equal(out$tables$coefficient_cv$mean_value, f$engine$glmnet.fit$cvm)
  expect_equal(out$tables$coefficient_cv$standard_error, f$engine$glmnet.fit$cvsd)
  expect_false(any(grepl("feature[0-9]", out$tables$coefficients$description)))
  for (p in sb_plot_pages(out$plots)) expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})
