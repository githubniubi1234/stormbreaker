test_that("GAM smooths match native prediction matrices and retain training rugs", {
  skip_if_not_installed("mgcv")
  d <- iris
  d$id <- seq_len(nrow(d))
  set.seed(27)
  d$Species <- sample(d$Species)
  f <- sb_fit_multiclass(sb_multiclass_task(d, "id", "Species", names(iris)[1:2]), "gam", 3, 1)
  out <- sb_multiclass_native_plots(f)
  expect_equal(nrow(out$tables$training_rug), nrow(d) * 4L)
  expect_setequal(out$tables$smooths$class, levels(d$Species)[2:3])
  s <- f$engine$smooth[[1]]
  p <- s$first.para:s$last.para
  tab <- out$tables$smooths[1:100, ]
  nd <- f$engine$model[rep(1, 100), , drop = FALSE]
  nd[[s$term]] <- tab$value
  X <- stats::predict(f$engine, nd, type = "lpmatrix")[, p, drop = FALSE]
  expect_equal(tab$estimate, as.numeric(X %*% stats::coef(f$engine)[p]), ignore_attr = TRUE)
  expect_equal(tab$se, sqrt(diag(X %*% f$engine$Vp[p, p] %*% t(X))), ignore_attr = TRUE)
  pages <- sb_plot_pages(out$plots)
  for (plot in pages) expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
})
