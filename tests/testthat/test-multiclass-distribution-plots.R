test_that("native distribution plots preserve densities and fitted means", {
  skip_if_not_installed("klaR")
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  for (method in c("kernel_naive_bayes", "rda")) {
    f <- sb_fit_multiclass(task, method, 3, 1)
    p <- sb_multiclass_native_plots(f)
    if (method == "kernel_naive_bayes") {
      expected <- unlist(lapply(f$engine$tables, function(t) unlist(lapply(t, `[[`, "y"))), use.names = FALSE)
      expect_equal(p$tables$class_densities$density, expected)
      expect_equal(nrow(p$tables$class_densities), 4L * 3L * 512L)
    } else {
      expect_equal(p$tables$class_means$value, as.vector(f$engine$means))
      co <- p$tables$class_correlations
      expect_equal(nrow(co), 4L * 4L * 3L)
      expect_equal(co$correlation[co$feature == co$other], rep(1, 12))
      expect_setequal(co$class, f$classes)
    }
    for (plot in p$plots) expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  }
})
