test_that("specialist figures preserve native means, spectra and loadings", {
  for (package in c("mclust", "HDclassif", "hda")) skip_if_not_installed(package)
  d <- withr::with_seed(701, data.frame(id = 1:90,
    x = rnorm(90), y = rnorm(90), z = rnorm(90), w = rnorm(90),
    label = rep(c("A", "B", "C"), 30)))
  task <- sb_multiclass_task(d, "id", "label", c("x", "y", "z", "w"))
  for (method in c("edda", "mclustda", "hdda", "hda")) {
    f <- sb_fit_multiclass(task, method, 3, 1)
    p <- sb_multiclass_native_plots(f)
    if (method %in% c("edda", "mclustda")) {
      expect_equal(p$tables$component_means$value,
        unlist(lapply(f$engine$models, function(x) as.vector(x$parameters$mean)), use.names = FALSE))
      expect_setequal(p$tables$component_means$class, f$classes)
    } else if (method == "hdda") {
      expect_equal(p$tables$eigenvalue_spectrum$value, as.vector(f$engine$ev))
      expect_equal(sum(p$tables$eigenvalue_spectrum$retained), sum(f$engine$d))
    } else {
      expect_equal(p$tables$projection_loadings$value, as.vector(f$engine$hda.loadings))
      expect_equal(sum(p$tables$projection_loadings$retained), length(f$features) * f$engine$reduced.dimension)
    }
    for (plot in p$plots) expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  }
})
