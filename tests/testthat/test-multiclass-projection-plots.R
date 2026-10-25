test_that("projection plots retain native coordinates and component identities", {
  for (package in c("mda", "ipred", "sda")) skip_if_not_installed(package)
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  for (method in c("fda", "mda", "pda", "slda", "sda")) {
    f <- sb_fit_multiclass(task, method, 3, 1)
    p <- sb_multiclass_native_plots(f)
    if (method %in% c("fda", "mda", "pda")) {
      expect_equal(p$tables$projected_centroids$value, as.vector(f$engine$means))
      expect_length(unique(p$tables$projected_centroids$label), nrow(f$engine$means))
    } else if (method == "slda") {
      expect_equal(p$tables$projection_loadings$value, as.vector(f$engine$scores))
      expect_equal(p$tables$projected_centroids$value, as.vector(f$engine$mylda$means))
    } else {
      expect_equal(p$tables$discriminant_coefficients$value, as.vector(t(f$engine$beta)))
      expect_equal(p$tables$shrinkage$intensity, unname(f$engine$regularization))
    }
    if (method == "mda") expect_equal(p$tables$subclass_weights$probability,
      unname(unlist(f$engine$sub.prior)))
    for (plot in p$plots) expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  }
})
