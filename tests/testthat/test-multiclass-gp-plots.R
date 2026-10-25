test_that("Gaussian-process figures preserve all native pairwise basis weights", {
  skip_if_not_installed("kernlab")
  d <- iris
  d$id <- seq_len(nrow(d))
  f <- sb_fit_multiclass(sb_multiclass_task(d, "id", "Species", names(iris)[1:4]), "gaussian_process", 3, 1)
  out <- sb_multiclass_native_plots(f)
  tab <- out$tables$kernel_weights
  expect_equal(tab$weight, as.numeric(unlist(f$engine@alpha)))
  expect_equal(tab$reference_row, unlist(f$engine@alphaindex, use.names = FALSE))
  pairs <- apply(utils::combn(f$classes, 2), 2, paste, collapse = " vs ")
  expect_equal(unique(tab$pair), pairs)
  expect_equal(nrow(tab), sum(lengths(f$engine@alpha)))
  expect_s3_class(ggplot2::ggplot_build(out$plots$kernel_weights), "ggplot_built")
})
