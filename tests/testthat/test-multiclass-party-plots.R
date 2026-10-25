test_that("multiclass conditional trees retain original nodes and memberships", {
  skip_if_not_installed("partykit")
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  f <- sb_fit_multiclass(task, "ctree", 3, 1)
  p <- sb_multiclass_native_plots(f)
  tab <- p$tables$tree_structure
  expect_setequal(tab$node, partykit::nodeids(f$engine))
  expect_equal(sum(tab$n[tab$terminal]), nrow(d))
  expect_equal(tab$n[is.na(tab$parent)], nrow(d))
  expect_length(grep("feature[0-9]", tab$incoming_rule), 0L)
  expect_s3_class(ggplot2::ggplot_build(p$plots$tree_structure), "ggplot_built")
})
