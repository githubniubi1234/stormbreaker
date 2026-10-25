test_that("tree package plotting retains the complete fitted tree", {
  skip_if_not_installed("tree")
  d <- iris
  d$id <- seq_len(nrow(d))
  fit <- sb_fit_multiclass(sb_multiclass_task(d, "id", "Species", names(iris)[1:4]), "tree", 3, 1)
  out <- sb_multiclass_native_plots(fit)
  tab <- out$tables$tree_structure
  frame <- fit$engine$frame
  expect_equal(tab$node, as.numeric(rownames(frame)))
  expect_equal(tab$n, frame$n)
  expect_equal(tab$prediction, fit$classes[match(as.character(frame$yval), fit$codes)])
  expect_equal(sum(tab$n[tab$terminal]), nrow(d))
  expect_true(all(tab$parent[-1] %in% tab$node))
  expect_false(any(grepl("feature[0-9]", tab$incoming_rule)))
  built <- ggplot2::ggplot_build(out$plots$tree_structure)
  expect_equal(nrow(built$data[[1]]), nrow(frame) - 1L)
  expect_equal(nrow(built$data[[2]]), nrow(frame))
})
