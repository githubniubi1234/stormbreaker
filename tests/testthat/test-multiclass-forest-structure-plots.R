test_that("forest structure figures count every original tree", {
  for (package in c("partykit", "ipred", "RLT", "aorsf", "wsrf")) skip_if_not_installed(package)
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  for (method in c("cforest", "bagging", "rlt", "aorsf", "wsrf")) {
    f <- sb_fit_multiclass(task, method, 3, 1)
    p <- sb_multiclass_native_plots(f)
    trees <- switch(method, cforest = length(f$engine$nodes), bagging = length(f$engine$mtrees),
      rlt = length(f$engine$FittedTrees), aorsf = f$engine$n_tree, wsrf = length(f$engine$trees))
    expect_equal(nrow(p$tables$tree_size), trees)
    expect_equal(sum(p$plots$tree_size$data$trees), trees)
    axes <- ggplot2::ggplot_build(p$plots$tree_size)$layout$panel_params[[1]]
    for (axis in c("x", "y")) {
      ticks <- stats::na.omit(axes[[axis]]$get_breaks())
      expect_equal(as.numeric(ticks), floor(as.numeric(ticks)))
    }
    if (method %in% c("cforest", "bagging", "rlt")) {
      expect_equal(sum(p$tables$split_usage$splits), sum(p$tables$tree_size$nodes - 1L))
      expect_setequal(as.character(p$tables$split_usage$feature), f$features)
    }
    if (method == "wsrf") expect_equal(p$tables$individual_tree_oob$error,
      as.numeric(f$engine$treeOOBErrorRates))
    for (plot in p$plots) expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  }
})
