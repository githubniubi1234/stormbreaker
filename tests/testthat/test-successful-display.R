test_that("default performance display omits failed labels while preserving source", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("binary",TRUE)$external
  x$evaluations[[1]]$status <- "failed"
  before <- x
  p <- sb_plot_performance(x,"roc_auc")
  expect_identical(p$scales$get_scales("x")$limits,x$evaluations[[2]]$pipeline_id)
  expect_identical(x,before)
  audit <- sb_plot_performance(x,"roc_auc",show_unavailable=TRUE)
  expect_identical(audit$scales$get_scales("x")$limits,rev(names(x$evaluations)))
})

test_that("default calibration and DCA omit failed facets but retain diagnostics", {
  skip_if_not_installed("ggplot2")
  for(type in c("binary","survival")) {
    x <- make_sb_curve_view_fixture(type,TRUE)$external
    x$evaluations[[1]]$status <- "failed"
    before <- x
    for(fun in list(sb_plot_calibration,sb_plot_dca)) {
      p <- fun(x)
      expect_identical(nrow(ggplot2::ggplot_build(p)$layout$layout),1L)
      expect_identical(unique(as.character(p$data$pipeline_id)),x$evaluations[[2]]$pipeline_id)
      expect_identical(nrow(attr(p,"stormbreaker_curve_status")),2L)
      audit <- fun(x,show_unavailable=TRUE)
      expect_identical(nrow(ggplot2::ggplot_build(audit)$layout$layout),2L)
    }
    expect_identical(x,before)
  }
})

test_that("all unavailable displays clearly report absence without pipeline axes", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("binary")$external
  x$status <- "failed"
  for(p in list(sb_plot_performance(x,"roc_auc"),sb_plot_calibration(x),sb_plot_dca(x))) {
    expect_match(ggplot2::ggplot_build(p)$data[[1]]$label,"No usable results",fixed=TRUE)
    expect_identical(nrow(ggplot2::ggplot_build(p)$layout$layout),1L)
  }
})
