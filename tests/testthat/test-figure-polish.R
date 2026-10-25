test_that("shared styling and OOB colours preserve recorded values", {
  skip_if_not_installed("ggplot2")
  t <- sb_viz_theme()
  expect_identical(t$text$family, "sans")
  expect_identical(t$axis.text$colour, "#333333")
  d <- data.frame(trees = 1:3, `0` = c(.4, .3, .2), `1` = c(.5, .4, .3),
    oob = c(.45, .35, .25), check.names = FALSE)
  p <- sb_selector_oob_plot(d)
  expect_equal(p$data$error, unlist(d[-1], use.names = FALSE))
  expect_equal(ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")$map(c("0", "1", "oob")),
    c("#2878A0", "#B86236", "#686868"))
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})

test_that("mixture correlation labels fit without changing correlations", {
  skip_if_not_installed("ggplot2")
  d <- expand.grid(feature = paste0("gene_", 1:12), feature2 = paste0("gene_", 1:12),
    class = c("0", "1"), component = 1)
  d$correlation <- as.numeric(d$feature == d$feature2)
  p <- sb_plot_mclust_da_diagnostic(list(type = "mclust_da_correlations", table = d))
  expect_identical(p$data, d)
  expect_equal(p$theme$axis.text.x$angle, 45)
  expect_equal(p$coordinates$ratio, 1)
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})
