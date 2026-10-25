test_that("empty plot collections remain empty", {
  expect_identical(sb_plot_pages(list()), list())
})

test_that("two-column feature facets paginate without losing support rows", {
  d <- expand.grid(feature = paste0("gene_", 1:12), value = 1:5)
  d$effect <- d$value / 2
  p <- ggplot2::ggplot(d, ggplot2::aes(value, effect)) + ggplot2::geom_line() +
    ggplot2::geom_rug(data = d) + ggplot2::facet_wrap(~feature, ncol = 2)
  pages <- sb_plot_pages(p)
  expect_length(pages, 2)
  expect_equal(sum(vapply(pages, function(z) nrow(z$data), integer(1))), nrow(d))
  expect_equal(sum(vapply(pages, function(z) nrow(z$layers[[2]]$data), integer(1))), nrow(d))
  expect_equal(vapply(pages, function(z) z$facet$params$ncol, numeric(1)), c(plot1 = 2, plot2 = 2), ignore_attr = TRUE)
})

test_that("query pagination retains every observation and ribbon row", {
  d <- data.frame(query = rep(1:15, each = 3), time = rep(1:3, 15), value = seq_len(45))
  ribbon <- transform(d, lower = value - 1, upper = value + 1)
  p <- ggplot2::ggplot(d, ggplot2::aes(time, value)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(data = ribbon, ggplot2::aes(ymin = lower, ymax = upper)) +
    ggplot2::facet_wrap(~query)
  pages <- sb_plot_pages(p)
  expect_length(pages, 3)
  expect_equal(do.call(rbind, lapply(pages, function(z) z$data))$value, d$value)
  expect_equal(sum(vapply(pages, function(z) nrow(z$layers[[2]]$data), integer(1))), nrow(ribbon))
  expect_equal(nrow(p$layers[[2]]$data), 45L)
  expect_equal(vapply(pages, function(z) length(unique(z$data$query)), integer(1)),
    c("plot / 1" = 6L, "plot / 2" = 6L, "plot / 3" = 3L))
  expect_length(lapply(pages, ggplot2::ggplot_build), 3)
})

test_that("labelled posterior plots paginate and retain separate notes", {
  d <- data.frame(query = 1:60, mean = seq_len(60) / 61)
  p <- ggplot2::ggplot(d, ggplot2::aes(factor(query), mean)) + ggplot2::geom_point() +
    ggplot2::labs(caption = "Pointwise 95% posterior credible intervals")
  p <- sb_compact_plot(p)
  pages <- sb_plot_pages(p)
  expect_length(pages, 3)
  expect_null(p$labels$caption)
  expect_identical(attr(pages[[1]], "stormbreaker_notes"), "Pointwise 95% posterior credible intervals")
  expect_equal(sum(vapply(pages, function(z) nrow(z$data), integer(1))), 60L)
})
