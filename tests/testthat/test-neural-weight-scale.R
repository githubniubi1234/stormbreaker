test_that("neural color mapping preserves extreme and small signed weights", {
  d <- data.frame(unit = 1:5, weight = c(-3000, -1, 0, 1, 2000))
  p <- ggplot2::ggplot(d, ggplot2::aes(unit, 1, fill = weight)) +
    ggplot2::geom_tile() + sb_neural_weight_scale()
  b <- ggplot2::ggplot_build(p)
  expect_identical(p$data, d)
  expect_equal(nrow(b$data[[1]]), nrow(d))
  expect_false(anyNA(b$data[[1]]$fill))
  expect_equal(length(unique(b$data[[1]]$fill)), nrow(d))
  scale <- b$plot$scales$get_scales("fill")
  expect_equal(scale$trans$inverse(scale$trans$transform(d$weight)), d$weight)
})
