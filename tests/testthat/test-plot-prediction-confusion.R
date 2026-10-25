test_that("direct confusion counts and cell orientation match supplied predictions", {
  y <- c(rep(0, 72), rep(1, 199))
  prob <- c(rep(.1, 60), rep(.9, 12), rep(.1, 8), rep(.9, 191))
  p <- sb_plot_prediction_confusion(y, prob, "apparent", cohort = "Train",
    class_labels = c("Control", "Treat"))
  expect_equal(p$data$count, c(60, 12, 8, 191))
  b <- ggplot2::ggplot_build(p)
  expect_equal(as.numeric(b$data[[1]]$x), c(1, 1, 2, 2))
  expect_equal(as.numeric(b$data[[1]]$y), c(2, 1, 2, 1))
  expect_equal(nrow(b$data[[2]]), 4L)
  expect_identical(anyDuplicated(b$data[[2]][c("x", "y")]), 0L)
  context <- attr(p, "stormbreaker_prediction_context")
  expect_identical(context$kind, "apparent")
  expect_equal(context$n, 271)
  expect_equal(context$counts$tp, 191)
  expect_null(p$labels$subtitle)
  expect_null(p$labels$caption)
  expect_identical(p$labels$title, "Train")
})

test_that("threshold ties and absent classes retain all cells", {
  p <- sb_plot_prediction_confusion(c(0, 1), c(.5, .5), "external")
  expect_equal(p$data$count, c(0, 1, 0, 1))
  q <- sb_plot_prediction_confusion(c(0, 0), c(.1, .9), "oof", normalize = "observed")
  expect_equal(q$data$value, c(.5, .5, NA, NA))
  expect_equal(nrow(ggplot2::ggplot_build(q)$data[[2]]), 4)
  expect_equal(sb_plot_prediction_confusion(c(0, 1), c(0, 1), "external", threshold = 1)$data$count,
    c(1, 0, 0, 1))
})

test_that("direct confusion rejects invalid inputs rather than dropping rows", {
  expect_error(sb_plot_prediction_confusion(c(0, 1), .5, "apparent"), "aligned")
  expect_error(sb_plot_prediction_confusion(c(0, NA), c(.1, .2), "apparent"), "aligned")
  expect_error(sb_plot_prediction_confusion(c(0, 1), c(.1, 2), "apparent"), "aligned")
  expect_error(sb_plot_prediction_confusion(c(0, 1), c(.1, .2), "apparent", threshold = NA_real_), "threshold")
  expect_error(sb_plot_prediction_confusion(c(0, 1), c(.1, .2), "apparent", class_labels = c("A", "A")), "distinct")
})
