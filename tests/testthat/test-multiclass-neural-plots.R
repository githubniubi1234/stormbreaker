test_that("multiclass network heatmaps preserve every native connection", {
  skip_if_not_installed("nnet")
  x <- scale(as.matrix(iris[1:4]))
  target <- stats::model.matrix(~ Species - 1, iris)
  m <- withr::with_seed(41, nnet::nnet(x, target, size = 3, softmax = TRUE,
    maxit = 500, decay = 0.1, trace = FALSE))
  fit <- list(method = "nnet", engine = m, features = names(iris)[1:4],
    classes = levels(iris$Species))
  result <- sb_multiclass_native_plots(fit)
  d <- result$tables$connection_weights
  expect_equal(d$weight, m$wts)
  expect_equal(nrow(d), length(m$wts))
  expect_setequal(as.character(d$destination), c(paste("Hidden", 1:3), fit$classes))
  expect_equal(sum(d$source == "Bias"), 6L)
  expect_setequal(d$layer, c("Input to hidden", "Hidden to output"))
  expect_s3_class(ggplot2::ggplot_build(result$plots$connection_weights), "ggplot_built")
})
