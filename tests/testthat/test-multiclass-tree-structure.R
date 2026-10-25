test_that("multiclass tree displays every original node and class", {
  skip_if_not_installed("rpart")
  d <- iris
  names(d)[1:4] <- paste0("feature", 1:4)
  d$Species <- factor(paste0("class", as.integer(d$Species)))
  engine <- rpart::rpart(Species ~ ., d, method = "class", control = rpart::rpart.control(xval = 0))
  fit <- list(engine = engine, feature_codes = paste0("feature", 1:4),
    features = names(iris)[1:4], codes = paste0("class", 1:3), classes = levels(iris$Species))
  tree <- sb_multiclass_tree_structure(fit)
  expect_equal(nrow(tree$table), nrow(engine$frame))
  expect_equal(tree$table$prediction, fit$classes[engine$frame$yval])
  expect_length(grep("feature[1-4]", tree$table$incoming_rule), 0)
  plot <- sb_plot_rpart_tree(tree)
  expect_equal(nrow(ggplot2::ggplot_build(plot)$data[[2]]), nrow(engine$frame))
})
