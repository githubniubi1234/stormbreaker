test_that("multinomial plots retain coefficients and native AIC history", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("brglm2")
  d <- iris
  d$Species <- withr::with_seed(91, sample(d$Species))
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  for (method in c("logistic", "stepglm", "brglm2")) {
    fit <- sb_fit_multiclass(task, method, 3, 1)
    plots <- sb_multiclass_native_plots(fit)
    expect_equal(plots$tables$coefficients$coefficient, as.vector(stats::coef(fit$engine)))
    expect_length(grep("feature[0-9]", plots$tables$coefficients$feature), 0L)
    expect_setequal(unique(plots$tables$coefficients$contrast),
      paste(fit$classes[-1], "vs", fit$classes[1]))
    expect_s3_class(ggplot2::ggplot_build(plots$plots$coefficients), "ggplot_built")
    if (method == "stepglm") expect_equal(plots$tables$selection_path$AIC, fit$engine$anova$AIC)
  }
})
