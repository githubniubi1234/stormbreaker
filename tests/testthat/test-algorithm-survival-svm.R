test_that("all twelve survival-SVM variants retain native scores and basis weights", {
  skip_if_not_installed("survivalsvm")
  task <- make_survival_svm_diagnostic_task()
  for (id in sb_survivalsvm_learner_ids()) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success", info = id)
    if (fit$status != "success") next
    before <- serialize(fit, NULL, version = 2L)
    model <- fit$diagnostic_artifacts$learner
    a <- model$audit
    x <- fit$diagnostic_artifacts$predictors
    z <- as.data.frame(scale(as.matrix(x), center = a$center, scale = a$scale))
    score <- sb_survival_svm_diagnostics(fit)
    native <- as.numeric(stats::predict(model$engine, newdata = z)$predicted)
    expect_equal(score$table$value, native, tolerance = 1e-10)
    q <- sb_survival_svm_diagnostics(fit, "basis_weights")
    m <- model$engine$model.fit
    w <- if (model$params$type == "regression") as.numeric(m$Beta) else {
      v <- as.numeric(t(m$DifMat$Mat) %*% m$Alpha)
      if (model$params$type == "hybrid") v <- v + m$Beta - m$Delta * m$Betarstar
      as.numeric(v)
    }
    expect_equal(q$table$value, w)
    expect_true(q$sensitive_training_data)
    constant <- q
    constant$table$value[] <- 0
    constant_plot <- ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(constant))
    expect_equal(constant_plot$data[[2L]]$y, 1)
    expect_true(score$sensitive_training_data)
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(score)), "ggplot_built")
    if (id %in% sb_survivalsvm_linear_learner_ids()) {
      beta <- sb_survival_svm_diagnostics(fit, "coefficients")
      expect_equal(as.numeric(as.matrix(x) %*% beta$table$value + beta$intercept), native, tolerance = 1e-10)
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(beta)), "ggplot_built")
    } else expect_error(sb_survival_svm_diagnostics(fit, "coefficients"), class = "stormbreaker_unavailable_diagnostic")
    expect_true(identical(serialize(fit, NULL, version = 2L), before))
    fit$diagnostic_artifacts <- NULL
    expect_error(sb_survival_svm_diagnostics(fit), class = "stormbreaker_missing_diagnostic_artifacts")
  }
})
