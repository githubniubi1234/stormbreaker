test_that("flexsurv diagnostics preserve original parameter scales and spline baseline", {
  skip_if_not_installed("flexsurv")
  d <- sb_example_survival(160L, seed = 8612L)
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
  for (id in c(sb_flexsurv_parametric_learner_ids(), sb_flexsurvspline_learner_ids())) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    engine <- fit$diagnostic_artifacts$learner$engine
    engine_features <- sprintf(".stormbreaker_x%04d", seq_along(features))
    numerical <- fit$learner_state$numerical_fitting
    original_beta <- unname(stats::coef(engine)[engine_features] / numerical$scale)
    location_shift <- sum(original_beta * numerical$center)
    view <- sb_flexsurv_diagnostics(fit)
    expect_identical(view$table$feature, features)
    expect_equal(view$table$coefficient, original_beta, tolerance = 0)
    spline <- id %in% sb_flexsurvspline_learner_ids()
    if (spline) {
      original_gamma <- stats::coef(engine)[names(view$gamma)]
      original_gamma[[1L]] <- original_gamma[[1L]] - location_shift
      expect_equal(view$gamma, original_gamma, tolerance = 0)
      expect_equal(view$knots, unname(engine$knots), tolerance = 0)
      baseline <- sb_flexsurv_diagnostics(fit, "spline_baseline")
      zero <- as.data.frame(as.list(stats::setNames(-numerical$center / numerical$scale, engine_features)))
      native <- get("summary.flexsurvreg", asNamespace("flexsurv"))(engine,
        newdata = zero, type = "survival", t = baseline$table$time, ci = FALSE)[[1]]
      expect_equal(baseline$table$survival, native$est, tolerance = 1e-12)
      expect_equal(range(log(baseline$table$time)), range(engine$knots), tolerance = 1e-12)
      expect_equal(sum(exp(engine$knots) %in% baseline$table$time), length(engine$knots))
      expect_identical(baseline$params, fit$learner_state$params)
      expect_identical(baseline$time_unit, "month")
      expect_identical(sb_plot_algorithm_diagnostic(baseline)$labels$x, "Time (month)")
    } else {
      location <- fit$learner_state$params$location_parameter
      expect_equal(view$intercept, unname(stats::coef(engine)[location]) - location_shift, tolerance = 0)
      expect_equal(view$ancillary, stats::setNames(engine$res[names(view$ancillary), "est"], names(view$ancillary)), tolerance = 0)
    }
    compact <- fit
    compact$diagnostic_artifacts <- NULL
    expect_identical(sb_flexsurv_diagnostics(compact), view)
    bundle <- sb_native_diagnostics(compact, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", if (spline) 2L else 1L))
    for (diagnostic in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(diagnostic)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_null(p$labels$subtitle)
      expect_null(p$labels$caption)
    }
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
