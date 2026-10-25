test_that("AFT diagnostic time ratios reproduce native conditional quantile ratios", {
  skip_if_not_installed("survival")
  d <- sb_example_survival(160L, seed = 8612L)
  features <- c("gene_c", "gene_a", "gene_b")
  d$gene_c <- d$gene_c * 2 + 4
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
  for (id in sb_survreg_learner_ids()) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    engine <- fit$diagnostic_artifacts$learner$engine
    view <- sb_aft_diagnostics(fit)
    expect_identical(view$table$feature, features)
    expect_equal(view$table$coefficient, unname(stats::coef(engine)[-1]), tolerance = 0)
    expect_equal(view$intercept, unname(stats::coef(engine)[1]), tolerance = 0)
    expect_equal(view$distribution_scale, unname(engine$scale), tolerance = 0)
    expect_identical(view$distribution, engine$dist)
    rows <- d[1:5, features, drop = FALSE]
    native_predict <- get("predict.survreg", asNamespace("survival"))
    base <- native_predict(engine, newdata = rows, type = "quantile", p = c(.2, .5, .8))
    for (j in seq_along(features)) {
      shifted <- rows
      shifted[[j]] <- shifted[[j]] + 1
      quantiles <- native_predict(engine, newdata = shifted, type = "quantile", p = c(.2, .5, .8))
      expect_equal(unname(quantiles / base), matrix(view$table$time_ratio[j], 5, 3), tolerance = 1e-12)
    }
    compact <- fit
    compact$diagnostic_artifacts <- NULL
    expect_identical(sb_aft_diagnostics(compact), view)
    bundle <- sb_native_diagnostics(compact, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", 2L))
    for (diagnostic in bundle$diagnostics) {
      expect_identical(diagnostic$table, view$table)
      p <- sb_plot_algorithm_diagnostic(diagnostic)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_null(p$labels$subtitle)
      expect_null(p$labels$caption)
    }
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
