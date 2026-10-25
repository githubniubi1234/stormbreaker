test_that("timereg views retain native cumulative processes and coefficient meanings", {
  skip_if_not_installed("timereg")
  d <- sb_example_survival(160L, seed = 8612L)
  names(d)[names(d) == "gene_c"] <- "gene c"
  names(d)[names(d) == "gene_a"] <- "gene a"
  features <- c("gene c", "gene a", "gene_b")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
  for (id in sb_timereg_survival_learner_ids()) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    aalen <- id == "aalen_additive_survival"
    native <- withr::with_seed(1L, if (aalen) timereg::aalen(
      survival::Surv(time, event) ~ `gene c` + `gene a` + gene_b, data = d,
      robust = 0, residuals = 0, n.sim = 0, weighted.test = 0, covariance = 0, resample.iid = 0, silent = 1) else
      timereg::prop.odds(timereg::Event(time, event) ~ `gene c` + `gene a` + gene_b, data = d,
        Nit = 20, detail = 0, n.sim = 0, weighted.test = 0, profile = 1, sym = 0, baselinevar = 1))
    effects <- sb_timereg_diagnostics(fit)
    baseline <- sb_timereg_diagnostics(fit, "baseline")
    expect_equal(baseline$table$time, unname(native$cum[, 1]), tolerance = 0)
    expect_equal(baseline$table$value, unname(native$cum[, 2]), tolerance = 0)
    expect_identical(effects$features, features)
    if (aalen) {
      expect_equal(effects$table$value, as.numeric(native$cum[, -(1:2)]), tolerance = 0)
      expect_equal(effects$table$time, rep(unname(native$cum[, 1]), length(features)), tolerance = 0)
      expect_equal(nrow(effects$table), nrow(native$cum) * length(features))
      expect_gt(sum(effects$table$value < 0), 0)
    } else expect_equal(effects$table$value, as.numeric(native$gamma[, 1]), tolerance = 0)
    compact <- fit
    compact$diagnostic_artifacts <- NULL
    expect_identical(sb_timereg_diagnostics(compact), effects)
    bundle <- sb_native_diagnostics(compact, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", 2L))
    for (view in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(view)
      built <- ggplot2::ggplot_build(p)
      expect_s3_class(built, "ggplot_built")
      expect_equal(if (view$type == "timereg_odds_effects") built$data[[2]]$x else built$data[[2]]$y,
        view$table$value, tolerance = 0)
      expect_null(p$labels$subtitle)
      expect_null(p$labels$caption)
    }
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
