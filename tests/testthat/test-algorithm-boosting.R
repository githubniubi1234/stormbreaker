test_that("all twenty mboost options retain original loss, scale, stopping and components", {
  skip_if_not_installed("mboost")
  skip_if_not_installed("partykit")
  ids <- sb_mboost_diagnostic_ids()
  expect_length(ids, 20L)
  for (id in ids) {
    binary <- grepl("binary", id, fixed = TRUE)
    d <- if (binary) make_sb_binary_fixture(n = 100L, seed = 108L) else make_sb_fixture(n = 100L, seed = 8612L)
    features <- c("gene_a", "gene_b", "gene_c")
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    expect_null(fit$diagnostic_artifacts)
    original <- serialize(fit, NULL)
    engine <- unserialize(fit$learner_state$engine_raw)
    q <- sb_boosting_diagnostics(fit)
    expect_equal(q$table$loss, unname(mboost::risk(engine)), tolerance = 0)
    expect_equal(q$table$iteration, seq.int(0L, mboost::mstop(engine)))
    expect_equal(q$stopping_iteration, mboost::mstop(engine))
    if (grepl("unoc", id, fixed = TRUE)) expect_equal(sum(q$table$loss < 0), nrow(q$table))
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
    }
    if (!startsWith(id, "blackboost_")) {
      components <- sb_boosting_diagnostics(fit, "component_usage")
      native_selected <- as.integer(mboost::selected(engine))
      expect_identical(components$selections$component_id, native_selected)
      expect_equal(components$table$updates, tabulate(native_selected, nrow(components$table)))
      expect_equal(sum(components$table$updates), mboost::mstop(engine))
      expect_setequal(setdiff(components$table$feature, "(Intercept)"), features)
      expect_equal(components$selections$feature,
        components$table$feature[components$selections$component_id])
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        p <- sb_plot_algorithm_diagnostic(components)
        expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
        expect_equal(nrow(p$data), nrow(components$table))
      }
    } else {
      err <- tryCatch(sb_boosting_diagnostics(fit, "component_usage"), error = identity)
      expect_s3_class(err, "stormbreaker_unavailable_diagnostic")
    }
    expect_identical(serialize(fit, NULL), original)
    expect_equal(length(grepRaw("B001", serialize(q, NULL), fixed = TRUE)), 0L)
  }
})
