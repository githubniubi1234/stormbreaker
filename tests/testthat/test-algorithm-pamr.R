test_that("PAMR views retain complete original and native shrunken centroids", {
  skip_if_not_installed("pamr")
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in c("pamr_binary", "pamr_survival")) {
    d <- if (id == "pamr_binary") sb_example_binary(160L, seed = 108L) else sb_example_survival(160L, seed = 8612L)
    task <- if (id == "pamr_binary") sb_binary_task(d, "patient_id", "outcome", features) else
      sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    prediction <- sb_predict_final(fit, d, "patient_id")
    e <- if (id == "pamr_binary") unserialize(fit$learner_state$engine_raw) else fit$learner_state$engine
    v <- sb_pamr_diagnostics(fit)
    original <- v$table[v$table$stage == "Before", ]
    after <- v$table[v$table$stage == "After", ]
    expect_equal(original$original_value, as.numeric(e$centroids), tolerance = 0)
    native <- pamr::pamr.predict(e, t(as.matrix(d[features])), threshold = 1, type = "centroid")
    expect_equal(after$original_value, as.numeric(native), tolerance = 0)
    delta <- sweep((e$centroids - e$centroid.overall) / e$sd, 2, e$threshold.scale * e$se.scale, "/")
    soft <- sign(delta) * pmax(abs(delta) - 1, 0)
    shrink <- sweep(soft, 2, e$threshold.scale * e$se.scale, "*")
    expect_equal(after$value, as.numeric(shrink), tolerance = 1e-12)
    expect_equal(after$original_value, as.numeric(e$centroid.overall + shrink * e$sd), tolerance = 1e-12)
    expect_equal(original$feature, rep(features, 2L))
    expect_equal(after$feature, rep(features, 2L))
    expect_equal(nrow(v$table), 12L)
    expect_gt(sum(after$value == 0), 0)
    expect_null(fit$diagnostic_artifacts)
    expect_identical(sb_native_diagnostics(fit, complete = TRUE)$status$status, "available")
    p <- sb_plot_algorithm_diagnostic(v)
    expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), nrow(v$table))
    expect_equal(p$data$value, v$table$value, tolerance = 0)
    expect_null(p$labels$caption)
    expect_null(p$labels$subtitle)
    expect_identical(sb_predict_final(fit, d, "patient_id"), prediction)
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
