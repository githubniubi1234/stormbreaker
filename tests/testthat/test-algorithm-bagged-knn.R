test_that("bagged KNN diagnostics preserve all original bags and native curves", {
  skip_if_not_installed("bnnSurvival")
  d <- sb_example_survival(100L, seed = 8612L)
  for (features in list(c("gene_c", "gene_a", "gene_b"), "gene_c")) {
    task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = "bagged_knn_survival"))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    e <- unserialize(fit$learner_state$engine_raw)
    bags <- methods::slot(e, "base_learners")
    z <- sb_bagged_knn_diagnostics(fit)
    expect_equal(nrow(z$table), 50L * nrow(d))
    for (i in seq_along(bags)) {
      index <- methods::slot(bags[[i]], "bootstrap_sample")
      expect_identical(z$table$multiplicity[z$table$bag == i],
        as.integer(table(factor(index, levels = seq_len(nrow(d))))))
    }
    query <- d[c(1L, 7L), features, drop = FALSE]
    query[[1L]] <- query[[1L]] + 0.17
    curves <- sb_bagged_knn_diagnostics(fit, "bag_curves", query)
    model <- sb_restore_bagged_knn_survival_model(fit$learner_state)
    original <- sb_predict_bagged_knn_survival_curves(model, query)
    for (i in seq_len(nrow(query))) {
      summary <- curves$summary[curves$summary$query == i, ]
      expect_equal(summary$time, original$time)
      expect_equal(summary$survival, unname(original$survival[i, ]), tolerance = 1e-12)
    }
    native_predict <- get("predict", asNamespace("bnnSurvival"))
    native <- native_predict(bags[[1L]], methods::slot(e, "train_data"),
      data.matrix(sb_bagged_knn_internal_data(query)$data), methods::slot(e, "timepoints"),
      methods::slot(e, "metric"), methods::slot(e, "weighting_function"), methods::slot(e, "k"))
    expect_equal(curves$table$survival[curves$table$bag == 1L], as.numeric(t(native)))
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(curves)), "ggplot_built")
    expect_identical(sb_native_diagnostics(fit, complete = TRUE)$status$status,
      c("available", "unavailable"))
    expect_identical(sb_native_diagnostics(fit, complete = TRUE, newdata = query)$status$status,
      c("available", "available"))
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
