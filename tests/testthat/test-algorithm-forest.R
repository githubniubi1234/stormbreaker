test_that("all seven ranger variants expose exact saved importance and tree counts", {
  skip_if_not_installed("ranger")
  for (id in sb_ranger_learner_ids()) {
    binary <- id %in% sb_ranger_binary_learner_ids()
    d <- if (binary) make_sb_binary_fixture(n = 100L, seed = 108L) else make_sb_fixture(n = 100L, seed = 8612L)
    features <- c("gene_a", "gene_b", "gene_c")
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    expect_null(fit$diagnostic_artifacts)
    original <- fit
    q <- sb_forest_diagnostics(fit)
    expect_identical(q$table, fit$learner_state$feature_importance)
    sizes <- sb_forest_diagnostics(fit, "tree_size")
    forest <- unserialize(fit$learner_state$model_raw)
    native_object <- structure(list(forest = forest), class = "ranger")
    native <- vapply(seq_len(forest$num.trees), function(i) {
      nrow(subset(ranger::treeInfo(native_object, tree = i), terminal))
    }, integer(1))
    expect_equal(sizes$table$terminal_nodes, native)
    expect_equal(nrow(sizes$table), fit$learner_state$num_trees)
    expect_identical(fit, original)
    err <- tryCatch(sb_forest_diagnostics(fit, "oob_error"), error = identity)
    expect_s3_class(err, "stormbreaker_unavailable_diagnostic")
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      for (z in list(q, sizes)) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
      expect_equal(sum(sb_plot_algorithm_diagnostic(sizes)$data$trees), forest$num.trees)
    }
  }
})

test_that("randomForest and randomForestSRC diagnostics preserve native OOB missingness", {
  skip_if_not_installed("randomForest")
  skip_if_not_installed("randomForestSRC")
  for (id in c("random_forest_binary", "rsf_survival")) {
    binary <- id == "random_forest_binary"
    d <- if (binary) make_sb_binary_fixture(n = 100L, seed = 108L) else make_sb_fixture(n = 100L, seed = 8612L)
    features <- c("gene_a", "gene_b", "gene_c")
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    original <- fit
    engine <- fit$diagnostic_artifacts$learner$engine
    sizes <- sb_forest_diagnostics(fit, "tree_size")
    native <- if (binary) (engine$forest$ndbigtree + 1L) / 2L else engine$leaf.count
    expect_equal(sizes$table$terminal_nodes, native)
    q <- sb_forest_diagnostics(fit, "oob_error")
    expect_equal(q$table$error, as.vector(engine$err.rate))
    expect_equal(q$table$available, is.finite(as.vector(engine$err.rate)))
    expect_identical(fit, original)
    if (binary) {
      importance <- sb_forest_diagnostics(fit)
      expect_equal(importance$table$importance, unname(engine$importance[, "MeanDecreaseGini"]))
    } else {
      err <- tryCatch(sb_forest_diagnostics(fit), error = identity)
      expect_s3_class(err, "stormbreaker_unavailable_diagnostic")
    }
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      for (z in list(q, sizes)) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
      expect_equal(nrow(sb_plot_algorithm_diagnostic(q)$data), sum(q$table$available))
    }
    fit$diagnostic_artifacts <- NULL
    expect_identical(sb_forest_diagnostics(fit, "tree_size"), sizes)
    err <- tryCatch(sb_forest_diagnostics(fit, "oob_error"), error = identity)
    expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
  }
})
