test_that("all LightGBM variants keep native gain and count actual rather than configured trees", {
  skip_if_not_installed("lightgbm")
  skip_if_not_installed("jsonlite")
  d <- make_sb_binary_fixture(n = 100L, seed = 108L)
  d$constant <- 0
  features <- c("gene_a", "gene_b", "gene_c", "constant")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  for (id in sb_lightgbm_binary_learner_ids()) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    original <- fit$learner_state
    importance <- sb_boosting_diagnostics(fit, "importance")
    expect_identical(importance$native_table, original$feature_importance)
    expect_equal(importance$table$value, original$feature_importance$gain)
    expect_equal(importance$table$value[importance$table$feature == "constant"], 0)
    sizes <- sb_boosting_diagnostics(fit, "tree_size")
    # Native model text independently exposes a num_leaves line for each tree.
    lines <- strsplit(original$model_string, "\n", fixed = TRUE)[[1]]
    native_leaves <- as.integer(sub("num_leaves=", "", grep("^num_leaves=", lines, value = TRUE), fixed = TRUE))
    expect_equal(sizes$table$terminal_nodes, native_leaves)
    expect_equal(sizes$table$native_tree_index, seq_along(native_leaves) - 1L)
    expect_equal(sizes$requested_iterations, original$nrounds)
    if (id == "lightgbm_binary_goss") expect_lt(nrow(sizes$table), sizes$requested_iterations)
    expect_identical(fit$learner_state, original)
    err <- tryCatch(sb_boosting_diagnostics(fit, "training_loss"), error = identity)
    expect_s3_class(err, "stormbreaker_unavailable_diagnostic")
    if (requireNamespace("ggplot2", quietly = TRUE)) for (q in list(importance, sizes)) expect_s3_class(
      ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
  }
})

test_that("native empty-pair failures are explicit and preserve their original cause", {
  skip_if_not_installed("catboost")
  d <- make_sb_binary_fixture(n = 100L, seed = 108L)
  fit <- sb_fit_final(
    sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b")),
    sb_pipeline(learner = "catboost_binary", catboost_depth = 1L)
  )
  expect_identical(fit$status, "success")
  original <- fit$learner_state
  native <- tryCatch(catboost::catboost.get_feature_importance(
    unserialize(original$model_raw), type = "Interaction", thread_count = 1L
  ), error = identity)
  observed <- tryCatch(sb_boosting_diagnostics(fit, "interaction"), error = identity)
  if (inherits(native, "error")) {
    expect_s3_class(observed, "stormbreaker_unavailable_diagnostic")
    expect_identical(conditionMessage(observed$parent), conditionMessage(native))
  } else {
    # Future CatBoost releases may fix the native empty-result wrapper.
    expect_s3_class(observed, "stormbreaker_algorithm_diagnostic")
    expect_equal(nrow(observed$table), 0L)
  }
  expect_identical(fit$learner_state, original)
})

test_that("all CatBoost task types map native interaction indices and limit only display", {
  skip_if_not_installed("catboost")
  for (id in sb_catboost_diagnostic_ids()) {
    binary <- id == "catboost_binary"
    d <- if (binary) make_sb_binary_fixture(n = 100L, seed = 108L) else make_sb_fixture(n = 100L, seed = 8612L)
    features <- c("gene_c", "gene_a", "gene_b")
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    original <- fit$learner_state
    importance <- sb_boosting_diagnostics(fit, "importance")
    expect_identical(importance$native_table, original$feature_importance)
    q <- sb_boosting_diagnostics(fit, "interaction", max_pairs = 1L)
    native <- catboost::catboost.get_feature_importance(unserialize(original$model_raw), type = "Interaction", thread_count = 1L)
    expect_identical(q$native_table, native)
    expect_equal(q$table$feature_1, features[native[, 1] + 1L])
    expect_equal(q$table$feature_2, features[native[, 2] + 1L])
    expect_equal(q$table$score, unname(native[, 3]))
    expect_identical(fit$learner_state, original)
    expect_equal(nrow(q$table), nrow(native))
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- sb_plot_algorithm_diagnostic(q)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_equal(nrow(p$data), 1L)
      expect_equal(p$data$score, max(native[, 3]))
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(importance)), "ggplot_built")
      empty <- q
      empty$table <- empty$table[0, ]
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(empty)), "ggplot_built")
    }
    err <- tryCatch(sb_boosting_diagnostics(fit, "interaction", max_pairs = 0L), error = identity)
    expect_s3_class(err, "error")
  }
})

test_that("one-feature CatBoost fits have no invented interaction pairs", {
  skip_if_not_installed("catboost")
  for (id in sb_catboost_diagnostic_ids()) {
  binary <- id == "catboost_binary"
  d <- if (binary) make_sb_binary_fixture(n = 100L, seed = 108L) else make_sb_fixture(n = 100L, seed = 8612L)
  task <- if (binary) sb_binary_task(d, "patient_id", "outcome", "gene_a") else sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month")
  fit <- sb_fit_final(task, sb_pipeline(learner = id))
  expect_identical(fit$status, "success")
  original <- fit$learner_state
  q <- sb_boosting_diagnostics(fit, "interaction")
  expect_equal(nrow(q$table), 0L)
  expect_equal(dim(q$native_table), c(0L, 3L))
  expect_identical(fit$learner_state, original)
  expect_match(q$interpretation, "extractor was not called", fixed = TRUE)
  if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(
    ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
  }
})
