test_that("bagging views preserve every original tree and native importance entry", {
  skip_if_not_installed("ipred")
  d <- sb_example_binary(100L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features), sb_pipeline(learner = "bagging_binary"))
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  prediction <- sb_predict_final(fit, d, "patient_id")
  e <- unserialize(fit$learner_state$engine_raw)
  sizes <- sb_bagging_diagnostics(fit, "tree_size")
  usage <- sb_bagging_diagnostics(fit)
  imp <- sb_bagging_diagnostics(fit, "importance")
  for (i in seq_along(e$mtrees)) {
    tree <- e$mtrees[[i]]$btree
    ids <- as.integer(rownames(tree$frame))
    expect_equal(sizes$table$terminal_nodes[i], sum(!((2L * ids) %in% ids)))
    expect_equal(sizes$table$nodes[i], length(ids))
    rows <- imp$table[imp$table$tree == i, ]
    expect_identical(rows$feature, features)
    native <- tree$variable.importance
    for (j in seq_along(features)) expect_equal(rows$importance[j],
      if (features[j] %in% names(native)) unname(native[features[j]]) else 0, tolerance = 0)
  }
  expect_equal(nrow(imp$table), 300L)
  expect_equal(sum(usage$table$splits), sum(sizes$table$nodes - sizes$table$terminal_nodes))
  vars <- unlist(lapply(e$mtrees, function(x) x$btree$frame$var))
  expect_equal(usage$table$splits, as.numeric(table(factor(vars, levels = features))))
  expect_null(fit$diagnostic_artifacts)
  bundle <- sb_native_diagnostics(fit, complete = TRUE)
  expect_identical(bundle$status$status, rep("available", 3L))
  for (v in bundle$diagnostics) {
    p <- sb_plot_algorithm_diagnostic(v)
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    expect_null(p$labels$caption)
    expect_null(p$labels$subtitle)
    if (v$type == "bagging_tree_size") expect_equal(sum(p$data$trees), 100L)
    if (v$type == "bagging_importance") expect_equal(p$data$value,
      as.numeric(tapply(imp$table$importance, imp$table$feature, mean)), tolerance = 0)
  }
  expect_identical(sb_predict_final(fit, d, "patient_id"), prediction)
  expect_identical(serialize(fit, NULL, version = 2L), before)
})
