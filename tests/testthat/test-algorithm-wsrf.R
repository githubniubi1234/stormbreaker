test_that("WSRF diagnostics never confuse compact placeholders with native metrics", {
  skip_if_not_installed("wsrf")
  d <- sb_example_binary(160L, seed = 108L)
  for (features in list(c("gene_c", "gene_a", "gene_b"), "gene_a")) {
    task <- sb_binary_task(d, "patient_id", "outcome", features)
    fit <- sb_fit_final(task, sb_pipeline(learner = "wsrf_binary"), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    prediction <- sb_predict_final(fit, d, "patient_id")
    x <- d[features]
    names(x) <- paste0(".stormbreaker_wsrf_x", seq_along(features))
    if (length(features) == 1L) x$.stormbreaker_wsrf_zero <- 0
    native <- withr::with_seed(1L, wsrf::wsrf(x,
      factor(d$outcome, levels = c(0, 1), labels = c("control", "case")),
      mtry = as.integer(floor(log2(length(features)) + 1)), ntree = 500L,
      weights = TRUE, parallel = FALSE, importance = FALSE, nodesize = 2L))
    importance <- sb_wsrf_diagnostics(fit)
    nodes <- sb_wsrf_diagnostics(fit, "tree_nodes")
    oob <- sb_wsrf_diagnostics(fit, "tree_oob")
    expect_equal(importance$table$importance, unname(native$importance[, "MeanDecreaseIGR"]), tolerance = 0)
    expect_equal(importance$table$feature[!importance$table$artificial_bridge], features)
    expect_equal(sum(importance$table$artificial_bridge), as.integer(length(features) == 1L))
    expect_equal(nodes$table$nodes, lengths(native$trees))
    expect_equal(oob$table$error, native$treeOOBErrorRates, tolerance = 0)
    expect_equal(oob$forest_oob, native$RFOOBErrorRate, tolerance = 0)
    expect_gt(sum(oob$table$error > 0), 0)
    expect_null(fit$learner_state$diagnostic_metrics)
    compact <- fit
    compact$diagnostic_artifacts <- NULL
    expect_equal(unserialize(compact$learner_state$engine_raw)$treeOOBErrorRates, rep(0, 500))
    expect_identical(sb_wsrf_diagnostics(compact, "tree_nodes"), nodes)
    expect_identical(sb_native_diagnostics(compact, complete = TRUE)$status$status,
      c("missing_capture", "available", "missing_capture"))
    bundle <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", 3))
    for (v in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(v)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_null(p$labels$caption)
      expect_null(p$labels$subtitle)
      if (v$type == "wsrf_tree_nodes") expect_equal(sum(p$data$trees), 500L)
      if (v$type == "wsrf_tree_oob") expect_equal(p$data$error, native$treeOOBErrorRates, tolerance = 0)
    }
    expect_identical(sb_predict_final(fit, d, "patient_id"), prediction)
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
