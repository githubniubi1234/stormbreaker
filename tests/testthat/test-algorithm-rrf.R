test_that("RRF forest views retain native importance, all trees and OOB histories", {
  skip_if_not_installed("RRF")
  d <- sb_example_binary(160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features),
    sb_pipeline(learner = "rrf_binary"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  prediction <- sb_predict_final(fit, d, "patient_id")
  e <- fit$diagnostic_artifacts$learner$engine
  importance <- sb_forest_diagnostics(fit)
  sizes <- sb_forest_diagnostics(fit, "tree_size")
  oob <- sb_forest_diagnostics(fit, "oob_error")
  expect_equal(importance$table$feature, rownames(e$importance))
  expect_equal(importance$table$importance, unname(e$importance[, "MeanDecreaseGini"]), tolerance = 0)
  expect_equal(nrow(importance$table), length(features))
  expect_gt(sum(importance$table$importance == 0), 0)
  native_sizes <- vapply(seq_len(e$ntree), function(k) {
    sum(RRF::getTree(e, k = k)[, "status"] == -1)
  }, integer(1))
  expect_equal(sizes$table$terminal_nodes, native_sizes, tolerance = 0)
  expect_equal(oob$table$error, as.numeric(e$err.rate), tolerance = 0)
  expect_equal(oob$table$available, is.finite(as.numeric(e$err.rate)))
  expect_equal(oob$table$trees, rep(seq_len(e$ntree), 3))
  compact <- fit
  compact$diagnostic_artifacts <- NULL
  expect_identical(sb_forest_diagnostics(compact), importance)
  expect_identical(sb_forest_diagnostics(compact, "tree_size"), sizes)
  expect_identical(sb_native_diagnostics(compact, complete = TRUE)$status$status,
    c("available", "available", "missing_capture"))
  bundle <- sb_native_diagnostics(fit, complete = TRUE)
  expect_identical(bundle$status$status, rep("available", 3L))
  for (v in bundle$diagnostics) expect_s3_class(
    ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(v)), "ggplot_built")
  expect_identical(sb_predict_final(fit, d, "patient_id"), prediction)
  expect_identical(serialize(fit, NULL, version = 2L), before)
})
