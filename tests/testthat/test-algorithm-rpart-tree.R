test_that("single-tree diagnostics preserve all nodes, exact splits and native predictions", {
  skip_if_not_installed("rpart")
  for (kind in c("binary", "survival")) {
    d <- if (kind == "binary") make_sb_binary_fixture(n = 150L, seed = 108L) else
      make_sb_fixture(n = 150L, seed = 8612L)
    features <- c("gene_a", "gene_b", "gene_c")
    task <- if (kind == "binary") sb_binary_task(d, "patient_id", "outcome", features) else
      sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = if (kind == "binary") "decision_tree_binary" else "rpart_survival"),
      capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    original <- fit
    q <- sb_tree_structure(fit)
    engine <- fit$diagnostic_artifacts$learner$engine
    expect_equal(q$table$node, as.numeric(rownames(engine$frame)))
    expect_identical(q$native_splits, engine$splits)
    expect_equal(q$table$n, engine$frame$n)
    expect_equal(q$table$parent[-1], floor(q$table$node[-1] / 2))
    expected <- if (kind == "binary") engine$frame$yval2[, 5] else engine$frame$yval
    expect_equal(q$table$prediction, unname(expected))
    expect_equal(sum(q$table$n[q$table$terminal]), nrow(d))
    expect_identical(fit, original)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- sb_plot_algorithm_diagnostic(q)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_equal(nrow(p$data), nrow(engine$frame))
    }
  }
})

test_that("single-node trees and missing captures are explicit", {
  skip_if_not_installed("rpart")
  d <- make_sb_binary_fixture(n = 60L, seed = 108L)
  d$gene_a[] <- 0
  d$gene_b[] <- 0
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "decision_tree_binary"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  q <- sb_tree_structure(fit)
  expect_equal(nrow(q$table), 1)
  expect_identical(q$table$incoming_rule, "root")
  expect_equal(q$table$prediction, mean(d$outcome))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
  }
  fit$diagnostic_artifacts <- NULL
  err <- tryCatch(sb_tree_structure(fit), error = identity)
  expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
})
