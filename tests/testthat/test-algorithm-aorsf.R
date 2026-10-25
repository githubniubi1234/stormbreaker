test_that("aorsf views preserve native leaf counts and final-only OOB scores", {
  skip_if_not_installed("aorsf")
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in c("aorsf_binary", "aorsf_survival")) {
    binary <- id == "aorsf_binary"
    d <- if (binary) sb_example_binary(100L, seed = 108L) else sb_example_survival(100L, seed = 8612L)
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    prediction <- sb_predict_final(fit, d, "patient_id")
    e <- unserialize(fit$learner_state$engine_raw)
    sizes <- sb_aorsf_diagnostics(fit)
    oob <- sb_aorsf_diagnostics(fit, "oob_score")
    expect_equal(mean(sizes$table$terminal_nodes), e$get_mean_leaves_per_tree(), tolerance = 1e-12)
    expect_equal(sizes$table$terminal_nodes,
      vapply(e$forest$leaf_summary, function(z) length(which(z != 0)), integer(1)))
    expect_equal(nrow(sizes$table), 500L)
    expect_equal(oob$table$score, as.numeric(e$eval_oobag$stat_values), tolerance = 0)
    expect_identical(oob$measure, if (binary) "AUC-ROC" else "Harrell's C-index")
    expect_equal(oob$table$trees, 500)
    expect_equal(e$oobag_eval_every, 500)
    expect_equal(nrow(oob$table), 1L)
    expect_null(fit$diagnostic_artifacts)
    expect_null(e$importance)
    bundle <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", 2L))
    for (v in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(v)
      b <- ggplot2::ggplot_build(p)
      expect_s3_class(b, "ggplot_built")
      expect_null(p$labels$caption)
      expect_null(p$labels$subtitle)
      if (v$type == "aorsf_tree_size") expect_equal(sum(p$data$trees), 500L) else {
        expect_length(b$data, 1L)
        expect_equal(b$data[[1]]$y, oob$table$score, tolerance = 0)
        expect_equal(nrow(b$data[[1]]), 1L)
      }
    }
    expect_identical(sb_predict_final(fit, d, "patient_id"), prediction)
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
