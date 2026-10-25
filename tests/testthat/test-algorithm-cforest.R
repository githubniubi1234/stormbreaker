test_that("cforest diagnostics map native split columns and retain every tree", {
  skip_if_not_installed("partykit")
  features <- c("gene_c", "gene_a", "gene_b")
  for (id in c("cforest_binary", "cforest_survival")) {
    binary <- id == "cforest_binary"
    d <- if (binary) sb_example_binary(100L, seed = 108L) else sb_example_survival(100L, seed = 8612L)
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    prediction <- sb_predict_final(fit, d, "patient_id")
    e <- unserialize(fit$learner_state$engine_raw)
    sizes <- sb_cforest_diagnostics(fit, "tree_size")
    usage <- sb_cforest_diagnostics(fit)
    expect_equal(sizes$table$terminal_nodes, vapply(e$nodes, function(z) length(partykit::nodeids(z, terminal = TRUE)), integer(1)))
    expect_equal(sizes$table$nodes, vapply(e$nodes, function(z) length(partykit::nodeids(z)), integer(1)))
    native_features <- unlist(lapply(e$nodes, function(z) {
      ids <- setdiff(partykit::nodeids(z), partykit::nodeids(z, terminal = TRUE))
      unlist(partykit::nodeapply(z, ids = ids, FUN = function(n) {
        partykit::character_split(partykit::split_node(n), data = e$data)$name
      }), use.names = FALSE)
    }), use.names = FALSE)
    expect_equal(usage$table$splits, as.numeric(table(factor(native_features, levels = features))))
    expect_equal(usage$table$feature, features)
    expect_equal(sum(usage$table$splits), sum(sizes$table$nodes - sizes$table$terminal_nodes))
    expect_equal(nrow(sizes$table), 500L)
    expect_null(fit$diagnostic_artifacts)
    bundle <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", 2L))
    for (v in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(v)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      expect_null(p$labels$caption)
      expect_null(p$labels$subtitle)
      if (v$type == "cforest_tree_size") expect_equal(sum(p$data$trees), 500L)
    }
    expect_identical(sb_predict_final(fit, d, "patient_id"), prediction)
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
