test_that("MCP and SCAD expose their actual paths for both task types", {
  skip_if_not_installed("ncvreg")
  for (kind in c("binary", "survival")) {
    d <- if (kind == "binary") make_sb_binary_fixture(n = 140L, seed = 108L) else
      make_sb_fixture(n = 160L, seed = 8612L)
    task <- if (kind == "binary") sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")) else
      sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
    for (penalty in c("mcp", "scad")) {
      learner <- paste0("ncvreg_", penalty, if (kind == "binary") "_binary_bic" else "_cox_bic")
      fit <- sb_fit_final(task, sb_pipeline(learner = learner), capture_diagnostics = TRUE)
      expect_identical(fit$status, "success")
      native <- fit$diagnostic_artifacts$learner$diagnostic_path
      q <- sb_algorithm_diagnostics(fit, "coefficient_path")
      expect_equal(q$table$coefficient, as.vector(as.matrix(native$beta)))
      expect_equal(unique(q$table$lambda), native$lambda)
      b <- sb_algorithm_diagnostics(fit, "bic_path")
      expect_equal(b$table$bic, as.numeric(native$bic))
      expect_equal(b$selected_lambda, b$table$lambda[which.min(b$table$bic)])
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
        expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(b)), "ggplot_built")
      }
    }
  }
})
