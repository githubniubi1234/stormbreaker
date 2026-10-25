test_that("C5.0 native feature summaries and trial sizes retain original values", {
  skip_if_not_installed("C50")
  d <- sb_example_binary(160L,seed=108L); d$gene_a <- d$gene_a+1.5*d$outcome
  task <- sb_binary_task(d,"patient_id","outcome",c("gene_c","gene_a","gene_b"))
  for(id in c("c50_binary","c50_rules_binary")) {
    fit <- sb_fit_final(task,sb_pipeline(learner=id),capture_diagnostics=TRUE)
    expect_identical(fit$status,"success")
    before <- serialize(fit,NULL,version=2L)
    prediction <- sb_predict_final(fit,d,"patient_id")
    e <- fit$diagnostic_artifacts$learner$engine
    for(type in c("usage","splits","size")) {
      z <- sb_c50_diagnostics(fit,type)
      if(type=="size") expect_equal(z$table$size,e$size) else
        expect_equal(z$table$value,C50::C5imp(e,metric=type)[fit$learner_state$features,"Overall"])
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)),"ggplot_built")
    }
    expect_identical(sb_native_diagnostics(fit,complete=TRUE)$status$status,rep("available",3))
    compact <- fit; compact$diagnostic_artifacts <- NULL
    expect_identical(sb_native_diagnostics(compact,complete=TRUE)$status$status,rep("missing_capture",3))
    expect_identical(sb_predict_final(compact,d,"patient_id"),prediction)
    expect_identical(serialize(fit,NULL,version=2L),before)
  }
})

test_that("unsplit C5.0 trials do not invent feature importance", {
  skip_if_not_installed("C50")
  d <- sb_example_binary(160L,seed=108L)
  fit <- sb_fit_final(sb_binary_task(d,"patient_id","outcome",c("gene_c","gene_a","gene_b")),
    sb_pipeline(learner="c50_binary"),capture_diagnostics=TRUE)
  expect_identical(fit$status,"success")
  expect_identical(sb_native_diagnostics(fit,complete=TRUE)$status$status,c("unavailable","unavailable","available"))
  expect_equal(sb_c50_diagnostics(fit,"size")$table$size,rep(1,10))
})
