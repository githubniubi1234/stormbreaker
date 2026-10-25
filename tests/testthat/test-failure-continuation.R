test_that("valid comparisons continue after a model-dependent failure below 100 inputs", {
  withr::local_seed(18)
  for(type in c("binary","survival")) {
    d <- if(type=="binary") sb_example_binary(45,seed=11) else sb_example_survival(45,seed=11)
    predictors <- paste0("input_",seq_len(35))
    for(name in predictors) d[[name]] <- rnorm(nrow(d))
    task <- if(type=="binary") sb_binary_task(d,"patient_id","outcome",predictors) else
      sb_surv_task(d,"patient_id","time","event",predictors,"month")
    pipelines <- if(type=="binary") list(sb_pipeline(learner="logistic"),sb_pipeline(learner="gaussian_naive_bayes")) else
      list(sb_pipeline(),sb_pipeline(selector="univ_cox_top_k",selector_n_features=1))
    x <- sb_compare_cv(task,pipelines,folds=3,seed=4)
    expect_identical(length(predictors),35L)
    expect_identical(unname(vapply(x$results,`[[`,character(1),"status")),c("failed","success"))
    expect_identical(names(x$results),vapply(pipelines,function(p) p$id,character(1)))
    issues <- sb_failure_table(x)
    failed_runs <- issues[issues$issue_type=="run",]
    expect_identical(failed_runs$fold,1:3)
    expect_identical(failed_runs$pipeline_id,rep(pipelines[[1]]$id,3L))
    expect_identical(all(!is.na(failed_runs$error_message)),TRUE)
  }
})
