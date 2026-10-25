test_that("issue explanations keep original reasons and distinguish stages", {
  local_mocked_bindings(sb_failure_table = function(x) sb_failure_empty())
  reasons <- c("Install missing_package", "cannot retain 1 or less variables",
    "Original path was not captured", "Predictor columns must match",
    "This nonlinear kernel has no single exact coefficient vector", "Unrecognized engine exception")
  r <- list(comparison = structure(list(), class = "stormbreaker_cv_comparison"),
    pipelines = list(), unavailable = list(), plot_errors = list(), details = list())
  for(i in seq_along(reasons))r$details[[paste0("model",i)]]<-list(
    fit=list(status="failed",error_message=reasons[i]),figures=NULL)
  tab <- sb_run_issues(r)
  expect_identical(tab$category, c("dependency","features","evidence","data","unsupported","other"))
  expect_identical(tab$reason, reasons)
  expect_identical(tab$stage, rep("final_fit",6))
  expect_identical(all(nzchar(tab$explanation)), TRUE)
  r$details <- list(ok=list(fit=list(status="success"),figures=list(unavailable=list(path=reasons[3]))))
  tab <- sb_run_issues(r)
  expect_identical(tab$stage,"plot path")
  expect_identical(r$details$ok$fit$status,"success")
})

test_that("survival doctor runs do not request ordinary binary ROC", {
  d <- sb_example_survival(90,19)
  r <- surv_run(d,time="time",event="event",id="patient_id",models="Cox",folds=3)
  expect_identical(r$comparison$status,"success")
  expect_null(r$plot_errors$roc)
  expect_null(r$plots$roc)
  expect_identical(r$issues,sb_run_issues(r))
  expect_s3_class(r$plots$performance,"ggplot")
})
