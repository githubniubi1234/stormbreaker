test_that("failure tables distinguish run and metric failures without mutation", {
  x <- make_sb_curve_view_fixture("binary")$external
  x$status <- "failed"
  x$error_stage <- "patient_match"
  x$error_class <- "saved_match_error"
  x$error_message <- "Saved original reason"
  x$metrics <- list()
  before <- x
  tab <- sb_failure_table(x,cohort_label="External A")
  expect_identical(tab$pipeline_id,x$pipeline_id)
  expect_identical(tab$issue_type,"pipeline")
  expect_identical(tab$stage,"external_evaluation")
  expect_identical(tab$source_error_stage,"patient_match")
  expect_identical(tab$error_message,"Saved original reason")
  expect_identical(tab$cohort_label,"External A")
  expect_identical(x,before)

  y <- make_sb_curve_view_fixture("survival")$external
  y$status <- "success"
  y$metrics <- list(harrell_c=y$metrics$harrell_c)
  y$metrics$harrell_c$status <- "not_evaluable"
  y$metrics$harrell_c$error_message <- "No comparable pairs saved"
  tab <- sb_failure_table(y)
  expect_identical(tab$issue_type,"metric")
  expect_identical(tab$metric,"harrell_c")
  expect_identical(tab$run_status,"success")
  expect_identical(tab$status,"not_evaluable")
  expect_identical(tab$prediction_horizon,6)
  expect_identical(tab$time_unit,"month")
})

test_that("fold records retain saved details and successful metrics are not failures", {
  x <- make_sb_curve_view_fixture("binary")$cv
  x$metrics <- list(roc_auc=x$metrics$roc_auc)
  expect_identical(nrow(sb_failure_table(x)),0L)
  x$status <- "failed"
  x$runs$status[1] <- "failed"
  x$runs$error_class[1] <- "saved_fit_error"
  x$runs$error_message[1] <- "Original fold reason"
  tab <- sb_failure_table(x)
  folds <- tab[tab$issue_type=="run",]
  expect_identical(folds$fold,as.integer(x$runs$fold[1]))
  expect_identical(folds$error_message,"Original fold reason")
  expect_identical(folds$stage,"cv_fold")
  expect_identical(folds$source_error_stage,NA_character_)
  expect_identical("patient_id" %in% names(tab),FALSE)
})

test_that("named cohorts remain distinct and unrequested metrics are not errors", {
  x <- make_sb_curve_view_fixture("binary")$external
  x$metrics <- list(log_loss=x$metrics$log_loss,unused=list(status="not_requested"))
  x$metrics$log_loss$estimate <- Inf
  tab <- sb_failure_table(list(A=x,B=x))
  expect_identical(tab$cohort_label,c("A","B"))
  expect_identical(tab$issue_type,rep("metric_nonfinite",2L))
  expect_identical(tab$status,rep("success",2L))
  expect_identical(tab$estimate,rep(Inf,2L))
})

test_that("failure CSV preserves errors, empty headers, and existing output", {
  x <- make_sb_curve_view_fixture("binary")$external
  x$metrics <- list(roc_auc=x$metrics$roc_auc)
  file <- tempfile(fileext=".csv")
  path <- sb_export_failures(x,file)
  expect_identical(path,normalizePath(file,winslash="/"))
  expect_identical(nrow(read.csv(file)),0L)
  saved <- readLines(file)
  expect_error(sb_export_failures(x,file),class="stormbreaker_output_exists")
  expect_identical(readLines(file),saved)
  x$status <- "failed"; x$error_message <- 'Reason, with "quotes" and newline\nsecond line'
  file <- tempfile(fileext=".csv")
  sb_export_failures(x,file)
  tab <- read.csv(file,stringsAsFactors=FALSE)
  expect_identical(tab$error_message[tab$issue_type=="pipeline"],x$error_message)
  expect_error(sb_failure_table(list(x)),class="stormbreaker_bad_failure_view")
})
