test_that("reports preserve boundaries and omit private source data", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("binary")
  x$external$scored_predictions$patient_id[1] <- "PRIVATE_PATIENT_SENTINEL"
  x$external$error_message <- "PRIVATE_ERROR_SENTINEL"
  before <- x
  path <- tempfile(fileext = ".html")
  result <- sb_render_report(list(Development = x$cv, Validation = x$external), path,
                             title = "My <report> & evidence")
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_identical(result, normalizePath(path, winslash = "/"))
  for (text in c("Development CV", "External evaluation", "Calibration", "Decision curves",
                 "data:image/svg+xml", "My &lt;report&gt; &amp; evidence", "roc_auc")) {
    expect_match(html, text, fixed = TRUE)
  }
  for (text in c("PRIVATE_PATIENT_SENTINEL", "PRIVATE_ERROR_SENTINEL", "<script", "<iframe")) {
    expect_identical(grepl(text, html, fixed = TRUE), FALSE)
  }
  expect_identical(x, before)
})

test_that("reports never overwrite an existing file or accept ambiguous sections", {
  x <- make_sb_curve_view_fixture("binary")$cv
  path <- tempfile(fileext = ".html")
  writeLines("PRESERVE_THIS_FILE", path)
  expect_error(sb_render_report(x, path), class = "stormbreaker_report_exists")
  expect_identical(readLines(path), "PRESERVE_THIS_FILE")
  for (bad in list(list(x), list(), structure(list(x,x), names=c("Same","Same")), 1)) {
    target <- tempfile(fileext = ".html")
    expect_error(sb_render_report(bad, target), class = "stormbreaker_bad_report_input")
    expect_identical(file.exists(target), FALSE)
  }
})

test_that("survival report retains saved horizon and failed curves", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("survival", TRUE)$external
  x$evaluations[[1]]$status <- "failed"
  x$evaluations[[1]]$metrics <- list()
  x$evaluations[[1]]$error_message <- "SECRET_FAILURE_DETAIL"
  path <- tempfile(fileext = ".html")
  sb_render_report(x, path)
  html <- paste(readLines(path, warn=FALSE),collapse="\n")
  for (text in c("harrell_c", "prediction_horizon", "month", "not_saved", "failed")) {
    expect_match(html,text,fixed=TRUE)
  }
  expect_identical(grepl("SECRET_FAILURE_DETAIL",html,fixed=TRUE),FALSE)
})

test_that("report rendering cannot fit or evaluate models", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("binary")$external
  local_mocked_bindings(
    sb_fit_final=function(...) stop("must not fit"),
    sb_run_cv=function(...) stop("must not run CV"),
    sb_evaluate_external=function(...) stop("must not evaluate"),
    sb_external_binary_auc_metric=function(...) stop("must not compute")
  )
  path <- tempfile(fileext=".html")
  sb_render_report(x,path)
  expect_identical(file.exists(path),TRUE)
})

test_that("survival CV without requested horizon has explicit missing curves", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_survival(60, seed=1)
  task <- sb_surv_task(d,"patient_id","time","event","gene_a","month")
  x <- sb_run_cv(task,folds=3,seed=1)
  path <- tempfile(fileext=".html")
  sb_render_report(x,path)
  html <- paste(readLines(path,warn=FALSE),collapse="\n")
  expect_match(html,"not_requested",fixed=TRUE)
  expect_match(html,"Development CV",fixed=TRUE)
  expect_match(html,"harrell_c",fixed=TRUE)
})

test_that("report tables faithfully format saved data without raw values hidden", {
  tab <- data.frame(pipeline_id="<unsafe&label>", estimate=c(0.123456789),
                    status="success", raw_patient="HIDDEN_PATIENT")
  html <- sb_report_table(tab,c("pipeline_id","estimate","status"),"Saved <metric>")
  expect_match(html,"0.12345679",fixed=TRUE)
  expect_match(html,"&lt;unsafe&amp;label&gt;",fixed=TRUE)
  expect_match(html,"Saved &lt;metric&gt;",fixed=TRUE)
  expect_identical(grepl("HIDDEN_PATIENT",html,fixed=TRUE),FALSE)
})

test_that("report writes separate failures and excludes failed pipeline names from HTML", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("binary",TRUE)$external
  x$evaluations[[1]]$status <- "failed"
  x$evaluations[[1]]$error_message <- 'PRIVATE reason, with "quotes"'
  before <- x
  path <- tempfile(fileext=".html")
  sb_render_report(x,path)
  csv <- sub("\\.html$","_failures.csv",path)
  html <- paste(readLines(path,warn=FALSE),collapse="\n")
  issues <- read.csv(csv,stringsAsFactors=FALSE)
  expect_identical(grepl(x$evaluations[[1]]$pipeline_id,html,fixed=TRUE),FALSE)
  expect_match(html,"Download failures CSV",fixed=TRUE)
  expect_identical(grepl("PRIVATE reason",html,fixed=TRUE),FALSE)
  expect_identical(issues$error_message[issues$issue_type=="pipeline"],x$evaluations[[1]]$error_message)
  expect_identical(nrow(issues),nrow(sb_failure_table(list(Result=x))))
  expect_identical(x,before)
})

test_that("existing companion CSV prevents either report output from being overwritten", {
  x <- make_sb_curve_view_fixture("binary")$external
  path <- tempfile(fileext=".html")
  csv <- sub("\\.html$","_failures.csv",path)
  writeLines("PRESERVE_CSV",csv)
  expect_error(sb_render_report(x,path),class="stormbreaker_output_exists")
  expect_identical(file.exists(path),FALSE)
  expect_identical(readLines(csv),"PRESERVE_CSV")
})

test_that("all-failed report explicitly states no usable results and exports failures", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("binary")$external
  x$status <- "failed"; x$metrics <- list()
  path <- tempfile(fileext=".html")
  sb_render_report(x,path)
  html <- paste(readLines(path,warn=FALSE),collapse="\n")
  expect_match(html,"No usable results in this section",fixed=TRUE)
  expect_match(html,"This is not a successful analysis",fixed=TRUE)
  expect_identical(grepl(x$pipeline_id,html,fixed=TRUE),FALSE)
  csv <- read.csv(sub("\\.html$","_failures.csv",path),stringsAsFactors=FALSE)
  expect_identical(csv$pipeline_id,x$pipeline_id)
  expect_identical(csv$status,"failed")
})
