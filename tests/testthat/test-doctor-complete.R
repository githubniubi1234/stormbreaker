test_that("completion reuses the chosen fit and isolates external failures", {
  set.seed(31)
  d <- data.frame(id = paste0("d", 1:100), a = rnorm(100), b = rnorm(100))
  d$y <- rbinom(100, 1, plogis(d$a))
  run <- bi_run(d, "y", c("a", "b"), "id", models = "Logistic", folds = 3)
  e <- d
  e$id <- paste0("e", 1:100)
  saved <- run$details$Logistic$fit
  local_mocked_bindings(sb_fit_final = function(...) stop("Unexpected refit"),
    sb_compare_cv = function(...) stop("Unexpected comparison"))
  bundle <- sb_complete_run(run, "Logistic", list(valid = e, broken = e["id"]), explain = "valid")
  expect_identical(bundle$fit, saved)
  expect_identical(bundle$comparison, run$comparison$comparison)
  expect_identical(bundle$evaluations$valid$status, "success")
  expect_equal(nrow(bundle$predictions$valid), 100L)
  expect_setequal(bundle$genes$gene, c("a", "b"))
  expect_type(bundle$issues$external$broken, "character")
  expect_s3_class(bundle$plots$external$valid, "ggplot")
  expect_equal(nrow(bundle$explanation$contributions), 100L)
  path <- tempfile("workflow-export-")
  sb_export_run(bundle, path)
  expect_equal(nrow(read.csv(file.path(path, "cohort-001-predictions.csv"))), 100L)
  expect_equal(read.csv(file.path(path, "genes-by-stage.csv"))$gene, bundle$genes$gene)
  expect_gt(file.info(file.path(path, "figure-001.pdf"))$size, 0)
  expect_equal(any(grepl("broken", read.csv(file.path(path, "issues.csv"))$scope)), TRUE)
  broken_explanation <- sb_complete_run(run, "Logistic", list(valid = e, broken = e["id"]), explain = "broken")
  expect_identical(broken_explanation$evaluations$valid$status, "success")
  expect_type(broken_explanation$issues$algorithm, "character")
})

test_that("survival completion fits only the selected pipeline when plots were skipped", {
  set.seed(32)
  d <- data.frame(id = paste0("s", 1:100), a = rnorm(100), b = rnorm(100),
    time = rexp(100), event = rep(0:1, 50))
  run <- surv_run(d, "time", "event", c("a", "b"), "id", models = "Cox", folds = 3, plots = FALSE)
  e <- d
  e$id <- paste0("v", 1:100)
  local_mocked_bindings(sb_compare_cv = function(...) stop("Unexpected comparison"))
  bundle <- sb_complete_run(run, "Cox", list(validation = e), plots = TRUE)
  expect_identical(bundle$fit$status, "success")
  expect_identical(bundle$evaluations$validation$status, "success")
  expect_equal(nrow(bundle$predictions$validation), 100L)
  expect_s3_class(bundle$plots$external$validation, "patchwork")
})
