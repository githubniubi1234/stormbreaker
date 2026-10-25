test_that("batch defaults defer final fits but preserve the comparison", {
  d <- sb_example_binary(70, 19)
  args <- list(data = d, outcome = "outcome", predictors = c("gene_a", "gene_b"),
    id = "patient_id", models = c("Logistic", "NaiveBayes"), folds = 3, seed = 11)
  eager <- do.call(bi_run, c(args, list(plots = TRUE)))
  fit_original <- sb_fit_final
  calls <- 0L
  local_mocked_bindings(sb_fit_final = function(...) {
    calls <<- calls + 1L
    fit_original(...)
  })
  deferred <- do.call(bi_run, args)
  expect_identical(calls, 0L)
  expect_length(deferred$details, 0L)
  expect_gt(length(deferred$plots), 0L)
  expect_equal(deferred$comparison$comparison, eager$comparison$comparison)
  local_mocked_bindings(sb_compare_cv = function(...) stop("Comparison must not repeat"))
  chosen <- sb_complete_run(deferred, "Logistic")
  expect_identical(calls, 1L)
  expect_equal(chosen$genes$gene, args$predictors)
  expect_equal(sb_predict_final(chosen$fit, d, "patient_id")$predictions,
    sb_predict_final(eager$details$Logistic$fit, d, "patient_id")$predictions)
  expect_setequal(names(chosen$plots$algorithm), names(eager$details$Logistic$figures$plots))
})

test_that("calculation entrypoints do not request browsers by default", {
  expect_identical(formals(sb_compare_models)$view, FALSE)
  expect_identical(formals(sb_compare_cohorts)$view, FALSE)
})

test_that("survival batch defaults preserve results and defer selection to completion", {
  d <- sb_example_survival(80, 21)
  args <- list(data = d, time = "time", event = "event", predictors = c("gene_a", "gene_b"),
    id = "patient_id", models = c("Cox", "survreg_weibull"), folds = 3, seed = 11)
  eager <- do.call(surv_run, c(args, list(plots = TRUE)))
  original <- sb_fit_final
  calls <- 0L
  local_mocked_bindings(sb_fit_final = function(...) {
    calls <<- calls + 1L
    original(...)
  })
  run <- do.call(surv_run, args)
  expect_identical(calls, 0L)
  expect_length(run$details, 0L)
  expect_equal(run$comparison$comparison, eager$comparison$comparison)
  expect_gt(length(run$plots), 0L)
  local_mocked_bindings(sb_compare_cv = function(...) stop("Unexpected comparison"))
  done <- sb_complete_run(run, "Cox")
  expect_identical(calls, 1L)
  expect_equal(sb_predict_final(done$fit, d, "patient_id")$predictions,
    sb_predict_final(eager$details$Cox$fit, d, "patient_id")$predictions)
  expect_setequal(done$genes$gene, args$predictors)
  expect_setequal(names(done$plots$algorithm), names(eager$details$Cox$figures$plots))
})

test_that("single pipeline still draws details and explicit FALSE skips them", {
  d <- sb_example_binary(60, 21)
  args <- list(data = d, outcome = "outcome", predictors = c("gene_a", "gene_b"),
    id = "patient_id", models = "Logistic", folds = 3)
  one <- do.call(bi_run, args)
  expect_identical(one$details$Logistic$fit$status, "success")
  expect_gt(length(one$details$Logistic$figures$plots), 0L)
  none <- do.call(bi_run, c(args, list(plots = FALSE)))
  expect_length(none$details, 0L)
  expect_length(none$plots, 0L)
  expect_equal(one$comparison$comparison, none$comparison$comparison)
})

test_that("multiclass batch comparisons avoid extra method fits", {
  local_mocked_bindings(sb_analyze_methods = function(...) stop("Unexpected final fits"))
  out <- multi_run(iris, "Species", models = c("lda", "qda"), folds = 3)
  expect_null(out$details)
  expect_identical(out$comparison$comparison$status, c("success", "success"))
})


test_that("default comparisons never open pages but explicit requests still do", {
  d <- sb_example_binary(60, 13)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  opened <- character()
  local_mocked_bindings(sb_open_result_page = function(path) opened <<- c(opened, path))
  x <- sb_compare_models(task, models = "logistic", folds = 3)
  d$patient_id <- paste0("external_", d$patient_id)
  external <- list(validation = sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b")))
  y <- sb_compare_cohorts(task, external, models = "logistic", folds = 3)
  expect_length(opened, 0L)
  expect_null(x$report_file)
  expect_null(y$report_file)
  visible <- sb_compare_cohorts(task, external, models = "logistic", folds = 3, view = TRUE)
  expect_length(opened, 1L)
  expect_equal(visible$development$comparison, y$development$comparison)
})

test_that("deferred screening retains selected genes and selector figures", {
  d <- sb_example_binary(200, 21)
  withr::local_seed(81)
  d$outcome <- rbinom(nrow(d), 1, plogis(2 * d$gene_a - 2 * d$gene_b))
  args <- list(data = d, outcome = "outcome", id = "patient_id",
    predictors = c("gene_a", "gene_b", "gene_c"),
    models = c("Logistic", "NaiveBayes"), selector = bi_selector("Lasso"), folds = 3)
  eager <- do.call(bi_run, c(args, list(plots = TRUE)))
  run <- do.call(bi_run, args)
  name <- names(run$pipelines)[1L]
  done <- sb_complete_run(run, name)
  reference <- sb_complete_run(eager, name)
  expect_identical(done$fit$status, "success")
  expect_equal(done$genes, reference$genes)
  expect_equal(sb_predict_final(done$fit, d, "patient_id")$predictions,
    sb_predict_final(reference$fit, d, "patient_id")$predictions)
  expect_setequal(names(done$plots$algorithm), names(reference$plots$algorithm))
})
