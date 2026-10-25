page_tasks <- function() {
  make <- function(prefix, seed, n) {
    d <- sb_example_binary(n = n, seed = seed)
    d$patient_id <- paste0(prefix, "_", d$patient_id)
    sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  }
  list(development = make("dev", 312, 48),
       validation = make("val", 313, 30),
       test = make("test", 314, 30))
}

test_that("cohort comparison reuses the development fits and saved metrics", {
  tasks <- page_tasks()
  p <- list(sb_pipeline(learner = "logistic"),
            sb_pipeline(learner = "gaussian_naive_bayes"))
  x <- sb_compare_cohorts(tasks$development, tasks[c("validation", "test")],
    roles = c(validation = "validation", test = "test"), pipelines = p,
    folds = 3, seed = 52, view = FALSE)
  direct <- sb_compare_cv(tasks$development, p, folds = 3, seed = 52)
  expect_equal(x$development$comparison, direct$comparison, tolerance = 0)
  expect_identical(x$cohort_roles, c(Development = "development_cv",
    validation = "validation", test = "test"))
  for (i in seq_along(p)) {
    fit <- sb_fit_final(tasks$development, p[[i]])
    expect_equal(x$fits[[i]], fit, tolerance = 0)
    pred <- sb_predict_final(fit, tasks$test$data, "patient_id")
    expect_equal(x$predictions$test[[i]], pred, tolerance = 0)
  }
  page <- sb_result_page_data(x)
  expect_identical(page$default_sort, "Development")
  expect_identical(vapply(page$datasets, `[[`, character(1), "role"),
    c("development_cv", "validation", "test"))
  expect_length(page$models, 2)
  auc <- Filter(function(cell) cell$metric == "roc_auc" && cell$dataset == "test", page$cells)
  expect_equal(vapply(auc, `[[`, numeric(1), "estimate"),
    x$cohorts$test$comparison$roc_auc, tolerance = 0)
})

test_that("held-out outcomes cannot alter the development models or predictions", {
  tasks <- page_tasks()
  args <- list(development = tasks$development, cohorts = list(test = tasks$test),
    roles = c(test = "test"), models = "logistic", folds = 3, view = FALSE)
  a <- do.call(sb_compare_cohorts, args)
  changed <- tasks$test$data
  changed$outcome <- 1L - changed$outcome
  args$cohorts$test <- sb_binary_task(changed, "patient_id", "outcome", c("gene_a", "gene_b"))
  b <- do.call(sb_compare_cohorts, args)
  expect_equal(a$fits, b$fits, tolerance = 0)
  expect_equal(a$predictions, b$predictions, tolerance = 0)
  expect_equal(a$development, b$development, tolerance = 0)
})

test_that("friendly labels describe algorithms without changing configuration identity", {
  forest <- sb_display_learner("ranger_binary")
  expect_match(forest, "ranger", fixed = TRUE)
  expect_false(grepl("fixed|_binary", forest))
  expect_identical(sb_display_selector("none"), "None")
  expect_identical(sb_display_learner("glmnet_logistic", list(alpha = 0)), "Ridge logistic")
  expect_identical(sb_display_learner("glmnet_logistic", list(alpha = 1)), "Lasso logistic")
  expect_false(identical(sb_display_learner("xgboost_binary"),
    sb_display_learner("xgboost_binary_dart")))
  catalog <- sb_model_catalog("binary")
  expect_identical(catalog$model_label[catalog$model_id == "ridge_logistic"], "Ridge logistic")
  expect_identical(catalog$model_label[catalog$model_id == "glmnet_logistic"], "Lasso logistic")
  expect_false(any(grepl("fixed|_binary|_survival", catalog$model_label)))
  expect_false(any(grepl("fixed|_binary|_survival", sb_model_catalog("survival")$model_label)))
})

test_that("result pages contain summary data with escaped labels and no patient rows", {
  tasks <- page_tasks()
  cv <- sb_compare_models(tasks$development, models = "logistic", folds = 3, view = FALSE)
  label <- "</script><img src=x onerror=alert(1)>"
  named <- stats::setNames(list(cv), label)
  path <- tempfile(fileext = ".html")
  expect_identical(sb_view_results(named, path = path, open = FALSE), normalizePath(path, winslash = "/"))
  html <- paste(readLines(path, encoding = "UTF-8"), collapse = "\n")
  expect_false(grepl(label, html, fixed = TRUE))
  expect_false(grepl(tasks$development$data$patient_id[1], html, fixed = TRUE))
  expect_match(html, "Logistic regression", fixed = TRUE)
  expect_error(sb_view_results(named, path = path, open = FALSE),
    class = "stormbreaker_result_page_exists")
  expect_identical(sb_view_results(named, path = path, open = FALSE, overwrite = TRUE),
    normalizePath(path, winslash = "/"))
})

test_that("report opening is optional and never discards a completed comparison", {
  tasks <- page_tasks()
  opened <- character()
  local_mocked_bindings(sb_open_result_page = function(path) opened <<- c(opened, path))
  x <- sb_compare_models(tasks$development, models = "logistic", folds = 3, view = TRUE)
  expect_length(opened, 1)
  expect_identical(x$report_file, opened[[1]])
  expect_true(file.exists(x$report_file))
  opened <- character()
  y <- sb_compare_models(tasks$development, models = "logistic", folds = 3, view = FALSE)
  expect_length(opened, 0)
  expect_null(y$report_file)
  expect_equal(x$comparison, y$comparison, tolerance = 0)
})

test_that("roles and cross-cohort identity checks are explicit", {
  tasks <- page_tasks()
  expect_error(sb_compare_cohorts(tasks$development, list(test = tasks$test),
    roles = c(test = "typo"), models = "logistic", view = FALSE),
    class = "stormbreaker_bad_cohort_roles")
  expect_error(sb_compare_cohorts(tasks$development, list(test = tasks$development),
    roles = c(test = "test"), models = "logistic", view = FALSE),
    class = "stormbreaker_internal_external_patient_overlap")
})

test_that("survival cells preserve horizon-specific estimates and unavailable values", {
  make <- function(prefix, seed) {
    d <- sb_example_survival(45, seed)
    d$patient_id <- paste0(prefix, d$patient_id)
    sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"),
      time_unit = "month", prediction_horizon = 6)
  }
  dev <- make("d", 88)
  ext <- make("e", 89)
  x <- sb_compare_cohorts(dev, list(Test = ext), roles = c(Test = "test"),
    models = c("coxph", "survreg_weibull"), folds = 3, view = FALSE)
  page <- sb_result_page_data(x)
  cells <- Filter(function(cell) cell$metric == "time_auc" && cell$dataset == "Test", page$cells)
  expect_equal(vapply(cells, function(cell) cell$estimate, numeric(1)),
    x$cohorts$Test$comparison$time_auc, tolerance = 0)
  expect_identical(vapply(cells, function(cell) cell$horizon, numeric(1)), c(6, 6))
  expect_identical(vapply(cells, function(cell) cell$time_unit, character(1)), c("month", "month"))
  x$cohorts$Test$evaluations[[1]]$metrics$time_auc$estimate <- Inf
  x$results$Test <- x$cohorts$Test
  missing <- Filter(function(cell) cell$metric == "time_auc" && cell$dataset == "Test",
    sb_result_page_data(x)$cells)
  expect_false(missing[[1]]$available)
  expect_match(sb_page_json(missing[[1]]), '"estimate":null', fixed = TRUE)
})

test_that("JSON encoding preserves strings and blocks script termination", {
  expect_identical(sb_page_json(list(text = '"\\\n</script>&')),
    '{"text":"\\"\\\\\\u000a\\u003c/script\\u003e\\u0026"}')
  expect_identical(sb_page_json(list(value = NA_real_, values = list(1, TRUE, NULL))),
    '{"value":null,"values":[1,true,null]}')
})

test_that("opening failure retains the saved page and completed results", {
  tasks <- page_tasks()
  local_mocked_bindings(sb_open_result_page = function(path) stop("browser unavailable"))
  expect_warning(
    x <- sb_compare_models(tasks$development, models = "logistic", folds = 3, view = TRUE),
    "page was saved"
  )
  expect_true(file.exists(x$report_file))
  expect_identical(x$comparison$status, "success")
})

test_that("a saved failed run remains a visible row", {
  tasks <- page_tasks()
  cv <- sb_compare_models(tasks$development, models = "logistic", folds = 3, view = FALSE)
  cv$results[[1]]$status <- "failed"
  page <- sb_result_page_data(cv)
  expect_length(page$models, 1)
  expect_false(any(vapply(page$cells, function(cell) cell$available, logical(1))))
  expect_identical(page$cells[[1]]$run_status, "failed")
})
