test_that("scalar tables copy all five supported metrics without changing input", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  x <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1)
  before <- x
  for (metric in c("roc_auc", "average_precision", "brier_score", "log_loss")) {
    tab <- sb_result_table(x, metric = metric)
    saved <- x$metrics[[metric]]
    expect_identical(nrow(tab), 1L)
    expect_identical(tab$pipeline_id, x$pipeline$id)
    expect_identical(tab$task_type, "binary")
    expect_identical(tab$metric, metric)
    expect_identical(tab$estimate, saved$estimate)
    expect_identical(tab$status, saved$status)
    expect_identical(tab$run_status, x$status)
    expect_identical(tab$scope, saved$scope)
    expect_identical(tab$aggregation, saved$aggregation)
  }
  expect_identical(x, before)

  s <- make_sb_fixture()
  st <- sb_surv_task(s, "patient_id", "time", "event", "gene_a", "month")
  sx <- sb_run_cv(st, folds = 3, seed = 1)
  tab <- sb_result_table(sx, metric = "harrell_c")
  expect_identical(tab$task_type, "survival")
  expect_identical(tab$estimate, sx$metrics$harrell_c$estimate)
  expect_identical(tab$scope, sx$metrics$harrell_c$scope)
  expect_identical(tab$aggregation, sx$metrics$harrell_c$aggregation)
})

test_that("prediction and run tables preserve saved columns and row order", {
  for (binary in c(TRUE, FALSE)) {
    if (binary) {
      d <- make_sb_binary_fixture()
      task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
      pipeline <- sb_pipeline(learner = "logistic")
    } else {
      d <- make_sb_fixture()
      task <- sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month")
      pipeline <- sb_pipeline()
    }
    x <- sb_run_cv(task, pipeline, folds = 3, seed = 1)
    before <- x
    for (component in c("predictions", "runs")) {
      tab <- sb_result_table(x, component)
      expect_identical(tab[names(x[[component]])], x[[component]])
      expect_identical(tab$task_type, rep(if (binary) "binary" else "survival",
                                         nrow(x[[component]])))
    }
    expect_identical(x, before)
    expect_identical("error_stage" %in% names(sb_result_table(x, "runs")), FALSE)
  }
})

test_that("comparison tables retain every supplied pipeline and fold", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  x <- sb_compare_cv(task, list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), folds = 3, seed = 1)
  before <- x
  tab <- sb_result_table(x, metric = "roc_auc")
  expect_identical(tab$pipeline_id, names(x$results))
  expect_identical(tab$estimate, unname(vapply(
    x$results, function(z) z$metrics$roc_auc$estimate, numeric(1)
  )))
  for (component in c("predictions", "runs")) {
    tab <- sb_result_table(x, component)
    expected <- do.call(rbind, lapply(unname(x$results), `[[`, component))
    rownames(expected) <- NULL
    expect_identical(tab[names(expected)], expected)
  }
  expect_identical(x, before)
})

test_that("incomplete and nonfinite estimates are retained rather than repaired", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  a <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1)
  b <- a
  b$pipeline$id <- "synthetic_failure"
  b$status <- "failed"
  b$metrics$log_loss$status <- "incomplete"
  b$metrics$log_loss$estimate <- NA_real_
  c <- a
  c$pipeline$id <- "synthetic_infinite_loss"
  c$metrics$log_loss$estimate <- Inf
  x <- structure(list(results = list(b, a, c)),
                 class = "stormbreaker_cv_comparison")
  tab <- sb_result_table(x, metric = "log_loss")
  expect_identical(tab$pipeline_id,
                   c("synthetic_failure", a$pipeline$id, "synthetic_infinite_loss"))
  expect_identical(tab$estimate, c(NA_real_, a$metrics$log_loss$estimate, Inf))
  expect_identical(tab$run_status, c("failed", "success", "success"))
  expect_identical(tab$status, c("incomplete", a$metrics$log_loss$status,
                               a$metrics$log_loss$status))
})

test_that("extraction never calls fitting or metric computation", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  x <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1)
  local_mocked_bindings(
    sb_run_cv = function(...) stop("must not run CV"),
    sb_fit_final = function(...) stop("must not fit"),
    sb_binary_auc_metric = function(...) stop("must not recompute")
  )
  expect_identical(sb_result_table(x, metric = "roc_auc")$estimate,
                   x$metrics$roc_auc$estimate)
  expect_identical(nrow(sb_result_table(x, "predictions")), nrow(d))
})

test_that("unsupported requests have explicit errors", {
  d <- make_sb_fixture()
  task <- sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month")
  x <- sb_run_cv(task, folds = 3, seed = 1)
  for (metric in list(NULL, NA_character_, character(), c("harrell_c", "roc_auc"),
                      "calibration", "unsupported_score", 1)) {
    expect_error(sb_result_table(x, metric = metric),
                 class = "stormbreaker_bad_result_metric")
  }
  expect_error(sb_result_table(x, metric = "roc_auc"),
               class = "stormbreaker_result_metric_not_saved")
  expect_error(sb_result_table(x, "predictions", metric = "harrell_c"),
               class = "stormbreaker_unused_result_metric")
  expect_error(sb_result_table(list(), metric = "harrell_c"),
               class = "stormbreaker_bad_result_view")
  expect_error(sb_result_table(structure(list(results = list()),
                                        class = "stormbreaker_cv_comparison"),
                               metric = "harrell_c"),
               class = "stormbreaker_bad_result_view")
})
