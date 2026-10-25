test_that("repeated tables keep repetitions without pooling or changing inputs", {
  d <- sb_example_binary(60, seed = 1)
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  x <- sb_compare_repeated_cv(task, list(sb_pipeline(learner = "logistic")),
                              folds = 3, seeds = c(1, 1))
  before <- serialize(x, NULL)
  tab <- sb_result_table(x, metric = "roc_auc")
  expect_equal(tab$repeat_id, 1:2)
  expect_equal(tab$estimate, x$comparison$roc_auc)
  expect_equal(nrow(sb_result_table(x, "predictions")), 120)
  expect_equal(nrow(sb_result_table(x, "runs")), 6)
  expect_identical(serialize(x, NULL), before)
})

test_that("native stack and repeated stack views retain member and stack values", {
  skip_if_not_installed("glmnet")
  for (type in c("binary", "survival")) {
    fixture <- make_sb_curve_view_fixture(type)
    pipelines <- if (type == "binary") list(sb_pipeline(learner = "logistic"),
      sb_pipeline(learner = "gaussian_naive_bayes")) else list(sb_pipeline(),
      sb_pipeline(selector = "univ_cox_top_k", selector_n_features = 1))
    stack <- sb_stack(pipelines, inner_folds = 3, inner_seed = 31)
    x <- sb_run_stack_cv(fixture$task, stack, folds = 3, seed = 37)
    before <- serialize(x, NULL)
    metrics <- if (type == "binary") c("roc_auc", "average_precision", "brier_score", "log_loss") else
      c("harrell_c", "time_auc", "uno_c", "ipcw_brier")
    for (metric in metrics) {
      tab <- sb_result_table(x, metric = metric)
      expect_equal(tab$model_role, c("member", "member", "stack"))
      expect_equal(tail(tab$estimate, 1), x$metrics[[metric]]$estimate)
      expect_equal(head(tab$estimate, 2), sb_result_table(x$members, metric = metric)$estimate)
      if (type == "survival" && metric != "harrell_c") {
        expect_equal(tab$time_unit, rep("month", 3))
        expect_equal(tab$prediction_horizon, rep(6, 3))
      }
    }
    predictions <- sb_result_table(x, "predictions")
    expect_equal(predictions[predictions$model_role == "stack", names(x$predictions)],
                 x$predictions, ignore_attr = TRUE)
    expect_equal(nrow(sb_result_table(x, "runs")), 9)
    repeated <- structure(list(repetitions = list(x, x)),
                          class = "stormbreaker_repeated_cv_comparison")
    tab <- sb_result_table(repeated, metric = metrics[1])
    expect_equal(tab$repeat_id, rep(1:2, each = 3))
    expect_equal(length(unique(tab$view_id)), 6)
    expect_identical(serialize(x, NULL), before)
  }
})

test_that("stack metadata uses saved prediction units only when unambiguous", {
  x <- make_sb_curve_view_fixture("survival")$cv
  stack <- structure(list(status = x$status, task_type = "survival",
    predictions = x$predictions, metrics = x$metrics, runs = x$runs,
    task = list(), stack = list(id = "saved-stack"),
    members = structure(list(results = list(x)), class = "stormbreaker_cv_comparison")),
    class = "stormbreaker_stack_cv_result")
  stack$predictions$time_unit[1] <- "day"
  tab <- sb_result_table(stack, metric = "time_auc")
  expect_true(is.na(tail(tab$time_unit, 1)))
  stack$predictions$time_unit <- NULL
  tab <- sb_result_table(stack, metric = "time_auc")
  expect_true(is.na(tail(tab$time_unit, 1)))
})

test_that("composite readers cannot fit or recalculate and preserve missing values", {
  x <- make_sb_curve_view_fixture("survival")$cv
  y <- x
  y$metrics$ipcw_brier$estimate <- 1.25
  y$metrics$ipcw_brier$prediction_horizon <- 12
  y$metrics$ipcw_brier$time_unit <- "day"
  repeated <- structure(list(repetitions = list(x, y)),
                        class = "stormbreaker_repeated_cv_comparison")
  before <- serialize(repeated, NULL)
  fail <- function(...) stop("Must not execute fitting or metrics")
  local_mocked_bindings(sb_run_cv = fail, sb_fit_final = fail,
                       sb_compare_external = fail, sb_evaluate_external = fail)
  tab <- sb_result_table(repeated, metric = "ipcw_brier")
  expect_equal(tab$prediction_horizon, c(6, 12))
  expect_equal(tab$estimate[2], 1.25)
  expect_equal(tab$time_unit[2], "day")
  expect_identical(serialize(repeated, NULL), before)
})

test_that("native IECV stacks preserve their actual saved evaluation roles", {
  skip_if_not_installed("glmnet")
  cohorts <- lapply(1:2, function(i) {
    d <- sb_example_binary(90, seed = i)
    d$patient_id <- paste0(i, "_", d$patient_id)
    sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  })
  names(cohorts) <- c("A", "B")
  stack <- sb_stack(list(sb_pipeline(learner = "logistic"),
                        sb_pipeline(learner = "gaussian_naive_bayes")), inner_folds = 3)
  x <- sb_run_internal_external_cv(cohorts, stack)
  before <- serialize(x, NULL)
  tab <- sb_result_table(x, metric = "roc_auc")
  expect_equal(tab$model_role, rep("stack", 2))
  expect_equal(tab$heldout_cohort_id, c("A", "B"))
  expect_equal(tab$estimate, x$results$roc_auc)
  expect_identical(serialize(x, NULL), before)
})

test_that("empty prediction leaves and raw evaluation failures keep their identity", {
  x <- make_sb_curve_view_fixture()$external
  failed <- x
  failed$status <- "failed"
  failed$scored_predictions <- x$scored_predictions[0, , drop = FALSE]
  failed$metrics <- list()
  failed$error_message <- "saved ID mismatch"
  iecv <- structure(list(folds = list(
    list(heldout_cohort_id = "A", comparison = x),
    list(heldout_cohort_id = "B", comparison = failed))),
    class = "stormbreaker_internal_external_cv")
  expect_equal(nrow(sb_result_table(iecv, "predictions")), nrow(x$scored_predictions))
  tab <- sb_result_table(iecv, metric = "roc_auc")
  expect_true(is.na(tab$estimate[2]))
  expect_equal(tab$status[2], "not_saved")
  failures <- sb_failure_table(iecv)
  expect_equal(failures$heldout_cohort_id, "B")
  expect_equal(failures$error_message, "saved ID mismatch")
})

test_that("IECV tables identify held-out cohorts and development validation", {
  cohorts <- lapply(1:2, function(i) {
    d <- sb_example_binary(60, seed = i)
    d$patient_id <- paste0(i, "_", d$patient_id)
    sb_binary_task(d, "patient_id", "outcome", "gene_a")
  })
  names(cohorts) <- c("site <A>", "site B")
  x <- sb_run_internal_external_cv(cohorts, list(sb_pipeline(learner = "logistic")))
  tab <- sb_result_table(x, metric = "roc_auc")
  expect_equal(tab$heldout_cohort_id, names(cohorts))
  expect_equal(tab$estimate, x$results$roc_auc)
  expect_true(all(tab$validation_scheme == "internal_external_cv"))
  expect_equal(nrow(sb_result_table(x, "predictions")), 120)
  expect_equal(nrow(sb_result_table(x, "runs")), 2)
})
