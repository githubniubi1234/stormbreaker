test_that("binary stacking learns only from inner OOF member predictions", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 60L, seed = 101L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c"),
    decision_thresholds = c(0.2, 0.5)
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 31L)

  result <- sb_run_stack_cv(task, stack, folds = 3L, seed = 37L)

  expect_s3_class(result, "stormbreaker_stack_cv_result")
  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(result$predictions$patient_id, data$patient_id)
  expect_true(all(is.finite(result$predictions$probability)))
  expect_identical(
    result$provenance$fold_assignment,
    result$members$fold_assignment
  )
  expect_length(result$inner_folds, 3L)
  for (audit in result$inner_folds) {
    analysis_ids <- result$folds[[audit$fold]]$analysis_ids
    expect_identical(rownames(audit$member_matrix), analysis_ids)
    expect_identical(names(audit$member_matrix), stack$pipeline_ids)
    for (pipeline_id in stack$pipeline_ids) {
      predictions <- audit$comparison$results[[pipeline_id]]$predictions
      expect_setequal(predictions$patient_id, analysis_ids)
      expect_identical(anyDuplicated(predictions$patient_id), 0L)
      expect_identical(all(predictions$status == "success"), TRUE)
    }
  }
  expect_identical(
    tail(result$comparison$pipeline_id, 1L),
    stack$id
  )
  expect_named(
    result$metrics,
    c(
      "roc_auc", "average_precision", "brier_score", "log_loss",
      "calibration", "moderate_calibration", "decision_curve"
    )
  )
})

test_that("outer assessment outcomes cannot alter their stacked predictions", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 48L, seed = 102L)
  assignment <- setNames(rep(1:3, length.out = nrow(data)), data$patient_id)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 41L)

  original <- sb_run_stack_cv(task, stack, folds = assignment, seed = 43L)
  changed <- data
  changed_rows <- assignment[data$patient_id] == 1L
  changed$outcome[changed_rows] <- 1L - changed$outcome[changed_rows]
  changed_task <- sb_binary_task(
    changed,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  perturbed <- sb_run_stack_cv(
    changed_task,
    stack,
    folds = assignment,
    seed = 43L
  )

  expect_equal(
    original$predictions$probability[changed_rows],
    perturbed$predictions$probability[changed_rows],
    tolerance = 0
  )
  expect_identical(
    original$inner_folds[[1L]]$member_matrix,
    perturbed$inner_folds[[1L]]$member_matrix
  )
  expect_identical(
    original$inner_folds[[1L]]$meta$state,
    perturbed$inner_folds[[1L]]$meta$state
  )
})

test_that("survival stacking returns risk and locked-horizon evaluation", {
  skip_if_not_installed("glmnet")

  data <- make_sb_fixture(n = 72L, seed = 103L)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 2,
    decision_thresholds = c(0.2, 0.4)
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "coxph"),
    sb_pipeline(learner = "survreg_weibull")
  ), inner_folds = 3L, inner_seed = 47L)

  result <- sb_run_stack_cv(task, stack, folds = 3L, seed = 53L)

  expect_s3_class(result, "stormbreaker_stack_cv_result")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(result$complete_survival_probability_oof, TRUE)
  expect_true(all(is.finite(result$predictions$risk)))
  expect_true(all(
    result$predictions$survival_probability >= 0 &
      result$predictions$survival_probability <= 1
  ))
  expect_named(
    result$metrics,
    c(
      "harrell_c", "ipcw_brier", "time_auc", "uno_c",
      "calibration", "weak_calibration", "decision_curve"
    )
  )
  expect_identical(
    vapply(
      result$inner_folds,
      function(audit) audit$meta$state$baseline_hazard$prediction_horizon,
      numeric(1)
    ),
    rep(2, 3L)
  )
})

test_that("stack folds fail visibly when the ridge dependency is unavailable", {
  data <- make_sb_binary_fixture(n = 48L, seed = 104L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L)
  local_mocked_bindings(sb_glmnet_available = function() FALSE)

  result <- sb_run_stack_cv(task, stack, folds = 3L, seed = 59L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$error_class), "stormbreaker_missing_glmnet")
  expect_true(all(is.na(result$predictions$probability)))
  expect_identical(names(result$members$results), stack$pipeline_ids)
  expect_identical(nrow(result$comparison), length(stack$pipelines) + 1L)
})
