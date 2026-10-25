test_that("member matrices realign binary OOF predictions by patient ID", {
  data <- make_sb_binary_fixture(n = 48L, seed = 91L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 17L)
  comparison <- sb_compare_cv(
    task,
    stack$pipelines,
    folds = 3L,
    seed = 17L
  )
  second_id <- stack$pipeline_ids[[2L]]
  comparison$results[[second_id]]$predictions <-
    comparison$results[[second_id]]$predictions[48:1, , drop = FALSE]

  matrix <- sb_stack_member_matrix(comparison, task, stack)

  expect_s3_class(matrix, "data.frame")
  expect_identical(names(matrix), stack$pipeline_ids)
  expect_identical(rownames(matrix), data$patient_id)
  for (pipeline_id in stack$pipeline_ids) {
    predictions <- comparison$results[[pipeline_id]]$predictions
    expected <- predictions$probability[
      match(data$patient_id, predictions$patient_id)
    ]
    expect_equal(matrix[[pipeline_id]], expected, tolerance = 0)
  }
})

test_that("member matrices reject incomplete OOF predictions", {
  data <- make_sb_binary_fixture(n = 48L, seed = 92L)
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
  comparison <- sb_compare_cv(task, stack$pipelines, folds = 3L, seed = 9L)
  pipeline_id <- stack$pipeline_ids[[1L]]

  duplicate <- comparison
  duplicate$results[[pipeline_id]]$predictions$patient_id[[1L]] <-
    duplicate$results[[pipeline_id]]$predictions$patient_id[[2L]]
  duplicate_error <- capture_sb_error(
    sb_stack_member_matrix(duplicate, task, stack)
  )
  expect_s3_class(
    duplicate_error,
    "stormbreaker_incomplete_stack_member_oof"
  )

  failed <- comparison
  failed$results[[pipeline_id]]$predictions$status[[1L]] <- "failed"
  failed_error <- capture_sb_error(
    sb_stack_member_matrix(failed, task, stack)
  )
  expect_s3_class(
    failed_error,
    "stormbreaker_incomplete_stack_member_oof"
  )

  nonfinite <- comparison
  nonfinite$results[[pipeline_id]]$predictions$probability[[1L]] <- Inf
  nonfinite_error <- capture_sb_error(
    sb_stack_member_matrix(nonfinite, task, stack)
  )
  expect_s3_class(
    nonfinite_error,
    "stormbreaker_incomplete_stack_member_oof"
  )
})

test_that("ridge logistic stack state matches direct coefficient prediction", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 60L, seed = 93L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 23L)
  comparison <- sb_compare_cv(
    task,
    stack$pipelines,
    folds = 3L,
    seed = stack$inner_seed
  )
  matrix <- sb_stack_member_matrix(comparison, task, stack)

  fitted <- sb_fit_stack_meta(task, stack, matrix)
  observed <- sb_predict_stack_meta(fitted$state, matrix)
  expected <- stats::plogis(
    fitted$state$intercept +
      as.vector(as.matrix(matrix) %*% fitted$state$coefficient)
  )

  expect_s3_class(fitted, "stormbreaker_stack_meta_fit")
  expect_identical(fitted$meta_learner, "ridge_logistic")
  expect_identical(fitted$tuning$alpha, 0)
  expect_identical(fitted$tuning$selection_rule, "lambda.min")
  expect_named(observed, "probability")
  expect_equal(observed$probability, expected, tolerance = 1e-14)
  expect_identical(fitted$state$features, stack$pipeline_ids)
})

test_that("ridge Cox stack state saves its OOF Breslow horizon", {
  skip_if_not_installed("glmnet")

  data <- make_sb_fixture(n = 60L, seed = 94L)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 2
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "coxph"),
    sb_pipeline(learner = "rpart_survival")
  ), inner_folds = 3L, inner_seed = 29L)
  comparison <- sb_compare_cv(
    task,
    stack$pipelines,
    folds = 3L,
    seed = stack$inner_seed
  )
  matrix <- sb_stack_member_matrix(comparison, task, stack)

  fitted <- sb_fit_stack_meta(task, stack, matrix)
  observed <- sb_predict_stack_meta(fitted$state, matrix)
  expected_risk <- as.vector(
    as.matrix(matrix) %*% fitted$state$coefficient
  )
  expected_probability <- sb_survival_probability_from_hazard(
    expected_risk,
    fitted$state$baseline_hazard$cumulative_baseline_hazard
  )

  expect_s3_class(fitted, "stormbreaker_stack_meta_fit")
  expect_identical(fitted$meta_learner, "ridge_cox")
  expect_identical(fitted$tuning$alpha, 0)
  expect_equal(observed$risk, expected_risk, tolerance = 1e-14)
  expect_equal(
    observed$survival_probability,
    expected_probability,
    tolerance = 1e-14
  )
  expect_identical(observed$prediction_horizon, 2)
  expect_identical(observed$time_unit, "month")
  expect_identical(
    fitted$state$baseline_hazard$source,
    "glmnet_fixed_beta_breslow"
  )
})
