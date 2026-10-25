test_that("a stack preserves distinct predeclared binary members", {
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  )

  stack <- sb_stack(pipelines, inner_folds = 3L, inner_seed = 41L)

  expect_s3_class(stack, "stormbreaker_stack")
  expect_identical(stack$task_type, "binary")
  expect_identical(stack$pipelines, pipelines)
  expect_identical(stack$pipeline_ids, vapply(
    pipelines,
    function(pipeline) pipeline$id,
    character(1)
  ))
  expect_identical(stack$meta_learner, "ridge_logistic")
  expect_identical(stack$inner_folds, 3L)
  expect_identical(stack$inner_seed, 41L)
  expect_identical(
    stack$id,
    paste0("ridge_stack__", paste(stack$pipeline_ids, collapse = "__"))
  )
  expect_identical(sb_validate_stack_contract(stack, "binary"), stack)
  expect_match(capture_output(print(stack)), "Members: 2")
})

test_that("a survival stack declares ridge Cox meta semantics", {
  stack <- sb_stack(list(
    sb_pipeline(learner = "coxph"),
    sb_pipeline(learner = "rpart_survival")
  ), inner_folds = 3L)

  expect_identical(stack$task_type, "survival")
  expect_identical(stack$meta_learner, "ridge_cox")
  expect_identical(stack$meta_alpha, 0)
  expect_identical(stack$meta_standardize, TRUE)
  expect_identical(stack$meta_selection_rule, "lambda.min")
  expect_identical(stack$meta_measure, "deviance")
})

test_that("stack construction rejects ambiguous member collections", {
  one_member <- capture_sb_error(sb_stack(list(
    sb_pipeline(learner = "logistic")
  )))
  expect_s3_class(one_member, "stormbreaker_bad_stack_pipelines")

  binary <- sb_pipeline(learner = "logistic")
  duplicate <- capture_sb_error(sb_stack(list(binary, binary)))
  expect_s3_class(
    duplicate,
    "stormbreaker_duplicate_stack_pipeline_id"
  )

  mixed <- capture_sb_error(sb_stack(list(
    binary,
    sb_pipeline(learner = "coxph")
  )))
  expect_s3_class(mixed, "stormbreaker_stack_task_mismatch")

  not_pipelines <- capture_sb_error(sb_stack(list(list(), list())))
  expect_s3_class(not_pipelines, "stormbreaker_bad_stack_pipelines")

  bad_folds <- capture_sb_error(sb_stack(
    list(
      sb_pipeline(learner = "logistic"),
      sb_pipeline(learner = "gaussian_naive_bayes")
    ),
    inner_folds = 2L
  ))
  expect_s3_class(bad_folds, "stormbreaker_bad_integer")
})

test_that("task slicing preserves binary roles and exact patient order", {
  data <- make_sb_binary_fixture(n = 40L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b"),
    decision_thresholds = c(0.2, 0.5)
  )
  rows <- sort(c(
    head(which(data$outcome == 0L), 8L),
    head(which(data$outcome == 1L), 8L)
  ))

  sliced <- sb_subset_task(task, rows)

  expect_s3_class(sliced, "stormbreaker_binary_task")
  expect_identical(sliced$data, data[rows, , drop = FALSE])
  expect_identical(sliced$predictors, task$predictors)
  expect_identical(sliced$decision_thresholds, task$decision_thresholds)
})

test_that("task slicing preserves survival horizon metadata", {
  data <- make_sb_fixture(n = 40L)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 2,
    decision_thresholds = c(0.2, 0.4)
  )
  rows <- sort(c(
    head(which(data$event == 0L), 8L),
    head(which(data$event == 1L), 8L)
  ))

  sliced <- sb_subset_task(task, rows)

  expect_s3_class(sliced, "stormbreaker_surv_task")
  expect_identical(sliced$data, data[rows, , drop = FALSE])
  expect_identical(sliced$time_unit, "month")
  expect_identical(sliced$prediction_horizon, 2)
  expect_identical(sliced$decision_thresholds, c(0.2, 0.4))
})
