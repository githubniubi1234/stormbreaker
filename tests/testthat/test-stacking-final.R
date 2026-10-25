test_that("binary stack final fit predicts without external outcomes", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 60L, seed = 111L)
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
  ), inner_folds = 3L, inner_seed = 61L)

  fit <- sb_fit_stack_final(task, stack)
  new_data <- data[1:8, c("patient_id", "gene_a", "gene_b", "gene_c")]
  prediction <- sb_predict_stack_final(fit, new_data, "patient_id")
  expected <- stats::plogis(
    fit$meta$state$intercept +
      as.vector(
        as.matrix(prediction$member_predictions) %*%
          fit$meta$state$coefficient
      )
  )

  expect_s3_class(fit, "stormbreaker_stack_final_fit")
  expect_identical(fit$status, "success")
  expect_identical(names(fit$member_fits), stack$pipeline_ids)
  expect_identical(fit$meta$state$features, stack$pipeline_ids)
  expect_identical(fit$required_features, task$predictors)
  expect_s3_class(prediction, "stormbreaker_final_prediction")
  expect_identical(prediction$status, "success")
  expect_identical(prediction$predictions$patient_id, new_data$patient_id)
  expect_equal(prediction$predictions$probability, expected, tolerance = 1e-14)
  expect_identical(prediction$decision_thresholds, c(0.2, 0.5))

  endpoint <- sb_binary_endpoint(
    data[1:8, c("patient_id", "outcome")],
    "patient_id",
    "outcome"
  )
  evaluation <- sb_evaluate_external(prediction, endpoint)
  expect_s3_class(evaluation, "stormbreaker_external_evaluation")
  expect_identical(evaluation$status, "success")
})

test_that("named member lists fit and predict like ordinary stack inputs", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 48L, seed = 117L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    first = sb_pipeline(learner = "logistic"),
    second = sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 73L)

  fit <- sb_fit_stack_final(task, stack)
  prediction <- sb_predict_stack_final(
    fit,
    data[1:6, c("patient_id", "gene_a", "gene_b")],
    "patient_id"
  )

  expect_null(names(stack$pipeline_ids))
  expect_identical(fit$status, "success")
  expect_identical(prediction$status, "success")
})

test_that("serialized stack fits reproduce binary predictions", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 54L, seed = 112L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 67L)
  fit <- sb_fit_stack_final(task, stack)
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(fit, path)
  restored <- readRDS(path)
  new_data <- data[1:6, c("patient_id", "gene_a", "gene_b")]

  original <- sb_predict_stack_final(fit, new_data, "patient_id")
  replayed <- sb_predict_stack_final(restored, new_data, "patient_id")

  expect_identical(replayed, original)
})

test_that("survival stack final prediction reuses its saved OOF baseline", {
  skip_if_not_installed("glmnet")

  data <- make_sb_fixture(n = 72L, seed = 113L)
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
  ), inner_folds = 3L, inner_seed = 71L)

  fit <- sb_fit_stack_final(task, stack)
  new_data <- data[, c(
    "patient_id", "gene_a", "gene_b", "gene_c"
  )]
  prediction <- sb_predict_stack_final(fit, new_data, "patient_id")
  expected_risk <- as.vector(
    as.matrix(prediction$member_predictions) %*%
      fit$meta$state$coefficient
  )
  expected_probability <- sb_survival_probability_from_hazard(
    expected_risk,
    fit$meta$state$baseline_hazard$cumulative_baseline_hazard
  )

  expect_identical(fit$status, "success")
  expect_identical(prediction$status, "success")
  expect_equal(prediction$predictions$risk, expected_risk, tolerance = 1e-14)
  expect_equal(
    prediction$predictions$survival_probability,
    expected_probability,
    tolerance = 1e-14
  )
  expect_identical(
    prediction$predictions$prediction_horizon,
    rep(2, nrow(new_data))
  )

  endpoint <- sb_surv_endpoint(
    data[, c("patient_id", "time", "event")],
    "patient_id",
    "time",
    "event",
    "month",
    prediction_horizon = 2
  )
  evaluation <- sb_evaluate_external(prediction, endpoint)
  expect_s3_class(evaluation, "stormbreaker_external_evaluation")
  expect_identical(evaluation$status, "success")
})

test_that("final stack fitting reports a missing ridge dependency", {
  data <- make_sb_binary_fixture(n = 48L, seed = 114L)
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

  fit <- sb_fit_stack_final(task, stack)

  expect_identical(fit$status, "failed")
  expect_identical(fit$error_stage, "meta_fit")
  expect_identical(fit$error_class, "stormbreaker_missing_glmnet")
  expect_length(fit$member_fits, 0L)
})

test_that("a failed member prediction stops the stack without fallback", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 48L, seed = 115L)
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
  fit <- sb_fit_stack_final(task, stack)
  calls <- 0L
  local_mocked_bindings(
    sb_predict_final = function(...) {
      calls <<- calls + 1L
      list(
        status = "failed",
        error_class = "stormbreaker_synthetic_member_failure",
        error_message = "synthetic member failure"
      )
    }
  )
  new_data <- data[1:4, c("patient_id", "gene_a", "gene_b")]

  prediction <- sb_predict_stack_final(fit, new_data, "patient_id")

  expect_identical(prediction$status, "failed")
  expect_identical(prediction$error_stage, "member_prediction")
  expect_identical(
    prediction$error_class,
    "stormbreaker_stack_member_prediction_failure"
  )
  expect_identical(calls, 1L)
})
