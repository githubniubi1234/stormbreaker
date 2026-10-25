test_that("shared selector fitting uses the declared endpoint", {
  survival_data <- make_sb_fixture(n = 60L, seed = 20260819L)
  survival_x <- survival_data[c("gene_a", "gene_b", "gene_c")]
  survival_endpoint <- list(
    task_type = "survival",
    time = survival_data$time,
    event = as.integer(survival_data$event)
  )
  survival_selector <- sb_univ_cox_top_k_selector_adapter(2L)
  survival_fit <- sb_fit_runtime_selector(
    survival_selector,
    survival_x,
    survival_endpoint
  )
  direct_survival <- survival_selector$fit(
    survival_x,
    survival_endpoint$time,
    survival_endpoint$event
  )

  expect_identical(
    survival_fit$selected_features,
    direct_survival$selected_features
  )
  expect_identical(survival_fit$scores, direct_survival$scores)

  binary_data <- sb_example_binary(n = 60L, seed = 20260819L)
  binary_x <- binary_data[c("gene_a", "gene_b", "gene_c")]
  binary_endpoint <- list(
    task_type = "binary",
    outcome = as.integer(binary_data$outcome)
  )
  binary_selector <- sb_univ_binary_auc_top_k_selector_adapter(2L)
  binary_fit <- sb_fit_runtime_selector(
    binary_selector,
    binary_x,
    binary_endpoint
  )
  direct_binary <- binary_selector$fit(
    binary_x,
    binary_endpoint$outcome
  )

  expect_identical(
    binary_fit$selected_features,
    direct_binary$selected_features
  )
  expect_identical(binary_fit$scores, direct_binary$scores)

  none_fit <- sb_fit_runtime_selector(
    sb_none_selector_adapter(),
    binary_x,
    binary_endpoint
  )
  expect_identical(none_fit$selected_features, names(binary_x))
})

test_that("shared learner preparation reruns glmnet inner folds on its input", {
  skip_if_not_installed("glmnet")

  survival_data <- make_sb_fixture(n = 60L, seed = 20260819L)
  survival_x <- survival_data[c("gene_a", "gene_b", "gene_c")]
  survival_ids <- as.character(survival_data$patient_id)
  survival_pipeline <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 101L
  )
  survival_runtime <- sb_resolve_pipeline_runtime(survival_pipeline)
  survival_prepared <- sb_prepare_runtime_learner(
    survival_runtime$learner,
    survival_pipeline,
    survival_x,
    list(
      task_type = "survival",
      time = survival_data$time,
      event = as.integer(survival_data$event)
    ),
    survival_ids
  )

  expect_identical(survival_prepared$learner$id, "glmnet_cox")
  expect_identical(
    names(survival_prepared$tuning$fold_id),
    survival_ids
  )
  expect_length(survival_prepared$tuning$fold_id, nrow(survival_data))

  binary_data <- sb_example_binary(n = 60L, seed = 20260819L)
  binary_x <- binary_data[c("gene_a", "gene_b", "gene_c")]
  binary_ids <- as.character(binary_data$patient_id)
  binary_pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 101L
  )
  binary_runtime <- sb_resolve_pipeline_runtime(binary_pipeline)
  binary_prepared <- sb_prepare_runtime_learner(
    binary_runtime$learner,
    binary_pipeline,
    binary_x,
    list(
      task_type = "binary",
      outcome = as.integer(binary_data$outcome)
    ),
    binary_ids
  )

  expect_identical(binary_prepared$learner$id, "glmnet_logistic")
  expect_identical(names(binary_prepared$tuning$fold_id), binary_ids)
  expect_length(binary_prepared$tuning$fold_id, nrow(binary_data))
})

test_that("shared learner fitting dispatches survival and binary endpoints", {
  survival_data <- make_sb_fixture(n = 40L, seed = 20260819L)
  survival_x <- survival_data[c("gene_a", "gene_b")]
  survival_model <- sb_fit_runtime_learner(
    sb_coxph_learner_adapter(),
    survival_x,
    list(
      task_type = "survival",
      time = survival_data$time,
      event = as.integer(survival_data$event)
    )
  )
  expect_s3_class(survival_model, "coxph")

  binary_data <- sb_example_binary(n = 60L, seed = 20260819L)
  binary_x <- binary_data[c("gene_a", "gene_b")]
  binary_model <- sb_fit_runtime_learner(
    sb_logistic_learner_adapter(),
    binary_x,
    list(
      task_type = "binary",
      outcome = as.integer(binary_data$outcome)
    )
  )
  expect_s3_class(binary_model, "glm")
})
