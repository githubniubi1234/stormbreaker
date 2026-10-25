test_that("the minimum pipeline has one explicit selector and learner", {
  pipeline <- sb_pipeline()

  expect_s3_class(pipeline, "stormbreaker_pipeline")
  expect_identical(pipeline$id, "none__coxph")
  expect_identical(pipeline$task_type, "survival")
  expect_identical(pipeline$selector, "none")
  expect_identical(pipeline$learner, "coxph")
  expect_identical(pipeline$risk_type, "linear_predictor")
  expect_identical(pipeline$risk_reference, "zero")
  expect_identical(pipeline$risk_formula, "x_beta")
  expect_identical(pipeline$risk_direction, "higher_is_higher_hazard")
})

test_that("the minimum binary pipeline predicts class-one probability", {
  pipeline <- sb_pipeline(learner = "logistic")

  expect_s3_class(pipeline, "stormbreaker_pipeline")
  expect_identical(pipeline$id, "none__logistic")
  expect_identical(pipeline$task_type, "binary")
  expect_identical(pipeline$selector, "none")
  expect_identical(pipeline$learner, "logistic")
  expect_identical(pipeline$prediction_type, "probability")
  expect_identical(pipeline$positive_class, 1L)
  expect_identical(
    pipeline$prediction_direction,
    "higher_is_more_likely_positive"
  )
  expect_identical(pipeline$fallback, "forbidden")
})

test_that("nested glmnet logistic declares fixed alpha and inner resampling", {
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 5L,
    inner_seed = 1701L
  )

  expect_identical(pipeline$task_type, "binary")
  expect_identical(pipeline$selector, "none")
  expect_identical(pipeline$learner, "glmnet_logistic")
  expect_identical(pipeline$prediction_type, "probability")
  expect_identical(pipeline$positive_class, 1L)
  expect_identical(
    pipeline$prediction_direction,
    "higher_is_more_likely_positive"
  )
  expect_identical(
    pipeline$learner_params,
    list(
      alpha = 1,
      standardize = TRUE,
      inner_folds = 5L,
      inner_seed = 1701L,
      selection_rule = "lambda.min",
      measure = "deviance"
    )
  )
  expect_identical(
    pipeline$id,
    "none__glmnet_logistic__alpha1__inner5__seed1701"
  )
})

test_that("glmnet logistic IDs distinguish Ridge and Elastic Net", {
  ridge <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 0,
    inner_folds = 3L,
    inner_seed = 101L
  )
  elastic <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 0.4,
    inner_folds = 3L,
    inner_seed = 202L
  )

  expect_identical(
    ridge$id,
    "none__glmnet_logistic__alpha0__inner3__seed101"
  )
  expect_identical(
    elastic$id,
    "none__glmnet_logistic__alpha0.4__inner3__seed202"
  )
  expect_identical(ridge$learner_params$alpha, 0)
  expect_identical(elastic$learner_params$alpha, 0.4)
})

test_that("glmnet standardization is explicit without changing old defaults", {
  default <- sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 101L
  )
  explicit_default <- sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_standardize = TRUE
  )
  disabled <- sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_standardize = FALSE
  )

  expect_identical(explicit_default, default)
  expect_identical(disabled$learner_params$standardize, FALSE)
  expect_identical(
    disabled$id,
    paste0(
      "none__glmnet_logistic__alpha1__standardize0",
      "__inner3__seed101"
    )
  )
  expect_identical(
    sb_validate_pipeline_contract(disabled, "binary"),
    disabled
  )

  for (value in list(NA, 0, "FALSE", c(TRUE, FALSE))) {
    error <- capture_sb_error(sb_pipeline(
      learner = "glmnet_logistic",
      inner_folds = 3L,
      inner_seed = 101L,
      glmnet_standardize = value
    ))
    expect_s3_class(error, "stormbreaker_bad_glmnet_standardize")
  }
  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "logistic",
    glmnet_standardize = FALSE
  ))
  expect_s3_class(irrelevant, "stormbreaker_irrelevant_glmnet_parameter")
})

test_that("glmnet lambda rule is predeclared without changing old defaults", {
  default <- sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L
  )
  explicit_default <- sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_lambda_rule = "lambda.min"
  )
  one_se <- sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_lambda_rule = "lambda.1se"
  )
  one_se_without_standardization <- sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_standardize = FALSE,
    glmnet_lambda_rule = "lambda.1se"
  )
  named_default <- sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_lambda_rule = c(rule = "lambda.min")
  )
  named_one_se <- sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_lambda_rule = c(rule = "lambda.1se")
  )

  expect_identical(explicit_default, default)
  expect_identical(named_default, default)
  expect_identical(named_one_se, one_se)
  expect_false(identical(named_default$id, one_se$id))
  expect_identical(one_se$learner_params$selection_rule, "lambda.1se")
  expect_identical(
    one_se$id,
    "none__glmnet_cox__alpha1__lambda1se__inner3__seed101"
  )
  expect_identical(
    one_se_without_standardization$id,
    paste0(
      "none__glmnet_cox__alpha1__standardize0__lambda1se",
      "__inner3__seed101"
    )
  )
  expect_identical(sb_validate_pipeline_contract(one_se, "survival"), one_se)

  for (value in list(NA_character_, "min", 1, c("lambda.min", "lambda.1se"))) {
    error <- capture_sb_error(sb_pipeline(
      learner = "glmnet_cox",
      inner_folds = 3L,
      inner_seed = 101L,
      glmnet_lambda_rule = value
    ))
    expect_s3_class(error, "stormbreaker_bad_glmnet_lambda_rule")
  }
  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "coxph",
    glmnet_lambda_rule = "lambda.1se"
  ))
  expect_s3_class(irrelevant, "stormbreaker_irrelevant_glmnet_parameter")
})

test_that("glmnet tuning measure is predeclared and task specific", {
  default <- sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 101L
  )
  explicit_default <- sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_tuning_measure = "deviance"
  )
  auc <- sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_tuning_measure = c(metric = "auc")
  )
  concordance <- sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_tuning_measure = "C"
  )

  expect_identical(explicit_default, default)
  expect_identical(auc$learner_params$measure, "auc")
  expect_identical(concordance$learner_params$measure, "C")
  expect_identical(
    auc$id,
    "none__glmnet_logistic__alpha1__measureauc__inner3__seed101"
  )
  expect_identical(
    concordance$id,
    "none__glmnet_cox__alpha1__measureC__inner3__seed101"
  )
  expect_identical(sb_validate_pipeline_contract(auc, "binary"), auc)
  expect_identical(
    sb_validate_pipeline_contract(concordance, "survival"),
    concordance
  )

  binary_error <- capture_sb_error(sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_tuning_measure = "C"
  ))
  survival_error <- capture_sb_error(sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_tuning_measure = "auc"
  ))
  malformed_error <- capture_sb_error(sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L,
    glmnet_tuning_measure = c("deviance", "C")
  ))
  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "coxph",
    glmnet_tuning_measure = "C"
  ))
  expect_s3_class(binary_error, "stormbreaker_bad_glmnet_tuning_measure")
  expect_s3_class(survival_error, "stormbreaker_bad_glmnet_tuning_measure")
  expect_s3_class(malformed_error, "stormbreaker_bad_glmnet_tuning_measure")
  expect_s3_class(irrelevant, "stormbreaker_irrelevant_glmnet_parameter")
})

test_that("nested glmnet logistic requires safe declared components", {
  missing_error <- capture_sb_error(
    sb_pipeline(learner = "glmnet_logistic")
  )
  expect_s3_class(missing_error, "stormbreaker_missing_inner_parameter")

  selector_error <- capture_sb_error(sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "glmnet_logistic",
    selector_n_features = 2L,
    inner_folds = 3L,
    inner_seed = 101L
  ))
  expect_s3_class(
    selector_error,
    "stormbreaker_incompatible_binary_selector"
  )
})

test_that("the binary minimum rejects unresolved selection and tuning", {
  selector_error <- capture_sb_error(sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "logistic",
    selector_n_features = 1L
  ))
  expect_s3_class(
    selector_error,
    "stormbreaker_incompatible_binary_selector"
  )

  inner_error <- capture_sb_error(sb_pipeline(
    learner = "logistic",
    inner_folds = 3L,
    inner_seed = 1L
  ))
  expect_s3_class(inner_error, "stormbreaker_irrelevant_inner_parameter")

  penalty_error <- capture_sb_error(sb_pipeline(
    learner = "logistic",
    glmnet_alpha = 1
  ))
  expect_s3_class(
    penalty_error,
    "stormbreaker_irrelevant_glmnet_parameter"
  )
})

test_that("unsupported components stop instead of being substituted", {
  selector_error <- capture_sb_error(sb_pipeline(selector = "lasso"))
  expect_s3_class(selector_error, "stormbreaker_unsupported_selector")

  learner_error <- capture_sb_error(sb_pipeline(learner = "random_forest"))
  expect_s3_class(learner_error, "stormbreaker_unsupported_learner")
})

test_that("a fixed Top-K univariable Cox pipeline is explicit", {
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "coxph",
    selector_n_features = 2L
  )

  expect_s3_class(pipeline, "stormbreaker_pipeline")
  expect_identical(pipeline$id, "univ_cox_top_k__coxph__k2")
  expect_identical(pipeline$selector, "univ_cox_top_k")
  expect_identical(pipeline$selector_params, list(n_features = 2L))
  expect_identical(pipeline$learner, "coxph")
  expect_identical(pipeline$fallback, "forbidden")
})

test_that("Top-K selector parameters are required and predeclared", {
  missing_error <- capture_sb_error(
    sb_pipeline(selector = "univ_cox_top_k")
  )
  expect_s3_class(missing_error, "stormbreaker_missing_selector_parameter")

  invalid_error <- capture_sb_error(
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 1.5
    )
  )
  expect_s3_class(invalid_error, "stormbreaker_bad_integer")

  irrelevant_error <- capture_sb_error(
    sb_pipeline(selector = "none", selector_n_features = 2L)
  )
  expect_s3_class(
    irrelevant_error,
    "stormbreaker_irrelevant_selector_parameter"
  )
})

test_that("a nested glmnet Cox pipeline declares its inner resampling", {
  pipeline <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    inner_folds = 5L,
    inner_seed = 1701L
  )

  expect_identical(pipeline$learner, "glmnet_cox")
  expect_identical(
    pipeline$learner_params,
    list(
      alpha = 1,
      standardize = TRUE,
      inner_folds = 5L,
      inner_seed = 1701L,
      selection_rule = "lambda.min",
      measure = "deviance"
    )
  )
  expect_identical(pipeline$fallback, "forbidden")
  expect_identical(
    pipeline$id,
    "none__glmnet_cox__alpha1__inner5__seed1701"
  )
})

test_that("glmnet pipeline IDs distinguish alpha and inner seeds", {
  first <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 101L
  )
  second <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0.9,
    inner_folds = 3L,
    inner_seed = 202L
  )

  expect_identical(
    first$id,
    "none__glmnet_cox__alpha0.5__inner3__seed101"
  )
  expect_identical(
    second$id,
    "none__glmnet_cox__alpha0.9__inner3__seed202"
  )
  expect_identical(first$learner_params$alpha, 0.5)
  expect_identical(second$learner_params$alpha, 0.9)
  expect_identical(identical(first$id, second$id), FALSE)
})

test_that("glmnet alpha is finite and bounded", {
  invalid_values <- list(-0.1, 1.1, NA_real_, Inf, "0.5", c(0.3, 0.7))

  for (value in invalid_values) {
    error <- capture_sb_error(sb_pipeline(
      learner = "glmnet_cox",
      glmnet_alpha = value,
      inner_folds = 3L,
      inner_seed = 101L
    ))
    expect_s3_class(error, "stormbreaker_bad_glmnet_alpha")
  }
})

test_that("alpha zero is an explicit Ridge Cox pipeline", {
  ridge <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0,
    inner_folds = 5L,
    inner_seed = 20260817L
  )

  expect_identical(ridge$learner_params$alpha, 0)
  expect_identical(
    ridge$id,
    "none__glmnet_cox__alpha0__inner5__seed20260817"
  )
  printed <- paste(capture.output(print(ridge)), collapse = "\n")
  expect_match(printed, "Penalty family: Ridge", fixed = TRUE)
})

test_that("glmnet alpha is rejected for ordinary Cox", {
  error <- capture_sb_error(sb_pipeline(glmnet_alpha = 0.5))

  expect_s3_class(error, "stormbreaker_irrelevant_glmnet_parameter")
})

test_that("nested glmnet Cox admits fold-local Top-K selection", {
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "glmnet_cox",
    selector_n_features = 2L,
    inner_folds = 5L,
    inner_seed = 1701L
  )

  expect_s3_class(pipeline, "stormbreaker_pipeline")
  expect_identical(
    pipeline$id,
    "univ_cox_top_k__glmnet_cox__k2__alpha1__inner5__seed1701"
  )

  too_small <- capture_sb_error(sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "glmnet_cox",
    selector_n_features = 1L,
    inner_folds = 5L,
    inner_seed = 1701L
  ))
  expect_s3_class(
    too_small,
    "stormbreaker_glmnet_too_few_selected_features"
  )
})

test_that("inner resampling arguments are required only for glmnet Cox", {
  missing_error <- capture_sb_error(sb_pipeline(learner = "glmnet_cox"))
  expect_s3_class(missing_error, "stormbreaker_missing_inner_parameter")

  irrelevant_error <- capture_sb_error(sb_pipeline(inner_folds = 5L))
  expect_s3_class(
    irrelevant_error,
    "stormbreaker_irrelevant_inner_parameter"
  )

  invalid_error <- capture_sb_error(sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 2L,
    inner_seed = 1701L
  ))
  expect_s3_class(invalid_error, "stormbreaker_bad_integer")
})

test_that("binary AUC Top-K is an explicit fold-local logistic pipeline", {
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 3L
  )

  expect_s3_class(pipeline, "stormbreaker_pipeline")
  expect_identical(
    pipeline$id,
    "univ_binary_auc_top_k__logistic__k3"
  )
  expect_identical(pipeline$task_type, "binary")
  expect_identical(pipeline$selector, "univ_binary_auc_top_k")
  expect_identical(pipeline$selector_params, list(n_features = 3L))
  expect_identical(pipeline$learner, "logistic")
  expect_identical(pipeline$prediction_type, "probability")
  expect_identical(pipeline$positive_class, 1L)
  expect_identical(pipeline$fallback, "forbidden")
})

test_that("binary AUC Top-K requires a predeclared positive K", {
  missing_error <- capture_sb_error(sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic"
  ))
  expect_s3_class(missing_error, "stormbreaker_missing_selector_parameter")

  invalid_error <- capture_sb_error(sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 0L
  ))
  expect_s3_class(invalid_error, "stormbreaker_bad_integer")
})

test_that("binary AUC Top-K is endpoint-specific and fold-local in glmnet", {
  survival_error <- capture_sb_error(sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "coxph",
    selector_n_features = 2L
  ))
  expect_s3_class(
    survival_error,
    "stormbreaker_incompatible_survival_selector"
  )

  nested <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "glmnet_logistic",
    selector_n_features = 2L,
    inner_folds = 3L,
    inner_seed = 101L
  )
  expect_identical(
    nested$id,
    paste0(
      "univ_binary_auc_top_k__glmnet_logistic__k2",
      "__alpha1__inner3__seed101"
    )
  )
})

test_that("the shared pipeline validator reconstructs every admitted family", {
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      learner = "coxph",
      selector_n_features = 2L
    ),
    sb_pipeline(
      learner = "glmnet_cox",
      glmnet_alpha = 0.5,
      inner_folds = 3L,
      inner_seed = 101L
    ),
    sb_pipeline(learner = "logistic"),
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "logistic",
      selector_n_features = 2L
    ),
    sb_pipeline(
      learner = "glmnet_logistic",
      glmnet_alpha = 0.5,
      inner_folds = 3L,
      inner_seed = 101L
    ),
    sb_pipeline(
      learner = "gaussian_naive_bayes",
      naive_bayes_var_smoothing = 1e-8
    ),
    sb_pipeline(
      learner = "lightgbm_binary",
      lightgbm_nrounds = 7L,
      lightgbm_learning_rate = 0.05,
      lightgbm_num_leaves = 9L,
      lightgbm_min_data_in_leaf = 3L,
      lightgbm_seed = 17L
    ),
    sb_pipeline(
      learner = "catboost_binary",
      catboost_iterations = 7L,
      catboost_learning_rate = 0.05,
      catboost_depth = 3L,
      catboost_l2_leaf_reg = 4,
      catboost_seed = 17L
    ),
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "ranger_binary",
      selector_n_features = 2L,
      ranger_num_trees = 7L,
      ranger_mtry = 1L,
      ranger_min_node_size = 10L,
      ranger_seed = 17L
    ),
    sb_pipeline(
      selector = "univ_cox_top_k",
      learner = "catboost_cox",
      selector_n_features = 2L,
      catboost_iterations = 7L,
      catboost_learning_rate = 0.05,
      catboost_depth = 3L,
      catboost_l2_leaf_reg = 4,
      catboost_seed = 17L
    ),
    sb_pipeline(
      selector = "univ_cox_top_k",
      learner = "ranger_survival",
      selector_n_features = 2L,
      ranger_num_trees = 7L,
      ranger_mtry = 1L,
      ranger_min_node_size = 3L,
      ranger_seed = 17L
    ),
    sb_pipeline(
      learner = "coxboost",
      coxboost_stepno = 7L,
      coxboost_penalty = 10
    ),
    sb_pipeline(
      learner = "plsrcox",
      plsrcox_components = 1L
    ),
    sb_pipeline(
      learner = "superpc",
      superpc_threshold = 0,
      superpc_components = 1L
    ),
    sb_pipeline(
      learner = "survivalsvm",
      survivalsvm_gamma_mu = 0.5
    )
  )

  for (pipeline in pipelines) {
    observed <- sb_validate_pipeline_contract(
      pipeline,
      pipeline$task_type
    )
    expect_identical(observed, pipeline)
  }
})

test_that("the shared pipeline validator rejects shortcuts and mutations", {
  list_error <- capture_sb_error(sb_validate_pipeline_contract(
    list(sb_pipeline()),
    "survival"
  ))
  expect_s3_class(list_error, "stormbreaker_bad_pipeline")

  mismatch_error <- capture_sb_error(sb_validate_pipeline_contract(
    sb_pipeline(learner = "logistic"),
    "survival"
  ))
  expect_s3_class(
    mismatch_error,
    "stormbreaker_pipeline_task_mismatch"
  )

  mutated <- sb_pipeline()
  mutated$risk_formula <- "different_formula"
  mutation_error <- capture_sb_error(sb_validate_pipeline_contract(
    mutated,
    "survival"
  ))
  expect_s3_class(mutation_error, "stormbreaker_unresolved_pipeline")
})

test_that("Gaussian Naive Bayes declares fixed variance smoothing", {
  pipeline <- sb_pipeline(learner = "gaussian_naive_bayes")

  expect_identical(pipeline$task_type, "binary")
  expect_identical(pipeline$selector, "none")
  expect_identical(pipeline$learner, "gaussian_naive_bayes")
  expect_identical(
    pipeline$id,
    "none__gaussian_naive_bayes__smooth1e-09"
  )
  expect_identical(
    pipeline$learner_params,
    list(
      var_smoothing = 1e-9,
      variance_estimator = "population",
      prior = "empirical"
    )
  )
  expect_identical(pipeline$prediction_type, "probability")
  expect_identical(pipeline$positive_class, 1L)
  expect_identical(pipeline$fallback, "forbidden")
  printed <- paste(capture.output(print(pipeline)), collapse = "\n")
  expect_match(printed, "Variance smoothing: 1e-09", fixed = TRUE)
})

test_that("Gaussian Naive Bayes IDs distinguish smoothing and selector", {
  no_smoothing <- sb_pipeline(
    learner = "gaussian_naive_bayes",
    naive_bayes_var_smoothing = 0
  )
  selected <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "gaussian_naive_bayes",
    selector_n_features = 3L,
    naive_bayes_var_smoothing = 1e-6
  )

  expect_identical(
    no_smoothing$id,
    "none__gaussian_naive_bayes__smooth0"
  )
  expect_identical(
    selected$id,
    "univ_binary_auc_top_k__gaussian_naive_bayes__k3__smooth1e-06"
  )
  expect_identical(selected$selector_params, list(n_features = 3L))
  expect_identical(selected$learner_params$var_smoothing, 1e-6)
})

test_that("Gaussian Naive Bayes IDs preserve adjacent smoothing values", {
  smoothing <- c(1.0000000000000001e-9, 1.0000000000000003e-9)
  pipelines <- lapply(smoothing, function(value) {
    sb_pipeline(
      learner = "gaussian_naive_bayes",
      naive_bayes_var_smoothing = value
    )
  })
  encoded <- sub("^.*__smooth", "", vapply(pipelines, `[[`, "", "id"))

  expect_false(identical(smoothing[[1L]], smoothing[[2L]]))
  expect_false(identical(pipelines[[1L]]$id, pipelines[[2L]]$id))
  expect_identical(as.numeric(encoded), smoothing)
})

test_that("Gaussian Naive Bayes parameters and selectors are explicit", {
  invalid_values <- list(-1e-9, NA_real_, Inf, "1e-9", c(1e-9, 1e-8))
  for (value in invalid_values) {
    error <- capture_sb_error(sb_pipeline(
      learner = "gaussian_naive_bayes",
      naive_bayes_var_smoothing = value
    ))
    expect_s3_class(
      error,
      "stormbreaker_bad_naive_bayes_var_smoothing"
    )
  }

  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "logistic",
    naive_bayes_var_smoothing = 1e-9
  ))
  expect_s3_class(
    irrelevant,
    "stormbreaker_irrelevant_naive_bayes_parameter"
  )

  survival_selector <- capture_sb_error(sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "gaussian_naive_bayes",
    selector_n_features = 1L
  ))
  expect_s3_class(
    survival_selector,
    "stormbreaker_incompatible_binary_selector"
  )
})

test_that("fixed LightGBM declares its complete binary configuration", {
  pipeline <- sb_pipeline(learner = "lightgbm_binary")

  expect_identical(pipeline$task_type, "binary")
  expect_identical(pipeline$selector, "none")
  expect_identical(pipeline$learner, "lightgbm_binary")
  expect_identical(
    pipeline$id,
    paste0(
      "none__lightgbm_binary__rounds100__lr0.1",
      "__leaves31__minleaf20__seed1"
    )
  )
  expect_identical(pipeline$learner_params$nrounds, 100L)
  expect_identical(pipeline$learner_params$learning_rate, 0.1)
  expect_identical(pipeline$learner_params$num_leaves, 31L)
  expect_identical(pipeline$learner_params$min_data_in_leaf, 20L)
  expect_identical(pipeline$learner_params$seed, 1L)
  expect_identical(pipeline$learner_params$objective, "binary")
  expect_identical(pipeline$learner_params$boosting, "gbdt")
  expect_identical(pipeline$learner_params$device_type, "cpu")
  expect_identical(pipeline$learner_params$num_threads, 1L)
  expect_identical(pipeline$learner_params$deterministic, TRUE)
  expect_identical(pipeline$learner_params$force_col_wise, TRUE)
  expect_identical(pipeline$learner_params$feature_fraction, 1)
  expect_identical(pipeline$learner_params$bagging_fraction, 1)
  expect_identical(pipeline$learner_params$bagging_freq, 0L)
  expect_identical(pipeline$prediction_type, "probability")
  expect_identical(pipeline$positive_class, 1L)
  expect_identical(pipeline$fallback, "forbidden")
  printed <- paste(capture.output(print(pipeline)), collapse = "\n")
  expect_match(printed, "Fixed LightGBM: 100 rounds", fixed = TRUE)
  expect_match(printed, "learning rate 0.1", fixed = TRUE)
})

test_that("fixed LightGBM IDs encode selector and declared parameters", {
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "lightgbm_binary",
    selector_n_features = 3L,
    lightgbm_nrounds = 12L,
    lightgbm_learning_rate = 0.025,
    lightgbm_num_leaves = 7L,
    lightgbm_min_data_in_leaf = 3L,
    lightgbm_seed = 99L
  )

  expect_identical(
    pipeline$id,
    paste0(
      "univ_binary_auc_top_k__lightgbm_binary__k3",
      "__rounds12__lr0.025__leaves7__minleaf3__seed99"
    )
  )
  expect_identical(pipeline$selector_params, list(n_features = 3L))
  expect_identical(pipeline$learner_params$nrounds, 12L)
  expect_identical(pipeline$learner_params$learning_rate, 0.025)
  expect_identical(pipeline$learner_params$num_leaves, 7L)
  expect_identical(pipeline$learner_params$min_data_in_leaf, 3L)
  expect_identical(pipeline$learner_params$seed, 99L)
})

test_that("fixed LightGBM IDs preserve adjacent learning rates", {
  rates <- c(0.10000000000000001, 0.10000000000000003)
  pipelines <- lapply(rates, function(rate) {
    sb_pipeline(
      learner = "lightgbm_binary",
      lightgbm_learning_rate = rate
    )
  })
  ids <- vapply(pipelines, `[[`, "", "id")
  encoded <- sub("^.*__lr", "", ids)
  encoded <- sub("__leaves.*$", "", encoded)

  expect_false(identical(rates[[1L]], rates[[2L]]))
  expect_false(identical(ids[[1L]], ids[[2L]]))
  expect_identical(as.numeric(encoded), rates)
})

test_that("fixed LightGBM parameters are validated and learner-specific", {
  integer_cases <- list(
    list(lightgbm_nrounds = 0L),
    list(lightgbm_num_leaves = 1L),
    list(lightgbm_min_data_in_leaf = -1L),
    list(lightgbm_seed = 0L),
    list(lightgbm_nrounds = 1.5)
  )
  for (arguments in integer_cases) {
    error <- capture_sb_error(do.call(
      sb_pipeline,
      c(list(learner = "lightgbm_binary"), arguments)
    ))
    expect_s3_class(error, "stormbreaker_bad_integer")
  }

  invalid_rates <- list(0, -0.1, NA_real_, Inf, "0.1", c(0.1, 0.2))
  for (rate in invalid_rates) {
    error <- capture_sb_error(sb_pipeline(
      learner = "lightgbm_binary",
      lightgbm_learning_rate = rate
    ))
    expect_s3_class(error, "stormbreaker_bad_lightgbm_learning_rate")
  }

  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "logistic",
    lightgbm_nrounds = 10L
  ))
  expect_s3_class(
    irrelevant,
    "stormbreaker_irrelevant_lightgbm_parameter"
  )

  survival_selector <- capture_sb_error(sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "lightgbm_binary",
    selector_n_features = 1L
  ))
  expect_s3_class(
    survival_selector,
    "stormbreaker_incompatible_binary_selector"
  )
})

test_that("fixed CatBoost declares its complete binary configuration", {
  pipeline <- sb_pipeline(learner = "catboost_binary")

  expect_identical(pipeline$task_type, "binary")
  expect_identical(pipeline$selector, "none")
  expect_identical(pipeline$learner, "catboost_binary")
  expect_identical(
    pipeline$id,
    paste0(
      "none__catboost_binary__iter100__lr0.1",
      "__depth6__l2reg3__seed1"
    )
  )
  expect_identical(pipeline$learner_params$iterations, 100L)
  expect_identical(pipeline$learner_params$learning_rate, 0.1)
  expect_identical(pipeline$learner_params$depth, 6L)
  expect_identical(pipeline$learner_params$l2_leaf_reg, 3)
  expect_identical(pipeline$learner_params$random_seed, 1L)
  expect_identical(pipeline$learner_params$loss_function, "Logloss")
  expect_identical(pipeline$learner_params$eval_metric, "Logloss")
  expect_identical(pipeline$learner_params$task_type, "CPU")
  expect_identical(pipeline$learner_params$thread_count, 1L)
  expect_identical(pipeline$learner_params$bootstrap_type, "Bayesian")
  expect_identical(pipeline$learner_params$bagging_temperature, 1)
  expect_identical(pipeline$learner_params$random_strength, 1)
  expect_identical(pipeline$learner_params$rsm, 1)
  expect_identical(pipeline$learner_params$grow_policy, "SymmetricTree")
  expect_identical(pipeline$learner_params$boosting_type, "Plain")
  expect_identical(pipeline$learner_params$use_best_model, FALSE)
  expect_identical(pipeline$learner_params$allow_writing_files, FALSE)
  expect_identical(pipeline$learner_params$save_snapshot, FALSE)
  expect_identical(pipeline$prediction_type, "probability")
  expect_identical(pipeline$positive_class, 1L)
  expect_identical(pipeline$fallback, "forbidden")
  printed <- paste(capture.output(print(pipeline)), collapse = "\n")
  expect_match(printed, "Fixed CatBoost: 100 iterations", fixed = TRUE)
  expect_match(printed, "learning rate 0.1", fixed = TRUE)
})

test_that("fixed CatBoost IDs encode selector and declared parameters", {
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "catboost_binary",
    selector_n_features = 3L,
    catboost_iterations = 12L,
    catboost_learning_rate = 0.025,
    catboost_depth = 4L,
    catboost_l2_leaf_reg = 2.5,
    catboost_seed = 99L
  )

  expect_identical(
    pipeline$id,
    paste0(
      "univ_binary_auc_top_k__catboost_binary__k3",
      "__iter12__lr0.025__depth4__l2reg2.5__seed99"
    )
  )
  expect_identical(pipeline$selector_params, list(n_features = 3L))
  expect_identical(pipeline$learner_params$iterations, 12L)
  expect_identical(pipeline$learner_params$learning_rate, 0.025)
  expect_identical(pipeline$learner_params$depth, 4L)
  expect_identical(pipeline$learner_params$l2_leaf_reg, 2.5)
  expect_identical(pipeline$learner_params$random_seed, 99L)
})

test_that("fixed CatBoost IDs preserve adjacent numeric settings", {
  rates <- c(0.10000000000000001, 0.10000000000000003)
  regularization <- c(3, 3.0000000000000004)
  pipelines <- Map(function(rate, l2) {
    sb_pipeline(
      learner = "catboost_binary",
      catboost_learning_rate = rate,
      catboost_l2_leaf_reg = l2
    )
  }, rates, regularization)
  ids <- vapply(pipelines, `[[`, "", "id")
  encoded_rates <- sub("^.*__lr", "", ids)
  encoded_rates <- sub("__depth.*$", "", encoded_rates)
  encoded_l2 <- sub("^.*__l2reg", "", ids)
  encoded_l2 <- sub("__seed.*$", "", encoded_l2)

  expect_false(identical(rates[[1L]], rates[[2L]]))
  expect_false(identical(regularization[[1L]], regularization[[2L]]))
  expect_false(identical(ids[[1L]], ids[[2L]]))
  expect_identical(as.numeric(encoded_rates), rates)
  expect_identical(as.numeric(encoded_l2), regularization)
})

test_that("fixed CatBoost parameters are validated and learner-specific", {
  integer_cases <- list(
    list(catboost_iterations = 0L),
    list(catboost_depth = 0L),
    list(catboost_seed = 0L),
    list(catboost_iterations = 1.5)
  )
  for (arguments in integer_cases) {
    error <- capture_sb_error(do.call(
      sb_pipeline,
      c(list(learner = "catboost_binary"), arguments)
    ))
    expect_s3_class(error, "stormbreaker_bad_integer")
  }

  depth_error <- capture_sb_error(sb_pipeline(
    learner = "catboost_binary",
    catboost_depth = 17L
  ))
  expect_s3_class(depth_error, "stormbreaker_bad_catboost_depth")

  invalid_rates <- list(0, -0.1, NA_real_, Inf, "0.1", c(0.1, 0.2))
  for (rate in invalid_rates) {
    error <- capture_sb_error(sb_pipeline(
      learner = "catboost_binary",
      catboost_learning_rate = rate
    ))
    expect_s3_class(error, "stormbreaker_bad_catboost_learning_rate")
  }

  invalid_l2 <- list(0, -1, NA_real_, Inf, "3", c(2, 3))
  for (l2 in invalid_l2) {
    error <- capture_sb_error(sb_pipeline(
      learner = "catboost_binary",
      catboost_l2_leaf_reg = l2
    ))
    expect_s3_class(error, "stormbreaker_bad_catboost_l2_leaf_reg")
  }

  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "logistic",
    catboost_iterations = 10L
  ))
  expect_s3_class(
    irrelevant,
    "stormbreaker_irrelevant_catboost_parameter"
  )

  survival_selector <- capture_sb_error(sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "catboost_binary",
    selector_n_features = 1L
  ))
  expect_s3_class(
    survival_selector,
    "stormbreaker_incompatible_binary_selector"
  )
})
