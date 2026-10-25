final_fit_contains_function <- function(value) {
  if (is.function(value) || is.environment(value) ||
    typeof(value) == "externalptr") {
    return(TRUE)
  }
  if (!is.list(value)) {
    return(FALSE)
  }
  any(vapply(unclass(value), final_fit_contains_function, logical(1)))
}

final_state_oracle <- function(state, x) {
  switch(state$learner,
    coxph = sb_predict_cox_learner_state(state, x),
    glmnet_cox = sb_predict_cox_learner_state(state, x),
    coxboost = sb_predict_cox_learner_state(state, x),
    plsrcox = sb_predict_cox_learner_state(state, x),
    superpc = sb_predict_cox_learner_state(state, x),
    survivalsvm = sb_predict_cox_learner_state(state, x),
    logistic = sb_predict_logistic_learner_state(state, x),
    glmnet_logistic = sb_predict_logistic_learner_state(state, x),
    gaussian_naive_bayes = sb_predict_naive_bayes_learner_state(state, x),
    lightgbm_binary = sb_predict_lightgbm_learner_state(state, x),
    catboost_binary = sb_predict_catboost_learner_state(state, x),
    catboost_cox = sb_predict_catboost_cox_learner_state(state, x),
    gbm_cox = sb_predict_gbm_cox_learner_state(state, x),
    ranger_binary = sb_predict_ranger_binary_learner_state(state, x),
    adaboost_binary =
      sb_predict_classical_binary_boosting_learner_state(state, x),
    logitboost_binary =
      sb_predict_classical_binary_boosting_learner_state(state, x),
    ranger_survival = sb_predict_ranger_survival_learner_state(state, x)
  )
}

expect_final_prediction_failure <- function(result, error_class, n_rows) {
  expect_s3_class(result, "stormbreaker_final_prediction")
  expect_identical(result$status, "failed")
  expect_identical(result$error_class, error_class)
  expect_identical(nrow(result$predictions), as.integer(n_rows))
  if (n_rows > 0L) {
    expect_true(all(result$predictions$status == "failed"))
    expect_true(all(result$predictions$error_class == error_class))
  }
}

test_that("ordinary Cox final fit uses every development patient", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_b", "gene_a"),
    "month"
  )
  pipeline <- sb_pipeline()
  fit <- sb_fit_final(task, pipeline)

  expect_s3_class(fit, "stormbreaker_final_fit")
  expect_identical(fit$status, "success")
  expect_identical(fit$task_type, "survival")
  expect_identical(fit$pipeline, pipeline)
  expect_identical(
    fit$development$patient_ids,
    as.character(d$patient_id)
  )
  expect_identical(fit$development$n, 60L)
  expect_identical(fit$development$independent_unit, "patient")
  expect_identical(
    fit$development$candidate_features,
    c("gene_b", "gene_a")
  )
  expect_identical(fit$development$endpoint$n_events, task$n_events)
  expect_identical(fit$development$endpoint$n_censored, task$n_censored)
  expect_identical(fit$selection$status, "success")
  expect_identical(
    fit$selection$selected_features,
    c("gene_b", "gene_a")
  )
  expect_s3_class(fit$selection$state, "stormbreaker_selector_fit")
  expect_identical(fit$tuning$status, "not_required")
  expect_identical(fit$learner$status, "success")
  expect_s3_class(fit$learner_state, "stormbreaker_learner_state")
  expect_identical(fit$required_features, c("gene_b", "gene_a"))
  expect_identical(fit$prediction$risk_formula, "x_beta")
  expect_identical(fit$prediction$risk_reference, "zero")
  expect_true(all(c(
    "package", "package_version", "r_version", "survival_version"
  ) %in% names(fit$provenance)))
  expect_false(any(c(
    "metrics", "folds", "comparison", "selection_stability",
    "predictions", "runs", "task", "data"
  ) %in% names(fit)))
  expect_false("data" %in% names(fit$development))
  expect_false(final_fit_contains_function(fit))

  printed <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(printed, "final development fit", fixed = TRUE)
  expect_match(printed, pipeline$id, fixed = TRUE)
  expect_match(printed, "survival", fixed = TRUE)
  expect_match(printed, "Required features: 2", fixed = TRUE)
  expect_false(grepl("AUC|Brier|concordance|validation performance", printed))
})

test_that("Top-K Cox final selection is fitted on all development rows", {
  d <- make_sb_fixture(n = 70L, seed = 20260820L)
  predictors <- c("gene_c", "gene_a", "gene_b")
  task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month"
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "coxph",
    selector_n_features = 2L
  )
  fit <- sb_fit_final(task, pipeline)
  direct <- sb_fit_univ_cox_top_k(
    d[predictors],
    d$time,
    as.integer(d$event),
    n_features = 2L
  )

  expect_identical(fit$status, "success")
  expect_identical(
    fit$selection$selected_features,
    direct$selected_features
  )
  expect_identical(fit$selection$scores, direct$scores)
  expect_identical(fit$selection$state, direct)
  expect_identical(fit$required_features, direct$selected_features)

  new_data <- d[
    1:4,
    c("patient_id", fit$required_features),
    drop = FALSE
  ]
  prediction <- sb_predict_final(fit, new_data, "patient_id")
  expect_identical(prediction$status, "success")
  expect_identical(nrow(prediction$predictions), 4L)
})

test_that("glmnet Cox final tuning covers the complete development cohort", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60L, seed = 20260821L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month"
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 1701L
  )
  fit <- sb_fit_final(task, pipeline)
  direct <- sb_prepare_runtime_learner(
    sb_resolve_pipeline_runtime(pipeline)$learner,
    pipeline,
    d[predictors],
    list(
      task_type = "survival",
      time = d$time,
      event = as.integer(d$event)
    ),
    d$patient_id
  )$tuning

  expect_identical(fit$status, "success")
  expect_identical(fit$tuning$status, "success")
  expect_identical(fit$tuning$lambda, direct$lambda)
  expect_identical(fit$tuning$curve, direct$curve)
  expect_identical(fit$tuning$fold_id, direct$fold_id)
  expect_identical(names(fit$tuning$fold_id), as.character(d$patient_id))
  expect_identical(length(fit$tuning$fold_id), nrow(d))
  expect_true(fit$tuning$lambda %in% fit$tuning$curve$lambda)
  expect_identical(
    fit$provenance$glmnet_version,
    as.character(utils::packageVersion("glmnet"))
  )
})

test_that("core binary learners produce portable final states", {
  d <- make_sb_binary_fixture(n = 80L, seed = 20260822L)
  predictors <- c("gene_b", "gene_a", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(
      learner = "gaussian_naive_bayes",
      naive_bayes_var_smoothing = 1e-8
    )
  )

  for (pipeline in pipelines) {
    fit <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success", info = pipeline$id)
    expect_identical(fit$task_type, "binary", info = pipeline$id)
    expect_identical(fit$pipeline, pipeline, info = pipeline$id)
    expect_identical(fit$selection$status, "success", info = pipeline$id)
    expect_identical(fit$tuning$status, "not_required", info = pipeline$id)
    expect_identical(fit$learner$status, "success", info = pipeline$id)
    expect_s3_class(
      fit$learner_state,
      "stormbreaker_learner_state"
    )
    expect_identical(fit$required_features, predictors, info = pipeline$id)
    expect_identical(fit$prediction$positive_class, 1L, info = pipeline$id)
    expect_false(final_fit_contains_function(fit), info = pipeline$id)
  }
})

test_that("glmnet logistic final tuning uses all development patients", {
  skip_if_not_installed("glmnet")

  d <- make_sb_binary_fixture(n = 80L, seed = 20260823L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 0.25,
    inner_folds = 4L,
    inner_seed = 1702L,
    glmnet_standardize = FALSE,
    glmnet_lambda_rule = "lambda.1se",
    glmnet_tuning_measure = "auc"
  )
  fit <- sb_fit_final(task, pipeline)
  direct <- sb_prepare_runtime_learner(
    sb_resolve_pipeline_runtime(pipeline)$learner,
    pipeline,
    d[predictors],
    list(task_type = "binary", outcome = as.integer(d$outcome)),
    d$patient_id
  )$tuning

  expect_identical(fit$status, "success")
  expect_identical(fit$tuning$status, "success")
  expect_identical(fit$tuning$lambda, direct$lambda)
  expect_identical(fit$tuning$curve, direct$curve)
  expect_identical(fit$tuning$fold_id, direct$fold_id)
  expect_identical(fit$tuning$standardize, FALSE)
  expect_identical(fit$tuning$selection_rule, "lambda.1se")
  expect_identical(fit$tuning$measure, "auc")
  expect_identical(fit$learner_state$params$standardize, FALSE)
  expect_identical(names(fit$tuning$fold_id), as.character(d$patient_id))
  expect_identical(length(fit$tuning$fold_id), nrow(d))
  expect_true(fit$tuning$lambda %in% fit$tuning$curve$lambda)
})

test_that("fixed tree learners produce portable final states", {
  d <- make_sb_binary_fixture(n = 80L, seed = 20260824L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)

  if (sb_lightgbm_available()) {
    pipeline <- sb_pipeline(
      learner = "lightgbm_binary",
      lightgbm_nrounds = 5L,
      lightgbm_num_leaves = 7L,
      lightgbm_min_data_in_leaf = 3L,
      lightgbm_seed = 101L
    )
    fit <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success")
    expect_type(fit$learner_state$model_string, "character")
    expect_false("booster" %in% names(fit$learner_state))
    expect_false(final_fit_contains_function(fit))
    expect_identical(
      fit$provenance$lightgbm_version,
      as.character(utils::packageVersion("lightgbm"))
    )
  }

  if (sb_catboost_available()) {
    pipeline <- sb_pipeline(
      learner = "catboost_binary",
      catboost_iterations = 5L,
      catboost_depth = 3L,
      catboost_seed = 101L
    )
    fit <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success")
    expect_type(fit$learner_state$model_raw, "raw")
    expect_false("booster" %in% names(fit$learner_state))
    expect_false(final_fit_contains_function(fit))
    expect_identical(
      fit$provenance$catboost_version,
      as.character(utils::packageVersion("catboost"))
    )
  }
})

test_that("final fitting requires one explicit task-matched pipeline", {
  d <- make_sb_fixture(n = 40L, seed = 20260825L)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month"
  )

  expect_s3_class(
    capture_sb_error(sb_fit_final(task, NULL)),
    "stormbreaker_bad_pipeline"
  )
  expect_s3_class(
    capture_sb_error(sb_fit_final(task, list(sb_pipeline()))),
    "stormbreaker_bad_pipeline"
  )
  expect_s3_class(
    capture_sb_error(sb_fit_final(
      task,
      structure(list(), class = "stormbreaker_cv_result")
    )),
    "stormbreaker_bad_pipeline"
  )
  expect_s3_class(
    capture_sb_error(sb_fit_final(
      task,
      sb_pipeline(learner = "logistic")
    )),
    "stormbreaker_pipeline_task_mismatch"
  )
})

test_that("final-fit failures retain completed earlier stages", {
  selector_data <- make_sb_fixture(n = 30L, seed = 20260826L)
  selector_data$gene_a <- 1
  selector_data$gene_b <- 1
  selector_task <- sb_surv_task(
    selector_data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month"
  )
  selector_fit <- sb_fit_final(
    selector_task,
    sb_pipeline(
      selector = "univ_cox_top_k",
      learner = "coxph",
      selector_n_features = 1L
    )
  )
  expect_identical(selector_fit$status, "failed")
  expect_identical(selector_fit$error_stage, "selection")
  expect_identical(selector_fit$selection$status, "failed")
  expect_identical(selector_fit$learner$status, "not_fitted")
  expect_null(selector_fit$learner_state)

  tuning_data <- make_sb_binary_fixture(n = 10L, seed = 20260827L)
  tuning_task <- sb_binary_task(
    tuning_data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  tuning_fit <- sb_fit_final(
    tuning_task,
    sb_pipeline(
      learner = "glmnet_logistic",
      inner_folds = 6L,
      inner_seed = 19L
    )
  )
  expect_identical(tuning_fit$status, "failed")
  expect_identical(tuning_fit$error_stage, "tuning")
  expect_identical(tuning_fit$selection$status, "success")
  expect_identical(tuning_fit$tuning$status, "failed")
  expect_identical(tuning_fit$learner$status, "not_fitted")
  expect_null(tuning_fit$learner_state)

  missing_pipeline <- sb_pipeline(
    learner = "lightgbm_binary",
    lightgbm_nrounds = 3L,
    lightgbm_num_leaves = 5L,
    lightgbm_min_data_in_leaf = 2L
  )
  testthat::local_mocked_bindings(
    sb_lightgbm_available = function() FALSE
  )
  missing_fit <- sb_fit_final(tuning_task, missing_pipeline)
  expect_identical(missing_fit$status, "failed")
  expect_identical(missing_fit$error_stage, "learner_fit")
  expect_identical(missing_fit$selection$status, "success")
  expect_identical(missing_fit$learner$status, "failed")
  expect_identical(
    missing_fit$error_class,
    "stormbreaker_missing_lightgbm"
  )
})

test_that("learner and state failures remain distinct", {
  d <- make_sb_fixture(n = 12L, seed = 20260828L)
  withr::with_seed(20260828L, {
    for (index in seq_len(12L)) {
      d[[sprintf("x_%02d", index)]] <- stats::rnorm(nrow(d))
    }
  })
  predictors <- sprintf("x_%02d", seq_len(12L))
  task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month"
  )
  learner_fit <- sb_fit_final(task, sb_pipeline())

  expect_identical(learner_fit$status, "failed")
  expect_identical(learner_fit$error_stage, "learner_fit")
  expect_identical(learner_fit$selection$status, "success")
  expect_identical(learner_fit$learner$status, "failed")
  expect_null(learner_fit$learner_state)

  ordinary_task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month"
  )
  pipeline <- sb_pipeline()
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  runtime$learner$extract_state <- function(model, context = list()) {
    sb_abort(
      "Synthetic state extraction failure.",
      "stormbreaker_test_state_extraction_failure"
    )
  }
  testthat::local_mocked_bindings(
    sb_resolve_pipeline_runtime = function(pipeline) runtime
  )
  state_fit <- sb_fit_final(ordinary_task, pipeline)

  expect_identical(state_fit$status, "failed")
  expect_identical(state_fit$error_stage, "state_extraction")
  expect_identical(state_fit$selection$status, "success")
  expect_identical(state_fit$learner$status, "success")
  expect_identical(state_fit$required_features, c("gene_a", "gene_b"))
  expect_null(state_fit$learner_state)
  expect_identical(
    state_fit$error_class,
    "stormbreaker_test_state_extraction_failure"
  )
})

test_that("final prediction has an endpoint-blind three-argument interface", {
  expect_identical(
    names(formals(sb_predict_final)),
    c("fit", "new_data", "id_col")
  )
})

test_that("survival final prediction uses only saved state", {
  d <- make_sb_fixture(n = 60L, seed = 20260829L)
  predictors <- c("gene_b", "gene_a")
  rows <- c(9L, 2L, 14L, 5L)

  risk_task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month"
  )
  pipelines <- list(sb_pipeline())
  if (sb_glmnet_available()) {
    pipelines[[length(pipelines) + 1L]] <- sb_pipeline(
      learner = "glmnet_cox",
      glmnet_alpha = 0.5,
      inner_folds = 3L,
      inner_seed = 2026L
    )
  }

  for (pipeline in pipelines) {
    fit <- sb_fit_final(risk_task, pipeline)
    new_data <- d[rows, c("patient_id", rev(predictors)), drop = FALSE]
    prediction <- sb_predict_final(fit, new_data, "patient_id")
    oracle <- final_state_oracle(
      fit$learner_state,
      new_data[fit$required_features]
    )

    expect_identical(prediction$status, "success", info = pipeline$id)
    expect_identical(
      prediction$predictions$patient_id,
      as.character(d$patient_id[rows]),
      info = pipeline$id
    )
    expect_true(all(
      prediction$predictions$pipeline_id == pipeline$id
    ), info = pipeline$id)
    expect_equal(
      prediction$predictions$risk,
      oracle$risk,
      tolerance = 0,
      info = pipeline$id
    )
    expect_false(any(c("time", "event", "outcome") %in%
      names(prediction$predictions)))

    restored <- unserialize(serialize(fit, NULL))
    restored_prediction <- sb_predict_final(
      restored,
      new_data,
      "patient_id"
    )
    expect_identical(restored_prediction, prediction, info = pipeline$id)
  }

  horizon_task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 12
  )
  horizon_fit <- sb_fit_final(horizon_task, sb_pipeline())
  horizon_data <- d[rows, c("patient_id", predictors), drop = FALSE]
  horizon_prediction <- sb_predict_final(
    horizon_fit,
    horizon_data,
    "patient_id"
  )
  risk <- unname(as.vector(
    as.matrix(horizon_data[predictors]) %*%
      horizon_fit$learner_state$coefficient
  ))
  probability <- exp(
    -horizon_fit$learner_state$baseline_hazard$
      cumulative_baseline_hazard * exp(risk)
  )

  expect_identical(horizon_prediction$status, "success")
  expect_equal(horizon_prediction$predictions$risk, risk, tolerance = 1e-12)
  expect_equal(
    horizon_prediction$predictions$survival_probability,
    probability,
    tolerance = 1e-12
  )
  expect_true(all(
    horizon_prediction$predictions$prediction_horizon == 12
  ))
  expect_true(all(horizon_prediction$predictions$time_unit == "month"))
})

test_that("declared survival outputs remain bound to their saved state", {
  d <- make_sb_fixture(n = 60L, seed = 20260907L)
  predictors <- c("gene_a", "gene_b")
  horizon_task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 12
  )
  horizon_fit <- sb_fit_final(horizon_task, sb_pipeline())
  new_data <- d[1:4, c("patient_id", predictors)]

  missing_baseline <- horizon_fit
  missing_baseline$learner_state$baseline_hazard <- NULL
  expect_final_prediction_failure(
    sb_predict_final(missing_baseline, new_data, "patient_id"),
    "stormbreaker_final_state_semantics_mismatch",
    4L
  )

  missing_unit <- horizon_fit
  missing_unit$learner_state$baseline_hazard$time_unit <- NULL
  expect_final_prediction_failure(
    sb_predict_final(missing_unit, new_data, "patient_id"),
    "stormbreaker_final_state_semantics_mismatch",
    4L
  )

  wrong_horizon <- horizon_fit
  wrong_horizon$learner_state$baseline_hazard$prediction_horizon <- 6
  expect_final_prediction_failure(
    sb_predict_final(wrong_horizon, new_data, "patient_id"),
    "stormbreaker_final_state_semantics_mismatch",
    4L
  )

  for (bad_hazard in list(-1, NA_real_, Inf)) {
    invalid_hazard <- horizon_fit
    invalid_hazard$learner_state$baseline_hazard$
      cumulative_baseline_hazard <- bad_hazard
    expect_final_prediction_failure(
      sb_predict_final(invalid_hazard, new_data, "patient_id"),
      "stormbreaker_final_state_semantics_mismatch",
      4L
    )
  }

  risk_task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month"
  )
  risk_fit <- sb_fit_final(risk_task, sb_pipeline())
  risk_fit$learner_state$baseline_hazard <-
    horizon_fit$learner_state$baseline_hazard
  expect_final_prediction_failure(
    sb_predict_final(risk_fit, new_data, "patient_id"),
    "stormbreaker_final_state_semantics_mismatch",
    4L
  )
})

test_that("every binary final learner predicts class-one probability", {
  d <- make_sb_binary_fixture(n = 80L, seed = 20260830L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(
      learner = "gaussian_naive_bayes",
      naive_bayes_var_smoothing = 1e-8
    )
  )
  if (sb_glmnet_available()) {
    pipelines[[length(pipelines) + 1L]] <- sb_pipeline(
      learner = "glmnet_logistic",
      glmnet_alpha = 0.5,
      inner_folds = 4L,
      inner_seed = 303L
    )
  }
  if (sb_lightgbm_available()) {
    pipelines[[length(pipelines) + 1L]] <- sb_pipeline(
      learner = "lightgbm_binary",
      lightgbm_nrounds = 5L,
      lightgbm_num_leaves = 7L,
      lightgbm_min_data_in_leaf = 3L,
      lightgbm_seed = 304L
    )
  }
  if (sb_catboost_available()) {
    pipelines[[length(pipelines) + 1L]] <- sb_pipeline(
      learner = "catboost_binary",
      catboost_iterations = 5L,
      catboost_depth = 3L,
      catboost_seed = 305L
    )
  }

  rows <- c(11L, 3L, 8L, 1L, 17L)
  new_data <- d[rows, c("patient_id", rev(predictors)), drop = FALSE]
  for (pipeline in pipelines) {
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(fit, new_data, "patient_id")
    oracle <- final_state_oracle(
      fit$learner_state,
      new_data[fit$required_features]
    )

    expect_identical(prediction$status, "success", info = pipeline$id)
    expect_identical(
      prediction$predictions$patient_id,
      as.character(d$patient_id[rows]),
      info = pipeline$id
    )
    expect_true(all(
      prediction$predictions$pipeline_id == pipeline$id
    ), info = pipeline$id)
    expect_equal(
      prediction$predictions$probability,
      oracle$probability,
      tolerance = 0,
      info = pipeline$id
    )
    expect_true(all(prediction$predictions$status == "success"))
    expect_false(any(c("outcome", "time", "event") %in%
      names(prediction$predictions)))

    restored <- unserialize(serialize(fit, NULL))
    restored_prediction <- sb_predict_final(
      restored,
      new_data,
      "patient_id"
    )
    expect_identical(restored_prediction, prediction, info = pipeline$id)
  }
})

test_that("binary decision thresholds propagate without changing probabilities", {
  d <- make_sb_binary_fixture(n = 80L, seed = 20260908L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  thresholds <- c(0.10, 0.25, 0.50)
  plain_task <- sb_binary_task(
    d, "patient_id", "outcome", predictors
  )
  threshold_task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    predictors,
    decision_thresholds = thresholds
  )
  pipeline <- sb_pipeline(learner = "logistic")
  plain_fit <- sb_fit_final(plain_task, pipeline)
  threshold_fit <- sb_fit_final(threshold_task, pipeline)
  new_data <- d[1:8, c("patient_id", predictors), drop = FALSE]
  plain_prediction <- sb_predict_final(
    plain_fit, new_data, "patient_id"
  )
  threshold_prediction <- sb_predict_final(
    threshold_fit, new_data, "patient_id"
  )

  expect_identical(
    threshold_fit$development$endpoint$decision_thresholds,
    thresholds
  )
  expect_identical(
    threshold_fit$prediction$decision_thresholds,
    thresholds
  )
  expect_identical(
    threshold_prediction$decision_thresholds,
    thresholds
  )
  expect_equal(
    threshold_prediction$predictions$probability,
    plain_prediction$predictions$probability,
    tolerance = 0
  )
  expect_false("decision_thresholds" %in% names(plain_prediction))
  expect_false(any(c(
    "predicted_class", "decision", "action", "treat", "selected"
  ) %in% names(threshold_prediction$predictions)))
})

test_that("binary final-fit threshold copies must remain valid and consistent", {
  d <- make_sb_binary_fixture(n = 60L, seed = 20260909L)
  predictors <- c("gene_a", "gene_b")
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    predictors,
    decision_thresholds = c(0.2, 0.5)
  )
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  new_data <- d[1:5, c("patient_id", predictors), drop = FALSE]

  missing_copy <- fit
  missing_copy$development$endpoint$decision_thresholds <- NULL
  expect_final_prediction_failure(
    sb_predict_final(missing_copy, new_data, "patient_id"),
    "stormbreaker_final_state_semantics_mismatch",
    5L
  )

  inconsistent_copy <- fit
  inconsistent_copy$development$endpoint$decision_thresholds <- c(0.2, 0.6)
  expect_final_prediction_failure(
    sb_predict_final(inconsistent_copy, new_data, "patient_id"),
    "stormbreaker_final_state_semantics_mismatch",
    5L
  )

  invalid_copy <- fit
  invalid_copy$prediction$decision_thresholds <- c(0.2, 1)
  expect_final_prediction_failure(
    sb_predict_final(invalid_copy, new_data, "patient_id"),
    "stormbreaker_bad_decision_thresholds",
    5L
  )
})

test_that("irrelevant endpoint columns are never inspected", {
  binary <- make_sb_binary_fixture(n = 60L, seed = 20260831L)
  binary_task <- sb_binary_task(
    binary,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  binary_fit <- sb_fit_final(
    binary_task,
    sb_pipeline(learner = "logistic")
  )
  binary_new <- binary[1:8, c("patient_id", "gene_a", "gene_b")]
  binary_new$outcome <- binary$outcome[1:8]
  binary_new$time <- seq_len(8L)
  binary_new$event <- rep(c(0L, 1L), 4L)
  original <- sb_predict_final(binary_fit, binary_new, "patient_id")
  changed <- binary_new
  changed$outcome <- c(NA, Inf, "changed", rep("irrelevant", 5L))
  changed$time <- list(NULL, NaN, Inf, -Inf, "x", NA, list(1), factor("a"))
  changed$event <- rev(changed$outcome)
  observed <- sb_predict_final(binary_fit, changed, "patient_id")

  expect_identical(observed, original)

  survival <- make_sb_fixture(n = 60L, seed = 20260901L)
  survival_task <- sb_surv_task(
    survival,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 12
  )
  survival_fit <- sb_fit_final(survival_task, sb_pipeline())
  survival_new <- survival[1:8, c(
    "patient_id", "gene_a", "gene_b", "time", "event"
  )]
  first <- sb_predict_final(survival_fit, survival_new, "patient_id")
  survival_new$time <- rep(NA_character_, 8L)
  survival_new$event <- rep(Inf, 8L)
  second <- sb_predict_final(survival_fit, survival_new, "patient_id")

  expect_identical(second, first)
})

test_that("prediction restores fitted feature order and keeps row order", {
  d <- make_sb_binary_fixture(n = 60L, seed = 20260902L)
  names(d)[names(d) == "gene_a"] <- "gene A"
  names(d)[names(d) == "gene_b"] <- "gene-B"
  d[["gene-B"]] <- as.integer(seq_len(nrow(d)) %% 7L)
  predictors <- c("gene A", "gene-B")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  rows <- c(12L, 1L, 7L, 4L)
  new_data <- d[rows, c("gene-B", "patient_id", "gene A"), drop = FALSE]
  prediction <- sb_predict_final(fit, new_data, "patient_id")
  oracle <- sb_predict_logistic_learner_state(
    fit$learner_state,
    new_data[predictors]
  )

  expect_identical(prediction$status, "success")
  expect_identical(
    prediction$predictions$patient_id,
    as.character(d$patient_id[rows])
  )
  expect_equal(
    prediction$predictions$probability,
    oracle$probability,
    tolerance = 0
  )
})

test_that("invalid prediction inputs fail visibly without fallback", {
  d <- make_sb_binary_fixture(n = 60L, seed = 20260903L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  valid <- d[1:6, c("patient_id", "gene_a", "gene_b")]

  missing <- valid[c("patient_id", "gene_a")]
  expect_final_prediction_failure(
    sb_predict_final(fit, missing, "patient_id"),
    "stormbreaker_missing_final_features",
    6L
  )

  duplicate_feature <- valid
  duplicate_feature$copy <- duplicate_feature$gene_a
  names(duplicate_feature)[names(duplicate_feature) == "copy"] <- "gene_a"
  expect_final_prediction_failure(
    sb_predict_final(fit, duplicate_feature, "patient_id"),
    "stormbreaker_duplicate_final_feature_columns",
    6L
  )

  duplicate_id_column <- valid
  duplicate_id_column$copy <- duplicate_id_column$patient_id
  names(duplicate_id_column)[names(duplicate_id_column) == "copy"] <-
    "patient_id"
  expect_final_prediction_failure(
    sb_predict_final(fit, duplicate_id_column, "patient_id"),
    "stormbreaker_duplicate_prediction_id_column",
    0L
  )

  duplicate_ids <- valid
  duplicate_ids$patient_id[2L] <- duplicate_ids$patient_id[1L]
  expect_final_prediction_failure(
    sb_predict_final(fit, duplicate_ids, "patient_id"),
    "stormbreaker_duplicate_prediction_id",
    0L
  )

  for (bad_id in list(NA_character_, "")) {
    bad <- valid
    bad$patient_id[2L] <- bad_id
    expect_final_prediction_failure(
      sb_predict_final(fit, bad, "patient_id"),
      "stormbreaker_bad_prediction_id",
      0L
    )
  }

  nonnumeric <- valid
  nonnumeric$gene_a <- as.character(nonnumeric$gene_a)
  expect_final_prediction_failure(
    sb_predict_final(fit, nonnumeric, "patient_id"),
    "stormbreaker_bad_adapter_predictor_type",
    6L
  )

  for (bad_value in list(NA_real_, Inf)) {
    invalid <- valid
    invalid$gene_a[1L] <- bad_value
    expect_final_prediction_failure(
      sb_predict_final(fit, invalid, "patient_id"),
      "stormbreaker_bad_adapter_predictor_values",
      6L
    )
  }
})

test_that("failed fits and malformed states remain failed predictions", {
  d <- make_sb_binary_fixture(n = 10L, seed = 20260904L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  failed_fit <- sb_fit_final(
    task,
    sb_pipeline(
      learner = "glmnet_logistic",
      inner_folds = 6L,
      inner_seed = 29L
    )
  )
  new_data <- d[1:4, c("patient_id", "gene_a", "gene_b")]
  expect_final_prediction_failure(
    sb_predict_final(failed_fit, new_data, "patient_id"),
    "stormbreaker_failed_final_fit",
    4L
  )

  good_fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  bad_state <- good_fit
  bad_state$learner_state$coefficient[1L] <- NA_real_
  expect_final_prediction_failure(
    sb_predict_final(bad_state, new_data, "patient_id"),
    "stormbreaker_bad_learner_state",
    4L
  )

  testthat::local_mocked_bindings(
    sb_predict_final_state = function(state, x) {
      list(probability = rep.int(2, nrow(x)))
    }
  )
  invalid_output <- sb_predict_final(good_fit, new_data, "patient_id")
  expect_final_prediction_failure(
    invalid_output,
    "stormbreaker_invalid_final_probability",
    4L
  )
})

test_that("a corrupt saved tree fails at state loading", {
  skip_if_not_installed("lightgbm")

  d <- make_sb_binary_fixture(n = 60L, seed = 20260905L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  fit <- sb_fit_final(
    task,
    sb_pipeline(
      learner = "lightgbm_binary",
      lightgbm_nrounds = 3L,
      lightgbm_num_leaves = 5L,
      lightgbm_min_data_in_leaf = 2L
    )
  )
  fit$learner_state$model_string <- "not a LightGBM model"
  new_data <- d[1:4, c("patient_id", "gene_a", "gene_b")]
  prediction <- sb_predict_final(fit, new_data, "patient_id")

  expect_final_prediction_failure(
    prediction,
    "stormbreaker_invalid_lightgbm_state",
    4L
  )
})

test_that("final prediction printing reports no validation performance", {
  d <- make_sb_binary_fixture(n = 40L, seed = 20260906L)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b")
  )
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  prediction <- sb_predict_final(
    fit,
    d[1:5, c("patient_id", "gene_a", "gene_b")],
    "patient_id"
  )
  printed <- paste(capture.output(print(prediction)), collapse = "\n")

  expect_match(printed, "final state-only prediction", fixed = TRUE)
  expect_match(printed, fit$pipeline$id, fixed = TRUE)
  expect_match(printed, "Patients: 5", fixed = TRUE)
  expect_false(grepl("AUC|Brier|concordance|validation performance", printed))
})

test_that("survival final prediction carries declared event-risk thresholds", {
  d <- make_sb_fixture(n = 60L, seed = 20260910L)
  predictors <- c("gene_a", "gene_b")
  thresholds <- c(0.20, 0.50)
  plain_task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month",
    prediction_horizon = 12
  )
  threshold_task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month",
    prediction_horizon = 12,
    decision_thresholds = thresholds
  )
  plain_fit <- sb_fit_final(plain_task, sb_pipeline())
  threshold_fit <- sb_fit_final(threshold_task, sb_pipeline())
  new_data <- d[1:8, c("patient_id", predictors), drop = FALSE]
  plain_prediction <- sb_predict_final(
    plain_fit,
    new_data,
    "patient_id"
  )
  threshold_prediction <- sb_predict_final(
    threshold_fit,
    new_data,
    "patient_id"
  )

  expect_identical(
    threshold_fit$development$endpoint$decision_thresholds,
    thresholds
  )
  expect_identical(threshold_fit$prediction$decision_thresholds, thresholds)
  expect_identical(threshold_prediction$decision_thresholds, thresholds)
  expect_equal(
    threshold_prediction$predictions$risk,
    plain_prediction$predictions$risk,
    tolerance = 0
  )
  expect_equal(
    threshold_prediction$predictions$survival_probability,
    plain_prediction$predictions$survival_probability,
    tolerance = 0
  )
  expect_false("decision_thresholds" %in% names(plain_prediction))
  expect_false(any(c(
    "predicted_class", "decision", "action", "treat", "selected"
  ) %in% names(threshold_prediction$predictions)))
})

test_that("survival final-fit threshold copies must remain consistent", {
  d <- make_sb_fixture(n = 60L, seed = 20260911L)
  predictors <- c("gene_a", "gene_b")
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 12,
    decision_thresholds = c(0.20, 0.50)
  )
  fit <- sb_fit_final(task, sb_pipeline())
  new_data <- d[1:5, c("patient_id", predictors), drop = FALSE]

  missing_copy <- fit
  missing_copy$development$endpoint$decision_thresholds <- NULL
  expect_final_prediction_failure(
    sb_predict_final(missing_copy, new_data, "patient_id"),
    "stormbreaker_final_state_semantics_mismatch",
    5L
  )

  inconsistent_copy <- fit
  inconsistent_copy$development$endpoint$decision_thresholds <- c(0.20, 0.60)
  expect_final_prediction_failure(
    sb_predict_final(inconsistent_copy, new_data, "patient_id"),
    "stormbreaker_final_state_semantics_mismatch",
    5L
  )

  invalid_copy <- fit
  invalid_copy$prediction$decision_thresholds <- c(0.20, 1)
  expect_final_prediction_failure(
    sb_predict_final(invalid_copy, new_data, "patient_id"),
    "stormbreaker_bad_decision_thresholds",
    5L
  )
})
