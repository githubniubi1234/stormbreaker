test_that("fixed XGBoost specifications are explicit", {
  binary <- sb_xgboost_params("xgboost_binary")
  survival <- sb_xgboost_params("xgboost_cox")

  expect_identical(binary$nrounds, 100L)
  expect_identical(binary$eta, 0.1)
  expect_identical(binary$max_depth, 3L)
  expect_identical(binary$min_child_weight, 1)
  expect_identical(binary$subsample, 1)
  expect_identical(binary$colsample_bytree, 1)
  expect_identical(binary$nthread, 1L)
  expect_identical(binary$seed, 1L)
  expect_identical(binary$objective, "binary:logistic")
  expect_identical(survival$objective, "survival:cox")
})

test_that("fixed XGBoost uses one thread for matrix construction", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  data <- make_sb_binary_fixture(n = 36L, seed = 7900L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  params <- sb_xgboost_params("xgboost_binary")
  params$nrounds <- 5L
  adapter <- sb_xgboost_learner_adapter(
    "xgboost_binary",
    params
  )
  original_dmatrix <- xgboost::xgb.DMatrix
  observed_threads <- integer()
  observed_matrix_input <- logical()

  testthat::with_mocked_bindings(
    xgb.DMatrix = function(data, nthread = NULL, ...) {
      observed_threads <<- c(
        observed_threads,
        if (is.null(nthread)) NA_integer_ else as.integer(nthread)
      )
      observed_matrix_input <<- c(observed_matrix_input, is.matrix(data))
      original_dmatrix(
        data = data,
        nthread = nthread,
        ...
      )
    },
    {
      model <- adapter$fit(x, data$outcome)
      adapter$predict(model, x[1:3, , drop = FALSE])
    },
    .package = "xgboost"
  )

  expect_identical(
    observed_threads[observed_matrix_input],
    c(1L, 1L)
  )
})

test_that("fixed XGBoost adapters preserve native predictions and state", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  binary <- make_sb_binary_fixture(n = 54L, seed = 7901L)
  binary_x <- binary[c("gene_a", "gene_b", "gene_c")]
  binary_params <- sb_xgboost_params("xgboost_binary")
  binary_params$nrounds <- 12L
  binary_adapter <- sb_xgboost_learner_adapter(
    "xgboost_binary",
    binary_params
  )
  binary_model <- binary_adapter$fit(binary_x, binary$outcome)
  binary_new <- binary_x[c(4L, 21L, 43L), , drop = FALSE]
  binary_probability <- binary_adapter$predict(binary_model, binary_new)
  binary_state <- binary_adapter$extract_state(binary_model)
  binary_replay <- binary_adapter$predict_state(
    unserialize(serialize(binary_state, NULL)),
    binary_new
  )

  expect_equal(
    binary_probability,
    as.numeric(stats::predict(binary_model$booster, as.matrix(binary_new))),
    tolerance = 0
  )
  expect_equal(binary_replay$probability, binary_probability, tolerance = 0)
  expect_type(binary_state$model_raw, "raw")
  expect_false("booster" %in% names(binary_state))

  survival <- make_sb_fixture(n = 60L, seed = 7902L)
  survival_x <- survival[c("gene_a", "gene_b", "gene_c")]
  survival_params <- sb_xgboost_params("xgboost_cox")
  survival_params$nrounds <- 12L
  survival_adapter <- sb_xgboost_learner_adapter(
    "xgboost_cox",
    survival_params
  )
  survival_model <- survival_adapter$fit(
    survival_x,
    survival$time,
    survival$event
  )
  survival_new <- survival_x[c(3L, 28L, 51L), , drop = FALSE]
  risk <- survival_adapter$predict(survival_model, survival_new)
  probability <- survival_adapter$predict_survival(
    survival_model,
    survival_x,
    survival$time,
    survival$event,
    survival_new,
    prediction_horizon = 7
  )
  survival_state <- survival_adapter$extract_state(survival_model, list(
    x = survival_x,
    time = survival$time,
    event = survival$event,
    prediction_horizon = 7,
    time_unit = "month"
  ))
  survival_replay <- survival_adapter$predict_state(
    unserialize(serialize(survival_state, NULL)),
    survival_new
  )

  native_margin <- as.numeric(stats::predict(
    survival_model$booster,
    as.matrix(survival_new),
    outputmargin = TRUE
  ))
  expect_equal(risk, native_margin, tolerance = 0)
  expect_equal(survival_replay$risk, risk, tolerance = 0)
  expect_equal(
    survival_replay$survival_probability,
    probability$probability,
    tolerance = 0
  )
  expect_identical(probability$source, "xgboost_cox_breslow_calibration")
  expect_identical(probability$ties, "breslow")
  expect_true(all(probability$probability >= 0 & probability$probability <= 1))
})

test_that("XGBoost routes run binary and survival workflows", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  binary <- make_sb_binary_fixture(n = 60L, seed = 7911L)
  binary_task <- sb_binary_task(
    binary,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  binary_pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "xgboost_binary",
    selector_n_features = 2L
  )
  binary_result <- sb_run_cv(
    binary_task,
    binary_pipeline,
    folds = 3L,
    seed = 7912L
  )
  binary_fit <- sb_fit_final(binary_task, binary_pipeline)
  binary_prediction <- sb_predict_final(
    binary_fit,
    binary[c(5L, 24L, 47L), c("patient_id", "gene_c", "gene_a", "gene_b")],
    "patient_id"
  )

  expect_identical(binary_result$status, "success")
  expect_true(binary_result$complete_oof)
  expect_identical(binary_fit$status, "success")
  expect_identical(binary_prediction$status, "success")
  expect_true(all(binary_prediction$predictions$probability >= 0 &
    binary_prediction$predictions$probability <= 1))

  survival <- make_sb_fixture(n = 60L, seed = 7913L)
  survival_task <- sb_surv_task(
    survival,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  survival_pipeline <- sb_pipeline(learner = "xgboost_cox")
  survival_result <- sb_run_cv(
    survival_task,
    survival_pipeline,
    folds = 3L,
    seed = 7914L
  )
  survival_fit <- sb_fit_final(survival_task, survival_pipeline)
  survival_prediction <- sb_predict_final(
    survival_fit,
    survival[c(7L, 30L, 55L), c("patient_id", "gene_b", "gene_c", "gene_a")],
    "patient_id"
  )

  expect_identical(survival_pipeline$risk_type, "cox_raw_margin")
  expect_identical(survival_pipeline$risk_reference, "xgboost_native")
  expect_identical(survival_result$status, "success")
  expect_true(survival_result$complete_oof)
  expect_true(survival_result$complete_survival_probability_oof)
  expect_identical(survival_fit$status, "success")
  expect_identical(survival_prediction$status, "success")
})
