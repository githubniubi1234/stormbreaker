test_that("XGBoost DART and linear boosters are distinct fixed learners", {
  binary_ids <- c(
    "xgboost_binary_dart",
    "xgboost_binary_linear"
  )
  cox_ids <- c(
    "xgboost_cox_dart",
    "xgboost_cox_linear"
  )

  expect_setequal(sb_xgboost_binary_learner_ids(), c(
    "xgboost_binary", binary_ids
  ))
  expect_setequal(sb_xgboost_cox_learner_ids(), c(
    "xgboost_cox", cox_ids
  ))

  dart <- sb_xgboost_params("xgboost_binary_dart")
  expect_identical(dart, list(
    nrounds = 100L,
    booster = "dart",
    eta = 0.1,
    max_depth = 3L,
    min_child_weight = 1,
    subsample = 1,
    colsample_bytree = 1,
    sample_type = "uniform",
    normalize_type = "tree",
    rate_drop = 0.1,
    skip_drop = 0,
    one_drop = FALSE,
    nthread = 1L,
    seed = 1L,
    objective = "binary:logistic"
  ))
  expect_identical(
    sb_xgboost_params("xgboost_cox_dart")$objective,
    "survival:cox"
  )

  linear <- sb_xgboost_params("xgboost_binary_linear")
  expect_identical(linear, list(
    nrounds = 100L,
    booster = "gblinear",
    eta = 0.1,
    lambda = 1,
    alpha = 0,
    updater = "coord_descent",
    feature_selector = "cyclic",
    nthread = 1L,
    seed = 1L,
    objective = "binary:logistic"
  ))
  expect_identical(
    sb_xgboost_params("xgboost_cox_linear")$objective,
    "survival:cox"
  )

  registry <- sb_pipeline_registry()
  expected_selectors <- list(
    binary = expected_fixed_learner_selectors("binary"),
    survival = expected_fixed_learner_selectors("survival")
  )
  for (learner in c(binary_ids, cox_ids)) {
    rows <- registry[registry$learner == learner, , drop = FALSE]
    task <- if (learner %in% binary_ids) "binary" else "survival"
    expect_identical(unique(rows$task_type), task)
    expect_setequal(rows$selector, expected_selectors[[task]])
    expect_true(all(grepl("xgboost", rows$optional_dependency, fixed = TRUE)))
  }
})

test_that("XGBoost booster variants preserve native predictions and state", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  binary <- make_sb_binary_fixture(n = 64L, seed = 8711L)
  binary_x <- binary[c("gene_a", "gene_b", "gene_c")]
  names(binary_x) <- c("gene A", "gene-B", "gene_c")
  binary_new <- binary_x[c(4L, 27L, 51L), , drop = FALSE]
  withr::local_seed(8710L)

  for (learner in c(
    "xgboost_binary_dart",
    "xgboost_binary_linear"
  )) {
    params <- sb_xgboost_params(learner)
    params$nrounds <- 12L
    adapter <- sb_xgboost_learner_adapter(learner, params)
    rng_before <- .Random.seed
    model <- adapter$fit(binary_x, binary$outcome)
    expect_identical(.Random.seed, rng_before)
    repeated <- adapter$fit(binary_x, binary$outcome)
    expect_identical(
      xgboost::xgb.save.raw(model$booster),
      xgboost::xgb.save.raw(repeated$booster)
    )

    probability <- adapter$predict(model, binary_new)
    native <- as.numeric(stats::predict(
      model$booster,
      as.matrix(binary_new)
    ))
    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    replay <- adapter$predict_state(state, binary_new)$probability
    config <- xgboost::xgb.config(model$booster)

    expect_equal(probability, native, tolerance = 0)
    expect_equal(replay, probability, tolerance = 0)
    expect_true(all(probability >= 0 & probability <= 1))
    if (is.list(config)) {
      expect_identical(config$learner$learner_train_param$booster, params$booster)
    } else {
      expect_true(grepl(
        paste0('"booster":"', params$booster, '"'),
        config,
        fixed = TRUE
      ))
    }
    expect_type(state$model_raw, "raw")
    expect_false("booster" %in% names(state))
  }

  survival <- make_sb_fixture(n = 70L, seed = 8712L)
  survival_x <- survival[c("gene_a", "gene_b", "gene_c")]
  survival_new <- survival_x[c(6L, 31L, 58L), , drop = FALSE]
  for (learner in c(
    "xgboost_cox_dart",
    "xgboost_cox_linear"
  )) {
    params <- sb_xgboost_params(learner)
    params$nrounds <- 12L
    adapter <- sb_xgboost_learner_adapter(learner, params)
    model <- adapter$fit(survival_x, survival$time, survival$event)
    risk <- adapter$predict(model, survival_new)
    native <- as.numeric(stats::predict(
      model$booster,
      as.matrix(survival_new),
      outputmargin = TRUE
    ))
    probability <- adapter$predict_survival(
      model,
      survival_x,
      survival$time,
      survival$event,
      survival_new,
      prediction_horizon = 7
    )
    state <- adapter$extract_state(model, list(
      x = survival_x,
      time = survival$time,
      event = survival$event,
      prediction_horizon = 7,
      time_unit = "month"
    ))
    replay <- adapter$predict_state(
      unserialize(serialize(state, NULL)),
      survival_new
    )

    expect_equal(risk, native, tolerance = 0)
    expect_equal(replay$risk, risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      probability$probability,
      tolerance = 0
    )
    expect_true(all(
      probability$probability >= 0 & probability$probability <= 1
    ))
    expect_identical(probability$source, "xgboost_cox_breslow_calibration")
  }
})

test_that("XGBoost forwards the declared DART and linear controls", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  data <- make_sb_binary_fixture(n = 40L, seed = 8713L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  original_train <- xgboost::xgb.train
  observed <- list()

  testthat::with_mocked_bindings(
    xgb.train = function(params = list(), data, nrounds, ...) {
      observed[[length(observed) + 1L]] <<- params
      original_train(
        params = params,
        data = data,
        nrounds = nrounds,
        ...
      )
    },
    {
      for (learner in c(
        "xgboost_binary_dart",
        "xgboost_binary_linear"
      )) {
        params <- sb_xgboost_params(learner)
        params$nrounds <- 3L
        sb_xgboost_learner_adapter(learner, params)$fit(
          x,
          data$outcome
        )
      }
    },
    .package = "xgboost"
  )

  expect_identical(observed[[1L]], list(
    objective = "binary:logistic",
    booster = "dart",
    eta = 0.1,
    max_depth = 3L,
    min_child_weight = 1,
    subsample = 1,
    colsample_bytree = 1,
    nthread = 1L,
    sample_type = "uniform",
    normalize_type = "tree",
    rate_drop = 0.1,
    skip_drop = 0,
    one_drop = FALSE
  ))
  expect_identical(observed[[2L]], list(
    objective = "binary:logistic",
    booster = "gblinear",
    eta = 0.1,
    lambda = 1,
    alpha = 0,
    updater = "coord_descent",
    feature_selector = "cyclic",
    nthread = 1L
  ))
})

test_that("XGBoost booster variants execute through CV and final prediction", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  binary <- make_sb_binary_fixture(n = 60L, seed = 8721L)
  binary_task <- sb_binary_task(
    binary,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  for (learner in c(
    "xgboost_binary_dart",
    "xgboost_binary_linear"
  )) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 2L
    )
    comparison <- sb_run_cv(binary_task, pipeline, folds = 3L, seed = 8722L)
    fit <- sb_fit_final(binary_task, pipeline)
    prediction <- sb_predict_final(
      fit,
      binary[1:5, c("patient_id", "gene_c", "gene_a", "gene_b")],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$probability, 5L)
  }

  survival <- make_sb_fixture(n = 66L, seed = 8723L)
  survival_task <- sb_surv_task(
    survival,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  for (learner in c(
    "xgboost_cox_dart",
    "xgboost_cox_linear"
  )) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 2L
    )
    comparison <- sb_run_cv(
      survival_task,
      pipeline,
      folds = 3L,
      seed = 8724L
    )
    fit <- sb_fit_final(survival_task, pipeline)
    prediction <- sb_predict_final(
      fit,
      survival[1:5, c("patient_id", "gene_b", "gene_c", "gene_a")],
      "patient_id"
    )

    expect_identical(pipeline$risk_type, "cox_raw_margin")
    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(
      comparison$complete_survival_probability_oof,
      TRUE
    )
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
  }
})
