test_that("XGBoost AFT boosters have one explicit native specification", {
  learners <- c(
    "xgboost_aft",
    "xgboost_aft_dart",
    "xgboost_aft_linear"
  )
  expect_identical(sb_xgboost_aft_learner_ids(), learners)

  for (learner in learners) {
    params <- sb_xgboost_params(learner)
    expect_identical(params$nrounds, 100L)
    expect_identical(params$nthread, 1L)
    expect_identical(params$seed, 1L)
    expect_identical(params$objective, "survival:aft")
    expect_identical(params$aft_loss_distribution, "normal")
    expect_identical(params$aft_loss_distribution_scale, 1)
  }
  expect_identical(
    sb_xgboost_booster("xgboost_aft"),
    "gbtree"
  )
  expect_identical(
    sb_xgboost_booster("xgboost_aft_dart"),
    "dart"
  )
  expect_identical(
    sb_xgboost_booster("xgboost_aft_linear"),
    "gblinear"
  )

  registry <- sb_pipeline_registry()
  selectors <- expected_fixed_learner_selectors("survival")
  for (learner in learners) {
    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_identical(unique(rows$task_type), "survival")
    expect_setequal(rows$selector, selectors)
    expect_true(all(grepl("xgboost", rows$optional_dependency, fixed = TRUE)))
  }
})

test_that("XGBoost AFT uses right-censoring bounds and declared controls", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  data <- make_sb_fixture(n = 56L, seed = 8831L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  original_train <- xgboost::xgb.train
  observed <- NULL

  testthat::with_mocked_bindings(
    xgb.train = function(params = list(), data, nrounds, ...) {
      observed <<- list(
        params = params,
        lower = xgboost::getinfo(data, "label_lower_bound"),
        upper = xgboost::getinfo(data, "label_upper_bound")
      )
      original_train(
        params = params,
        data = data,
        nrounds = nrounds,
        ...
      )
    },
    {
      params <- sb_xgboost_params("xgboost_aft")
      params$nrounds <- 4L
      sb_xgboost_learner_adapter(
        "xgboost_aft",
        params
      )$fit(x, data$time, data$event)
    },
    .package = "xgboost"
  )

  expect_equal(observed$lower, data$time, tolerance = 1e-5)
  expect_equal(
    observed$upper[data$event == 1L],
    data$time[data$event == 1L],
    tolerance = 1e-5
  )
  expect_true(all(is.infinite(observed$upper[data$event == 0L])))
  expect_identical(observed$params$objective, "survival:aft")
  expect_identical(observed$params$aft_loss_distribution, "normal")
  expect_identical(observed$params$aft_loss_distribution_scale, 1)
})

test_that("XGBoost AFT risk, native survival, and raw state replay exactly", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  data <- make_sb_fixture(n = 72L, seed = 8832L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("gene one", "2gene", "gene-three")
  new_x <- x[c(4L, 29L, 61L), , drop = FALSE]
  withr::local_seed(8833L)

  for (learner in sb_xgboost_aft_learner_ids()) {
    params <- sb_xgboost_params(learner)
    params$nrounds <- 12L
    adapter <- sb_xgboost_learner_adapter(learner, params)
    rng_before <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, rng_before)
    repeated <- adapter$fit(x, data$time, data$event)
    expect_identical(
      xgboost::xgb.save.raw(model$booster),
      xgboost::xgb.save.raw(repeated$booster)
    )

    native_input <- xgboost::xgb.DMatrix(as.matrix(new_x), nthread = 1L)
    native_margin <- as.numeric(stats::predict(
      model$booster,
      native_input,
      outputmargin = TRUE
    ))
    native_time <- as.numeric(stats::predict(
      model$booster,
      native_input
    ))
    risk <- adapter$predict(model, new_x)
    probability <- adapter$predict_survival(
      model,
      x,
      data$time,
      data$event,
      new_x,
      prediction_horizon = 7
    )
    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = 7,
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, new_x)

    expect_equal(native_time, exp(native_margin), tolerance = 1e-6)
    expect_equal(risk, -native_margin, tolerance = 0)
    expect_equal(
      probability$probability,
      stats::pnorm(native_margin - log(7)),
      tolerance = 0
    )
    expect_identical(probability$source, "xgboost_aft_normal_native")
    expect_identical(probability$ties, "not_applicable")
    expect_equal(replay$risk, risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      probability$probability,
      tolerance = 0
    )
    expect_identical(state$risk_reference, "xgboost_aft_normal_native")
    expect_identical(
      state$risk_formula,
      "negative_xgboost_aft_log_time_margin"
    )
    expect_type(state$model_raw, "raw")
    expect_false("booster" %in% names(state))
    expect_false(any(c("time", "event", "outcome") %in% names(state)))
  }
})

test_that("XGBoost AFT supports one feature and monotone horizon survival", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  data <- make_sb_fixture(n = 64L, seed = 8834L)
  x <- data["gene_a"]
  names(x) <- "1 gene-A"
  params <- sb_xgboost_params("xgboost_aft")
  params$nrounds <- 10L
  adapter <- sb_xgboost_learner_adapter(
    "xgboost_aft",
    params
  )
  model <- adapter$fit(x, data$time, data$event)
  horizons <- c(4, 7, 10)
  probability <- vapply(
    horizons,
    function(horizon) {
      adapter$predict_survival(
        model,
        x,
        data$time,
        data$event,
        x[1:6, , drop = FALSE],
        prediction_horizon = horizon
      )$probability
    },
    numeric(6L)
  )

  expect_true(all(probability >= 0 & probability <= 1))
  expect_true(all(probability[, 1L] >= probability[, 2L]))
  expect_true(all(probability[, 2L] >= probability[, 3L]))
})

test_that("XGBoost AFT boosters execute through CV and final prediction", {
  skip_if_not_installed("xgboost", minimum_version = "1.7.0")

  data <- make_sb_fixture(n = 72L, seed = 8835L)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  for (learner in sb_xgboost_aft_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 2L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 8836L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      data[1:5, c("patient_id", "gene_c", "gene_a", "gene_b")],
      "patient_id"
    )

    expect_identical(pipeline$risk_type, "aft_log_time")
    expect_identical(
      pipeline$risk_reference,
      "xgboost_aft_normal_native"
    )
    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(
      comparison$complete_survival_probability_oof,
      TRUE
    )
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_true(all(
      prediction$predictions$survival_probability >= 0 &
        prediction$predictions$survival_probability <= 1
    ))
  }
})
