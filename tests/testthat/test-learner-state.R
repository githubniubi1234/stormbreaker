sb_state_contains_live_object <- function(value) {
  forbidden_class <- inherits(
    value,
    c(
      "coxph", "glmnetfit", "lognet", "coxnet",
      "lgb.Booster", "catboost.Model"
    )
  )
  if (is.function(value) || is.environment(value) ||
      typeof(value) == "externalptr" || forbidden_class) {
    return(TRUE)
  }
  if (!is.list(value)) {
    return(FALSE)
  }
  any(vapply(unclass(value), sb_state_contains_live_object, logical(1)))
}

test_that("ordinary Cox state reconstructs risk and horizon probability", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  x <- d[c("gene_b", "gene_a")]
  adapter <- sb_coxph_learner_adapter()
  model <- adapter$fit(x, d$time, d$event)
  direct_risk <- adapter$predict(model, x)
  direct_probability <- adapter$predict_survival(
    model,
    x,
    d$time,
    d$event,
    x,
    prediction_horizon = 12
  )$probability
  state <- adapter$extract_state(model, list(
    x = x,
    time = d$time,
    event = as.integer(d$event),
    prediction_horizon = 12,
    time_unit = "month"
  ))
  rm(model)
  observed <- adapter$predict_state(state, x)
  reconstructed_risk <- unname(as.vector(
    as.matrix(x[state$features]) %*% state$coefficient
  ))
  reconstructed_probability <- exp(
    -state$baseline_hazard$cumulative_baseline_hazard *
      exp(reconstructed_risk)
  )

  expect_s3_class(state, "stormbreaker_learner_state")
  expect_identical(state$learner, "coxph")
  expect_identical(inherits(state, "coxph"), FALSE)
  expect_false(sb_state_contains_live_object(state))
  expect_equal(observed$risk, direct_risk, tolerance = 1e-12)
  expect_equal(observed$risk, reconstructed_risk, tolerance = 1e-12)
  expect_equal(
    observed$survival_probability,
    direct_probability,
    tolerance = 1e-12
  )
  expect_equal(
    observed$survival_probability,
    reconstructed_probability,
    tolerance = 1e-12
  )
  expect_identical(
    adapter$predict_state(unserialize(serialize(state, NULL)), x),
    observed
  )
})

test_that("glmnet Cox state predicts without retaining a glmnet fit", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_glmnet_cox_learner_adapter(
    alpha = 0.5,
    lambda = 0.05,
    standardize = TRUE
  )
  model <- adapter$fit(x, d$time, d$event)
  direct_risk <- adapter$predict(model, x)
  direct_probability <- adapter$predict_survival(
    model,
    x,
    d$time,
    d$event,
    x,
    prediction_horizon = 12
  )$probability
  state <- adapter$extract_state(model, list(
    x = x,
    time = d$time,
    event = as.integer(d$event),
    prediction_horizon = 12,
    time_unit = "month"
  ))
  rm(model)
  observed <- adapter$predict_state(state, x)
  reconstructed_risk <- unname(as.vector(
    as.matrix(x[state$features]) %*% state$coefficient
  ))
  reconstructed_probability <- exp(
    -state$baseline_hazard$cumulative_baseline_hazard *
      exp(reconstructed_risk)
  )

  expect_identical(state$learner, "glmnet_cox")
  expect_identical(inherits(state, "glmnetfit"), FALSE)
  expect_false(sb_state_contains_live_object(state))
  expect_equal(observed$risk, direct_risk, tolerance = 1e-12)
  expect_equal(observed$risk, reconstructed_risk, tolerance = 1e-12)
  expect_equal(
    observed$survival_probability,
    direct_probability,
    tolerance = 1e-12
  )
  expect_equal(
    observed$survival_probability,
    reconstructed_probability,
    tolerance = 1e-12
  )
  expect_identical(
    adapter$predict_state(unserialize(serialize(state, NULL)), x),
    observed
  )
})

test_that("logistic states reconstruct class-one probabilities", {
  d <- sb_example_binary(n = 80L, seed = 20260819L)
  x <- d[c("gene_b", "gene_a")]
  adapter <- sb_logistic_learner_adapter()
  model <- adapter$fit(x, d$outcome)
  direct <- adapter$predict(model, x)
  state <- adapter$extract_state(model, list())
  rm(model)
  observed <- adapter$predict_state(state, x)
  linear_predictor <- state$intercept + unname(as.vector(
    as.matrix(x[state$features]) %*% state$coefficient
  ))

  expect_identical(state$learner, "logistic")
  expect_identical(inherits(state, "glm"), FALSE)
  expect_false(sb_state_contains_live_object(state))
  expect_equal(observed$probability, direct, tolerance = 1e-12)
  expect_equal(
    observed$probability,
    stats::plogis(linear_predictor),
    tolerance = 1e-12
  )
})

test_that("glmnet logistic state predicts from saved coefficients only", {
  skip_if_not_installed("glmnet")

  d <- sb_example_binary(n = 80L, seed = 20260819L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_glmnet_logistic_learner_adapter(
    alpha = 0.5,
    lambda = 0.05,
    standardize = TRUE
  )
  model <- adapter$fit(x, d$outcome)
  direct <- adapter$predict(model, x)
  state <- adapter$extract_state(model, list())
  rm(model)
  observed <- adapter$predict_state(state, x)

  expect_identical(state$learner, "glmnet_logistic")
  expect_identical(inherits(state, "glmnetfit"), FALSE)
  expect_false(sb_state_contains_live_object(state))
  expect_equal(observed$probability, direct, tolerance = 1e-12)
  expect_identical(
    adapter$predict_state(unserialize(serialize(state, NULL)), x),
    observed
  )
})

test_that("Gaussian Naive Bayes state retains only portable distributions", {
  d <- sb_example_binary(n = 80L, seed = 20260819L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_gaussian_naive_bayes_learner_adapter(1e-9)
  model <- adapter$fit(x, d$outcome)
  direct <- adapter$predict(model, x)
  state <- adapter$extract_state(model, list())
  rm(model)
  observed <- adapter$predict_state(state, x)

  expect_identical(state$learner, "gaussian_naive_bayes")
  expect_false(sb_state_contains_live_object(state))
  expect_equal(observed$probability, direct, tolerance = 1e-14)
  expect_identical(
    adapter$predict_state(unserialize(serialize(state, NULL)), x),
    observed
  )
})

test_that("LightGBM state reloads its saved model string", {
  skip_if_not_installed("lightgbm")

  d <- sb_example_binary(n = 80L, seed = 20260819L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_lightgbm_binary_learner_adapter(
    nrounds = 5L,
    learning_rate = 0.1,
    num_leaves = 7L,
    min_data_in_leaf = 3L,
    seed = 101L
  )
  model <- adapter$fit(x, d$outcome)
  direct <- adapter$predict(model, x)
  state <- adapter$extract_state(model, list())
  rm(model)
  observed <- adapter$predict_state(state, x)

  expect_identical(state$learner, "lightgbm_binary")
  expect_type(state$model_string, "character")
  expect_identical("booster" %in% names(state), FALSE)
  expect_false(sb_state_contains_live_object(state))
  expect_equal(observed$probability, direct, tolerance = 1e-12)
  expect_equal(
    adapter$predict_state(unserialize(serialize(state, NULL)), x)$probability,
    observed$probability,
    tolerance = 1e-12
  )
})

test_that("CatBoost state reloads its saved raw model", {
  skip_if_not_installed("catboost")

  d <- sb_example_binary(n = 80L, seed = 20260819L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_catboost_binary_learner_adapter(
    iterations = 5L,
    learning_rate = 0.1,
    depth = 3L,
    l2_leaf_reg = 3,
    seed = 101L
  )
  model <- adapter$fit(x, d$outcome)
  direct <- adapter$predict(model, x)
  state <- adapter$extract_state(model, list())
  rm(model)
  observed <- adapter$predict_state(state, x)

  expect_identical(state$learner, "catboost_binary")
  expect_type(state$model_raw, "raw")
  expect_identical("booster" %in% names(state), FALSE)
  expect_false(sb_state_contains_live_object(state))
  expect_equal(observed$probability, direct, tolerance = 1e-12)
  expect_equal(
    adapter$predict_state(unserialize(serialize(state, NULL)), x)$probability,
    observed$probability,
    tolerance = 1e-12
  )
})

test_that("invalid portable states fail visibly without fallback", {
  d <- make_sb_fixture(n = 40L, seed = 20260819L)
  x <- d[c("gene_a", "gene_b")]
  adapter <- sb_coxph_learner_adapter()
  model <- adapter$fit(x, d$time, d$event)
  state <- adapter$extract_state(model, list())
  state$coefficient[1L] <- NA_real_
  coefficient_error <- capture_sb_error(adapter$predict_state(state, x))
  expect_s3_class(coefficient_error, "stormbreaker_bad_learner_state")

  if (sb_lightgbm_available()) {
    binary <- sb_example_binary(n = 60L, seed = 20260819L)
    binary_x <- binary[c("gene_a", "gene_b")]
    lightgbm <- sb_lightgbm_binary_learner_adapter(
      nrounds = 3L,
      num_leaves = 5L,
      min_data_in_leaf = 2L
    )
    lightgbm_model <- lightgbm$fit(binary_x, binary$outcome)
    lightgbm_state <- lightgbm$extract_state(lightgbm_model, list())
    lightgbm_state$model_string <- "not a LightGBM model"
    lightgbm_error <- capture_sb_error(
      lightgbm$predict_state(lightgbm_state, binary_x)
    )
    expect_s3_class(
      lightgbm_error,
      "stormbreaker_invalid_lightgbm_state"
    )
  }

  if (sb_catboost_available()) {
    binary <- sb_example_binary(n = 60L, seed = 20260819L)
    binary_x <- binary[c("gene_a", "gene_b")]
    catboost <- sb_catboost_binary_learner_adapter(
      iterations = 3L,
      depth = 2L
    )
    catboost_model <- catboost$fit(binary_x, binary$outcome)
    catboost_state <- catboost$extract_state(catboost_model, list())
    catboost_state$model_raw <- charToRaw("not a serialized model")
    catboost_error <- capture_sb_error(
      catboost$predict_state(catboost_state, binary_x)
    )
    expect_s3_class(catboost_error, "stormbreaker_invalid_catboost_state")
  }
})

test_that("portable tree states still require their declared dependencies", {
  skip_if_not_installed("lightgbm")

  d <- sb_example_binary(n = 60L, seed = 20260819L)
  x <- d[c("gene_a", "gene_b")]
  adapter <- sb_lightgbm_binary_learner_adapter(
    nrounds = 3L,
    num_leaves = 5L,
    min_data_in_leaf = 2L
  )
  model <- adapter$fit(x, d$outcome)
  state <- adapter$extract_state(model, list())
  testthat::local_mocked_bindings(
    sb_lightgbm_available = function() FALSE
  )

  error <- capture_sb_error(adapter$predict_state(state, x))
  expect_s3_class(error, "stormbreaker_missing_lightgbm")
})

test_that("portable CatBoost state still requires CatBoost", {
  skip_if_not_installed("catboost")

  d <- sb_example_binary(n = 60L, seed = 20260819L)
  x <- d[c("gene_a", "gene_b")]
  adapter <- sb_catboost_binary_learner_adapter(
    iterations = 3L,
    depth = 2L
  )
  model <- adapter$fit(x, d$outcome)
  state <- adapter$extract_state(model, list())
  testthat::local_mocked_bindings(
    sb_catboost_available = function() FALSE
  )

  error <- capture_sb_error(adapter$predict_state(state, x))
  expect_s3_class(error, "stormbreaker_missing_catboost")
})
