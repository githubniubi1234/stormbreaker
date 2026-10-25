test_that("fixed SurvReg learners declare four distribution assumptions", {
  expected <- list(
    survreg_weibull = list(
      distribution = "weibull",
      scale = "estimated"
    ),
    survreg_exponential = list(
      distribution = "exponential",
      scale = "fixed_1"
    ),
    survreg_lognormal = list(
      distribution = "lognormal",
      scale = "estimated"
    ),
    survreg_loglogistic = list(
      distribution = "loglogistic",
      scale = "estimated"
    )
  )

  expect_identical(sb_survreg_learner_ids(), names(expected))
  expect_true(all(names(expected) %in%
    sb_additional_survival_learner_ids()))
  for (learner in names(expected)) {
    expect_identical(
      sb_additional_survival_params(learner),
      expected[[learner]]
    )
    expect_identical(
      sb_additional_survival_dependency(learner),
      list(package = "survival", packages = "survival", available = TRUE)
    )
  }
})

test_that("fixed SurvReg risk and survival equal native distribution results", {
  data <- make_sb_fixture(n = 72L, seed = 8521L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:7, , drop = FALSE]
  horizon <- 7

  for (learner in sb_survreg_learner_ids()) {
    params <- sb_additional_survival_params(learner)
    adapter <- sb_additional_survival_learner_adapter(learner, params)
    model <- adapter$fit(x, data$time, data$event)
    observed_risk <- adapter$predict(model, assessment)
    observed_survival <- adapter$predict_survival(
      model,
      x,
      data$time,
      data$event,
      assessment,
      prediction_horizon = horizon
    )
    native_location <- stats::predict(
      model$engine,
      newdata = assessment,
      type = "lp"
    )
    native_survival <- 1 - survival::psurvreg(
      horizon,
      mean = native_location,
      scale = model$engine$scale,
      distribution = params$distribution
    )
    state <- adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = horizon,
      time_unit = "month"
    ))
    replay <- adapter$predict_state(
      unserialize(serialize(state, NULL)),
      assessment
    )

    expect_equal(observed_risk, -unname(native_location), tolerance = 1e-14)
    expect_equal(
      observed_survival$probability,
      unname(native_survival),
      tolerance = 1e-14
    )
    expect_equal(replay$risk, observed_risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      observed_survival$probability,
      tolerance = 0
    )
    expect_false("engine_raw" %in% names(state))
    expect_identical(
      state$risk_reference,
      paste0(params$distribution, "_aft_native")
    )
    expect_identical(
      state$survival_probability$source,
      paste0("survreg_", params$distribution, "_distribution_direct")
    )
  }
})

test_that("fixed SurvReg preserves quoted predictor names and native predictions", {
  data <- make_sb_fixture(n = 72L, seed = 8521L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("NKX3-1", "2 gene", "if")
  assessment <- x[1:7, , drop = FALSE]

  for (learner in sb_survreg_learner_ids()) {
    params <- sb_additional_survival_params(learner)
    formula_data <- sb_additional_survival_formula_data(x, data$time, data$event)
    native <- survival::survreg(
      formula_data$formula, data = formula_data$data,
      dist = params$distribution, model = FALSE, x = FALSE, y = FALSE
    )
    adapter <- sb_additional_survival_learner_adapter(learner, params)
    model <- adapter$fit(x, data$time, data$event)
    location <- unname(stats::predict(native, newdata = assessment, type = "lp"))
    probability <- unname(1 - survival::psurvreg(
      7, mean = location, scale = native$scale,
      distribution = params$distribution
    ))
    state <- adapter$extract_state(model, list(
      x = x, time = data$time, event = data$event,
      prediction_horizon = 7, time_unit = "month"
    ))
    replay <- adapter$predict_state(unserialize(serialize(state, NULL)), assessment)

    expect_identical(model$features, names(x))
    expect_identical(state$features, names(x))
    expect_identical(stats::coef(model$engine), stats::coef(native))
    expect_identical(model$engine$scale, native$scale)
    expect_equal(adapter$predict(model, assessment), -location, tolerance = 1e-14)
    expect_equal(replay$risk, -location, tolerance = 1e-14)
    expect_equal(replay$survival_probability, probability, tolerance = 1e-14)
  }
})

test_that("fixed SurvReg still rejects misaligned and nonfinite coefficients", {
  data <- make_sb_fixture(n = 72L, seed = 8521L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  formula_data <- sb_additional_survival_formula_data(x, data$time, data$event)
  native <- survival::survreg(
    formula_data$formula, data = formula_data$data, dist = "weibull"
  )
  adapter <- sb_additional_survival_learner_adapter("survreg_weibull")
  bad <- native
  names(bad$coefficients)[2:3] <- rev(names(bad$coefficients)[2:3])
  local_mocked_bindings(sb_fit_survreg = function(...) bad)
  expect_error(adapter$fit(x, data$time, data$event),
    class = "stormbreaker_incomplete_survreg_fit")
  bad <- native
  bad$coefficients[[2L]] <- NA_real_
  expect_error(adapter$fit(x, data$time, data$event),
    class = "stormbreaker_incomplete_survreg_fit")
})
