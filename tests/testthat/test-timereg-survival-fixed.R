make_timereg_survival_fixture <- function(n = 90L, seed = 20260831L) {
  withr::with_seed(seed, {
    gene_a <- stats::runif(n)
    gene_b <- stats::runif(n)
    gene_c <- stats::runif(n)
    event_time <- stats::rexp(
      n,
      rate = 0.08 + 0.08 * gene_a + 0.05 * gene_b + 0.03 * gene_c
    )
    censor_time <- stats::rexp(n, rate = 0.025)
    data.frame(
      patient_id = sprintf("TR%03d", seq_len(n)),
      time = pmax(pmin(event_time, censor_time), 0.01),
      event = as.integer(event_time <= censor_time),
      gene_a,
      gene_b,
      gene_c,
      check.names = FALSE
    )
  })
}

test_that("timereg contributes two distinct fixed survival learners", {
  expect_identical(
    sb_timereg_survival_learner_ids(),
    c("aalen_additive_survival", "prop_odds_survival")
  )

  registry <- sb_pipeline_registry()
  rows <- registry[
    registry$learner %in% sb_timereg_survival_learner_ids(), ,
    drop = FALSE
  ]
  expect_identical(nrow(rows), 2L * length(expected_fixed_learner_selectors("survival")))
  expect_setequal(
    unique(rows$selector),
    expected_fixed_learner_selectors("survival")
  )
  expect_true(all(grepl("timereg", rows$optional_dependency, fixed = TRUE)))

  aalen <- sb_pipeline(learner = "aalen_additive_survival")
  prop_odds <- sb_pipeline(learner = "prop_odds_survival")
  expect_identical(
    aalen$risk_formula,
    "cumulative_hazard_at_last_analysis_event"
  )
  expect_identical(aalen$learner_params$risk_time_rule, "last_analysis_event")
  expect_identical(prop_odds$risk_formula, "x_beta")
  expect_identical(prop_odds$risk_reference, "zero")
})

test_that("timereg portable predictions reproduce native model formulas", {
  skip_if_not_installed("timereg", minimum_version = "2.0.7")
  d <- make_timereg_survival_fixture(n = 120L, seed = 9431L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:9, , drop = FALSE]

  set.seed(9432L)
  caller_seed <- .Random.seed
  aalen_adapter <- sb_timereg_survival_learner_adapter(
    "aalen_additive_survival"
  )
  aalen <- aalen_adapter$fit(x, d$time, d$event)
  expect_identical(.Random.seed, caller_seed)

  aalen_risk <- aalen_adapter$predict(aalen, assessment)
  expected_aalen_risk <- unname(
    aalen$cumulative_intercept[[length(aalen$time)]] +
      as.vector(
        as.matrix(assessment) %*%
          aalen$cumulative_coefficient[nrow(aalen$cumulative_coefficient), ]
      )
  )
  expect_equal(aalen_risk, expected_aalen_risk, tolerance = 0)
  native_aalen <- withr::with_seed(
    1L,
    timereg::aalen(
      survival::Surv(time, event) ~ gene_a + gene_b + gene_c,
      data = d,
      robust = 0L,
      residuals = 0L,
      n.sim = 0L,
      weighted.test = 0L,
      covariance = 0L,
      resample.iid = 0L,
      silent = 1L
    )
  )
  native_aalen_survival <- stats::predict(
    native_aalen,
    X = cbind(1, as.matrix(assessment)),
    times = aalen$risk_reference_time,
    se = FALSE,
    uniform = FALSE
  )$S0
  expect_equal(
    aalen_risk,
    -log(unname(as.numeric(native_aalen_survival))),
    tolerance = 1e-12
  )
  native_aalen_horizon <- stats::predict(
    native_aalen,
    X = cbind(1, as.matrix(assessment)),
    times = 5,
    se = FALSE,
    uniform = FALSE
  )$S0
  portable_aalen_horizon <- aalen_adapter$predict_survival(
    aalen,
    x,
    d$time,
    d$event,
    assessment,
    5
  )$probability
  expect_equal(
    portable_aalen_horizon,
    unname(as.numeric(native_aalen_horizon)),
    tolerance = 1e-12
  )

  prop_adapter <- sb_timereg_survival_learner_adapter(
    "prop_odds_survival"
  )
  prop_odds <- prop_adapter$fit(x, d$time, d$event)
  prop_risk <- prop_adapter$predict(prop_odds, assessment)
  expect_equal(
    prop_risk,
    unname(as.vector(as.matrix(assessment) %*% prop_odds$coefficient)),
    tolerance = 0
  )

  horizon <- 5
  native_prop <- withr::with_seed(
    1L,
    timereg::prop.odds(
      timereg::Event(time, event) ~ gene_a + gene_b + gene_c,
      data = d,
      Nit = 20L,
      detail = 0L,
      n.sim = 0L,
      weighted.test = 0L,
      profile = 1L,
      sym = 0L,
      baselinevar = 1L
    )
  )
  native_survival <- stats::predict(
    native_prop,
    X = matrix(1, nrow(assessment), 1L),
    Z = as.matrix(assessment),
    times = horizon,
    se = FALSE,
    uniform = FALSE
  )$S0
  portable_survival <- prop_adapter$predict_survival(
    prop_odds,
    x,
    d$time,
    d$event,
    assessment,
    horizon
  )$probability
  expect_equal(
    portable_survival,
    unname(as.numeric(native_survival)),
    tolerance = 1e-12
  )
})

test_that("Aalen native survival is never silently clipped", {
  model <- structure(
    list(
      learner = "aalen_additive_survival",
      features = "gene_a",
      params = sb_timereg_survival_params(
        "aalen_additive_survival"
      ),
      time = c(0, 1),
      cumulative_intercept = c(0, -0.2),
      cumulative_coefficient = matrix(c(0, 0), ncol = 1L),
      risk_reference_time = 1
    ),
    class = "stormbreaker_timereg_survival_fixed"
  )

  expect_error(
    sb_timereg_direct_survival_probability(
      model,
      data.frame(gene_a = 0),
      horizon = 1
    ),
    class = "stormbreaker_invalid_timereg_survival_prediction"
  )
})

test_that("Aalen horizon state validation is invariant to patient row order", {
  model <- structure(
    list(
      learner = "aalen_additive_survival",
      features = "gene_a",
      params = sb_timereg_survival_params(
        "aalen_additive_survival"
      ),
      time = c(0, 1),
      cumulative_intercept = c(0, 0.2),
      cumulative_coefficient = matrix(c(0, 1), ncol = 1L),
      risk_reference_time = 1,
      baseline_cumulative_odds = NULL,
      coefficient = NULL
    ),
    class = "stormbreaker_timereg_survival_fixed"
  )
  context <- list(
    x = data.frame(gene_a = c(0, -1)),
    time = c(1, 2),
    event = c(1L, 0L),
    prediction_horizon = 1,
    time_unit = "month"
  )

  expect_error(
    sb_extract_timereg_survival_state(model, context),
    class = "stormbreaker_invalid_timereg_survival_prediction"
  )
  reversed <- c(2L, 1L)
  context$x <- context$x[reversed, , drop = FALSE]
  context$time <- context$time[reversed]
  context$event <- context$event[reversed]
  expect_error(
    sb_extract_timereg_survival_state(model, context),
    class = "stormbreaker_invalid_timereg_survival_prediction"
  )
})

test_that("timereg learners execute CV and portable final prediction", {
  skip_if_not_installed("timereg", minimum_version = "2.0.7")
  d <- make_timereg_survival_fixture(n = 90L, seed = 9433L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  risk_task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month"
  )

  for (learner in sb_timereg_survival_learner_ids()) {
    pipeline <- sb_pipeline(learner = learner)
    cv <- sb_run_cv(risk_task, pipeline, folds = 3L, seed = 9434L)
    fit <- sb_fit_final(risk_task, pipeline)
    state <- unserialize(serialize(fit$learner_state, NULL))
    fit$learner_state <- state
    new_data <- d[1:6, c("patient_id", "gene_c", "gene_a", "gene_b")]
    with_endpoints <- new_data
    with_endpoints$time <- -999
    with_endpoints$event <- 1L
    prediction <- sb_predict_final(fit, new_data, "patient_id")
    endpoint_prediction <- sb_predict_final(
      fit,
      with_endpoints,
      "patient_id"
    )

    expect_identical(cv$status, "success", info = learner)
    expect_identical(cv$complete_oof, TRUE, info = learner)
    expect_identical(fit$status, "success", info = learner)
    expect_identical(prediction$status, "success", info = learner)
    expect_equal(
      prediction$predictions$risk,
      endpoint_prediction$predictions$risk,
      tolerance = 0,
      info = learner
    )
    expect_identical(
      fit$provenance$timereg_version,
      as.character(utils::packageVersion("timereg")),
      info = learner
    )
    expect_false("engine_raw" %in% names(state), info = learner)
    expect_false(any(c("patient_ids", "time_outcome", "event") %in%
      names(state)), info = learner)
  }

  horizon_task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 5
  )
  fit <- sb_fit_final(
    horizon_task,
    sb_pipeline(learner = "prop_odds_survival")
  )
  prediction <- sb_predict_final(
    fit,
    d[1:6, c("patient_id", predictors)],
    "patient_id"
  )
  expect_identical(fit$status, "success")
  expect_identical(prediction$status, "success")
  expect_true(all(prediction$predictions$survival_probability >= 0 &
    prediction$predictions$survival_probability <= 1))
})
