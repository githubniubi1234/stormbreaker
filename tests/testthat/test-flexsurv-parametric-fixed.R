make_flexsurv_parametric_fixture <- function(
  n = 84L, p = 6L, seed = 9901L, tied = FALSE
) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n)
    true_risk <- 0.8 * x[, 1L]
    if (p >= 2L) {
      true_risk <- true_risk - 0.45 * x[, 2L]
    }
    event_time <- flexsurv::rgengamma(
      n,
      mu = 2.3 - true_risk,
      sigma = 0.65,
      Q = 0.55
    )
    censor_time <- stats::rexp(n, rate = 1 / 20)
    time <- pmin(event_time, censor_time)
    if (tied) {
      time <- pmax(round(time, 1), 0.1)
    }
    feature_names <- c(
      "gene one",
      if (p >= 2L) "2gene" else NULL,
      if (p >= 3L) paste0("gene-", seq_len(p - 2L)) else NULL
    )
    data.frame(
      patient_id = sprintf("FSP%03d", seq_len(n)),
      time = time,
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    )
  })
}

flexsurv_parametric_expected_params <- function() {
  list(
    flexsurv_gamma = list(
      distribution = "gamma",
      parameterization = "shape_rate",
      location_parameter = "rate",
      risk_transform = "log_location"
    ),
    flexsurv_gompertz = list(
      distribution = "gompertz",
      parameterization = "shape_rate",
      location_parameter = "rate",
      risk_transform = "log_location"
    ),
    flexsurv_gengamma = list(
      distribution = "gengamma",
      parameterization = "stable_mu_sigma_Q",
      location_parameter = "mu",
      risk_transform = "negative_location"
    )
  )
}

test_that("flexsurv admits three nonduplicate parametric survival families", {
  expected <- flexsurv_parametric_expected_params()

  expect_identical(sb_flexsurv_parametric_learner_ids(), names(expected))
  for (learner in names(expected)) {
    expect_identical(
      sb_flexsurv_parametric_params(learner),
      expected[[learner]]
    )
  }
  expect_false(any(c(
    "flexsurv_genf_fixed", "flexsurv_gengamma_orig_fixed",
    "flexsurv_genf_orig_fixed", "flexsurv_weibull_fixed",
    "flexsurv_exponential_fixed", "flexsurv_lognormal_fixed",
    "flexsurv_loglogistic_fixed"
  ) %in% sb_flexsurv_parametric_learner_ids()))

  rows <- sb_pipeline_registry()
  rows <- rows[
    rows$learner %in% sb_flexsurv_parametric_learner_ids(), ,
    drop = FALSE
  ]
  expect_identical(nrow(rows), 3L * length(expected_fixed_learner_selectors("survival")))
  expect_setequal(
    unique(rows$selector),
    expected_fixed_learner_selectors("survival")
  )
  expect_true(all(grepl("flexsurv", rows$optional_dependency, fixed = TRUE)))
})

test_that("flexsurv risk and horizon survival equal native predictions", {
  skip_if_not_installed("flexsurv", minimum_version = "2.3.2")

  data <- make_flexsurv_parametric_fixture()
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[c(3L, 29L, 71L), , drop = FALSE]
  horizon <- as.numeric(stats::quantile(data$time, 0.45))

  for (learner in sb_flexsurv_parametric_learner_ids()) {
    adapter <- sb_flexsurv_parametric_learner_adapter(learner)
    set.seed(9902L)
    caller_rng <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, caller_rng)

    engine_assessment <- sb_flexsurv_numerical_predictors(assessment, model$numerical_fitting)
    names(engine_assessment) <- model$engine_features
    native_location <- as.numeric(stats::predict(
      model$engine,
      newdata = engine_assessment,
      type = "link",
      conf.int = FALSE
    )$.pred_link)
    native_risk <- if (model$params$risk_transform == "log_location") {
      log(native_location)
    } else {
      -native_location
    }
    native_survival <- as.numeric(stats::predict(
      model$engine,
      newdata = engine_assessment,
      type = "survival",
      times = horizon,
      conf.int = FALSE
    )$.pred_survival)

    risk <- adapter$predict(model, assessment)
    survival <- adapter$predict_survival(
      model,
      x,
      data$time,
      data$event,
      assessment,
      horizon
    )
    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = horizon,
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, assessment)

    expect_equal(risk, native_risk, tolerance = 1e-12)
    expect_equal(survival$probability, native_survival, tolerance = 1e-12)
    expect_equal(replay$risk, risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      survival$probability,
      tolerance = 0
    )
    expect_identical(survival$ties, "not_applicable")
    expect_false("engine_raw" %in% names(state))
    expect_false(any(c("time", "event", "outcome") %in% names(state)))
    expect_identical(
      state$survival_probability$source,
      paste0(
        "flexsurv_", model$params$distribution,
        "_distribution_direct"
      )
    )
  }
})

test_that("flexsurv parametric families support one feature and tied times", {
  skip_if_not_installed("flexsurv", minimum_version = "2.3.2")

  cases <- list(
    one_feature = make_flexsurv_parametric_fixture(
      n = 64L, p = 1L, seed = 9903L
    ),
    tied_times = make_flexsurv_parametric_fixture(
      n = 84L, p = 5L, seed = 9904L, tied = TRUE
    )
  )
  for (data in cases) {
    features <- setdiff(names(data), c("patient_id", "time", "event"))
    x <- data[features]
    horizons <- as.numeric(stats::quantile(data$time, c(0.3, 0.7)))
    for (learner in sb_flexsurv_parametric_learner_ids()) {
      adapter <- sb_flexsurv_parametric_learner_adapter(learner)
      model <- adapter$fit(x, data$time, data$event)
      early <- adapter$predict_survival(
        model, x, data$time, data$event,
        x[1:5, , drop = FALSE], horizons[[1L]]
      )$probability
      late <- adapter$predict_survival(
        model, x, data$time, data$event,
        x[1:5, , drop = FALSE], horizons[[2L]]
      )$probability

      expect_true(all(is.finite(early)))
      expect_true(all(is.finite(late)))
      expect_true(all(early >= 0 & early <= 1))
      expect_true(all(late >= 0 & late <= 1))
      expect_true(all(early >= late))
    }
  }
})

test_that("flexsurv reports a nonfinite Hessian as a failed fit", {
  skip_if_not_installed("flexsurv", minimum_version = "2.3.2")

  data <- withr::with_seed(2105L, {
    n <- 44L
    p <- 52L
    x <- matrix(stats::rnorm(n * p), nrow = n)
    true_risk <- 0.8 * x[, 1L] - 0.45 * x[, 2L]
    event_time <- flexsurv::rgompertz(
      n,
      shape = 0.035,
      rate = exp(true_risk) / 14
    )
    censor_time <- stats::rexp(n, rate = 1 / 20)
    x <- as.data.frame(x)
    names(x) <- paste0("gene ", seq_len(p))
    list(
      x = x,
      time = pmin(event_time, censor_time),
      event = as.integer(event_time <= censor_time)
    )
  })

  expect_error(
    sb_flexsurv_parametric_learner_adapter(
      "flexsurv_gompertz"
    )$fit(data$x, data$time, data$event),
    class = "stormbreaker_flexsurv_nonconvergence"
  )
})

test_that("flexsurv parametric methods execute through CV and final prediction", {
  skip_if_not_installed("flexsurv", minimum_version = "2.3.2")

  data <- make_flexsurv_parametric_fixture(n = 78L, seed = 9906L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    features,
    "month",
    prediction_horizon = as.numeric(stats::quantile(data$time, 0.45))
  )
  new_data <- data[1:6, c("patient_id", rev(features)), drop = FALSE]

  for (learner in sb_flexsurv_parametric_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 3L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 9907L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      new_data,
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_true(comparison$complete_oof)
    expect_true(comparison$complete_survival_probability_oof)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(prediction$predictions$survival_probability >= 0 &
      prediction$predictions$survival_probability <= 1))
  }
})
