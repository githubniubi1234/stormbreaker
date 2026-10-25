make_flexsurvspline_fixture <- function(
  n = 120L, p = 4L, seed = 4111L, tied = FALSE
) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    true_risk <- 0.75 * x[, 1L]
    if (p >= 2L) {
      true_risk <- true_risk - 0.45 * x[, 2L]
    }
    if (p >= 3L) {
      true_risk <- true_risk + 0.20 * x[, 3L]
    }
    event_time <- (
      -log(stats::runif(n)) / exp(true_risk)
    )^(1 / 1.35) * 8
    censor_time <- stats::rexp(n, rate = 1 / 16)
    time <- pmin(event_time, censor_time)
    if (tied) {
      time <- pmax(round(time * 2) / 2, 0.25)
    }
    feature_names <- c(
      "gene one",
      if (p >= 2L) "2gene" else NULL,
      if (p >= 3L) paste0("gene-", seq_len(p - 2L)) else NULL
    )
    data.frame(
      patient_id = sprintf("FSS%03d", seq_len(n)),
      time = time,
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    )
  })
}

flexsurvspline_expected_params <- function() {
  list(
    flexsurvspline_hazard_k1 = list(
      scale = "hazard",
      k = 1L,
      timescale = "log",
      spline = "rp",
      knot_rule = "analysis_event_log_time_quantiles"
    ),
    flexsurvspline_odds_k1 = list(
      scale = "odds",
      k = 1L,
      timescale = "log",
      spline = "rp",
      knot_rule = "analysis_event_log_time_quantiles"
    ),
    flexsurvspline_normal_k1 = list(
      scale = "normal",
      k = 1L,
      timescale = "log",
      spline = "rp",
      knot_rule = "analysis_event_log_time_quantiles"
    )
  )
}

test_that("three nonduplicate one-knot Royston-Parmar scales are admitted", {
  expected <- flexsurvspline_expected_params()

  expect_identical(sb_flexsurvspline_learner_ids(), names(expected))
  for (learner in names(expected)) {
    expect_identical(
      sb_flexsurvspline_params(learner),
      expected[[learner]]
    )
  }
  expect_false(any(c(
    "flexsurvspline_hazard_k0_fixed",
    "flexsurvspline_odds_k0_fixed",
    "flexsurvspline_normal_k0_fixed",
    "flexsurvspline_hazard_k2_fixed",
    "flexsurvspline_hazard_splines2ns_k1_fixed"
  ) %in% sb_flexsurvspline_learner_ids()))

  rows <- sb_pipeline_registry()
  rows <- rows[
    rows$learner %in% sb_flexsurvspline_learner_ids(), ,
    drop = FALSE
  ]
  expect_identical(nrow(rows), 3L * length(expected_fixed_learner_selectors("survival")))
  expect_setequal(
    unique(rows$selector),
    expected_fixed_learner_selectors("survival")
  )
  expect_true(all(grepl("flexsurv", rows$optional_dependency, fixed = TRUE)))
})

test_that("spline risk, knots, and direct survival equal native predictions", {
  skip_if_not_installed("flexsurv", minimum_version = "2.3.2")

  data <- make_flexsurvspline_fixture()
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[c(2L, 37L, 91L), , drop = FALSE]
  horizons <- as.numeric(stats::quantile(data$time, c(0.35, 0.65)))
  expected_knots <- c(
    min(log(data$time[data$event == 1L])),
    as.numeric(stats::quantile(log(data$time[data$event == 1L]), 0.5)),
    max(log(data$time[data$event == 1L]))
  )

  for (learner in sb_flexsurvspline_learner_ids()) {
    adapter <- sb_flexsurvspline_learner_adapter(learner)
    set.seed(4112L)
    caller_rng <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, caller_rng)

    engine_assessment <- sb_flexsurv_numerical_predictors(assessment, model$numerical_fitting)
    names(engine_assessment) <- model$engine_features
    native_risk <- as.numeric(stats::predict(
      model$engine,
      newdata = engine_assessment,
      type = "link",
      conf.int = FALSE
    )$.pred_link)
    native_early <- as.numeric(stats::predict(
      model$engine,
      newdata = engine_assessment,
      type = "survival",
      times = horizons[[1L]],
      conf.int = FALSE
    )$.pred_survival)
    native_late <- as.numeric(stats::predict(
      model$engine,
      newdata = engine_assessment,
      type = "survival",
      times = horizons[[2L]],
      conf.int = FALSE
    )$.pred_survival)

    risk <- adapter$predict(model, assessment)
    early <- adapter$predict_survival(
      model, x, data$time, data$event, assessment, horizons[[1L]]
    )
    late <- adapter$predict_survival(
      model, x, data$time, data$event, assessment, horizons[[2L]]
    )
    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = horizons[[1L]],
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, assessment)

    expect_equal(risk, native_risk, tolerance = 1e-12)
    expect_equal(early$probability, native_early, tolerance = 1e-12)
    expect_equal(late$probability, native_late, tolerance = 1e-12)
    expect_true(all(early$probability >= late$probability))
    expect_equal(model$prediction$knots, expected_knots, tolerance = 1e-12)
    expect_equal(state$knots, expected_knots, tolerance = 1e-12)
    expect_equal(replay$risk, risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      early$probability,
      tolerance = 0
    )
    expect_false(any(c(
      "engine", "engine_raw", "time", "event", "outcome", "data"
    ) %in% names(state)))
    expect_identical(early$ties, "not_applicable")
    expect_identical(
      state$survival_probability$source,
      paste0(
        "flexsurvspline_", model$params$scale,
        "_k1_rp_distribution_direct"
      )
    )
  }
})

test_that("spline knots are learned only from the supplied analysis fold", {
  skip_if_not_installed("flexsurv", minimum_version = "2.3.2")

  data <- make_flexsurvspline_fixture(n = 121L, seed = 4120L)
  data$time[[121L]] <- max(data$time) * 50
  data$event[[121L]] <- 1L
  analysis <- data[seq_len(120L), , drop = FALSE]
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  expected_analysis_knots <- c(
    min(log(analysis$time[analysis$event == 1L])),
    as.numeric(stats::quantile(
      log(analysis$time[analysis$event == 1L]), 0.5
    )),
    max(log(analysis$time[analysis$event == 1L]))
  )
  full_boundary <- max(log(data$time[data$event == 1L]))

  for (learner in sb_flexsurvspline_learner_ids()) {
    model <- sb_flexsurvspline_learner_adapter(learner)$fit(
      analysis[features], analysis$time, analysis$event
    )
    expect_equal(
      model$prediction$knots,
      expected_analysis_knots,
      tolerance = 1e-12
    )
    expect_false(isTRUE(all.equal(
      model$prediction$knots[[3L]], full_boundary,
      tolerance = 1e-12
    )))
  }
})

test_that("one-feature and tied-time spline fits remain portable", {
  skip_if_not_installed("flexsurv", minimum_version = "2.3.2")

  cases <- list(
    one_feature = make_flexsurvspline_fixture(
      n = 96L, p = 1L, seed = 4201L
    ),
    tied_times = make_flexsurvspline_fixture(
      n = 120L, p = 5L, seed = 4202L, tied = TRUE
    )
  )
  for (data in cases) {
    features <- setdiff(names(data), c("patient_id", "time", "event"))
    x <- data[features]
    horizons <- as.numeric(stats::quantile(data$time, c(0.3, 0.7)))
    for (learner in sb_flexsurvspline_learner_ids()) {
      adapter <- sb_flexsurvspline_learner_adapter(learner)
      model <- adapter$fit(x, data$time, data$event)
      early <- adapter$predict_survival(
        model, x, data$time, data$event,
        x[1:5, , drop = FALSE], horizons[[1L]]
      )$probability
      late <- adapter$predict_survival(
        model, x, data$time, data$event,
        x[1:5, , drop = FALSE], horizons[[2L]]
      )$probability
      expect_true(all(is.finite(adapter$predict(
        model, x[1:5, , drop = FALSE]
      ))))
      expect_true(all(early >= late))
    }
  }
})

test_that("one-knot spline methods execute through CV and final prediction", {
  skip_if_not_installed("flexsurv", minimum_version = "2.3.2")

  data <- make_flexsurvspline_fixture(n = 120L, seed = 4111L)
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

  for (learner in sb_flexsurvspline_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 3L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 4113L)
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
