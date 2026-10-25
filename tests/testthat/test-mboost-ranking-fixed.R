make_mboost_ranking_fixture <- function(n = 84L, p = 6L, seed = 9871L) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n)
    signal <- -0.9 * x[, 1L]
    if (p >= 2L) {
      signal <- signal + 0.6 * x[, 2L]
    }
    if (p >= 3L) {
      signal <- signal - 0.4 * x[, 3L]
    }
    event_time <- exp(2.2 + signal + stats::rnorm(n, sd = 0.65))
    censor_time <- exp(stats::rnorm(n, mean = 2.65, sd = 0.8))
    feature_names <- c(
      "gene one",
      if (p >= 2L) "2gene" else NULL,
      if (p >= 3L) paste0("gene-", seq_len(p - 2L)) else NULL
    )
    colnames(x) <- feature_names
    data.frame(
      patient_id = sprintf("MBR%03d", seq_len(n)),
      time = pmin(event_time, censor_time),
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    )
  })
}

mboost_ranking_expected_specifications <- function() {
  data.frame(
    learner = c(
      "gamboost_gehan_aft",
      "gamboost_unoc",
      "glmboost_gehan_aft",
      "glmboost_unoc",
      "blackboost_gehan_aft",
      "blackboost_unoc"
    ),
    algorithm = rep(c("gamboost", "glmboost", "blackboost"), each = 2L),
    family = rep(c("Gehan", "Cindex"), times = 3L),
    objective = rep(c("gehan_rank_aft", "uno_ipcw_concordance"), times = 3L),
    stringsAsFactors = FALSE
  )
}

mboost_ranking_manual_residual_survival <- function(
  analysis_time, analysis_event, analysis_link, assessment_link, horizon
) {
  fit <- survival::survfit(
    survival::Surv(log(analysis_time) - analysis_link, analysis_event) ~ 1
  )
  index <- findInterval(log(horizon) - assessment_link, fit$time)
  probability <- rep(1, length(index))
  use <- index > 0L
  probability[use] <- fit$surv[index[use]]
  unname(probability)
}

test_that("installed mboost ranking losses use longer-survival link direction", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")

  time <- c(1, 2, 3, 4)
  response <- survival::Surv(time, rep(1L, length(time)))
  weights <- rep(1, length(time))
  expect_lt(
    mboost::Gehan()@risk(response, log(time), weights),
    mboost::Gehan()@risk(response, -log(time), weights)
  )
  expect_lt(
    mboost::Cindex(sigma = 0.1, ipcw = 1)@risk(response, time, weights),
    mboost::Cindex(sigma = 0.1, ipcw = 1)@risk(response, -time, weights)
  )
})

test_that("mboost ranking declares six objective-by-base-learner methods", {
  expected <- mboost_ranking_expected_specifications()
  expect_identical(sb_mboost_ranking_learner_ids(), expected$learner)

  for (i in seq_len(nrow(expected))) {
    learner <- expected$learner[[i]]
    params <- sb_mboost_ranking_params(learner)
    expect_identical(params$algorithm, expected$algorithm[[i]])
    expect_identical(params$family, expected$family[[i]])
    expect_identical(params$objective, expected$objective[[i]])
    expect_identical(params$mstop, 100L)
    expect_identical(params$nu, 0.1)
    expect_identical(params$risk, "inbag")
    if (identical(params$family, "Cindex")) {
      expect_identical(params$sigma, 0.1)
      expect_identical(params$ipcw, 1)
      expect_identical(params$calibration_ties, "breslow")
    }
    if (identical(params$algorithm, "blackboost")) {
      expect_identical(params$seed, 1L)
    } else {
      expect_false("seed" %in% names(params))
    }
  }

  registry <- sb_pipeline_registry()
  rows <- registry[
    registry$learner %in% expected$learner, ,
    drop = FALSE
  ]
  expect_identical(nrow(rows), 6L * length(expected_fixed_learner_selectors("survival")))
  expect_identical(anyDuplicated(rows[c("selector", "learner")]), 0L)
  expect_setequal(
    unique(rows$selector),
    expected_fixed_learner_selectors("survival")
  )
  expect_true(all(grepl("mboost", rows$optional_dependency, fixed = TRUE)))
  black_rows <- rows[grepl("^blackboost", rows$learner), , drop = FALSE]
  expect_true(all(grepl("partykit", black_rows$optional_dependency, fixed = TRUE)))
})

test_that("mboost ranking risk and horizon survival match declared methods", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  data <- make_mboost_ranking_fixture(seed = 9872L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[c(4L, 29L, 71L), , drop = FALSE]
  horizon <- as.numeric(stats::quantile(data$time, 0.45))

  for (learner in sb_mboost_ranking_learner_ids()) {
    adapter <- sb_mboost_ranking_learner_adapter(learner)
    set.seed(9873L)
    caller_rng <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, caller_rng)

    engine_analysis <- x
    engine_assessment <- assessment
    names(engine_analysis) <- model$engine_features
    names(engine_assessment) <- model$engine_features
    native_analysis_link <- as.numeric(stats::predict(
      model$engine,
      newdata = engine_analysis,
      type = "link"
    ))
    native_assessment_link <- as.numeric(stats::predict(
      model$engine,
      newdata = engine_assessment,
      type = "link"
    ))
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

    expect_s3_class(model$engine, model$params$algorithm)
    expect_identical(as.integer(mboost::mstop(model$engine)), 100L)
    expect_equal(risk, -native_assessment_link, tolerance = 0)
    expect_equal(replay$risk, risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      survival$probability,
      tolerance = 0
    )
    expect_type(state$engine_raw, "raw")
    expect_identical(
      state$training_data_retention,
      "required_by_mboost_prediction_state"
    )
    expect_false(any(c("time", "event", "outcome") %in% names(state)))

    if (identical(model$params$family, "Gehan")) {
      expected_probability <- mboost_ranking_manual_residual_survival(
        data$time,
        data$event,
        native_analysis_link,
        native_assessment_link,
        horizon
      )
      expect_equal(survival$probability, expected_probability, tolerance = 0)
      expect_identical(
        survival$source,
        sb_mboost_ranking_direct_survival_source(learner)
      )
      expect_null(state$baseline_hazard)
      expect_true(is.list(state$residual_survival))
    } else {
      calibration <- survival::coxph(
        survival::Surv(data$time, data$event) ~ raw_risk,
        data = data.frame(
          time = data$time,
          event = data$event,
          raw_risk = -native_analysis_link
        ),
        ties = "breslow",
        singular.ok = FALSE
      )
      coefficient <- unname(as.numeric(stats::coef(calibration)))
      calibrated_analysis <- coefficient * -native_analysis_link
      calibrated_assessment <- coefficient * -native_assessment_link
      baseline <- sb_breslow_baseline_hazard(
        data$time,
        data$event,
        calibrated_analysis,
        horizon
      )
      expected_probability <- sb_survival_probability_from_hazard(
        calibrated_assessment,
        baseline
      )
      expect_equal(
        survival$calibration_coefficient,
        coefficient,
        tolerance = 0
      )
      expect_equal(survival$probability, expected_probability, tolerance = 0)
      expect_null(state$residual_survival)
      expect_equal(state$calibration_coefficient, coefficient, tolerance = 0)
      expect_true(is.list(state$baseline_hazard))
    }

    if (identical(model$params$algorithm, "blackboost")) {
      repeated <- adapter$fit(x, data$time, data$event)
      expect_identical(
        serialize(model$engine, NULL),
        serialize(repeated$engine, NULL)
      )
    }
  }
})

test_that("mboost ranking supports one feature and p greater than n", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  cases <- list(
    one_feature = make_mboost_ranking_fixture(n = 54L, p = 1L, seed = 9874L),
    high_dimensional = make_mboost_ranking_fixture(
      n = 44L,
      p = 52L,
      seed = 9875L
    )
  )
  for (data in cases) {
    features <- setdiff(names(data), c("patient_id", "time", "event"))
    x <- data[features]
    horizons <- as.numeric(stats::quantile(data$time, c(0.3, 0.7)))
    for (learner in sb_mboost_ranking_learner_ids()) {
      adapter <- sb_mboost_ranking_learner_adapter(learner)
      model <- adapter$fit(x, data$time, data$event)
      risk <- adapter$predict(model, x[1:5, , drop = FALSE])
      early <- adapter$predict_survival(
        model, x, data$time, data$event,
        x[1:5, , drop = FALSE], horizons[[1L]]
      )$probability
      late <- adapter$predict_survival(
        model, x, data$time, data$event,
        x[1:5, , drop = FALSE], horizons[[2L]]
      )$probability

      expect_true(all(is.finite(risk)))
      expect_true(all(is.finite(early)))
      expect_true(all(is.finite(late)))
      expect_true(all(early >= 0 & early <= 1))
      expect_true(all(late >= 0 & late <= 1))
      expect_true(all(early >= late))
    }
  }
})

test_that("mboost ranking methods execute through CV and final prediction", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  data <- make_mboost_ranking_fixture(n = 72L, seed = 9876L)
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

  for (learner in sb_mboost_ranking_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 3L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 9877L)
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
