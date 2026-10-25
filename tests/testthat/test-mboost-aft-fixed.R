make_mboost_aft_fixture <- function(n = 84L, p = 6L, seed = 9811L) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n)
    signal <- -0.8 * x[, 1L]
    if (p >= 2L) {
      signal <- signal + 0.5 * x[, 2L]
    }
    event_time <- exp(2.2 + signal + stats::rnorm(n, sd = 0.65))
    censor_time <- exp(stats::rnorm(n, mean = 2.55, sd = 0.8))
    feature_names <- c(
      "gene one",
      if (p >= 2L) "2gene" else NULL,
      if (p >= 3L) paste0("gene-", seq_len(p - 2L)) else NULL
    )
    colnames(x) <- feature_names
    data.frame(
      patient_id = sprintf("MBA%03d", seq_len(n)),
      time = pmin(event_time, censor_time),
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    )
  })
}

mboost_aft_expected_specifications <- function() {
  data.frame(
    learner = c(
      "gamboost_weibull_aft",
      "gamboost_loglogistic_aft",
      "gamboost_lognormal_aft",
      "glmboost_weibull_aft",
      "glmboost_loglogistic_aft",
      "glmboost_lognormal_aft",
      "blackboost_weibull_aft",
      "blackboost_loglogistic_aft",
      "blackboost_lognormal_aft"
    ),
    algorithm = rep(c("gamboost", "glmboost", "blackboost"), each = 3L),
    family = rep(c("Weibull", "Loglog", "Lognormal"), times = 3L),
    distribution = rep(
      c("weibull", "loglogistic", "lognormal"),
      times = 3L
    ),
    stringsAsFactors = FALSE
  )
}

mboost_aft_native_survival <- function(distribution, link, scale, horizon) {
  transformed <- (log(horizon) - link) / scale
  switch(distribution,
    weibull = exp(-exp(transformed)),
    loglogistic = stats::plogis(-transformed),
    lognormal = stats::pnorm(transformed, lower.tail = FALSE)
  )
}

test_that("mboost AFT declares nine distribution-by-base-learner methods", {
  expected <- mboost_aft_expected_specifications()
  expect_identical(sb_mboost_aft_learner_ids(), expected$learner)

  for (i in seq_len(nrow(expected))) {
    learner <- expected$learner[[i]]
    params <- sb_mboost_aft_params(learner)
    expect_identical(params$algorithm, expected$algorithm[[i]])
    expect_identical(params$family, expected$family[[i]])
    expect_identical(params$distribution, expected$distribution[[i]])
    expect_identical(params$mstop, 100L)
    expect_identical(params$nu, 0.1)
    expect_identical(params$risk, "inbag")
    expect_identical(params$nuisance_range, c(0.1, 100))
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
  expect_identical(nrow(rows), 9L * length(expected_fixed_learner_selectors("survival")))
  expect_identical(anyDuplicated(rows[c("selector", "learner")]), 0L)
  expect_setequal(
    unique(rows$selector),
    expected_fixed_learner_selectors("survival")
  )
  expect_true(all(grepl("mboost", rows$optional_dependency, fixed = TRUE)))
  black_rows <- rows[grepl("^blackboost", rows$learner), , drop = FALSE]
  expect_true(all(grepl("partykit", black_rows$optional_dependency, fixed = TRUE)))
})

test_that("mboost AFT risk, direct survival, RNG, and state match native fits", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  data <- make_mboost_aft_fixture(seed = 9812L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[c(3L, 27L, 68L), , drop = FALSE]
  horizon <- as.numeric(stats::quantile(data$time, 0.45))

  for (learner in sb_mboost_aft_learner_ids()) {
    adapter <- sb_mboost_aft_learner_adapter(learner)
    set.seed(9813L)
    caller_rng <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, caller_rng)

    engine_assessment <- assessment
    names(engine_assessment) <- model$engine_features
    native_link <- as.numeric(stats::predict(
      model$engine,
      newdata = engine_assessment,
      type = "link"
    ))
    native_survival <- mboost_aft_native_survival(
      model$params$distribution,
      native_link,
      model$nuisance_scale,
      horizon
    )
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
    expect_equal(
      model$nuisance_scale,
      as.numeric(mboost::nuisance(model$engine)),
      tolerance = 0
    )
    expect_equal(risk, -native_link, tolerance = 0)
    expect_equal(survival$probability, native_survival, tolerance = 0)
    expect_identical(survival$ties, "not_applicable")
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

    if (identical(model$params$algorithm, "blackboost")) {
      repeated <- adapter$fit(x, data$time, data$event)
      expect_identical(
        serialize(model$engine, NULL),
        serialize(repeated$engine, NULL)
      )
    }
  }
})

test_that("mboost AFT supports one feature, p greater than n, and horizons", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  cases <- list(
    one_feature = make_mboost_aft_fixture(n = 54L, p = 1L, seed = 9814L),
    high_dimensional = make_mboost_aft_fixture(n = 44L, p = 52L, seed = 9815L)
  )
  for (data in cases) {
    features <- setdiff(names(data), c("patient_id", "time", "event"))
    x <- data[features]
    horizons <- as.numeric(stats::quantile(data$time, c(0.3, 0.7)))
    for (learner in sb_mboost_aft_learner_ids()) {
      adapter <- sb_mboost_aft_learner_adapter(learner)
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

test_that("mboost AFT methods execute through CV and final prediction", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  data <- make_mboost_aft_fixture(n = 72L, seed = 9816L)
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

  for (learner in sb_mboost_aft_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 3L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 9817L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      new_data,
      "patient_id"
    )

    expect_identical(pipeline$risk_type, "aft_log_time")
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
