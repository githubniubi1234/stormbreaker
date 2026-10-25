make_pamr_survival_fixture <- function(n = 90L, p = 12L, seed = 1721L) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n)
    colnames(x) <- paste0("gene_", seq_len(p))
    linear_predictor <- 0.8 * x[, 1L] - 0.55 * x[, 2L]
    event_time <- stats::rexp(n, rate = exp(linear_predictor) / 10)
    censor_time <- stats::rexp(n, rate = 1 / 18)
    data.frame(
      patient_id = sprintf("PS%03d", seq_len(n)),
      time = pmin(event_time, censor_time) + 0.05,
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    )
  })
}

test_that("PAMR survival canonicalizes equivalent learner identifiers", {
  plain <- sb_pipeline("none", "pamr_survival")
  named <- c(requested = "pamr_survival")
  classed <- structure(
    "pamr_survival",
    class = c("stormbreaker_test_learner", "character")
  )

  expect_identical(sb_pipeline("none", named), plain)
  expect_identical(sb_pipeline("none", classed), plain)
})

test_that("PAMR survival reports its native one-input boundary", {
  skip_if_not_installed("pamr", minimum_version = "1.56")
  data <- make_pamr_survival_fixture(n = 72L, p = 2L, seed = 1720L)

  expect_error(
    sb_pamr_survival_learner_adapter()$fit(
      data["gene_1"], data$time, data$event
    ),
    class = "stormbreaker_pamr_survival_too_few_features"
  )
})

test_that("PAMR survival risk is the native early-survival posterior", {
  skip_if_not_installed("pamr", minimum_version = "1.56")

  data <- make_pamr_survival_fixture(n = 60L, p = 90L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  adapter <- sb_pamr_survival_learner_adapter()
  params <- adapter$params
  set.seed(2721)
  rng_before <- .Random.seed
  model <- adapter$fit(x, data$time, data$event)
  expect_identical(.Random.seed, rng_before)

  native_data <- list(
    x = t(data.matrix(x)),
    survival.time = data$time,
    censoring.status = data$event,
    geneid = features,
    genenames = features
  )
  native <- NULL
  invisible(utils::capture.output(
    native <- pamr::pamr.train(
      native_data,
      threshold = params$threshold,
      scale.sd = params$scale_sd,
      offset.percent = params$offset_percent,
      prior = NULL,
      remove.zeros = params$remove_zero_thresholds,
      sign.contrast = "both",
      ngroup.survival = params$survival_groups
    ),
    type = "output"
  ))
  assessment <- x[1:8, , drop = FALSE]
  native_posterior <- pamr::pamr.predict(
    native,
    t(data.matrix(assessment)),
    threshold = params$threshold,
    type = "posterior"
  )
  risk <- adapter$predict(model, assessment)
  state <- unserialize(serialize(adapter$extract_state(model), NULL))

  expect_equal(risk, unname(native_posterior[, "1"]), tolerance = 0)
  expect_equal(
    adapter$predict_state(state, assessment)$risk,
    risk,
    tolerance = 0
  )
  expect_equal(
    model$audit$survival_cutoff,
    unname(native$cutoffs.survival),
    tolerance = 0
  )
  expect_setequal(
    names(model$engine),
    c(
      "centroids", "centroid.overall", "sd", "threshold.scale",
      "se.scale", "sign.contrast", "prior", "norm.cent"
    )
  )
  expect_identical(state$training_data_retention, "none")
})

test_that("PAMR survival horizon output uses analysis-fold Cox calibration", {
  skip_if_not_installed("pamr", minimum_version = "1.56")

  data <- make_pamr_survival_fixture(seed = 1722L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[1:7, , drop = FALSE]
  horizon <- unname(stats::quantile(data$time, 0.4))
  adapter <- sb_pamr_survival_learner_adapter()
  model <- adapter$fit(x, data$time, data$event)
  analysis_risk <- adapter$predict(model, x)
  calibration <- survival::coxph(
    survival::Surv(data$time, data$event) ~ analysis_risk,
    ties = "breslow",
    singular.ok = FALSE
  )
  coefficient <- unname(as.numeric(stats::coef(calibration)))
  cumulative_hazard <- sb_breslow_baseline_hazard(
    data$time,
    data$event,
    coefficient * analysis_risk,
    horizon
  )
  observed <- adapter$predict_survival(
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

  expect_equal(observed$calibration_coefficient, coefficient, tolerance = 0)
  expect_equal(
    observed$probability,
    sb_survival_probability_from_hazard(
      coefficient * adapter$predict(model, assessment),
      cumulative_hazard
    ),
    tolerance = 0
  )
  expect_equal(replay$risk, adapter$predict(model, assessment), tolerance = 0)
  expect_equal(
    replay$survival_probability,
    observed$probability,
    tolerance = 0
  )
  expect_identical(observed$source, sb_pamr_survival_source())
})

test_that("PAMR survival executes through CV and final prediction", {
  skip_if_not_installed("pamr", minimum_version = "1.56")

  data <- make_pamr_survival_fixture(seed = 1723L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  horizon <- unname(stats::quantile(data$time, 0.35))
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    features,
    "month",
    prediction_horizon = horizon
  )
  pipeline <- sb_pipeline(learner = "pamr_survival")
  comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 1724L)
  fit <- sb_fit_final(task, pipeline)
  prediction <- sb_predict_final(
    unserialize(serialize(fit, NULL)),
    data[1:6, c("patient_id", rev(features)), drop = FALSE],
    "patient_id"
  )

  expect_identical(comparison$status, "success")
  expect_identical(comparison$complete_oof, TRUE)
  expect_identical(comparison$complete_survival_probability_oof, TRUE)
  expect_identical(fit$status, "success")
  expect_identical(prediction$status, "success")
  expect_length(prediction$predictions$risk, 6L)
  expect_length(prediction$predictions$survival_probability, 6L)
  expect_identical(
    fit$provenance$pamr_version,
    as.character(utils::packageVersion("pamr"))
  )
})
