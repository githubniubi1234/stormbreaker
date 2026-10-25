make_sb_survival_weak_calibration_fixture <- function() {
  transformed_prediction <- c(0.8, -0.4, 1.2, -0.9, 0.3, 1.0, -1.1, 0.1)
  data.frame(
    patient_id = paste0("weak_", seq_along(transformed_prediction)),
    fold = rep(1:4, each = 2L),
    time = c(0.7, 1.3, 2.1, 3.4, 4.0, 5.2, 6.1, 7.3),
    event = c(1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L),
    survival_probability = exp(-exp(transformed_prediction)),
    stringsAsFactors = FALSE
  )
}

manual_sb_event_probability_at <- function(time, event, horizon) {
  event_times <- sort(unique(time[time <= horizon & event == 1L]))
  survival_probability <- 1
  for (event_time in event_times) {
    n_risk <- sum(time >= event_time)
    n_event <- sum(time == event_time & event == 1L)
    survival_probability <- survival_probability * (1 - n_event / n_risk)
  }
  1 - survival_probability
}

manual_sb_no_tie_cox_slope <- function(time, event, predictor, horizon) {
  administrative_time <- pmin(time, horizon)
  event_by_horizon <- as.integer(event == 1L & time <= horizon)
  event_rows <- which(event_by_horizon == 1L)
  partial_log_likelihood <- function(coefficient) {
    sum(vapply(event_rows, function(row) {
      risk_set <- administrative_time >= administrative_time[[row]]
      coefficient * predictor[[row]] - log(sum(exp(
        coefficient * predictor[risk_set]
      )))
    }, numeric(1)))
  }
  stats::optimize(
    function(coefficient) -partial_log_likelihood(coefficient),
    interval = c(-10, 10),
    tol = 1e-12
  )$minimum
}

manual_sb_efron_cox_slope <- function(time, event, predictor, horizon) {
  administrative_time <- pmin(time, horizon)
  event_by_horizon <- as.integer(event == 1L & time <= horizon)
  event_times <- sort(unique(administrative_time[event_by_horizon == 1L]))
  partial_log_likelihood <- function(coefficient) {
    total <- 0
    for (event_time in event_times) {
      event_rows <- which(
        administrative_time == event_time & event_by_horizon == 1L
      )
      risk_rows <- which(administrative_time >= event_time)
      n_tied_events <- length(event_rows)
      risk_sum <- sum(exp(coefficient * predictor[risk_rows]))
      tied_event_sum <- sum(exp(coefficient * predictor[event_rows]))
      denominators <- vapply(
        0:(n_tied_events - 1L),
        function(tie_index) {
          risk_sum - tie_index / n_tied_events * tied_event_sum
        },
        numeric(1)
      )
      total <- total + coefficient * sum(predictor[event_rows]) -
        sum(log(denominators))
    }
    total
  }
  stats::optimize(
    function(coefficient) -partial_log_likelihood(coefficient),
    interval = c(-10, 10),
    tol = 1e-12
  )$minimum
}

call_sb_survival_weak_calibration <- function(
    data,
    horizon = 4,
    fold = data$fold,
    prediction_status = rep("success", nrow(data))) {
  sb_survival_weak_calibration_metric(
    patient_id = data$patient_id,
    time = data$time,
    event = data$event,
    survival_probability = data$survival_probability,
    prediction_status = prediction_status,
    prediction_horizon = horizon,
    time_unit = "month",
    scope = "test_fixed_horizon",
    fold = fold
  )
}

test_that("survival weak calibration matches independent scalar oracles", {
  data <- make_sb_survival_weak_calibration_fixture()
  horizon <- 4

  observed <- call_sb_survival_weak_calibration(data, horizon)
  transformed_prediction <- log(-log(data$survival_probability))
  expected_event_probability <- manual_sb_event_probability_at(
    data$time,
    data$event,
    horizon
  )
  expected_mean_prediction <- mean(-expm1(log(
    data$survival_probability
  )))
  expected_slope <- manual_sb_no_tie_cox_slope(
    data$time,
    data$event,
    transformed_prediction,
    horizon
  )

  expect_identical(observed$status, "success")
  expect_identical(observed$input_status, "success")
  expect_identical(observed$mean_calibration$status, "success")
  expect_identical(observed$slope_fit$status, "success")
  expect_equal(
    observed$mean_calibration$observed_event_probability,
    expected_event_probability,
    tolerance = 1e-12
  )
  expect_equal(
    observed$mean_calibration$mean_predicted_event_probability,
    expected_mean_prediction,
    tolerance = 1e-15
  )
  expect_equal(
    observed$mean_calibration_ratio,
    expected_event_probability / expected_mean_prediction,
    tolerance = 1e-12
  )
  expect_equal(observed$calibration_slope, expected_slope, tolerance = 5e-6)
  expect_identical(observed$mean_calibration$ratio_direction, "observed_over_expected")
  expect_identical(observed$mean_calibration_target, 1)
  expect_identical(observed$slope_target, 1)
  expect_identical(observed$slope_fit$ties, "efron")
  expect_identical(
    observed$slope_fit$administrative_censoring_at_horizon,
    TRUE
  )
})

test_that("survival weak calibration uses the Efron tied-event likelihood", {
  transformed_prediction <- c(-1.1, 0.2, 0.8, -0.4, 1.5, 0.5)
  data <- data.frame(
    patient_id = paste0("tie_", seq_along(transformed_prediction)),
    fold = rep(1L, length(transformed_prediction)),
    time = c(1, 1, 2, 3, 4, 5),
    event = c(1L, 1L, 1L, 0L, 1L, 0L),
    survival_probability = exp(-exp(transformed_prediction)),
    stringsAsFactors = FALSE
  )

  observed <- call_sb_survival_weak_calibration(data, horizon = 5)
  expected_efron <- manual_sb_efron_cox_slope(
    data$time,
    data$event,
    transformed_prediction,
    5
  )
  breslow <- unname(stats::coef(survival::coxph(
    survival::Surv(time, event) ~ transformed_prediction,
    data = data,
    ties = "breslow"
  ))[[1L]])

  expect_identical(observed$slope_fit$status, "success")
  expect_identical(observed$slope_fit$ties, "efron")
  expect_equal(observed$calibration_slope, expected_efron, tolerance = 5e-6)
  expect_gt(abs(observed$calibration_slope - breslow), 0.05)
})

test_that("survival weak calibration preserves row evidence and row invariance", {
  data <- make_sb_survival_weak_calibration_fixture()
  original <- call_sb_survival_weak_calibration(data)
  order <- c(8, 2, 6, 1, 7, 3, 5, 4)
  shuffled <- call_sb_survival_weak_calibration(data[order, ])

  expect_equal(
    shuffled$mean_calibration_ratio,
    original$mean_calibration_ratio,
    tolerance = 1e-12
  )
  expect_equal(
    shuffled$calibration_slope,
    original$calibration_slope,
    tolerance = 1e-12
  )
  expect_identical(original$patient_predictions$row_instance, seq_len(nrow(data)))
  expect_identical(original$patient_predictions$patient_id, data$patient_id)
  expect_equal(
    original$patient_predictions$predicted_event_probability,
    -expm1(log(data$survival_probability)),
    tolerance = 0
  )
  expect_equal(
    original$patient_predictions$transformed_prediction,
    log(-log(data$survival_probability)),
    tolerance = 0
  )
  expect_identical(
    original$patient_predictions$administrative_time,
    pmin(data$time, 4)
  )
  expect_identical(
    original$patient_predictions$event_by_horizon,
    as.integer(data$event == 1L & data$time <= 4)
  )

  sampled <- data[c(1, 1, 3, 4, 6, 6, 8, 2), ]
  repeated <- call_sb_survival_weak_calibration(sampled)
  expect_identical(repeated$patient_predictions$row_instance, seq_len(nrow(sampled)))
  expect_identical(repeated$patient_predictions$patient_id, sampled$patient_id)
})

test_that("survival weak calibration ignores outcomes strictly after horizon", {
  data <- make_sb_survival_weak_calibration_fixture()
  changed <- data
  later <- changed$time > 4
  changed$time[later] <- changed$time[later] + 100
  changed$event[later] <- 1L - changed$event[later]

  original <- call_sb_survival_weak_calibration(data)
  modified <- call_sb_survival_weak_calibration(changed)

  expect_equal(
    modified$mean_calibration_ratio,
    original$mean_calibration_ratio,
    tolerance = 0
  )
  expect_equal(
    modified$calibration_slope,
    original$calibration_slope,
    tolerance = 0
  )
})

test_that("survival weak calibration keeps exact-horizon event semantics", {
  data <- data.frame(
    patient_id = paste0("horizon_", 1:7),
    fold = rep(1L, 7L),
    time = c(1, 2, 4, 4, 5, 6, 7),
    event = c(1L, 0L, 1L, 0L, 1L, 0L, 1L),
    survival_probability = exp(-exp(c(-0.5, 0.7, -0.1, 1.1, 0.4, -0.8, 0.2))),
    stringsAsFactors = FALSE
  )

  observed <- call_sb_survival_weak_calibration(data, horizon = 4)

  expect_identical(observed$patient_predictions$event_by_horizon, c(
    1L, 0L, 1L, 0L, 0L, 0L, 0L
  ))
  expect_identical(observed$patient_predictions$administrative_time, c(
    1, 2, 4, 4, 4, 4, 4
  ))
  expect_equal(
    observed$mean_calibration$observed_event_probability,
    manual_sb_event_probability_at(data$time, data$event, 4),
    tolerance = 0
  )
})

test_that("mean calibration survives boundary and constant predictions", {
  data <- make_sb_survival_weak_calibration_fixture()
  boundary <- data
  boundary$survival_probability[1:2] <- c(0, 1)
  boundary_result <- call_sb_survival_weak_calibration(boundary)

  expect_identical(boundary_result$mean_calibration$status, "success")
  expect_true(is.finite(boundary_result$mean_calibration_ratio))
  expect_identical(boundary_result$slope_fit$status, "not_evaluable")
  expect_identical(
    boundary_result$slope_fit$error_class,
    "stormbreaker_survival_weak_calibration_boundary_probability"
  )
  expect_identical(is.na(boundary_result$calibration_slope), TRUE)
  expect_identical(boundary_result$status, "incomplete")
  expect_identical(boundary_result$probability_clipping, "none")

  constant <- data
  constant$survival_probability <- rep(0.7, nrow(constant))
  constant_result <- call_sb_survival_weak_calibration(constant)

  expect_identical(constant_result$mean_calibration$status, "success")
  expect_true(is.finite(constant_result$mean_calibration_ratio))
  expect_identical(constant_result$slope_fit$status, "not_evaluable")
  expect_identical(
    constant_result$slope_fit$error_class,
    "stormbreaker_survival_weak_calibration_constant_prediction"
  )
})

test_that("events after horizon cannot rescue the calibration slope", {
  data <- data.frame(
    patient_id = paste0("later_", 1:6),
    fold = rep(1L, 6L),
    time = c(5, 6, 7, 8, 9, 10),
    event = c(1L, 0L, 1L, 0L, 1L, 0L),
    survival_probability = exp(-exp(c(-1, -0.5, 0, 0.5, 1, 1.5))),
    stringsAsFactors = FALSE
  )

  observed <- call_sb_survival_weak_calibration(data, horizon = 4)

  expect_identical(observed$mean_calibration$status, "success")
  expect_identical(observed$mean_calibration$observed_event_probability, 0)
  expect_identical(observed$mean_calibration_ratio, 0)
  expect_identical(observed$slope_fit$status, "not_evaluable")
  expect_identical(
    observed$slope_fit$error_class,
    "stormbreaker_survival_weak_calibration_no_event"
  )
  expect_identical(is.na(observed$calibration_slope), TRUE)
})

test_that("mean horizon support cannot erase an independently estimable slope", {
  transformed_prediction <- c(-0.8, 0.4, 1.1, -0.2, 0.7, -1.0, 0.1, 1.4)
  data <- data.frame(
    patient_id = paste0("short_", seq_along(transformed_prediction)),
    fold = rep(1L, length(transformed_prediction)),
    time = c(0.5, 0.9, 1.4, 2.0, 2.6, 3.1, 3.7, 4.4),
    event = c(1L, 0L, 1L, 1L, 0L, 1L, 0L, 0L),
    survival_probability = exp(-exp(transformed_prediction)),
    stringsAsFactors = FALSE
  )

  observed <- call_sb_survival_weak_calibration(data, horizon = 10)
  expected_slope <- manual_sb_no_tie_cox_slope(
    data$time,
    data$event,
    transformed_prediction,
    10
  )

  expect_identical(observed$mean_calibration$status, "not_evaluable")
  expect_identical(
    observed$mean_calibration$error_class,
    "stormbreaker_survival_weak_calibration_horizon_unsupported"
  )
  expect_identical(is.na(observed$mean_calibration_ratio), TRUE)
  expect_identical(observed$slope_fit$status, "success")
  expect_equal(observed$calibration_slope, expected_slope, tolerance = 5e-6)
  expect_identical(observed$status, "incomplete")

  absorbing <- data.frame(
    patient_id = paste0("absorbing_", 1:4),
    fold = rep(1L, 4L),
    time = c(1, 2, 3, 4),
    event = rep(1L, 4L),
    survival_probability = exp(-exp(c(-1, -0.2, 0.4, 1.1))),
    stringsAsFactors = FALSE
  )
  absorbing_result <- call_sb_survival_weak_calibration(
    absorbing,
    horizon = 10
  )
  expect_identical(absorbing_result$mean_calibration$status, "success")
  expect_identical(
    absorbing_result$mean_calibration$horizon_support,
    "survival_zero_before_horizon"
  )
  expect_identical(
    absorbing_result$mean_calibration$observed_event_probability,
    1
  )
})

test_that("survival weak calibration fit conditions stay component-local", {
  data <- make_sb_survival_weak_calibration_fixture()
  baseline <- call_sb_survival_weak_calibration(data)
  original_fit <- sb_survival_weak_calibration_fit_slope

  testthat::local_mocked_bindings(
    sb_survival_weak_calibration_fit_slope = function(fit_data) {
      warning("secondary slope warning")
      original_fit(fit_data)
    }
  )
  warned <- call_sb_survival_weak_calibration(data)

  expect_identical(warned$mean_calibration, baseline$mean_calibration)
  expect_identical(warned$slope_fit$status, "failed")
  expect_identical(
    warned$slope_fit$error_class,
    "stormbreaker_survival_weak_calibration_cox_warning"
  )
  expect_match(warned$slope_fit$error_message, "secondary slope warning", fixed = TRUE)

  testthat::local_mocked_bindings(
    sb_survival_weak_calibration_fit_slope = function(fit_data) {
      condition <- simpleError("secondary slope error")
      class(condition) <- c("test_weak_slope_error", class(condition))
      stop(condition)
    }
  )
  failed <- call_sb_survival_weak_calibration(data)

  expect_identical(failed$mean_calibration, baseline$mean_calibration)
  expect_identical(failed$slope_fit$status, "failed")
  expect_identical(failed$slope_fit$error_class, "test_weak_slope_error")
  expect_identical(failed$slope_fit$error_message, "secondary slope error")
})

test_that("survival weak calibration exposes input and unexpected failures", {
  data <- make_sb_survival_weak_calibration_fixture()
  failed_status <- rep("success", nrow(data))
  failed_status[[2L]] <- "failed"
  failed <- call_sb_survival_weak_calibration(
    data,
    prediction_status = failed_status
  )

  expect_identical(failed$status, "incomplete")
  expect_identical(failed$input_status, "prediction_failed")
  expect_identical(failed$mean_calibration$status, "not_evaluated")
  expect_identical(failed$slope_fit$status, "not_fitted")

  condition <- simpleError("outer weak calibration error")
  class(condition) <- c("test_weak_calibration_error", class(condition))
  outer <- sb_survival_weak_calibration_error_metric(
    error = condition,
    prediction_horizon = 4,
    time_unit = "month",
    scope = "test_fixed_horizon",
    n = nrow(data),
    fold = data$fold
  )

  expect_identical(outer$status, "metric_failed")
  expect_identical(outer$input_status, "unexpected_error")
  expect_identical(outer$error_class, "test_weak_calibration_error")
  expect_identical(outer$error_message, "outer weak calibration error")
  expect_identical(outer$n_folds, 4L)
})
