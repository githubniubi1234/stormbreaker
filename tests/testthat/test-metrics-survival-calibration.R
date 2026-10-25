make_sb_survival_calibration_fixture <- function(
    n = 160L,
    seed = 20260820L,
    horizon = 2) {
  withr::with_seed(seed, {
    prognostic_index <- stats::rnorm(n)
    event_time <- stats::rexp(
      n,
      rate = 0.15 * exp(0.45 * prognostic_index)
    )
    censoring_time <- stats::rexp(n, rate = 0.05)
    time <- pmin(event_time, censoring_time)
    event <- as.integer(event_time <= censoring_time)
    survival_probability <- exp(
      -0.15 * horizon * exp(0.45 * prognostic_index)
    )
    data.frame(
      patient_id = sprintf("C%03d", seq_len(n)),
      time = time,
      event = event,
      survival_probability = survival_probability,
      status = rep("success", n),
      fold = rep(seq_len(5L), length.out = n),
      stringsAsFactors = FALSE
    )
  })
}

manual_three_knot_rcs_basis <- function(x, knots) {
  knot1 <- knots[[1L]]
  knot2 <- knots[[2L]]
  knot3 <- knots[[3L]]
  scale <- (knot3 - knot1)^(2 / 3)
  nonlinear <- pmax((x - knot1) / scale, 0)^3 +
    (
      (knot2 - knot1) * pmax((x - knot3) / scale, 0)^3 -
        (knot3 - knot1) * pmax((x - knot2) / scale, 0)^3
    ) / (knot3 - knot2)
  cbind(rcs_linear = x, rcs_nonlinear = nonlinear)
}

manual_survival_calibration_oracle <- function(data, horizon) {
  x <- log(-log(data$survival_probability))
  knots <- unname(stats::quantile(
    x,
    probs = c(0.10, 0.50, 0.90),
    names = FALSE,
    type = 7
  ))
  basis <- manual_three_knot_rcs_basis(x, knots)
  fit_data <- data.frame(
    time = data$time,
    event = data$event,
    rcs_linear = basis[, 1L],
    rcs_nonlinear = basis[, 2L]
  )
  fit <- survival::coxph(
    survival::Surv(time, event) ~ rcs_linear + rcs_nonlinear,
    data = fit_data,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
  new_data <- fit_data
  new_data$time <- rep(horizon, nrow(new_data))
  new_data$event <- rep(1L, nrow(new_data))
  cumulative_hazard_predict <- unname(stats::predict(
    fit,
    newdata = new_data,
    type = "expected"
  ))

  baseline <- survival::basehaz(fit, centered = TRUE)
  eligible <- baseline$time <= horizon
  baseline_horizon <- if (any(eligible)) {
    unname(utils::tail(baseline$hazard[eligible], 1L))
  } else {
    0
  }
  centered_basis <- sweep(basis, 2L, fit$means, FUN = "-")
  centered_lp <- drop(centered_basis %*% stats::coef(fit))
  cumulative_hazard_centered <- if (baseline_horizon == 0) {
    rep(0, nrow(data))
  } else {
    exp(log(baseline_horizon) + centered_lp)
  }
  calibrated_survival <- exp(-cumulative_hazard_predict)
  absolute_error <- abs(
    calibrated_survival - data$survival_probability
  )
  list(
    knots = knots,
    fit = fit,
    cumulative_hazard_predict = cumulative_hazard_predict,
    cumulative_hazard_centered = cumulative_hazard_centered,
    calibrated_survival_probability = calibrated_survival,
    ici = mean(absolute_error),
    e50 = unname(stats::quantile(absolute_error, 0.50, type = 7)),
    e90 = unname(stats::quantile(absolute_error, 0.90, type = 7))
  )
}

call_sb_survival_calibration <- function(data, horizon = 2, fold = data$fold) {
  sb_survival_calibration_metric(
    patient_id = data$patient_id,
    time = data$time,
    event = data$event,
    survival_probability = data$survival_probability,
    prediction_status = data$status,
    prediction_horizon = horizon,
    time_unit = "year",
    scope = "pooled_patient_level_out_of_fold",
    fold = fold
  )
}

test_that("survival calibration matches an independent three-knot RCS oracle", {
  data <- make_sb_survival_calibration_fixture()
  expected <- manual_survival_calibration_oracle(data, horizon = 2)
  observed <- call_sb_survival_calibration(data)

  expect_identical(observed$metric, "survival_calibration_curve")
  expect_identical(observed$calibration_level, "moderate")
  expect_identical(
    observed$method,
    "cox_ph_three_knot_restricted_natural_cubic_spline"
  )
  expect_identical(observed$scope, "pooled_patient_level_out_of_fold")
  expect_identical(observed$purpose, "descriptive_evaluation_only")
  expect_identical(observed$event_probability_transform, "complementary_log_log")
  expect_identical(observed$probability_clipping, "none")
  expect_identical(observed$ties, "efron")
  expect_identical(observed$baseline_reference, "coxph_mean_basis_centered")
  expect_identical(observed$status, "success")
  expect_identical(observed$input_status, "success")
  expect_equal(observed$knots$quantile_probability, c(0.10, 0.50, 0.90))
  expect_equal(observed$knots$transformed_value, expected$knots, tolerance = 0)
  expect_equal(
    expected$cumulative_hazard_predict,
    expected$cumulative_hazard_centered,
    tolerance = 1e-11
  )
  expect_equal(
    observed$patient_calibration$calibrated_survival_probability,
    expected$calibrated_survival_probability,
    tolerance = 1e-10
  )
  expect_equal(observed$ici, expected$ici, tolerance = 1e-11)
  expect_equal(observed$e50, expected$e50, tolerance = 1e-11)
  expect_equal(observed$e90, expected$e90, tolerance = 1e-11)
  expect_equal(
    observed$fit$means,
    unname(observed$fit$means),
    ignore_attr = TRUE
  )
  expect_identical(
    observed$patient_calibration$row_instance,
    seq_len(nrow(data))
  )
  expect_identical(
    observed$patient_calibration$patient_id,
    data$patient_id
  )
  expect_identical(observed$n_folds, 5L)
  expect_identical(observed$n_evaluable_folds, 5L)
})

test_that("survival calibration keeps extreme interior probabilities finite", {
  data <- make_sb_survival_calibration_fixture(n = 180L, seed = 17L)
  data$survival_probability[1:3] <- c(
    1e-20,
    1e-300,
    1 - .Machine$double.eps
  )
  observed <- call_sb_survival_calibration(data)

  expect_identical(observed$status, "success")
  expect_identical(observed$n_boundary_probability, 0L)
  expect_true(all(is.finite(
    observed$patient_calibration$transformed_prediction[1:3]
  )))
  expect_equal(
    observed$patient_calibration$transformed_prediction[1:3],
    log(-log(data$survival_probability[1:3])),
    tolerance = 0
  )
  expect_equal(
    observed$patient_calibration$predicted_event_probability[1:3],
    -expm1(log(data$survival_probability[1:3])),
    tolerance = 0
  )
  expect_true(all(is.finite(
    observed$patient_calibration$absolute_error[1:3]
  )))
})

test_that("survival calibration exposes exact boundaries without clipping", {
  data <- make_sb_survival_calibration_fixture()
  data$survival_probability[1:2] <- c(0, 1)
  observed <- call_sb_survival_calibration(data)

  expect_identical(observed$status, "incomplete")
  expect_identical(observed$input_status, "not_evaluable")
  expect_identical(
    observed$error_class,
    "stormbreaker_survival_calibration_boundary_probability"
  )
  expect_identical(observed$n_boundary_probability, 2L)
  expect_identical(nrow(observed$patient_calibration), nrow(data))
  expect_equal(
    observed$patient_calibration$survival_probability[1:2],
    c(0, 1),
    tolerance = 0
  )
  expect_true(all(is.na(c(observed$ici, observed$e50, observed$e90))))
})

test_that("survival calibration does not move duplicated knots", {
  data <- make_sb_survival_calibration_fixture()
  data$survival_probability <- rep(0.8, nrow(data))
  data$survival_probability[seq_len(5L)] <- 0.7
  observed <- call_sb_survival_calibration(data)

  expect_identical(observed$status, "incomplete")
  expect_identical(observed$input_status, "not_evaluable")
  expect_identical(
    observed$error_class,
    "stormbreaker_survival_calibration_duplicate_knots"
  )
  expect_identical(observed$n_unique_survival_probability, 2L)
  expect_identical(nrow(observed$knots), 3L)
  expect_true(anyDuplicated(observed$knots$transformed_value) > 0L)
})

test_that("survival calibration distinguishes horizon support cases", {
  data <- make_sb_survival_calibration_fixture(n = 180L, seed = 29L)
  after_all_followup <- call_sb_survival_calibration(
    data,
    horizon = max(data$time) + 1
  )
  expect_identical(after_all_followup$status, "incomplete")
  expect_identical(
    after_all_followup$error_class,
    "stormbreaker_survival_calibration_horizon_unsupported"
  )
  expect_identical(after_all_followup$n_followed_to_horizon, 0L)

  early <- data
  early$time <- early$time + 1
  before_first_event <- call_sb_survival_calibration(early, horizon = 0.5)
  expect_identical(before_first_event$status, "success")
  expect_equal(
    before_first_event$fit$baseline_cumulative_hazard_horizon,
    0,
    tolerance = 0
  )
  expect_true(all(
    before_first_event$patient_calibration$calibrated_survival_probability == 1
  ))

  after_last_event <- data
  after_last_event$event[after_last_event$time > 3] <- 0L
  after_last_event$time[[1L]] <- max(after_last_event$time) + 5
  after_last_event$event[[1L]] <- 0L
  supported <- call_sb_survival_calibration(after_last_event, horizon = 4)
  expect_identical(supported$status, "success")
  expect_gt(supported$n_followed_to_horizon, 0L)
  expect_true(max(after_last_event$time[after_last_event$event == 1L]) < 4)
})

test_that("survival calibration keeps complete follow-up after the horizon", {
  data <- make_sb_survival_calibration_fixture(n = 180L, seed = 31L)
  horizon <- 1
  changed <- data
  after <- changed$time > horizon
  changed$event[after] <- 1L - changed$event[after]
  if (sum(changed$event) == 0L) {
    changed$event[which(after)[1L]] <- 1L
  }
  original <- call_sb_survival_calibration(data, horizon = horizon)
  modified <- call_sb_survival_calibration(changed, horizon = horizon)

  expect_identical(original$status, "success")
  expect_identical(modified$status, "success")
  expect_false(isTRUE(all.equal(
    original$patient_calibration$calibrated_survival_probability,
    modified$patient_calibration$calibrated_survival_probability,
    tolerance = 1e-12
  )))
})

test_that("survival calibration retains bootstrap row instances and time zero", {
  data <- make_sb_survival_calibration_fixture(n = 180L, seed = 41L)
  data$time[1:2] <- 0
  data$event[1:2] <- c(1L, 0L)
  sampled <- data[c(1:170, 10, 10, 20, 20, 30, 30, 40, 40, 50, 50), ]
  rownames(sampled) <- NULL
  observed <- call_sb_survival_calibration(sampled)

  expect_identical(observed$status, "success")
  expect_true(anyDuplicated(sampled$patient_id) > 0L)
  expect_identical(
    observed$patient_calibration$row_instance,
    seq_len(nrow(sampled))
  )
  expect_identical(
    observed$patient_calibration$patient_id,
    sampled$patient_id
  )
  expect_identical(nrow(observed$patient_calibration), nrow(sampled))
})

test_that("survival calibration is row-order invariant after row alignment", {
  data <- make_sb_survival_calibration_fixture(n = 180L, seed = 43L)
  order <- c(seq(2L, nrow(data), by = 2L), seq(1L, nrow(data), by = 2L))
  observed <- call_sb_survival_calibration(data)
  shuffled <- call_sb_survival_calibration(data[order, ])
  aligned <- match(
    observed$patient_calibration$patient_id,
    shuffled$patient_calibration$patient_id
  )

  expect_identical(observed$status, "success")
  expect_identical(shuffled$status, "success")
  expect_equal(observed$ici, shuffled$ici, tolerance = 1e-12)
  expect_equal(observed$e50, shuffled$e50, tolerance = 1e-12)
  expect_equal(observed$e90, shuffled$e90, tolerance = 1e-12)
  expect_equal(
    observed$patient_calibration$calibrated_survival_probability,
    shuffled$patient_calibration$calibrated_survival_probability[aligned],
    tolerance = 1e-11
  )
})

test_that("survival calibration fit warnings and errors remain metric-local", {
  data <- make_sb_survival_calibration_fixture()
  original <- sb_survival_calibration_fit_cox
  testthat::local_mocked_bindings(
    sb_survival_calibration_fit_cox = function(fit_data) {
      warning("forced calibration warning")
      original(fit_data)
    }
  )
  warned <- call_sb_survival_calibration(data)
  expect_identical(warned$status, "fit_failed")
  expect_match(warned$fit$warnings, "forced calibration warning")
  expect_true(all(is.na(c(warned$ici, warned$e50, warned$e90))))
})

test_that("survival calibration error audits preserve the original condition", {
  error <- structure(
    list(message = "forced optional calibration error", call = NULL),
    class = c("forced_calibration_error", "error", "condition")
  )
  observed <- sb_survival_calibration_error_metric(
    error = error,
    prediction_horizon = 2,
    time_unit = "year",
    scope = "patient_level_external_cohort_at_locked_horizon",
    n = 25L
  )

  expect_identical(observed$status, "metric_failed")
  expect_identical(observed$input_status, "unexpected_error")
  expect_identical(observed$error_class, "forced_calibration_error")
  expect_identical(
    observed$error_message,
    "forced optional calibration error"
  )
  expect_identical(observed$n, 25L)
})
