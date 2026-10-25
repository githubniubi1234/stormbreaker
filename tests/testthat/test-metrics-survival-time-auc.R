make_sb_time_auc_curve <- function(time, event) {
  sb_fit_censoring_km(time, event)
}

make_sb_time_auc_stratum <- function(time, event, marker, horizon, curve,
                                     prediction_ok = NULL, fold = 1L) {
  if (is.null(prediction_ok)) {
    prediction_ok <- rep(TRUE, length(time))
  }
  sb_time_auc_stratum(
    time = time,
    event = event,
    survival_probability = 1 - marker,
    prediction_ok = prediction_ok,
    censoring_curve = curve,
    prediction_horizon = horizon,
    fold = fold
  )
}

test_that("time AUC matches a literal unequal-IPCW pair oracle", {
  horizon <- 3
  time <- c(1.5, 2.5, 4, 5)
  event <- c(1L, 1L, 0L, 1L)
  marker <- c(0.9, 0.5, 0.5, 0.2)
  curve <- make_sb_time_auc_curve(
    c(1, 2, 3, 4, 5),
    c(0L, 0L, 1L, 1L, 1L)
  )
  cases <- event == 1L & time <= horizon
  controls <- time > horizon
  case_weights <- 1 / vapply(
    time[cases],
    function(value) sb_censoring_survival_at(
      curve,
      value,
      left_limit = TRUE
    ),
    numeric(1)
  )
  pair_credit <- outer(
    marker[cases],
    marker[controls],
    function(case, control) {
      (case > control) + 0.5 * (case == control)
    }
  )
  expected_concordant <- sum(case_weights * rowSums(
    outer(marker[cases], marker[controls], `>`)
  ))
  expected_tied <- sum(case_weights * rowSums(
    outer(marker[cases], marker[controls], `==`)
  ))
  expected_numerator <- sum(case_weights * rowSums(pair_credit))
  expected_denominator <- sum(case_weights) * sum(controls)

  observed <- make_sb_time_auc_stratum(
    time,
    event,
    marker,
    horizon,
    curve
  )

  expect_identical(observed$status, "success")
  expect_identical(observed$n_cases, 2L)
  expect_identical(observed$n_controls, 2L)
  expect_identical(observed$n_excluded_censored, 0L)
  expect_equal(observed$case_weight_sum, sum(case_weights), tolerance = 1e-15)
  expect_equal(
    observed$weighted_concordant_pairs,
    expected_concordant,
    tolerance = 1e-15
  )
  expect_equal(
    observed$weighted_tied_pairs,
    expected_tied,
    tolerance = 1e-15
  )
  expect_equal(
    observed$weighted_pair_credit,
    expected_numerator,
    tolerance = 1e-15
  )
  expect_equal(
    observed$pair_denominator,
    expected_denominator,
    tolerance = 1e-15
  )
  expect_equal(
    observed$estimate,
    expected_numerator / expected_denominator,
    tolerance = 1e-15
  )
})

test_that("time AUC fixes the discrete horizon case and control boundary", {
  horizon <- 2
  time <- c(2, 2, 3, 1.5)
  event <- c(1L, 0L, 1L, 0L)
  marker <- c(0.8, 0.9, 0.2, 0.1)
  curve <- make_sb_time_auc_curve(1:4, rep(1L, 4L))

  observed <- make_sb_time_auc_stratum(
    time,
    event,
    marker,
    horizon,
    curve
  )

  expect_identical(observed$status, "success")
  expect_identical(observed$n_cases, 1L)
  expect_identical(observed$n_controls, 1L)
  expect_identical(observed$n_excluded_censored, 2L)
  expect_equal(observed$estimate, 1, tolerance = 0)

  before_only <- make_sb_time_auc_stratum(
    time,
    event,
    marker,
    horizon - .Machine$double.eps^0.5,
    curve
  )
  expect_identical(before_only$status, "not_evaluable")
  expect_identical(before_only$n_cases, 0L)
})

test_that("time AUC never requires the cancelled control G at the horizon", {
  horizon <- 2
  curve <- make_sb_time_auc_curve(
    c(1, 2, 2),
    c(1L, 1L, 0L)
  )
  expect_equal(
    sb_censoring_survival_at(curve, horizon, left_limit = TRUE),
    1,
    tolerance = 0
  )
  expect_equal(
    sb_censoring_survival_at(curve, horizon, left_limit = FALSE),
    0,
    tolerance = 0
  )

  observed <- make_sb_time_auc_stratum(
    time = c(2, 3),
    event = c(1L, 1L),
    marker = c(0.8, 0.2),
    horizon = horizon,
    curve = curve
  )

  expect_identical(observed$status, "success")
  expect_equal(observed$case_weight_sum, 1, tolerance = 0)
  expect_equal(observed$pair_denominator, 1, tolerance = 0)
  expect_equal(observed$estimate, 1, tolerance = 0)
})

test_that("time AUC exposes target-specific unevaluable and failure states", {
  curve <- make_sb_time_auc_curve(1:4, rep(1L, 4L))
  no_case <- make_sb_time_auc_stratum(
    c(3, 4), c(1L, 0L), c(0.8, 0.2), 2, curve
  )
  no_control <- make_sb_time_auc_stratum(
    c(1, 2), c(1L, 0L), c(0.8, 0.2), 2, curve
  )
  bad_prediction <- make_sb_time_auc_stratum(
    c(1, 3), c(1L, 1L), c(0.8, 0.2), 2, curve,
    prediction_ok = c(TRUE, FALSE)
  )
  zero_curve <- make_sb_time_auc_curve(c(1, 1), c(1L, 0L))
  zero_case_support <- make_sb_time_auc_stratum(
    c(2, 3), c(1L, 1L), c(0.8, 0.2), 2, zero_curve
  )
  empty_curve <- make_sb_time_auc_stratum(
    c(1, 3), c(1L, 1L), c(0.8, 0.2), 2,
    sb_empty_censoring_curve()
  )

  expect_identical(no_case$status, "not_evaluable")
  expect_identical(no_case$error_class, "stormbreaker_no_time_auc_case")
  expect_identical(no_control$status, "not_evaluable")
  expect_identical(no_control$error_class, "stormbreaker_no_time_auc_control")
  expect_identical(bad_prediction$status, "prediction_failed")
  expect_identical(
    bad_prediction$error_class,
    "stormbreaker_invalid_survival_prediction"
  )
  expect_identical(zero_case_support$status, "censoring_failed")
  expect_identical(
    zero_case_support$error_class,
    "stormbreaker_zero_censoring_survival"
  )
  expect_identical(empty_curve$status, "censoring_failed")
  expect_identical(
    empty_curve$error_class,
    "stormbreaker_missing_censoring_curve"
  )
  expect_true(all(is.na(c(
    no_case$estimate,
    no_control$estimate,
    bad_prediction$estimate,
    zero_case_support$estimate,
    empty_curve$estimate
  ))))
})

test_that("time AUC is order invariant and obeys strict ranking transforms", {
  time <- c(1, 1.5, 3, 4)
  event <- rep(1L, 4L)
  marker <- c(0.8, 0.6, 0.7, 0.5)
  curve <- make_sb_time_auc_curve(1:4, rep(1L, 4L))
  observed <- make_sb_time_auc_stratum(time, event, marker, 2, curve)
  shuffled <- make_sb_time_auc_stratum(
    time[c(4, 2, 1, 3)],
    event[c(4, 2, 1, 3)],
    marker[c(4, 2, 1, 3)],
    2,
    curve
  )
  transformed <- make_sb_time_auc_stratum(
    time,
    event,
    marker^2,
    2,
    curve
  )
  reversed <- make_sb_time_auc_stratum(
    time,
    event,
    1 - marker,
    2,
    curve
  )

  expect_equal(observed$estimate, 0.75, tolerance = 0)
  expect_equal(shuffled$estimate, observed$estimate, tolerance = 0)
  expect_equal(transformed$estimate, observed$estimate, tolerance = 0)
  expect_equal(reversed$estimate, 1 - observed$estimate, tolerance = 0)
})

test_that("CV time AUC sums only within-fold weighted pair quantities", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:5),
    fold = c(1L, 1L, 2L, 2L, 2L),
    time = c(1, 3, 1, 3, 4),
    event = rep(1L, 5L),
    survival_probability = 1 - c(0.2, 0.1, 0.9, 0.8, 0.7),
    survival_probability_status = rep("success", 5L),
    status = rep("success", 5L),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(
    list(fold = 1L, time = 1:4, event = rep(1L, 4L)),
    list(fold = 2L, time = 1:4, event = rep(1L, 4L))
  )
  curves <- lapply(analysis_inputs, function(input) {
    make_sb_time_auc_curve(input$time, input$event)
  })
  names(curves) <- c("1", "2")

  observed <- sb_survival_time_auc_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 2,
    censoring_curves = curves
  )

  expect_identical(observed$status, "success")
  expect_equal(observed$estimate, 1, tolerance = 0)
  expect_equal(observed$weighted_pair_credit, 3, tolerance = 0)
  expect_equal(observed$pair_denominator, 3, tolerance = 0)
  expect_identical(observed$n_folds, 2L)
  expect_identical(observed$n_evaluable_folds, 2L)
  expect_equal(observed$per_fold$pair_denominator, c(1, 2), tolerance = 0)

  marker <- 1 - predictions$survival_probability
  pooled_credit <- outer(
    marker[predictions$time <= 2],
    marker[predictions$time > 2],
    function(case, control) (case > control) + 0.5 * (case == control)
  )
  expect_equal(mean(pooled_credit), 4 / 6, tolerance = 0)
  expect_gt(abs(observed$estimate - mean(pooled_credit)), 0)
  expect_equal(
    observed$pair_denominator,
    sum(observed$per_fold$pair_denominator),
    tolerance = 0
  )
})

test_that("CV time AUC requires exact analysis and censoring fold mapping", {
  predictions <- data.frame(
    patient_id = c("case", "control"),
    fold = rep(1L, 2L),
    time = c(1, 3),
    event = rep(1L, 2L),
    survival_probability = c(0.2, 0.8),
    survival_probability_status = rep("success", 2L),
    status = rep("success", 2L),
    stringsAsFactors = FALSE
  )
  input <- list(fold = 1L, time = 1:3, event = rep(1L, 3L))
  curve <- make_sb_time_auc_curve(input$time, input$event)
  duplicated_input <- sb_survival_time_auc_metric(
    predictions,
    list(input, input),
    2,
    list(`1` = curve)
  )
  missing_curve <- sb_survival_time_auc_metric(
    predictions,
    list(input),
    2,
    list(`2` = curve)
  )
  duplicate_curve <- sb_survival_time_auc_metric(
    predictions,
    list(input),
    2,
    structure(list(curve, curve), names = c("1", "1"))
  )

  expect_identical(duplicated_input$status, "incomplete")
  expect_identical(
    duplicated_input$per_fold$error_class,
    "stormbreaker_bad_time_auc_analysis_input"
  )
  expect_identical(missing_curve$status, "incomplete")
  expect_identical(
    missing_curve$per_fold$error_class,
    "stormbreaker_missing_censoring_curve"
  )
  expect_identical(duplicate_curve$status, "incomplete")
  expect_identical(
    duplicate_curve$per_fold$error_class,
    "stormbreaker_bad_time_auc_censoring_curve"
  )
})

test_that("an incomplete time AUC retains successful fold evidence", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:4),
    fold = c(1L, 1L, 2L, 2L),
    time = c(1, 3, 1, 1.5),
    event = rep(1L, 4L),
    survival_probability = c(0.2, 0.8, 0.2, 0.8),
    survival_probability_status = rep("success", 4L),
    status = rep("success", 4L),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(
    list(fold = 1L, time = 1:3, event = rep(1L, 3L)),
    list(fold = 2L, time = 1:3, event = rep(1L, 3L))
  )
  curves <- lapply(analysis_inputs, function(input) {
    make_sb_time_auc_curve(input$time, input$event)
  })
  names(curves) <- c("1", "2")

  observed <- sb_survival_time_auc_metric(
    predictions,
    analysis_inputs,
    2,
    curves
  )

  expect_identical(observed$status, "incomplete")
  expect_true(is.na(observed$estimate))
  expect_identical(observed$n_evaluable_folds, 1L)
  expect_equal(observed$weighted_pair_credit, 1, tolerance = 0)
  expect_equal(observed$pair_denominator, 1, tolerance = 0)
  expect_identical(observed$per_fold$status, c("success", "not_evaluable"))
})
