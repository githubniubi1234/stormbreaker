uno_test_curve <- function(time = c(1, 2, 3, 4, 5),
                           survival_after = c(1, 0.8, 0.8, 0.4, 0.4)) {
  data.frame(
    time = time,
    n_risk = rep.int(1L, length(time)),
    n_event = rep.int(0L, length(time)),
    n_censoring = rep.int(0L, length(time)),
    n_risk_for_censoring = rep.int(1L, length(time)),
    survival_before = c(1, head(survival_after, -1L)),
    survival_after = survival_after,
    stringsAsFactors = FALSE
  )
}

literal_uno_pairs <- function(time, event, risk, curve, horizon) {
  numerator <- denominator <- concordant <- tied_risk <- 0
  comparable_pairs <- tied_time_pairs <- 0L
  eligible <- which(event == 1L & time < horizon)
  for (index in eligible) {
    g_left <- if (any(curve$time < time[[index]])) {
      utils::tail(curve$survival_after[curve$time < time[[index]]], 1L)
    } else {
      1
    }
    weight <- 1 / g_left^2
    comparators <- which(
      time > time[[index]] |
        (time == time[[index]] & event == 0L)
    )
    comparable_pairs <- comparable_pairs + length(comparators)
    tied_time_pairs <- tied_time_pairs + sum(
      time[comparators] == time[[index]]
    )
    concordant <- concordant + weight * sum(
      risk[[index]] > risk[comparators]
    )
    tied_risk <- tied_risk + weight * sum(
      risk[[index]] == risk[comparators]
    )
    denominator <- denominator + weight * length(comparators)
  }
  numerator <- concordant + 0.5 * tied_risk
  list(
    estimate = numerator / denominator,
    numerator = numerator,
    denominator = denominator,
    concordant = concordant,
    tied_risk = tied_risk,
    comparable_pairs = as.integer(comparable_pairs),
    tied_time_pairs = as.integer(tied_time_pairs)
  )
}

test_that("Uno C uses squared left-limit IPCW and literal pair credit", {
  curve <- uno_test_curve()
  time <- c(1.5, 2.5, 3.5, 3.5, 5.5)
  event <- c(1L, 1L, 1L, 0L, 0L)
  risk <- c(3, 1, 2, 2, 0)
  expected <- literal_uno_pairs(time, event, risk, curve, 5)

  result <- sb_uno_c_stratum(
    time,
    event,
    risk,
    rep(TRUE, length(time)),
    curve,
    prediction_horizon = 5,
    fold = 1L
  )

  expect_identical(result$status, "success")
  expect_equal(result$weighted_concordant_pairs, 7.125)
  expect_equal(result$weighted_tied_risk_pairs, 1.5625)
  expect_equal(result$weighted_pair_credit, 7.90625)
  expect_equal(result$weighted_pair_denominator, 11.8125)
  expect_equal(result$estimate, expected$estimate)
  expect_identical(result$n_comparable_pairs, 9L)
  expect_identical(result$n_tied_time_comparable_pairs, 1L)
  expect_equal(
    result$estimate,
    result$weighted_pair_credit / result$weighted_pair_denominator
  )
  expect_false(isTRUE(all.equal(result$estimate, 7.25 / 10.25)))
})

test_that("Uno C locks strict horizon and observed-time tie conventions", {
  curve <- uno_test_curve(c(0.5, 1.5, 2.5), c(1, 1, 1))
  boundary <- sb_uno_c_stratum(
    time = c(1.999999, 2, 3),
    event = c(1L, 1L, 0L),
    risk = c(2, 100, 0),
    prediction_ok = rep(TRUE, 3L),
    censoring_curve = curve,
    prediction_horizon = 2,
    fold = 1L
  )
  expect_identical(boundary$status, "success")
  expect_identical(boundary$n_eligible_events, 1L)
  expect_identical(boundary$n_events_at_horizon_excluded, 1L)
  expect_equal(boundary$estimate, 0.5)

  tied_time <- sb_uno_c_stratum(
    time = c(1, 1, 1, 2),
    event = c(1L, 0L, 1L, 0L),
    risk = c(3, 1, 100, 0),
    prediction_ok = rep(TRUE, 4L),
    censoring_curve = curve,
    prediction_horizon = 2.5,
    fold = 1L
  )
  expect_identical(tied_time$status, "success")
  expect_identical(tied_time$n_comparable_pairs, 4L)
  expect_identical(tied_time$n_tied_time_comparable_pairs, 2L)
  expect_equal(tied_time$estimate, 1)
})

test_that("Uno C gives half credit only to exactly equal risks", {
  result <- sb_uno_c_stratum(
    time = c(1, 2, 3),
    event = c(1L, 0L, 0L),
    risk = c(1, 1, 1 + 1e-12),
    prediction_ok = rep(TRUE, 3L),
    censoring_curve = uno_test_curve(c(0.5, 2.5), c(1, 1)),
    prediction_horizon = 4,
    fold = 1L
  )
  expect_identical(result$status, "success")
  expect_equal(result$weighted_tied_risk_pairs, 1)
  expect_equal(result$estimate, 0.25)
})

test_that("Uno C reports unavailable inputs without fallback", {
  curve <- uno_test_curve(c(0.5, 1.5), c(1, 1))
  no_case <- sb_uno_c_stratum(
    c(1, 2), c(0L, 0L), c(2, 1), c(TRUE, TRUE), curve, 3, 1L
  )
  expect_identical(no_case$status, "not_evaluable")
  expect_identical(no_case$error_class, "stormbreaker_no_uno_c_case")
  expect_true(is.na(no_case$estimate))

  no_pair <- sb_uno_c_stratum(
    c(2, 1), c(1L, 0L), c(2, 1), c(TRUE, TRUE), curve, 3, 1L
  )
  expect_identical(no_pair$status, "not_evaluable")
  expect_identical(
    no_pair$error_class,
    "stormbreaker_no_uno_c_comparable_pair"
  )

  bad_risk <- sb_uno_c_stratum(
    c(1, 2), c(1L, 0L), c(NA_real_, 1), c(TRUE, TRUE), curve, 3, 1L
  )
  expect_identical(bad_risk$status, "prediction_failed")
  expect_identical(
    bad_risk$error_class,
    "stormbreaker_invalid_uno_c_risk"
  )

  bad_status <- sb_uno_c_stratum(
    c(1, 2), c(1L, 0L), c(2, 1), c(TRUE, FALSE), curve, 3, 1L
  )
  expect_identical(bad_status$status, "prediction_failed")

  zero_horizon_support <- sb_uno_c_stratum(
    c(0.5, 3),
    c(1L, 0L),
    c(2, 1),
    c(TRUE, TRUE),
    uno_test_curve(c(1, 2), c(0.5, 0)),
    3,
    1L
  )
  expect_identical(zero_horizon_support$status, "censoring_failed")
  expect_identical(
    zero_horizon_support$error_class,
    "stormbreaker_zero_uno_c_horizon_censoring_survival"
  )

  missing_curve <- sb_uno_c_stratum(
    c(1, 2),
    c(1L, 0L),
    c(2, 1),
    c(TRUE, TRUE),
    sb_empty_censoring_curve(),
    3,
    1L
  )
  expect_identical(missing_curve$status, "censoring_failed")
  expect_identical(
    missing_curve$error_class,
    "stormbreaker_missing_uno_c_censoring_curve"
  )
})

test_that("Uno C is invariant to row order and strict risk transformations", {
  curve <- uno_test_curve()
  time <- c(1.5, 2.5, 3.5, 5.5)
  event <- c(1L, 1L, 1L, 0L)
  risk <- c(3, 1, 2, 0)
  base <- sb_uno_c_stratum(
    time, event, risk, rep(TRUE, 4L), curve, 5, 1L
  )
  order <- c(4, 2, 1, 3)
  reordered <- sb_uno_c_stratum(
    time[order], event[order], risk[order], rep(TRUE, 4L), curve, 5, 1L
  )
  transformed <- sb_uno_c_stratum(
    time, event, exp(risk), rep(TRUE, 4L), curve, 5, 1L
  )
  reversed <- sb_uno_c_stratum(
    time, event, -risk, rep(TRUE, 4L), curve, 5, 1L
  )

  expect_equal(reordered$estimate, base$estimate)
  expect_equal(transformed$estimate, base$estimate)
  expect_equal(reversed$estimate, 1 - base$estimate)
})

test_that("Uno C aggregates weighted pair totals without cross-fold ranking", {
  curve <- uno_test_curve(c(0.25, 6), c(1, 1))
  predictions <- data.frame(
    patient_id = paste0("p", 1:6),
    fold = c(1L, 1L, 2L, 2L, 2L, 2L),
    time = c(0.5, 4, 2, 3, 4, 5),
    event = c(1L, 0L, 1L, 0L, 0L, 0L),
    risk = c(10, 0, 0, 10, 20, 30),
    status = "success",
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(
    list(fold = 1L, time = c(1, 2), event = c(1L, 1L)),
    list(fold = 2L, time = c(1, 2), event = c(1L, 1L))
  )
  curves <- list(`1` = curve, `2` = curve)
  result <- sb_survival_uno_c_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 6,
    censoring_curves = curves
  )

  expect_identical(result$status, "success")
  expect_equal(result$weighted_pair_credit, 1)
  expect_equal(result$weighted_pair_denominator, 4)
  expect_equal(result$estimate, 0.25)
  expect_false(isTRUE(all.equal(result$estimate, mean(c(1, 0)))))

  pooled <- literal_uno_pairs(
    predictions$time,
    predictions$event,
    predictions$risk,
    curve,
    6
  )
  expect_false(isTRUE(all.equal(result$estimate, pooled$estimate)))

  shifted <- predictions
  shifted$risk[shifted$fold == 2L] <-
    shifted$risk[shifted$fold == 2L] + 1000
  shifted_result <- sb_survival_uno_c_metric(
    shifted,
    analysis_inputs,
    prediction_horizon = 6,
    censoring_curves = curves
  )
  expect_equal(shifted_result$estimate, result$estimate)
})

test_that("Uno C requires exactly one analysis input and curve per fold", {
  predictions <- data.frame(
    patient_id = c("p1", "p2"),
    fold = c(1L, 1L),
    time = c(1, 2),
    event = c(1L, 0L),
    risk = c(2, 1),
    status = "success",
    stringsAsFactors = FALSE
  )
  curve <- uno_test_curve(c(0.5, 3), c(1, 1))

  missing_input <- sb_survival_uno_c_metric(
    predictions, list(), 3, list(`1` = curve)
  )
  expect_identical(
    missing_input$per_fold$error_class,
    "stormbreaker_bad_uno_c_analysis_input"
  )

  duplicate_input <- sb_survival_uno_c_metric(
    predictions,
    list(
      list(fold = 1L, time = 1, event = 1L),
      list(fold = 1L, time = 2, event = 0L)
    ),
    3,
    list(`1` = curve)
  )
  expect_identical(
    duplicate_input$per_fold$error_class,
    "stormbreaker_bad_uno_c_analysis_input"
  )

  missing_curve <- sb_survival_uno_c_metric(
    predictions,
    list(list(fold = 1L, time = c(1, 2), event = c(1L, 0L))),
    3,
    list()
  )
  expect_identical(
    missing_curve$per_fold$error_class,
    "stormbreaker_missing_uno_c_censoring_curve"
  )

  duplicate_curves <- list(curve, curve)
  names(duplicate_curves) <- c("1", "1")
  duplicate_curve <- sb_survival_uno_c_metric(
    predictions,
    list(list(fold = 1L, time = c(1, 2), event = c(1L, 0L))),
    3,
    duplicate_curves
  )
  expect_identical(
    duplicate_curve$per_fold$error_class,
    "stormbreaker_bad_uno_c_censoring_curve"
  )
})

test_that("unexpected Uno errors become fold-specific metric evidence", {
  predictions <- data.frame(
    patient_id = c("p1", "p2"),
    fold = c(1L, 1L),
    time = c(1, 2),
    event = c(1L, 0L),
    risk = c(2, 1),
    status = "success",
    stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    sb_uno_c_stratum = function(...) {
      condition <- simpleError("injected fold Uno failure")
      class(condition) <- c("test_fold_uno_error", class(condition))
      stop(condition)
    }
  )

  result <- sb_survival_uno_c_metric(
    predictions,
    list(list(fold = 1L, time = c(1, 2), event = c(1L, 0L))),
    3,
    list(`1` = uno_test_curve(c(0.5, 3), c(1, 1)))
  )

  expect_identical(result$status, "incomplete")
  expect_true(is.na(result$estimate))
  expect_identical(result$per_fold$status, "metric_failed")
  expect_identical(result$per_fold$error_class, "test_fold_uno_error")
  expect_identical(
    result$per_fold$error_message,
    "injected fold Uno failure"
  )
})
