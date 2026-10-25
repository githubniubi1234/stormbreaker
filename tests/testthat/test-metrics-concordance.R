test_that("Harrell C matches a three-patient hand calculation", {
  predictions <- data.frame(
    fold = rep(1L, 3L),
    time = c(2, 5, 6),
    event = c(1L, 1L, 0L),
    risk = c(0.8, 0.3, 0.5),
    status = rep("success", 3L)
  )

  metric <- sb_harrell_c_metric(predictions)

  required_metric_fields <- c(
    "metric", "scope", "aggregation", "time_weight",
    "risk_direction", "status", "estimate", "n_folds",
    "n_evaluable_folds", "concordant", "discordant",
    "tied_risk", "tied_time", "tied_both", "comparable_pairs",
    "per_fold"
  )
  expect_length(setdiff(required_metric_fields, names(metric)), 0L)
  expect_identical(metric$metric, "harrell_c")
  expect_identical(metric$scope, "within_outer_assessment_fold")
  expect_identical(metric$aggregation, "summed_comparable_pairs")
  expect_identical(metric$time_weight, "n")
  expect_identical(
    metric$risk_direction,
    "higher_is_higher_hazard"
  )
  expect_identical(metric$status, "success")
  expect_equal(metric$estimate, 2 / 3, tolerance = 1e-15)
  expect_identical(metric$n_folds, 1L)
  expect_identical(metric$n_evaluable_folds, 1L)
  expect_equal(metric$concordant, 2, tolerance = 0)
  expect_equal(metric$discordant, 1, tolerance = 0)
  expect_equal(metric$tied_risk, 0, tolerance = 0)
  expect_equal(metric$comparable_pairs, 3, tolerance = 0)
  required_fold_fields <- c(
    "fold", "status", "n", "n_events", "estimate",
    "concordant", "discordant", "tied_risk", "tied_time",
    "tied_both", "comparable_pairs", "error_class",
    "error_message"
  )
  expect_length(setdiff(required_fold_fields, names(metric$per_fold)), 0L)
})

test_that("Harrell C stratum exposes survival concordance counts", {
  audit <- sb_harrell_c_stratum(
    time = c(2, 5, 6),
    event = c(1L, 1L, 0L),
    risk = c(0.8, 0.3, 0.5)
  )

  expect_identical(audit$status, "success")
  expect_equal(audit$estimate, 2 / 3, tolerance = 1e-15)
  expect_equal(audit$concordant, 2, tolerance = 0)
  expect_equal(audit$discordant, 1, tolerance = 0)
  expect_equal(audit$tied_risk, 0, tolerance = 0)
  expect_equal(audit$comparable_pairs, 3, tolerance = 0)
})

test_that("Harrell C stratum keeps tied-risk half credit", {
  audit <- sb_harrell_c_stratum(
    time = c(1, 2),
    event = c(1L, 0L),
    risk = c(0.5, 0.5)
  )

  expect_identical(audit$status, "success")
  expect_equal(audit$estimate, 0.5, tolerance = 0)
  expect_equal(audit$tied_risk, 1, tolerance = 0)
})

test_that("Harrell C stratum and fold wrapper distinguish all tie counts", {
  count_fields <- c(
    "concordant", "discordant", "tied_risk", "tied_time",
    "tied_both", "comparable_pairs"
  )
  cases <- list(
    list(risk = c(0.9, 0.8, 0.2), expected = c(2, 0, 0, 1, 0, 2)),
    list(risk = c(0.5, 0.5, 0.2), expected = c(2, 0, 0, 0, 1, 2))
  )

  for (case in cases) {
    audit <- sb_harrell_c_stratum(
      time = c(1, 1, 2),
      event = c(1L, 1L, 0L),
      risk = case$risk
    )
    predictions <- data.frame(
      fold = rep(1L, 3L),
      time = c(1, 1, 2),
      event = c(1L, 1L, 0L),
      risk = case$risk,
      status = rep("success", 3L)
    )
    metric <- sb_harrell_c_metric(predictions)

    expect_equal(
      unname(unlist(audit[count_fields], use.names = FALSE)),
      case$expected,
      tolerance = 0
    )
    expect_equal(
      unname(unlist(metric$per_fold[1L, count_fields], use.names = FALSE)),
      case$expected,
      tolerance = 0
    )
  }
})

test_that("Harrell C compares only patients predicted by the same fold model", {
  predictions <- data.frame(
    fold = rep(1:2, each = 3L),
    time = c(2, 5, 6, 1, 4, 7),
    event = c(1L, 1L, 0L, 1L, 0L, 1L),
    risk = c(0.8, 0.3, 0.5, 0.9, 0.4, 0.2),
    status = rep("success", 6L)
  )

  observed <- sb_harrell_c_metric(predictions)
  shifted <- predictions
  shifted$risk[shifted$fold == 2L] <-
    shifted$risk[shifted$fold == 2L] + 1000
  shifted_observed <- sb_harrell_c_metric(shifted)
  reversed <- predictions
  reversed$risk <- -reversed$risk
  reversed_observed <- sb_harrell_c_metric(reversed)
  official <- survival::concordance(
    survival::Surv(time, event) ~ risk + survival::strata(fold),
    data = predictions,
    reverse = TRUE,
    timewt = "n"
  )

  expect_equal(observed$estimate, 4 / 5, tolerance = 1e-15)
  expect_equal(shifted_observed$estimate, observed$estimate, tolerance = 0)
  expect_equal(reversed_observed$estimate, 1 - observed$estimate,
               tolerance = 1e-15)
  expect_equal(observed$estimate, unname(official$concordance),
               tolerance = 1e-15)
  expect_equal(observed$per_fold$estimate, c(2 / 3, 1),
               tolerance = 1e-15)
  expect_equal(observed$per_fold$comparable_pairs, c(3, 2),
               tolerance = 0)
  expect_equal(
    observed$estimate,
    sum(observed$per_fold$concordant) /
      sum(observed$per_fold$comparable_pairs),
    tolerance = 1e-15
  )
})

test_that("risk-score ties receive half credit", {
  predictions <- data.frame(
    fold = c(1L, 1L),
    time = c(1, 2),
    event = c(1L, 0L),
    risk = c(0.5, 0.5),
    status = c("success", "success")
  )

  metric <- sb_harrell_c_metric(predictions)

  expect_identical(metric$status, "success")
  expect_equal(metric$estimate, 0.5, tolerance = 0)
  expect_equal(metric$tied_risk, 1, tolerance = 0)
  expect_equal(metric$comparable_pairs, 1, tolerance = 0)
})

test_that("a fold without comparable pairs is visibly not evaluable", {
  predictions <- data.frame(
    fold = rep(1L, 3L),
    time = c(5, 1, 2),
    event = c(1L, 0L, 0L),
    risk = c(0.8, 0.3, 0.5),
    status = rep("success", 3L)
  )

  metric <- sb_harrell_c_metric(predictions)

  expect_identical(metric$status, "incomplete")
  expect_identical(metric$n_evaluable_folds, 0L)
  expect_identical(is.na(metric$estimate), TRUE)
  expect_identical(metric$per_fold$status, "not_evaluable")
  expect_match(metric$per_fold$error_message, "No comparable patient pairs")
})

test_that("failed OOF predictions cannot yield a complete concordance", {
  predictions <- data.frame(
    fold = rep(1:2, each = 3L),
    time = c(2, 5, 6, 1, 4, 7),
    event = c(1L, 1L, 0L, 1L, 0L, 1L),
    risk = c(0.8, 0.3, 0.5, NA, NA, NA),
    status = c(rep("success", 3L), rep("failed", 3L))
  )

  metric <- sb_harrell_c_metric(predictions)

  expect_identical(metric$status, "incomplete")
  expect_identical(is.na(metric$estimate), TRUE)
  expect_identical(metric$n_folds, 2L)
  expect_identical(metric$n_evaluable_folds, 1L)
  expect_identical(
    metric$per_fold$status,
    c("success", "prediction_failed")
  )
  expect_match(
    metric$per_fold$error_message[2],
    "failed or non-finite"
  )
})
