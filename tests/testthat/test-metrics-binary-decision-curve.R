test_that("binary decision curve matches hand-counted net benefit", {
  scored <- data.frame(
    outcome = c(1L, 1L, 0L, 0L),
    probability = c(0.90, 0.50, 0.70, 0.20),
    status = rep("success", 4L)
  )

  metric <- stormbreaker:::sb_binary_decision_curve_metric(
    scored,
    c(0.30, 0.50, 0.80),
    "hand_calculation"
  )

  expect_identical(metric$status, "success")
  expect_identical(
    metric$decision_rule,
    "probability_greater_than_or_equal_to_threshold"
  )
  expect_identical(metric$curve$true_positive, c(2L, 2L, 1L))
  expect_identical(metric$curve$false_positive, c(1L, 1L, 0L))
  expect_identical(metric$curve$true_negative, c(1L, 1L, 2L))
  expect_identical(metric$curve$false_negative, c(0L, 0L, 1L))
  expect_identical(metric$curve$n_selected, c(3L, 3L, 1L))
  expect_equal(metric$curve$sensitivity, c(1, 1, 0.5), tolerance = 0)
  expect_equal(metric$curve$specificity, c(0.5, 0.5, 1), tolerance = 0)
  expect_equal(
    metric$curve$positive_predictive_value,
    c(2 / 3, 2 / 3, 1),
    tolerance = 0
  )
  expect_equal(
    metric$curve$negative_predictive_value,
    c(1, 1, 2 / 3),
    tolerance = 0
  )
  expect_equal(
    metric$curve$model_net_benefit,
    c(0.5 - 0.25 * 0.3 / 0.7, 0.5 - 0.25, 0.25),
    tolerance = 1e-15
  )
  expect_equal(
    metric$curve$treat_all_net_benefit,
    0.5 - 0.5 * c(0.3 / 0.7, 1, 4),
    tolerance = 1e-15
  )
  expect_identical(
    metric$curve$treat_none_net_benefit,
    c(0, 0, 0)
  )
})

test_that("binary decision curve is invariant to patient order", {
  scored <- data.frame(
    outcome = c(1L, 0L, 1L, 0L),
    probability = c(0.8, 0.6, 0.4, 0.1),
    status = rep("success", 4L)
  )

  direct <- stormbreaker:::sb_binary_decision_curve_metric(
    scored,
    c(0.2, 0.5),
    "test"
  )
  reversed <- stormbreaker:::sb_binary_decision_curve_metric(
    scored[4:1, ],
    c(0.2, 0.5),
    "test"
  )

  expect_identical(reversed$curve, direct$curve)
})

test_that("binary decision curve never promotes a partial patient subset", {
  scored <- data.frame(
    outcome = c(1L, 0L, 1L),
    probability = c(0.8, NA_real_, 0.2),
    status = c("success", "failed", "success")
  )

  metric <- stormbreaker:::sb_binary_decision_curve_metric(
    scored,
    c(0.2, 0.5),
    "test"
  )

  expect_identical(metric$status, "incomplete")
  expect_identical(metric$n, 3L)
  expect_identical(metric$n_scored, 0L)
  expect_true(all(is.na(metric$curve$model_net_benefit)))
  expect_match(metric$error_message, "incomplete or invalid")
})

test_that("binary decision curve handles boundaries and one class", {
  scored <- data.frame(
    outcome = c(0L, 0L, 0L),
    probability = c(0, 0.5, 1),
    status = rep("success", 3L)
  )

  metric <- stormbreaker:::sb_binary_decision_curve_metric(
    scored,
    c(0.25, 0.75),
    "test"
  )

  expect_identical(metric$status, "success")
  expect_identical(metric$n_positive, 0L)
  expect_identical(metric$n_negative, 3L)
  expect_identical(metric$curve$true_positive, c(0L, 0L))
  expect_identical(metric$curve$false_positive, c(2L, 1L))
  expect_true(all(is.na(metric$curve$sensitivity)))
  expect_equal(metric$curve$specificity, c(1 / 3, 2 / 3), tolerance = 0)
  expect_identical(metric$curve$positive_predictive_value, c(0, 0))
  expect_identical(metric$curve$negative_predictive_value, c(1, 1))
  expect_equal(
    metric$curve$treat_all_net_benefit,
    -c(0.25 / 0.75, 0.75 / 0.25),
    tolerance = 1e-15
  )

  scored$outcome <- rep.int(1L, nrow(scored))
  positive_only <- stormbreaker:::sb_binary_decision_curve_metric(
    scored,
    c(0.25, 0.75),
    "test"
  )
  expect_equal(
    positive_only$curve$sensitivity,
    c(2 / 3, 1 / 3),
    tolerance = 0
  )
  expect_true(all(is.na(positive_only$curve$specificity)))
  expect_identical(positive_only$curve$positive_predictive_value, c(1, 1))
  expect_identical(positive_only$curve$negative_predictive_value, c(0, 0))
})

test_that("undefined predictive-value denominators stay explicit", {
  scored <- data.frame(
    outcome = c(1L, 0L),
    probability = c(0.9, 0.8),
    status = rep("success", 2L)
  )

  metric <- stormbreaker:::sb_binary_decision_curve_metric(
    scored,
    c(0.5, 0.95),
    "test"
  )

  expect_identical(metric$curve$sensitivity, c(1, 0))
  expect_identical(metric$curve$specificity, c(0, 1))
  expect_equal(
    metric$curve$positive_predictive_value,
    c(0.5, NA_real_),
    tolerance = 0
  )
  expect_equal(
    metric$curve$negative_predictive_value,
    c(NA_real_, 0.5),
    tolerance = 0
  )
})

test_that("binary decision curve exposes reconstructible curve columns", {
  scored <- data.frame(
    outcome = c(1L, 0L),
    probability = c(0.9, 0.1),
    status = rep("success", 2L)
  )

  metric <- stormbreaker:::sb_binary_decision_curve_metric(
    scored,
    0.5,
    "test"
  )

  expect_named(
    metric$curve,
    c(
      "threshold_probability",
      "threshold_odds",
      "n",
      "n_positive",
      "n_negative",
      "n_selected",
      "true_positive",
      "false_positive",
      "true_negative",
      "false_negative",
      "sensitivity",
      "specificity",
      "positive_predictive_value",
      "negative_predictive_value",
      "model_net_benefit",
      "treat_all_net_benefit",
      "treat_none_net_benefit"
    )
  )
  expect_identical(metric$n_thresholds, 1L)
  expect_identical(metric$threshold_source, "study_declared_binary_task")
})
