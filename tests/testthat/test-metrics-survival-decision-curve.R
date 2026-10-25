make_sb_survival_dca_predictions <- function(
    time,
    event,
    event_probability,
    status = rep("success", length(time)),
    probability_status = rep("success", length(time))) {
  data.frame(
    patient_id = paste0("s", seq_along(time)),
    time = time,
    event = as.integer(event),
    survival_probability = 1 - event_probability,
    status = status,
    survival_probability_status = probability_status,
    stringsAsFactors = FALSE
  )
}

test_that("survival DCA reconstructs the published Kaplan-Meier formula", {
  predictions <- make_sb_survival_dca_predictions(
    time = 1:6,
    event = c(1L, 0L, 1L, 0L, 1L, 0L),
    event_probability = c(0.80, 0.70, 0.60, 0.40, 0.55, 0.10)
  )

  metric <- sb_survival_decision_curve_metric(
    predictions,
    decision_thresholds = 0.50,
    prediction_horizon = 4,
    time_unit = "month",
    scope = "test"
  )
  row <- metric$curve[1L, ]

  expect_identical(metric$status, "success")
  expect_identical(row$n_selected, 4L)
  expect_equal(row$selected_fraction, 4 / 6, tolerance = 0)
  expect_equal(row$selected_survival_probability, 3 / 8, tolerance = 0)
  expect_equal(row$selected_event_probability, 5 / 8, tolerance = 0)
  expect_equal(row$estimated_true_positive_rate, 5 / 12, tolerance = 1e-15)
  expect_equal(row$estimated_false_positive_rate, 1 / 4, tolerance = 0)
  expect_equal(row$model_net_benefit, 1 / 6, tolerance = 1e-15)
  expect_equal(metric$overall_survival_probability, 5 / 8, tolerance = 0)
  expect_equal(row$treat_all_net_benefit, -1 / 4, tolerance = 0)
  expect_identical(row$treat_none_net_benefit, 0)
})

test_that("survival DCA applies exact action and empty-selection boundaries", {
  predictions <- make_sb_survival_dca_predictions(
    time = c(1, 2, 3, 4),
    event = c(1L, 0L, 1L, 0L),
    event_probability = c(0.80, 0.50, 0.30, 0.10)
  )
  metric <- sb_survival_decision_curve_metric(
    predictions,
    decision_thresholds = c(0.05, 0.50, 0.90),
    prediction_horizon = 3,
    time_unit = "month",
    scope = "test"
  )

  all_row <- metric$curve[metric$curve$threshold_probability == 0.05, ]
  tied_row <- metric$curve[metric$curve$threshold_probability == 0.50, ]
  empty_row <- metric$curve[metric$curve$threshold_probability == 0.90, ]
  expect_identical(all_row$n_selected, 4L)
  expect_equal(
    all_row$model_net_benefit,
    all_row$treat_all_net_benefit,
    tolerance = 0
  )
  expect_identical(tied_row$n_selected, 2L)
  expect_identical(empty_row$n_selected, 0L)
  expect_identical(empty_row$horizon_support, "empty_selection")
  expect_identical(empty_row$selected_survival_probability, 1)
  expect_identical(empty_row$selected_event_probability, 0)
  expect_identical(empty_row$model_net_benefit, 0)
})

test_that("event-censor ties process the event before censoring", {
  predictions <- make_sb_survival_dca_predictions(
    time = c(1, 1),
    event = c(1L, 0L),
    event_probability = c(0.80, 0.80)
  )
  metric <- sb_survival_decision_curve_metric(
    predictions,
    decision_thresholds = 0.50,
    prediction_horizon = 1,
    time_unit = "month",
    scope = "test"
  )

  expect_identical(metric$status, "success")
  expect_equal(
    metric$curve$selected_survival_probability,
    0.5,
    tolerance = 0
  )
  expect_equal(metric$overall_survival_probability, 0.5, tolerance = 0)
})

test_that("zero survival is estimable beyond the last selected follow-up", {
  predictions <- make_sb_survival_dca_predictions(
    time = c(1, 2, 4),
    event = c(1L, 1L, 0L),
    event_probability = c(0.80, 0.80, 0.10)
  )
  metric <- sb_survival_decision_curve_metric(
    predictions,
    decision_thresholds = 0.50,
    prediction_horizon = 3,
    time_unit = "month",
    scope = "test"
  )

  expect_identical(metric$status, "success")
  expect_identical(
    metric$curve$horizon_support,
    "survival_zero_before_horizon"
  )
  expect_identical(metric$curve$selected_survival_probability, 0)
  expect_equal(metric$curve$model_net_benefit, 2 / 3, tolerance = 0)
  expect_equal(
    metric$curve$treat_all_net_benefit,
    1 / 3,
    tolerance = 1e-15
  )
})

test_that("positive KM survival is not extrapolated past observed support", {
  predictions <- make_sb_survival_dca_predictions(
    time = c(1, 2, 4),
    event = c(1L, 0L, 0L),
    event_probability = c(0.80, 0.80, 0.10)
  )
  metric <- sb_survival_decision_curve_metric(
    predictions,
    decision_thresholds = 0.50,
    prediction_horizon = 3,
    time_unit = "month",
    scope = "test"
  )

  expect_identical(metric$status, "incomplete")
  expect_identical(metric$n_evaluable_thresholds, 0L)
  expect_identical(metric$curve$status, "not_evaluable")
  expect_identical(metric$curve$horizon_support, "insufficient_follow_up")
  expect_true(is.na(metric$curve$model_net_benefit))
  expect_true(is.finite(metric$curve$treat_all_net_benefit))
  expect_identical(
    metric$curve$error_class,
    "stormbreaker_insufficient_survival_dca_follow_up"
  )
})

test_that("the treat-all reference uses the same zero-versus-positive rule", {
  zero <- make_sb_survival_dca_predictions(
    time = c(1, 2),
    event = c(1L, 1L),
    event_probability = c(0.10, 0.10)
  )
  positive <- make_sb_survival_dca_predictions(
    time = c(1, 2),
    event = c(1L, 0L),
    event_probability = c(0.10, 0.10)
  )
  zero_metric <- sb_survival_decision_curve_metric(
    zero, 0.50, 3, "month", "test"
  )
  positive_metric <- sb_survival_decision_curve_metric(
    positive, 0.50, 3, "month", "test"
  )

  expect_identical(zero_metric$overall_horizon_support,
                   "survival_zero_before_horizon")
  expect_identical(zero_metric$status, "success")
  expect_identical(zero_metric$curve$model_net_benefit, 0)
  expect_identical(zero_metric$curve$treat_all_net_benefit, 1)
  expect_identical(positive_metric$overall_horizon_support,
                   "insufficient_follow_up")
  expect_identical(positive_metric$status, "incomplete")
  expect_true(is.na(positive_metric$curve$treat_all_net_benefit))
})

test_that("invalid or failed horizon predictions leave DCA incomplete", {
  predictions <- make_sb_survival_dca_predictions(
    time = c(1, 2, 3),
    event = c(1L, 0L, 1L),
    event_probability = c(0, 0.5, 1),
    probability_status = c("success", "failed", "success")
  )
  metric <- sb_survival_decision_curve_metric(
    predictions, 0.50, 2, "month", "test"
  )

  expect_identical(metric$status, "incomplete")
  expect_identical(metric$curve$status, "prediction_failed")
  expect_true(is.na(metric$curve$model_net_benefit))
  expect_identical(
    metric$curve$error_class,
    "stormbreaker_invalid_survival_prediction"
  )
})
