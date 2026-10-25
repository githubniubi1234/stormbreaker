test_that("event KM is right-continuous and ignores censor-only jumps", {
  time <- c(1, 2, 2, 3, 4)
  event <- c(0L, 1L, 0L, 0L, 1L)
  curve <- sb_fit_event_km(time, event)
  independent <- survival::survfit(
    survival::Surv(time, event) ~ 1,
    se.fit = FALSE
  )

  expect_identical(
    names(curve),
    c(
      "time", "n_risk", "n_event", "n_censoring",
      "survival_before", "survival_after"
    )
  )
  expect_identical(curve$n_risk, c(5L, 4L, 2L, 1L))
  expect_identical(curve$n_event, c(0L, 1L, 0L, 1L))
  expect_identical(curve$n_censoring, c(1L, 1L, 1L, 0L))
  expect_equal(curve$survival_before, c(1, 1, 0.75, 0.75))
  expect_equal(curve$survival_after, c(1, 0.75, 0.75, 0))
  expect_equal(
    curve$survival_after,
    summary(independent, times = curve$time, extend = TRUE)$surv,
    tolerance = 1e-15
  )
  expect_equal(sb_event_survival_at(curve, 1), 1)
  expect_equal(sb_event_survival_at(curve, 2), 0.75)
  expect_equal(sb_event_survival_at(curve, 3), 0.75)
  expect_equal(sb_event_survival_at(curve, 4), 0)
})

test_that("IPCW Brier reduces to patient-level squared error without censoring", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:5),
    fold = c(1L, 1L, 2L, 2L, 2L),
    time = c(1, 3, 1, 3, 4),
    event = rep(1L, 5L),
    survival_probability = c(0.8, 0.4, 0.7, 0.5, 0.2),
    survival_probability_status = rep("success", 5L),
    status = rep("success", 5L),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(
    list(fold = 1L, time = c(1, 2, 3, 4), event = rep(1L, 4L)),
    list(fold = 2L, time = c(1, 2, 3, 4), event = rep(1L, 4L))
  )

  observed <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 2
  )
  expected_loss <- c(0.8^2, (1 - 0.4)^2, 0.7^2,
                     (1 - 0.5)^2, (1 - 0.2)^2)

  expect_identical(observed$status, "success")
  expect_identical(observed$n, 5L)
  expect_identical(observed$n_scored, 5L)
  expect_equal(observed$weight_sum, 5, tolerance = 0)
  expect_equal(
    observed$weighted_loss_sum,
    sum(expected_loss),
    tolerance = 1e-15
  )
  expect_equal(
    observed$estimate,
    sum(expected_loss) / 5,
    tolerance = 1e-15
  )
  expect_equal(
    observed$per_fold$estimate,
    c(sum(expected_loss[1:2]) / 2, sum(expected_loss[3:5]) / 3),
    tolerance = 0
  )
  expect_gt(
    abs(observed$estimate - mean(observed$per_fold$estimate)),
    0
  )
  expect_identical(observed$reference$status, "success")
  expect_equal(
    observed$reference$per_fold$reference_survival_horizon,
    c(0.5, 0.5),
    tolerance = 0
  )
  expect_equal(observed$reference$weighted_loss_sum, 1.25, tolerance = 0)
  expect_equal(observed$reference$estimate, 0.25, tolerance = 0)
  expect_identical(observed$skill_score$status, "success")
  expect_equal(
    observed$skill_score$estimate,
    1 - sum(expected_loss) / 1.25,
    tolerance = 1e-15
  )
})

test_that("survival Brier reference and skill aggregate patient losses", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:5),
    fold = c(1L, 1L, 2L, 2L, 2L),
    time = c(1, 3, 1, 3, 4),
    event = rep(1L, 5L),
    survival_probability = c(0.8, 0.4, 0.7, 0.5, 0.2),
    survival_probability_status = rep("success", 5L),
    status = rep("success", 5L),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(
    list(fold = 1L, time = c(1, 2, 3, 4), event = rep(1L, 4L)),
    list(
      fold = 2L,
      time = c(1, 2, 3, 4),
      event = c(1L, 0L, 1L, 1L)
    )
  )
  model_loss <- c(0.8^2, (1 - 0.4)^2, 0.7^2,
                  (1 - 0.5)^2, (1 - 0.2)^2)
  reference_loss <- c(0.5^2, (1 - 0.5)^2, 0.75^2,
                      (1 - 0.75)^2, (1 - 0.75)^2)

  observed <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 2
  )
  reference <- observed$reference

  expect_identical(reference$type, "outer_analysis_kaplan_meier")
  expect_identical(reference$source, "outer_analysis_event_only")
  expect_identical(reference$status, "success")
  expect_identical(reference$n, 5L)
  expect_identical(reference$n_scored, 5L)
  expect_equal(
    reference$per_fold$reference_survival_horizon,
    c(0.5, 0.75),
    tolerance = 0
  )
  expect_equal(
    reference$per_fold$weighted_loss_sum,
    c(sum(reference_loss[1:2]), sum(reference_loss[3:5])),
    tolerance = 0
  )
  expect_equal(reference$weighted_loss_sum, sum(reference_loss), tolerance = 0)
  expect_equal(reference$estimate, sum(reference_loss) / 5, tolerance = 0)
  expect_gt(
    abs(reference$estimate - mean(reference$per_fold$estimate)),
    0
  )
  expect_identical(observed$skill_score$metric, "ipcw_brier_skill_score")
  expect_identical(observed$skill_score$status, "success")
  expect_equal(
    observed$skill_score$estimate,
    1 - sum(model_loss) / sum(reference_loss),
    tolerance = 1e-15
  )
})

test_that("IPCW Brier uses analysis-fold G(T-) and G(t)", {
  predictions <- data.frame(
    patient_id = c("event", "survivor", "censored"),
    fold = rep(1L, 3L),
    time = c(2, 4, 2.5),
    event = c(1L, 1L, 0L),
    survival_probability = c(0.6, 0.7, 0.5),
    survival_probability_status = rep("success", 3L),
    status = rep("success", 3L),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(list(
    fold = 1L,
    time = c(1, 2, 2, 3, 4),
    event = c(0L, 1L, 0L, 0L, 1L)
  ))

  observed <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 3
  )
  curve <- observed$censoring_curves[[1L]]
  expected_event_loss <- (1 / 0.8) * 0.6^2
  expected_survivor_loss <- (1 / (0.8 * 2 / 3)) * (1 - 0.7)^2

  expect_identical(
    names(curve),
    c(
      "time", "n_risk", "n_event", "n_censoring",
      "n_risk_for_censoring", "survival_before", "survival_after"
    )
  )
  expect_equal(curve$survival_before, c(1, 0.8, 0.8 * 2 / 3, 0.8 / 3))
  expect_equal(curve$survival_after, c(0.8, 0.8 * 2 / 3, 0.8 / 3, 0.8 / 3))
  expect_identical(observed$per_fold$n_events_by_horizon, 1L)
  expect_identical(observed$per_fold$n_surviving_past_horizon, 1L)
  expect_identical(observed$per_fold$n_censored_by_horizon, 1L)
  expect_equal(
    observed$per_fold$censoring_survival_horizon,
    0.8 * 2 / 3
  )
  expect_equal(observed$weight_sum, 1 / 0.8 + 1 / (0.8 * 2 / 3))
  expect_equal(
    observed$weighted_loss_sum,
    expected_event_loss + expected_survivor_loss
  )
  expect_equal(
    observed$estimate,
    (expected_event_loss + expected_survivor_loss) / 3
  )
})

test_that("reverse KM orders tied events before censoring", {
  analysis_time <- c(1, 2, 2, 3, 4)
  analysis_event <- c(0L, 1L, 0L, 0L, 1L)

  curve <- sb_fit_censoring_km(analysis_time, analysis_event)
  censor_rows <- curve$n_censoring > 0L
  shifted_time <- analysis_time + ifelse(analysis_event == 0L, 0.25, 0)
  independent <- survival::survfit(
    survival::Surv(shifted_time, 1L - analysis_event) ~ 1,
    se.fit = FALSE
  )
  independent_after <- summary(
    independent,
    times = c(1.25, 2.25, 3.25),
    extend = TRUE
  )$surv

  expect_identical(curve$n_event, c(0L, 1L, 0L, 1L))
  expect_identical(curve$n_risk_for_censoring, c(5L, 3L, 2L, 0L))
  expect_equal(
    curve$survival_after[censor_rows],
    independent_after,
    tolerance = 1e-15
  )
  expect_equal(
    sb_censoring_survival_at(curve, 2, left_limit = TRUE),
    0.8,
    tolerance = 1e-15
  )
  expect_equal(
    sb_censoring_survival_at(curve, 3, left_limit = TRUE),
    0.8 * 2 / 3,
    tolerance = 1e-15
  )
})

test_that("censoring at the horizon has a known survivor outcome", {
  predictions <- data.frame(
    patient_id = c("event_at", "censored_at", "survivor"),
    fold = rep(1L, 3L),
    time = c(3, 3, 4),
    event = c(1L, 0L, 1L),
    survival_probability = c(0.6, 0.7, 0.8),
    survival_probability_status = rep("success", 3L),
    status = rep("success", 3L),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(list(
    fold = 1L,
    time = c(1, 2, 2, 3, 4),
    event = c(0L, 1L, 0L, 0L, 1L)
  ))

  observed <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 3
  )
  censoring_survival <- 0.8 * 2 / 3
  expected_loss <- (
    0.6^2 + (1 - 0.7)^2 + (1 - 0.8)^2
  ) / censoring_survival

  expect_identical(observed$status, "success")
  expect_identical(observed$n_scored, 3L)
  expect_equal(observed$weight_sum, 3 / censoring_survival)
  expect_equal(observed$weighted_loss_sum, expected_loss)
  expect_equal(observed$estimate, expected_loss / 3)
})

test_that("left-continuous G retains support at a censoring horizon", {
  predictions <- data.frame(
    patient_id = c("P1", "P2"),
    fold = c(1L, 1L),
    time = c(1, 3),
    event = c(1L, 1L),
    survival_probability = c(0.6, 0.7),
    survival_probability_status = c("success", "success"),
    status = c("success", "success"),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(list(
    fold = 1L,
    time = c(1, 2, 3),
    event = c(1L, 0L, 0L)
  ))

  observed <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 3
  )

  expect_identical(observed$status, "success")
  expect_equal(observed$per_fold$censoring_survival_horizon, 0.5)
  expect_equal(observed$weight_sum, 3)
  expect_equal(observed$weighted_loss_sum, 0.6^2 + 2 * 0.7^2)
  expect_equal(observed$estimate, (0.6^2 + 2 * 0.7^2) / 2)
  expect_length(observed$censoring_curves, 1L)
})

test_that("invalid censoring survival withholds model, reference, and skill", {
  predictions <- data.frame(
    patient_id = c("P1", "P2"),
    fold = rep(1L, 2L),
    time = c(1, 3),
    event = c(1L, 1L),
    survival_probability = c(0.6, 0.7),
    survival_probability_status = rep("success", 2L),
    status = rep("success", 2L),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(list(
    fold = 1L,
    time = c(1, 2),
    event = c(1L, 0L)
  ))
  forced_survival <- 0
  local_mocked_bindings(
    sb_fit_censoring_km = function(...) {
      data.frame(
        time = c(1, 2),
        n_risk = c(2L, 1L),
        n_event = c(1L, 0L),
        n_censoring = c(0L, 1L),
        n_risk_for_censoring = c(1L, 1L),
        survival_before = c(1, forced_survival),
        survival_after = c(forced_survival, forced_survival)
      )
    }
  )

  zero <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 2
  )
  forced_survival <- Inf
  non_finite <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 2
  )

  for (observed in list(zero, non_finite)) {
    expect_identical(observed$status, "incomplete")
    expect_identical(observed$per_fold$status, "censoring_failed")
    expect_identical(observed$reference$status, "incomplete")
    expect_identical(
      observed$reference$per_fold$status,
      "censoring_failed"
    )
    expect_identical(observed$skill_score$status, "incomplete")
    expect_identical(is.na(observed$skill_score$estimate), TRUE)
  }
})

test_that("failed probability folds retain partial IPCW audit evidence", {
  predictions <- data.frame(
    patient_id = c("P1", "P2", "P3", "P4"),
    fold = c(1L, 1L, 2L, 2L),
    time = c(1, 3, 1, 3),
    event = rep(1L, 4L),
    survival_probability = c(0.6, 0.7, NA, NA),
    survival_probability_status = c(
      "success", "success", "failed", "failed"
    ),
    status = rep("success", 4L),
    stringsAsFactors = FALSE
  )
  analysis_inputs <- list(
    list(fold = 1L, time = c(1, 2, 3), event = c(1L, 0L, 1L)),
    list(fold = 2L, time = c(1, 2, 3), event = c(1L, 0L, 1L))
  )

  observed <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    prediction_horizon = 2
  )

  expect_identical(observed$status, "incomplete")
  expect_true(is.na(observed$estimate))
  expect_identical(
    observed$per_fold$status,
    c("success", "prediction_failed")
  )
  expect_gt(observed$weighted_loss_sum, 0)
  expect_identical(observed$reference$status, "success")
  expect_identical(observed$reference$n_scored, 4L)
  expect_identical(observed$skill_score$status, "incomplete")
  expect_identical(is.na(observed$skill_score$estimate), TRUE)
})

test_that("boundary reference probabilities are valid but zero loss has no skill", {
  survivor_predictions <- data.frame(
    patient_id = c("P1", "P2"),
    fold = rep(1L, 2L),
    time = c(3, 4),
    event = c(1L, 0L),
    survival_probability = c(0.8, 0.7),
    survival_probability_status = rep("success", 2L),
    status = rep("success", 2L),
    stringsAsFactors = FALSE
  )
  survivor_analysis <- list(list(
    fold = 1L,
    time = c(3, 4),
    event = c(1L, 0L)
  ))

  reference_one <- sb_survival_ipcw_brier_metric(
    survivor_predictions,
    survivor_analysis,
    prediction_horizon = 2
  )

  expect_identical(reference_one$reference$status, "success")
  expect_identical(
    reference_one$reference$per_fold$reference_survival_horizon,
    1
  )
  expect_identical(reference_one$reference$weighted_loss_sum, 0)
  expect_identical(reference_one$reference$estimate, 0)
  expect_identical(reference_one$skill_score$status, "incomplete")
  expect_identical(is.na(reference_one$skill_score$estimate), TRUE)

  event_predictions <- survivor_predictions[1L, , drop = FALSE]
  event_predictions$time <- 3
  event_predictions$event <- 1L
  event_predictions$survival_probability <- 0.5
  event_analysis <- list(list(
    fold = 1L,
    time = c(1, 2),
    event = c(1L, 1L)
  ))
  reference_zero <- sb_survival_ipcw_brier_metric(
    event_predictions,
    event_analysis,
    prediction_horizon = 2
  )

  expect_identical(reference_zero$reference$status, "success")
  expect_identical(
    reference_zero$reference$per_fold$reference_survival_horizon,
    0
  )
  expect_identical(reference_zero$reference$estimate, 1)
  expect_identical(reference_zero$skill_score$status, "success")
  expect_identical(reference_zero$skill_score$estimate, 0.75)
})
