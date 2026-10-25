test_that("fold-stratified ROC AUC aggregates only within-model pairs", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:8),
    fold = rep(1:2, each = 4),
    outcome = c(1L, 1L, 0L, 0L, 1L, 1L, 0L, 0L),
    probability = c(0.9, 0.5, 0.4, 0.1, 0.6, 0.2, 0.6, 0.3),
    status = "success",
    stringsAsFactors = FALSE
  )

  metric <- sb_binary_auc_metric(predictions)

  expect_identical(metric$status, "success")
  expect_identical(metric$scope, "within_outer_assessment_fold")
  expect_identical(metric$comparable_pairs, 8)
  expect_identical(metric$correctly_ordered_pairs, 5)
  expect_identical(metric$tied_pairs, 1)
  expect_equal(metric$estimate, 5.5 / 8, tolerance = 0)
  expect_equal(metric$per_fold$estimate, c(1, 0.375), tolerance = 0)
})

test_that("failed binary predictions cannot produce a headline AUC", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:4),
    fold = rep(1:2, each = 2),
    outcome = c(1L, 0L, 1L, 0L),
    probability = c(0.8, 0.2, NA_real_, NA_real_),
    status = c("success", "success", "failed", "failed"),
    stringsAsFactors = FALSE
  )

  metric <- sb_binary_auc_metric(predictions)

  expect_identical(metric$status, "incomplete")
  expect_identical(is.na(metric$estimate), TRUE)
  expect_identical(
    metric$per_fold$status,
    c("success", "prediction_failed")
  )
})

test_that("binary AUC stratum exposes exact pair accounting", {
  audit <- sb_binary_auc_stratum(
    outcome = c(1L, 1L, 0L, 0L),
    probability = c(0.6, 0.2, 0.6, 0.3)
  )

  expect_identical(audit$status, "success")
  expect_identical(audit$n, 4L)
  expect_identical(audit$n_positive, 2L)
  expect_identical(audit$n_negative, 2L)
  expect_equal(audit$correctly_ordered_pairs, 1, tolerance = 0)
  expect_equal(audit$tied_pairs, 1, tolerance = 0)
  expect_equal(audit$comparable_pairs, 4, tolerance = 0)
  expect_equal(audit$estimate, 1.5 / 4, tolerance = 0)
})

test_that("binary AUC stratum is explicitly not evaluable with one class", {
  audit <- sb_binary_auc_stratum(
    outcome = c(0L, 0L),
    probability = c(0.1, 0.2)
  )

  expect_identical(audit$status, "not_evaluable")
  expect_identical(is.na(audit$estimate), TRUE)
  expect_equal(audit$comparable_pairs, 0, tolerance = 0)
})

test_that("average precision is tie-invariant and positive-weighted by fold", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:10),
    fold = c(rep(1L, 4L), rep(2L, 6L)),
    outcome = c(1L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L),
    probability = c(0.9, 0.8, 0.8, 0.1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4),
    status = "success",
    stringsAsFactors = FALSE
  )

  metric <- sb_binary_average_precision_metric(predictions)

  expect_identical(metric$metric, "average_precision")
  expect_identical(metric$scope, "within_outer_assessment_fold")
  expect_identical(
    metric$aggregation,
    "positive_weighted_mean_of_fold_average_precision"
  )
  expect_identical(metric$interpolation, "none")
  expect_identical(metric$status, "success")
  expect_identical(metric$n, 10L)
  expect_identical(metric$n_scored, 10L)
  expect_identical(metric$n_positive, 5L)
  expect_identical(metric$n_negative, 5L)
  expect_equal(
    metric$per_fold$estimate,
    c(5 / 6, 8 / 15),
    tolerance = 1e-15
  )
  expect_equal(
    metric$per_fold$positive_weighted_sum,
    c(5 / 3, 8 / 5),
    tolerance = 1e-15
  )
  expect_equal(metric$positive_weighted_sum, 49 / 15, tolerance = 1e-15)
  expect_equal(metric$estimate, 49 / 75, tolerance = 1e-15)
  expect_false(isTRUE(all.equal(
    metric$estimate,
    mean(metric$per_fold$estimate),
    tolerance = 0
  )))
  expect_length(metric$curves, 2L)
  expect_equal(metric$curves[[1L]]$threshold, c(0.9, 0.8, 0.1))
  expect_identical(metric$curves[[1L]]$n_at_threshold, c(1L, 2L, 1L))
  expect_identical(metric$curves[[1L]]$positive_at_threshold, c(1L, 1L, 0L))
  expect_identical(metric$curves[[1L]]$negative_at_threshold, c(0L, 1L, 1L))
  expect_equal(metric$curves[[1L]]$precision, c(1, 2 / 3, 1 / 2))
  expect_equal(metric$curves[[1L]]$recall, c(1 / 2, 1, 1))
  expect_equal(metric$curves[[1L]]$recall_increment, c(1 / 2, 1 / 2, 0))
  expect_equal(metric$curves[[1L]]$contribution, c(1 / 2, 1 / 3, 0))

  permuted <- predictions[c(4L, 3L, 2L, 1L, 10L, 6L, 9L, 5L, 8L, 7L), ]
  permuted_metric <- sb_binary_average_precision_metric(permuted)
  expect_identical(permuted_metric$per_fold, metric$per_fold)
  expect_identical(permuted_metric$curves, metric$curves)
  expect_identical(permuted_metric$estimate, metric$estimate)
})

test_that("average precision never ranks probabilities from different folds", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:4),
    fold = rep(1:2, each = 2L),
    outcome = rep(c(1L, 0L), 2L),
    probability = c(0.6, 0.5, 0.2, 0.1),
    status = "success",
    stringsAsFactors = FALSE
  )

  metric <- sb_binary_average_precision_metric(predictions)

  expect_identical(metric$status, "success")
  expect_equal(metric$per_fold$estimate, c(1, 1), tolerance = 0)
  expect_equal(metric$estimate, 1, tolerance = 0)

  pooled_curve <- sb_binary_precision_recall_curve(
    predictions$outcome,
    predictions$probability,
    fold_id = 1L
  )
  pooled_average_precision <- sum(pooled_curve$contribution)
  expect_equal(pooled_average_precision, 5 / 6, tolerance = 1e-15)
  expect_false(isTRUE(all.equal(
    metric$estimate,
    pooled_average_precision,
    tolerance = 0
  )))
})

test_that("incomplete folds cannot produce headline average precision", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:8),
    fold = rep(1:2, each = 4L),
    outcome = c(1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L),
    probability = c(0.9, 0.8, 0.4, 0.1, NA, NA, NA, NA),
    status = c(rep("success", 4L), rep("failed", 4L)),
    stringsAsFactors = FALSE
  )

  metric <- sb_binary_average_precision_metric(predictions)

  expect_identical(metric$status, "incomplete")
  expect_identical(is.na(metric$estimate), TRUE)
  expect_identical(metric$n, 8L)
  expect_identical(metric$n_scored, 4L)
  expect_identical(metric$n_evaluable_folds, 1L)
  expect_identical(
    metric$per_fold$status,
    c("success", "prediction_failed")
  )
  expect_equal(metric$positive_weighted_sum, 5 / 3, tolerance = 1e-15)
  expect_identical(nrow(metric$curves[[2L]]), 0L)

  out_of_range <- predictions
  out_of_range$status <- "success"
  out_of_range$probability[5:8] <- c(1.1, 0.7, 0.4, 0.2)
  invalid_metric <- sb_binary_average_precision_metric(out_of_range)
  expect_identical(invalid_metric$status, "incomplete")
  expect_identical(
    invalid_metric$per_fold$status,
    c("success", "prediction_failed")
  )

  one_class <- predictions
  one_class$status <- "success"
  one_class$probability[5:8] <- c(0.9, 0.7, 0.4, 0.2)
  one_class$outcome[5:8] <- 1L
  one_class_metric <- sb_binary_average_precision_metric(one_class)
  expect_identical(one_class_metric$status, "incomplete")
  expect_identical(
    one_class_metric$per_fold$status,
    c("success", "not_evaluable")
  )
})

test_that("OOF Brier score is patient-weighted across unequal folds", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:5),
    fold = c(1L, 1L, 2L, 2L, 2L),
    outcome = c(1L, 0L, 1L, 0L, 1L),
    probability = c(0.8, 0.3, 0.6, 0.4, 0.9),
    status = "success",
    stringsAsFactors = FALSE
  )
  squared_error <- (predictions$outcome - predictions$probability)^2

  metric <- sb_binary_brier_metric(predictions)

  expect_identical(metric$metric, "brier_score")
  expect_identical(metric$scope, "patient_level_out_of_fold")
  expect_identical(
    metric$aggregation,
    "summed_squared_error_over_patients"
  )
  expect_identical(metric$status, "success")
  expect_identical(metric$n, 5L)
  expect_identical(metric$n_scored, 5L)
  expect_identical(metric$n_folds, 2L)
  expect_identical(metric$n_evaluable_folds, 2L)
  expect_equal(metric$squared_error_sum, sum(squared_error), tolerance = 0)
  expect_equal(metric$estimate, mean(squared_error), tolerance = 0)
  expect_identical(metric$per_fold$n, c(2L, 3L))
  expect_equal(
    metric$per_fold$squared_error_sum,
    c(sum(squared_error[1:2]), sum(squared_error[3:5])),
    tolerance = 0
  )
  expect_equal(
    metric$per_fold$estimate,
    c(mean(squared_error[1:2]), mean(squared_error[3:5])),
    tolerance = 0
  )
  expect_false(isTRUE(all.equal(
    metric$estimate,
    mean(metric$per_fold$estimate),
    tolerance = 0
  )))
})

test_that("incomplete OOF probabilities cannot produce a headline Brier score", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:4),
    fold = rep(1:2, each = 2),
    outcome = c(1L, 0L, 1L, 0L),
    probability = c(0.8, 0.2, NA_real_, NA_real_),
    status = c("success", "success", "failed", "failed"),
    stringsAsFactors = FALSE
  )

  metric <- sb_binary_brier_metric(predictions)

  expect_identical(metric$status, "incomplete")
  expect_identical(is.na(metric$estimate), TRUE)
  expect_identical(metric$n, 4L)
  expect_identical(metric$n_scored, 2L)
  expect_equal(metric$squared_error_sum, 0.08, tolerance = 1e-15)
  expect_identical(
    metric$per_fold$status,
    c("success", "prediction_failed")
  )

  out_of_range <- predictions
  out_of_range$status <- "success"
  out_of_range$probability[3:4] <- c(1.1, 0.4)
  invalid_metric <- sb_binary_brier_metric(out_of_range)
  expect_identical(invalid_metric$status, "incomplete")
  expect_identical(is.na(invalid_metric$estimate), TRUE)
  expect_identical(
    invalid_metric$per_fold$status,
    c("success", "prediction_failed")
  )
})

test_that("fold-local prevalence reference yields patient-weighted Brier skill", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:5),
    fold = c(1L, 1L, 2L, 2L, 2L),
    outcome = c(1L, 0L, 1L, 0L, 1L),
    probability = c(0.8, 0.3, 0.6, 0.4, 0.9),
    status = "success",
    stringsAsFactors = FALSE
  )
  reference_input <- data.frame(
    fold = 1:2,
    analysis_prevalence = c(0.25, 0.75)
  )
  reference_probability <- c(0.25, 0.25, 0.75, 0.75, 0.75)
  model_loss <- (predictions$outcome - predictions$probability)^2
  reference_loss <- (predictions$outcome - reference_probability)^2

  metric <- sb_binary_brier_metric(predictions, reference_input)
  reference <- metric$reference
  skill <- metric$skill_score

  expect_identical(reference$type, "outer_analysis_prevalence")
  expect_identical(reference$source, "outer_analysis_outcome_only")
  expect_identical(reference$status, "success")
  expect_identical(reference$n, 5L)
  expect_identical(reference$n_scored, 5L)
  expect_equal(
    reference$squared_error_sum,
    sum(reference_loss),
    tolerance = 0
  )
  expect_equal(reference$estimate, mean(reference_loss), tolerance = 0)
  expect_identical(
    reference$per_fold$analysis_prevalence,
    c(0.25, 0.75)
  )
  expect_equal(
    reference$per_fold$squared_error_sum,
    c(sum(reference_loss[1:2]), sum(reference_loss[3:5])),
    tolerance = 0
  )
  expect_identical(skill$metric, "brier_skill_score")
  expect_identical(skill$status, "success")
  expect_equal(
    skill$estimate,
    1 - sum(model_loss) / sum(reference_loss),
    tolerance = 0
  )
})

test_that("Brier skill is unavailable when model OOF is incomplete", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:4),
    fold = rep(1:2, each = 2),
    outcome = c(1L, 0L, 1L, 0L),
    probability = c(0.8, 0.2, NA_real_, NA_real_),
    status = c("success", "success", "failed", "failed"),
    stringsAsFactors = FALSE
  )
  reference_input <- data.frame(
    fold = 1:2,
    analysis_prevalence = c(0.4, 0.6)
  )

  metric <- sb_binary_brier_metric(predictions, reference_input)

  expect_identical(metric$status, "incomplete")
  expect_identical(metric$reference$status, "success")
  expect_identical(metric$reference$n_scored, 4L)
  expect_identical(metric$skill_score$status, "incomplete")
  expect_identical(is.na(metric$skill_score$estimate), TRUE)
})

test_that("OOF log loss is patient-weighted with a fold-local reference", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:5),
    fold = c(1L, 1L, 2L, 2L, 2L),
    outcome = c(1L, 0L, 1L, 0L, 1L),
    probability = c(0.8, 0.25, 0.6, 0.4, 0.9),
    status = "success",
    stringsAsFactors = FALSE
  )
  reference_input <- data.frame(
    fold = 1:2,
    analysis_prevalence = c(0.5, 0.4)
  )
  case_loss <- c(
    -log(0.8), -log(0.75), -log(0.6), -log(0.6), -log(0.9)
  )
  reference_case_loss <- c(
    -log(0.5), -log(0.5), -log(0.4), -log(0.6), -log(0.4)
  )

  metric <- sb_binary_log_loss_metric(predictions, reference_input)

  expect_identical(metric$metric, "log_loss")
  expect_identical(metric$scope, "patient_level_out_of_fold")
  expect_identical(metric$aggregation, "summed_case_loss_over_patients")
  expect_identical(metric$probability_clipping, "none")
  expect_identical(metric$status, "success")
  expect_identical(metric$n, 5L)
  expect_identical(metric$n_scored, 5L)
  expect_identical(metric$n_finite_loss, 5L)
  expect_identical(metric$n_infinite_loss, 0L)
  expect_equal(metric$loss_sum, sum(case_loss), tolerance = 1e-15)
  expect_equal(metric$estimate, mean(case_loss), tolerance = 1e-15)
  expect_identical(metric$per_fold$n, c(2L, 3L))
  expect_equal(
    metric$per_fold$loss_sum,
    c(sum(case_loss[1:2]), sum(case_loss[3:5])),
    tolerance = 1e-15
  )
  expect_equal(
    metric$per_fold$estimate,
    c(mean(case_loss[1:2]), mean(case_loss[3:5])),
    tolerance = 1e-15
  )
  expect_false(isTRUE(all.equal(
    metric$estimate,
    mean(metric$per_fold$estimate),
    tolerance = 0
  )))
  expect_identical(metric$reference$status, "success")
  expect_equal(
    metric$reference$loss_sum,
    sum(reference_case_loss),
    tolerance = 1e-15
  )
  expect_equal(
    metric$reference$estimate,
    mean(reference_case_loss),
    tolerance = 1e-15
  )
  expect_equal(
    metric$reference$per_fold$analysis_prevalence,
    c(0.5, 0.4),
    tolerance = 0
  )
  expect_identical(metric$skill_score$status, "success")
  expect_equal(
    metric$skill_score$estimate,
    1 - sum(case_loss) / sum(reference_case_loss),
    tolerance = 1e-15
  )
})

test_that("log loss does not clip exact probability boundaries", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:4),
    fold = 1L,
    outcome = c(1L, 0L, 1L, 0L),
    probability = c(1, 0, 0, 1),
    status = "success",
    stringsAsFactors = FALSE
  )
  reference_input <- data.frame(
    fold = 1L,
    analysis_prevalence = 0.5
  )

  metric <- sb_binary_log_loss_metric(predictions, reference_input)

  expect_identical(metric$status, "success")
  expect_identical(metric$n_finite_loss, 2L)
  expect_identical(metric$n_infinite_loss, 2L)
  expect_identical(metric$per_fold$n_finite_loss, 2L)
  expect_identical(metric$per_fold$n_infinite_loss, 2L)
  expect_identical(metric$loss_sum, Inf)
  expect_identical(metric$estimate, Inf)
  expect_identical(metric$reference$status, "success")
  expect_identical(is.finite(metric$reference$estimate), TRUE)
  expect_identical(metric$skill_score$status, "success")
  expect_identical(metric$skill_score$estimate, -Inf)
})

test_that("failed log-loss folds retain evidence without a headline", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:4),
    fold = rep(1:2, each = 2L),
    outcome = rep(c(1L, 0L), 2L),
    probability = c(0.8, 0.2, NA_real_, NA_real_),
    status = c("success", "success", "failed", "failed"),
    stringsAsFactors = FALSE
  )
  reference_input <- data.frame(
    fold = 1:2,
    analysis_prevalence = c(0.5, 0.5)
  )

  metric <- sb_binary_log_loss_metric(predictions, reference_input)

  expect_identical(metric$status, "incomplete")
  expect_identical(is.na(metric$estimate), TRUE)
  expect_identical(metric$n_scored, 2L)
  expect_equal(metric$loss_sum, -2 * log(0.8), tolerance = 1e-15)
  expect_identical(
    metric$per_fold$status,
    c("success", "prediction_failed")
  )
  expect_identical(metric$reference$status, "success")
  expect_identical(metric$skill_score$status, "incomplete")
  expect_identical(is.na(metric$skill_score$estimate), TRUE)

  out_of_range <- predictions
  out_of_range$status <- "success"
  out_of_range$probability[3:4] <- c(1.1, 0.2)
  invalid_metric <- sb_binary_log_loss_metric(out_of_range, reference_input)
  expect_identical(invalid_metric$status, "incomplete")
  expect_identical(
    invalid_metric$per_fold$status,
    c("success", "prediction_failed")
  )

  non_finite <- out_of_range
  non_finite$probability[3:4] <- c(Inf, 0.2)
  non_finite_metric <- sb_binary_log_loss_metric(
    non_finite,
    reference_input
  )
  expect_identical(non_finite_metric$status, "incomplete")
  expect_identical(
    non_finite_metric$per_fold$status,
    c("success", "prediction_failed")
  )

  complete_predictions <- predictions
  complete_predictions$status <- "success"
  complete_predictions$probability[3:4] <- c(0.7, 0.3)
  missing_reference <- sb_binary_log_loss_metric(
    complete_predictions,
    reference_input[1L, , drop = FALSE]
  )
  expect_identical(missing_reference$status, "success")
  expect_identical(missing_reference$reference$status, "incomplete")
  expect_identical(missing_reference$skill_score$status, "incomplete")

  no_reference <- sb_binary_log_loss_metric(complete_predictions)
  expect_identical(no_reference$status, "success")
  expect_identical(no_reference$reference$status, "not_requested")
  expect_identical(no_reference$skill_score$status, "incomplete")
})

test_that("OOF calibration intercept and slope match direct logistic fits", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:20),
    fold = rep(1:5, each = 4L),
    outcome = c(
      1L, 1L, 0L, 1L,
      0L, 0L, 1L, 1L,
      1L, 0L, 1L, 0L,
      0L, 1L, 1L, 0L,
      1L, 0L, 1L, 0L
    ),
    probability = seq(0.08, 0.92, length.out = 20L),
    status = "success",
    stringsAsFactors = FALSE
  )
  logit_probability <- stats::qlogis(predictions$probability)
  direct_intercept <- stats::glm(
    predictions$outcome ~ 1,
    family = stats::binomial(),
    offset = logit_probability
  )
  direct_slope <- stats::glm(
    predictions$outcome ~ logit_probability,
    family = stats::binomial()
  )

  metric <- sb_binary_calibration_metric(predictions)

  expect_identical(metric$metric, "calibration_intercept_slope")
  expect_identical(metric$scope, "pooled_patient_level_out_of_fold")
  expect_identical(metric$purpose, "descriptive_evaluation_only")
  expect_identical(metric$probability_clipping, "none")
  expect_identical(metric$status, "success")
  expect_identical(metric$input_status, "success")
  expect_identical(metric$n, 20L)
  expect_identical(metric$n_scored, 20L)
  expect_identical(metric$n_positive, 11L)
  expect_identical(metric$n_negative, 9L)
  expect_identical(metric$n_boundary_probability, 0L)
  expect_identical(metric$intercept_target, 0)
  expect_identical(metric$slope_target, 1)
  expect_equal(
    metric$calibration_intercept,
    unname(stats::coef(direct_intercept)[[1L]]),
    tolerance = 1e-12
  )
  expect_equal(
    metric$calibration_slope,
    unname(stats::coef(direct_slope)[[2L]]),
    tolerance = 1e-12
  )
  expect_identical(metric$intercept_fit$status, "success")
  expect_identical(metric$slope_fit$status, "success")
  expect_identical(metric$intercept_fit$converged, TRUE)
  expect_identical(metric$slope_fit$converged, TRUE)
  expect_equal(
    metric$slope_fit$coefficients$estimate,
    unname(stats::coef(direct_slope)),
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(
    metric$calibration_intercept,
    metric$slope_fit$coefficients$estimate[[1L]],
    tolerance = 1e-6
  )))
  expect_identical(metric$per_fold$n, rep(4L, 5L))
  expect_identical(metric$per_fold$n_positive, c(3L, rep(2L, 4L)))
  expect_identical(metric$per_fold$n_negative, c(1L, rep(2L, 4L)))

  permuted <- predictions[c(20:1), ]
  permuted_metric <- sb_binary_calibration_metric(permuted)
  expect_equal(
    permuted_metric$calibration_intercept,
    metric$calibration_intercept,
    tolerance = 1e-12
  )
  expect_equal(
    permuted_metric$calibration_slope,
    metric$calibration_slope,
    tolerance = 1e-12
  )
  relabeled <- predictions
  relabeled$fold <- 6L - relabeled$fold
  relabeled_metric <- sb_binary_calibration_metric(relabeled)
  expect_equal(
    relabeled_metric$calibration_intercept,
    metric$calibration_intercept,
    tolerance = 1e-12
  )
  expect_equal(
    relabeled_metric$calibration_slope,
    metric$calibration_slope,
    tolerance = 1e-12
  )
})

test_that("empirically ideal probability groups hit calibration targets", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:30),
    fold = rep(1:5, each = 6L),
    outcome = c(
      rep(1L, 2L), rep(0L, 8L),
      rep(1L, 5L), rep(0L, 5L),
      rep(1L, 8L), rep(0L, 2L)
    ),
    probability = rep(c(0.2, 0.5, 0.8), each = 10L),
    status = "success",
    stringsAsFactors = FALSE
  )

  metric <- sb_binary_calibration_metric(predictions)

  expect_identical(metric$status, "success")
  expect_equal(metric$calibration_intercept, 0, tolerance = 1e-10)
  expect_equal(metric$calibration_slope, 1, tolerance = 1e-10)

  logit_probability <- stats::qlogis(predictions$probability)
  recalibration_intercept <-
    metric$slope_fit$coefficients$estimate[[1L]]
  recalibrated <- stats::plogis(
    recalibration_intercept + metric$calibration_slope * logit_probability
  )
  expect_equal(
    sum(predictions$outcome - recalibrated),
    0,
    tolerance = 1e-10
  )
  expect_equal(
    sum((predictions$outcome - recalibrated) * logit_probability),
    0,
    tolerance = 1e-10
  )
})

test_that("OOF calibration keeps non-estimable inputs explicit", {
  predictions <- data.frame(
    patient_id = paste0("P", 1:8),
    fold = rep(1:2, each = 4L),
    outcome = rep(c(0L, 1L), 4L),
    probability = seq(0.15, 0.85, length.out = 8L),
    status = "success",
    stringsAsFactors = FALSE
  )

  failed <- predictions
  failed$status[5:8] <- "failed"
  failed$probability[5:8] <- NA_real_
  failed_metric <- sb_binary_calibration_metric(failed)
  expect_identical(failed_metric$status, "incomplete")
  expect_identical(failed_metric$input_status, "prediction_failed")
  expect_identical(failed_metric$n_scored, 4L)
  expect_identical(
    failed_metric$per_fold$status,
    c("success", "prediction_failed")
  )
  expect_identical(failed_metric$intercept_fit$status, "not_fitted")
  expect_identical(failed_metric$slope_fit$status, "not_fitted")
  expect_identical(is.na(failed_metric$calibration_intercept), TRUE)
  expect_identical(is.na(failed_metric$calibration_slope), TRUE)

  boundary <- predictions
  boundary$probability[[1L]] <- 0
  boundary_metric <- sb_binary_calibration_metric(boundary)
  expect_identical(boundary_metric$status, "incomplete")
  expect_identical(boundary_metric$input_status, "not_evaluable")
  expect_identical(boundary_metric$n_scored, 8L)
  expect_identical(boundary_metric$n_boundary_probability, 1L)
  expect_identical(
    boundary_metric$error_class,
    "stormbreaker_binary_calibration_boundary_probability"
  )
  expect_match(
    boundary_metric$error_message,
    "without probability clipping",
    fixed = TRUE
  )

  out_of_range <- predictions
  out_of_range$probability[[1L]] <- 1.1
  invalid_metric <- sb_binary_calibration_metric(out_of_range)
  expect_identical(invalid_metric$status, "incomplete")
  expect_identical(invalid_metric$input_status, "prediction_failed")

  one_class <- predictions
  one_class$outcome <- 0L
  one_class_metric <- sb_binary_calibration_metric(one_class)
  expect_identical(one_class_metric$status, "incomplete")
  expect_identical(one_class_metric$input_status, "not_evaluable")
  expect_match(
    one_class_metric$error_message,
    "Both outcome classes",
    fixed = TRUE
  )

  constant_probability <- predictions
  constant_probability$probability <- 0.5
  constant_metric <- sb_binary_calibration_metric(constant_probability)
  expect_identical(constant_metric$status, "incomplete")
  expect_identical(constant_metric$input_status, "success")
  expect_identical(constant_metric$intercept_fit$status, "success")
  expect_identical(constant_metric$slope_fit$status, "failed")
  expect_identical(
    constant_metric$slope_fit$error_class,
    "stormbreaker_binary_calibration_nonfinite_coefficients"
  )
  expect_identical(
    constant_metric$error_class,
    "stormbreaker_binary_calibration_nonfinite_coefficients"
  )
  expect_match(
    constant_metric$slope_fit$error_message,
    "non-finite coefficients",
    fixed = TRUE
  )

  separated <- predictions
  separated$outcome <- c(rep(0L, 4L), rep(1L, 4L))
  separated$probability <- c(rep(0.2, 4L), rep(0.8, 4L))
  separated_metric <- sb_binary_calibration_metric(separated)
  expect_identical(separated_metric$status, "incomplete")
  expect_identical(separated_metric$input_status, "success")
  expect_identical(separated_metric$intercept_fit$status, "success")
  expect_identical(separated_metric$slope_fit$status, "failed")
  expect_identical(
    separated_metric$slope_fit$error_class,
    "stormbreaker_calibration_separation"
  )
  expect_match(
    separated_metric$slope_fit$error_message,
    "no finite calibration slope",
    fixed = TRUE
  )
  expect_identical(is.na(separated_metric$calibration_slope), TRUE)

  reversed_separation <- separated
  reversed_separation$probability <- 1 - reversed_separation$probability
  reversed_metric <- sb_binary_calibration_metric(reversed_separation)
  expect_identical(reversed_metric$status, "incomplete")
  expect_identical(
    reversed_metric$slope_fit$error_class,
    "stormbreaker_calibration_separation"
  )

  quasi_separation <- separated
  quasi_separation$probability <- c(0.2, 0.3, 0.4, 0.5, 0.5, 0.6, 0.7, 0.8)
  quasi_metric <- sb_binary_calibration_metric(quasi_separation)
  expect_identical(quasi_metric$status, "incomplete")
  expect_identical(
    quasi_metric$slope_fit$error_class,
    "stormbreaker_calibration_separation"
  )
})
