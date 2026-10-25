stable_survival_fixture <- function() {
  d <- sb_example_survival(100, seed = 20260819)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month", prediction_horizon = 6)
  fit <- sb_fit_final(task, sb_pipeline())
  external <- sb_example_survival(60, seed = 8921)
  external$patient_id <- paste0("external_", seq_len(nrow(external)))
  external$gene_a <- seq(10, 15, length.out = nrow(external)) /
    fit$learner_state$coefficient[1]
  external$gene_b <- external$gene_c <- 0
  new_data <- external[c("patient_id", "gene_a", "gene_b", "gene_c")]
  prediction <- sb_predict_final(fit, new_data, "patient_id")
  endpoint <- sb_surv_endpoint(external, "patient_id", "time", "event",
    "month", prediction_horizon = 6)
  list(fit = fit, external = external, new_data = new_data,
    prediction = prediction, endpoint = endpoint)
}

test_that("saturated Cox predictions preserve their model-defined finite scale", {
  x <- stable_survival_fixture()
  rows <- x$prediction$predictions
  expect_identical(x$prediction$status, "success")
  expect_true(all(rows$survival_probability == 0))
  expect_length(rows$log_cumulative_hazard, nrow(rows))
  expect_equal(rows$log_cumulative_hazard,
    log(x$fit$learner_state$baseline_hazard$cumulative_baseline_hazard) + rows$risk)
  fit <- x$fit
  d <- sb_example_survival(40, seed = 99)
  ordinary <- sb_predict_final(fit, d, "patient_id")
  expect_false("log_cumulative_hazard" %in% names(ordinary$predictions))
  expect_identical(fit, x$fit)
  # Genuine baseline-zero probabilities must not acquire a fabricated scale.
  fit$learner_state$baseline_hazard$cumulative_baseline_hazard <- 0
  expect_null(sb_final_cox_log_cumulative_hazard(fit,
    list(risk = c(1, 2), survival_probability = c(1, 1))))
  fit$learner_state$baseline_hazard$cumulative_baseline_hazard <- 0.5
  expect_equal(sb_final_cox_log_cumulative_hazard(fit,
    list(risk = c(-50, -60), survival_probability = c(1, 1))),
    log(0.5) + c(-50, -60))
  fit$learner_state$baseline_hazard$cumulative_baseline_hazard <- 1e-320
  expect_null(sb_final_cox_log_cumulative_hazard(fit,
    list(risk = c(710, 711), survival_probability = c(0, 0))))
  expect_error(sb_validate_stable_survival_scale(c(0, 1), c(0.3, 0.3)), "disagree")
  expect_error(sb_validate_stable_survival_scale(c(0, Inf), c(0.3, 0.3)), "finite")
})

test_that("external calibration and AUC use the scale, not clipped probabilities", {
  x <- stable_survival_fixture()
  evaluation <- sb_evaluate_external(x$prediction, x$endpoint)
  rows <- evaluation$scored_predictions
  expect_equal(rows$log_cumulative_hazard, x$prediction$predictions$log_cumulative_hazard)
  manual <- survival::coxph(survival::Surv(pmin(time, 6), event == 1 & time <= 6) ~
    log_cumulative_hazard, data = rows, ties = "efron")
  expect_equal(evaluation$metrics$weak_calibration$calibration_slope,
    unname(stats::coef(manual)), tolerance = 1e-10)
  expect_equal(evaluation$metrics$calibration$patient_calibration$transformed_prediction,
    rows$log_cumulative_hazard)
  expect_identical(evaluation$metrics$time_auc$marker, "log_cumulative_hazard")
  expect_identical(evaluation$metrics$calibration$status, "success")
  old <- x$prediction
  old$predictions$log_cumulative_hazard <- NULL
  legacy <- sb_evaluate_external(old, x$endpoint)
  expect_identical(legacy$metrics$weak_calibration$slope_fit$error_class,
    "stormbreaker_survival_weak_calibration_boundary_probability")
  expect_identical(legacy$metrics$calibration$error_class,
    "stormbreaker_survival_calibration_boundary_probability")
  expect_identical(evaluation$metrics$ipcw_brier, legacy$metrics$ipcw_brier)
  expect_identical(evaluation$metrics$harrell_c, legacy$metrics$harrell_c)
  reordered <- x$prediction
  reordered$predictions <- reordered$predictions[rev(seq_len(nrow(rows))), ]
  expect_equal(sb_evaluate_external(reordered, x$endpoint)$scored_predictions,
    evaluation$scored_predictions)
  malformed <- x$prediction
  malformed$predictions$log_cumulative_hazard[1] <- Inf
  expect_identical(sb_evaluate_external(malformed, x$endpoint)$status, "failed")
  contract <- sb_external_bootstrap_contract(evaluation)
  expect_equal(contract$scored$log_cumulative_hazard, rows$log_cumulative_hazard)
  boot <- sb_bootstrap_external(evaluation, reps = 12L, seed = 71L)
  expect_false(identical(boot$status, "failed"))
  expect_true(all(boot$intervals$n_finite_successful[
    boot$intervals$metric %in% c("time_auc", "calibration_slope", "calibration_ici")] > 0L))
})
