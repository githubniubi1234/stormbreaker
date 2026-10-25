# Original-direction metric tests explicitly opt out of direction correction.
evaluate_external_original <- function(prediction, endpoint) {
  sb_evaluate_external(prediction, endpoint, reverse = FALSE)
}

make_sb_binary_final_prediction <- function(
    patient_ids,
    probability,
    pipeline_id = "none__logistic",
    learner = "logistic") {
  n <- length(patient_ids)
  result <- list(
    status = "success",
    task_type = "binary",
    pipeline_id = pipeline_id,
    learner = learner,
    id_col = "sample_id",
    predictions = data.frame(
      patient_id = patient_ids,
      pipeline_id = rep(pipeline_id, n),
      learner = rep(learner, n),
      positive_class = rep.int(1L, n),
      probability = probability,
      status = rep("success", n),
      error_class = rep(NA_character_, n),
      error_message = rep(NA_character_, n),
      stringsAsFactors = FALSE
    ),
    error_stage = NA_character_,
    error_class = NA_character_,
    error_message = NA_character_
  )
  class(result) <- "stormbreaker_final_prediction"
  result
}

make_sb_surv_final_prediction <- function(
    patient_ids,
    risk,
    survival_probability = NULL,
    prediction_horizon = NULL,
    time_unit = NULL,
    pipeline_id = "none__coxph",
    learner = "coxph",
    decision_thresholds = NULL) {
  n <- length(patient_ids)
  rows <- data.frame(
    patient_id = patient_ids,
    pipeline_id = rep(pipeline_id, n),
    learner = rep(learner, n),
    risk = risk,
    stringsAsFactors = FALSE
  )
  if (!is.null(prediction_horizon)) {
    rows$prediction_horizon <- rep(prediction_horizon, n)
    rows$time_unit <- rep(time_unit, n)
    rows$survival_probability <- survival_probability
  }
  rows$status <- rep("success", n)
  rows$error_class <- rep(NA_character_, n)
  rows$error_message <- rep(NA_character_, n)
  result <- list(
    status = "success",
    task_type = "survival",
    pipeline_id = pipeline_id,
    learner = learner,
    id_col = "sample_id",
    predictions = rows,
    error_stage = NA_character_,
    error_class = NA_character_,
    error_message = NA_character_
  )
  if (!is.null(decision_thresholds)) {
    result$decision_thresholds <- decision_thresholds
  }
  class(result) <- "stormbreaker_final_prediction"
  result
}

expect_external_failure <- function(result, stage, class_name) {
  expect_s3_class(result, "stormbreaker_external_evaluation")
  expect_identical(result$status, "failed")
  expect_identical(result$error_stage, stage)
  expect_identical(result$error_class, class_name)
  expect_identical(nrow(result$scored_predictions), 0L)
  expect_identical(result$metrics, list())
}

test_that("external evaluation keeps prediction-first arguments and explicit direction control", {
  expect_identical(
    names(formals(sb_evaluate_external)),
    c("prediction", "endpoint", "reverse")
  )
})

test_that("binary outcomes align by ID into endpoint order", {
  prediction <- make_sb_binary_final_prediction(
    c("p3", "p1", "p2"),
    c(0.9, 0.1, 0.8)
  )
  endpoint <- sb_binary_endpoint(
    data.frame(
      subject = c("p1", "p2", "p3"),
      truth = c(0L, 1L, 1L)
    ),
    "subject",
    "truth"
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_s3_class(result, "stormbreaker_external_evaluation")
  expect_identical(result$status, "success")
  expect_identical(result$scored_predictions$patient_id, c("p1", "p2", "p3"))
  expect_identical(result$scored_predictions$outcome, c(0L, 1L, 1L))
  expect_identical(result$scored_predictions$probability, c(0.1, 0.8, 0.9))
  expect_identical(result$patient_match$status, "success")
  expect_identical(result$patient_match$exact_set, TRUE)
  expect_identical(result$patient_match$output_order, "external_endpoint")
  expect_identical(
    names(result$scored_predictions),
    c("patient_id", "outcome", "probability", "status")
  )
})

test_that("binary external mean calibration is the observed-to-expected ratio", {
  prediction <- make_sb_binary_final_prediction(
    c("p3", "p1", "p2"),
    c(0.7, 0.1, 0.4)
  )
  endpoint <- sb_binary_endpoint(
    data.frame(
      subject = c("p1", "p2", "p3"),
      truth = c(0L, 1L, 1L)
    ),
    "subject",
    "truth"
  )

  result <- evaluate_external_original(prediction, endpoint)
  calibration <- result$metrics$mean_calibration

  expect_identical(calibration$status, "success")
  expect_identical(calibration$ratio_direction, "observed_over_expected")
  expect_identical(calibration$observed_events, 2L)
  expect_equal(calibration$expected_events, 1.2, tolerance = 1e-15)
  expect_equal(calibration$observed_event_probability, 2 / 3)
  expect_equal(calibration$mean_predicted_event_probability, 0.4)
  expect_equal(calibration$mean_calibration_ratio, 2 / 1.2)
  expect_identical(calibration$target, 1)
  expect_match(
    paste(capture.output(print(result)), collapse = "\n"),
    "External mean calibration O:E:",
    fixed = TRUE
  )
})

test_that("binary mean calibration does not invent an O:E denominator", {
  prediction <- make_sb_binary_final_prediction(
    c("p1", "p2"),
    c(0, 0)
  )
  endpoint <- sb_binary_endpoint(
    data.frame(id = c("p1", "p2"), outcome = c(0L, 1L)),
    "id",
    "outcome"
  )

  calibration <- evaluate_external_original(prediction, endpoint)$metrics$
    mean_calibration

  expect_identical(calibration$status, "not_evaluable")
  expect_true(is.na(calibration$mean_calibration_ratio))
  expect_identical(
    calibration$error_class,
    "stormbreaker_binary_mean_calibration_zero_expected"
  )
})

test_that("binary mean calibration does not label an overflowed O:E successful", {
  calibration <- sb_external_binary_mean_calibration_metric(data.frame(
    outcome = 1L,
    probability = 1e-309
  ))

  expect_identical(calibration$status, "not_evaluable")
  expect_identical(is.na(calibration$estimate), TRUE)
  expect_identical(is.na(calibration$mean_calibration_ratio), TRUE)
  expect_identical(calibration$expected_events, 1e-309)
  expect_identical(
    calibration$error_class,
    "stormbreaker_binary_mean_calibration_nonfinite_ratio"
  )
})

test_that("survival outcomes align by ID rather than row position", {
  prediction <- make_sb_surv_final_prediction(
    c("p2", "p3", "p1"),
    c(2, 3, 1)
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      subject = c("p1", "p2", "p3"),
      follow_up = c(1, 2, 3),
      event = c(1L, 1L, 0L)
    ),
    "subject",
    "follow_up",
    "event",
    "month"
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_identical(result$scored_predictions$patient_id, c("p1", "p2", "p3"))
  expect_identical(result$scored_predictions$risk, c(1, 2, 3))
  expect_identical(result$scored_predictions$time, c(1, 2, 3))
  expect_identical(result$scored_predictions$event, c(1L, 1L, 0L))
})

test_that("external evaluation rejects any nonidentical patient set", {
  prediction <- make_sb_binary_final_prediction(
    c("p1", "p2", "p3"),
    c(0.1, 0.5, 0.9)
  )
  endpoint_inputs <- list(
    data.frame(id = c("p1", "p2"), y = c(0L, 1L)),
    data.frame(id = c("p1", "p2", "p3", "p4"), y = c(0L, 1L, 0L, 1L)),
    data.frame(id = c("p1", "p2", "p4"), y = c(0L, 1L, 0L))
  )

  for (data in endpoint_inputs) {
    endpoint <- sb_binary_endpoint(data, "id", "y")
    result <- evaluate_external_original(prediction, endpoint)
    expect_external_failure(
      result,
      "patient_match",
      "stormbreaker_external_patient_mismatch"
    )
    expect_identical(result$patient_match$status, "failed")
    expect_identical(result$patient_match$exact_set, FALSE)
  }
})

test_that("invalid prediction identity and provenance fail before scoring", {
  endpoint <- sb_binary_endpoint(
    data.frame(id = c("p1", "p2"), y = c(0L, 1L)),
    "id",
    "y"
  )
  base <- make_sb_binary_final_prediction(c("p1", "p2"), c(0.2, 0.8))
  duplicate_id <- base
  duplicate_id$predictions$patient_id[2L] <- "p1"
  blank_id <- base
  blank_id$predictions$patient_id[1L] <- ""
  wrong_class <- base
  wrong_class$predictions$positive_class[1L] <- 0L
  wrong_pipeline <- base
  wrong_pipeline$predictions$pipeline_id[1L] <- "other"
  wrong_learner <- base
  wrong_learner$predictions$learner[1L] <- "other"
  failed_row <- base
  failed_row$predictions$status[1L] <- "failed"
  invalid <- list(
    duplicate_id,
    blank_id,
    wrong_class,
    wrong_pipeline,
    wrong_learner,
    failed_row
  )

  for (prediction in invalid) {
    class(prediction) <- "stormbreaker_final_prediction"
    result <- evaluate_external_original(prediction, endpoint)
    expect_external_failure(
      result,
      "prediction_contract",
      "stormbreaker_bad_external_prediction"
    )
  }

  failed <- base
  failed$status <- "failed"
  result <- evaluate_external_original(failed, endpoint)
  expect_external_failure(
    result,
    "prediction_contract",
    "stormbreaker_failed_external_prediction"
  )
})

test_that("endpoint mutation and task-type mismatch remain explicit", {
  binary_prediction <- make_sb_binary_final_prediction(
    c("p1", "p2"),
    c(0.2, 0.8)
  )
  binary_endpoint <- sb_binary_endpoint(
    data.frame(id = c("p1", "p2"), y = c(0L, 1L)),
    "id",
    "y"
  )
  binary_endpoint$n_positive <- 2L
  mutated <- evaluate_external_original(binary_prediction, binary_endpoint)
  expect_external_failure(
    mutated,
    "endpoint_contract",
    "stormbreaker_mutated_endpoint"
  )

  survival_endpoint <- sb_surv_endpoint(
    data.frame(
      id = c("p1", "p2"),
      time = c(1, 2),
      event = c(1L, 0L)
    ),
    "id",
    "time",
    "event",
    "month"
  )
  mismatch <- evaluate_external_original(binary_prediction, survival_endpoint)
  expect_external_failure(
    mismatch,
    "endpoint_contract",
    "stormbreaker_external_task_mismatch"
  )

  modeling_data <- data.frame(
    id = c("p1", "p2"),
    y = c(0L, 1L),
    x = c(1, 2)
  )
  modeling_task <- sb_binary_task(
    modeling_data,
    "id",
    "y",
    "x"
  )
  task_result <- evaluate_external_original(binary_prediction, modeling_task)
  expect_external_failure(
    task_result,
    "endpoint_contract",
    "stormbreaker_bad_endpoint"
  )
})

test_that("an existing prediction is evaluated without reprediction or a fit", {
  development <- sb_example_binary(n = 30L, seed = 701L)
  task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  external <- sb_example_binary(n = 12L, seed = 702L)
  external$patient_id <- paste0("external_", seq_len(nrow(external)))
  prediction <- sb_predict_final(
    fit,
    external[c("patient_id", "gene_a", "gene_b")],
    "patient_id"
  )
  endpoint <- sb_binary_endpoint(
    external,
    "patient_id",
    "outcome"
  )
  rm(fit)
  testthat::local_mocked_bindings(
    sb_predict_final = function(...) stop("must not repredict")
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_identical(result$status, "success")
  expect_identical(
    result$evaluation_boundary$model_access,
    "none"
  )
  expect_identical(
    result$evaluation_boundary$outcome_use,
    "metric_evaluation_only"
  )
  expect_identical(
    result$scored_predictions$probability,
    prediction$predictions$probability
  )
  expect_null(result$fit)
  expect_null(result$model)
  expect_null(result[["task"]])
})

test_that("saved survival predictions integrate with endpoint-only evaluation", {
  development <- sb_example_survival(n = 30L, seed = 703L)
  task <- sb_surv_task(
    development,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 6
  )
  fit <- sb_fit_final(task, sb_pipeline())
  external <- sb_example_survival(n = 12L, seed = 704L)
  external$patient_id <- paste0("external_", seq_len(nrow(external)))
  prediction <- sb_predict_final(
    fit,
    external[c("patient_id", "gene_a", "gene_b")],
    "patient_id"
  )
  endpoint <- sb_surv_endpoint(
    external,
    "patient_id",
    "time",
    "event",
    "month"
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_false(identical(result$status, "failed"))
  expect_identical(result$patient_match$exact_set, TRUE)
  expect_identical(
    result$scored_predictions$risk,
    prediction$predictions$risk
  )
  expect_identical(
    result$scored_predictions$survival_probability,
    prediction$predictions$survival_probability
  )
  expect_identical(
    unique(result$scored_predictions$prediction_horizon),
    6
  )
  calibration <- result$metrics$calibration
  direct_calibration <- sb_survival_calibration_metric(
    patient_id = result$scored_predictions$patient_id,
    time = result$scored_predictions$time,
    event = result$scored_predictions$event,
    survival_probability = result$scored_predictions$survival_probability,
    prediction_status = result$scored_predictions$status,
    prediction_horizon = 6,
    time_unit = "month",
    scope = "patient_level_external_cohort_at_locked_horizon"
  )
  expect_identical(calibration$status, "success")
  expect_identical(
    calibration$scope,
    "patient_level_external_cohort_at_locked_horizon"
  )
  expect_identical(calibration$knots, direct_calibration$knots)
  expect_identical(calibration$curve, direct_calibration$curve)
  expect_identical(
    calibration$patient_calibration,
    direct_calibration$patient_calibration
  )
  expect_identical(
    calibration[c("ici", "e50", "e90")],
    direct_calibration[c("ici", "e50", "e90")]
  )
  weak_calibration <- result$metrics$weak_calibration
  direct_weak_calibration <- sb_survival_weak_calibration_metric(
    patient_id = result$scored_predictions$patient_id,
    time = result$scored_predictions$time,
    event = result$scored_predictions$event,
    survival_probability = result$scored_predictions$survival_probability,
    prediction_status = result$scored_predictions$status,
    prediction_horizon = 6,
    time_unit = "month",
    scope = "patient_level_external_cohort_at_locked_horizon"
  )
  expect_identical(weak_calibration$status, "success")
  expect_identical(
    weak_calibration$patient_predictions$patient_id,
    result$scored_predictions$patient_id
  )
  expect_identical(weak_calibration$per_fold, data.frame(
    fold = integer(),
    status = character(),
    n = integer(),
    n_scored = integer(),
    n_events = integer(),
    n_censored = integer(),
    n_boundary_probability = integer(),
    n_followed_to_horizon = integer(),
    error_class = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  ))
  expect_equal(
    weak_calibration$mean_calibration_ratio,
    direct_weak_calibration$mean_calibration_ratio,
    tolerance = 0
  )
  expect_equal(
    weak_calibration$calibration_slope,
    direct_weak_calibration$calibration_slope,
    tolerance = 0
  )

  shuffled_external <- external[rev(seq_len(nrow(external))), , drop = FALSE]
  shuffled <- evaluate_external_original(
    prediction,
    sb_surv_endpoint(
      shuffled_external,
      "patient_id",
      "time",
      "event",
      "month"
    )
  )
  expect_identical(
    shuffled$metrics$calibration$patient_calibration$patient_id,
    shuffled_external$patient_id
  )
  expect_equal(
    shuffled$metrics$calibration[c("ici", "e50", "e90")],
    calibration[c("ici", "e50", "e90")],
    tolerance = 1e-14
  )
  expect_identical(shuffled$metrics$calibration$knots, calibration$knots)
  expect_equal(
    shuffled$metrics$calibration$curve,
    calibration$curve,
    tolerance = 1e-14
  )
  expect_identical(
    shuffled$metrics$weak_calibration$patient_predictions$patient_id,
    shuffled_external$patient_id
  )
  expect_equal(
    shuffled$metrics$weak_calibration[c(
      "mean_calibration_ratio", "calibration_slope"
    )],
    weak_calibration[c("mean_calibration_ratio", "calibration_slope")],
    tolerance = 1e-13
  )
  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(printed, "External survival calibration ICI:", fixed = TRUE)
  expect_match(printed, "External mean calibration O:E:", fixed = TRUE)
  expect_match(printed, "External survival calibration slope:", fixed = TRUE)
})

test_that("external survival weak calibration is optional and isolated", {
  ids <- paste0("weak_external_", seq_len(12L))
  prediction <- make_sb_surv_final_prediction(
    ids,
    risk = seq(1.1, -0.7, length.out = 12L),
    survival_probability = exp(-exp(seq(-1.2, 0.9, length.out = 12L))),
    prediction_horizon = 6,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
      event = rep(c(1L, 0L), 6L)
    ),
    "id",
    "time",
    "event",
    "month"
  )
  baseline <- evaluate_external_original(prediction, endpoint)
  testthat::local_mocked_bindings(
    sb_survival_weak_calibration_metric = function(...) {
      condition <- simpleError("external weak calibration exploded")
      class(condition) <- c(
        "test_external_weak_calibration_error",
        class(condition)
      )
      stop(condition)
    }
  )

  observed <- evaluate_external_original(prediction, endpoint)

  expect_identical(observed$status, baseline$status)
  expect_identical(observed$scored_predictions, baseline$scored_predictions)
  expect_identical(
    observed$metrics[c(
      "harrell_c", "ipcw_brier", "time_auc", "uno_c", "calibration"
    )],
    baseline$metrics[c(
      "harrell_c", "ipcw_brier", "time_auc", "uno_c", "calibration"
    )]
  )
  expect_identical(
    observed$metrics$weak_calibration$status,
    "metric_failed"
  )
  expect_identical(
    observed$metrics$weak_calibration$error_class,
    "test_external_weak_calibration_error"
  )
})

test_that("external survival calibration is optional and isolated", {
  development <- sb_example_survival(n = 30L, seed = 703L)
  task <- sb_surv_task(
    development,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 6
  )
  fit <- sb_fit_final(task, sb_pipeline())
  external <- sb_example_survival(n = 12L, seed = 704L)
  external$patient_id <- paste0("external_", seq_len(nrow(external)))
  prediction <- sb_predict_final(
    fit,
    external[c("patient_id", "gene_a", "gene_b")],
    "patient_id"
  )
  endpoint <- sb_surv_endpoint(
    external,
    "patient_id",
    "time",
    "event",
    "month"
  )
  baseline <- evaluate_external_original(prediction, endpoint)
  testthat::local_mocked_bindings(
    sb_survival_calibration_metric = function(...) {
      condition <- simpleError("external calibration exploded")
      class(condition) <- c("test_external_calibration_error", class(condition))
      stop(condition)
    }
  )

  observed <- evaluate_external_original(prediction, endpoint)

  expect_identical(observed$status, baseline$status)
  expect_identical(observed$scored_predictions, baseline$scored_predictions)
  expect_identical(
    observed$metrics[c("harrell_c", "ipcw_brier", "time_auc", "uno_c")],
    baseline$metrics[c("harrell_c", "ipcw_brier", "time_auc", "uno_c")]
  )
  expect_identical(observed$metrics$calibration$status, "metric_failed")
  expect_identical(
    observed$metrics$calibration$error_class,
    "test_external_calibration_error"
  )
  expect_identical(
    observed$metrics$calibration$error_message,
    "external calibration exploded"
  )
  expect_match(
    paste(capture.output(print(observed)), collapse = "\n"),
    "External survival calibration: unavailable",
    fixed = TRUE
  )
})

test_that("boundary survival probability only disables external calibration", {
  ids <- paste0("p", seq_len(12L))
  prediction <- make_sb_surv_final_prediction(
    ids,
    risk = seq(1.2, 0.1, length.out = 12L),
    survival_probability = c(0, seq(0.15, 0.85, length.out = 11L)),
    prediction_horizon = 3,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 2, 2.5, 4, 5, 6, 1.5, 3.5, 7, 8, 9, 10),
      event = c(1L, 1L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L, 0L)
    ),
    "id",
    "time",
    "event",
    "month",
    prediction_horizon = 3
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_identical(result$status, "success")
  expect_identical(result$metrics$harrell_c$status, "success")
  expect_identical(result$metrics$ipcw_brier$status, "success")
  expect_identical(result$metrics$time_auc$status, "success")
  expect_identical(result$metrics$uno_c$status, "success")
  expect_identical(result$metrics$calibration$status, "incomplete")
  expect_identical(result$metrics$weak_calibration$status, "incomplete")
  expect_identical(
    result$metrics$weak_calibration$mean_calibration$status,
    "success"
  )
  expect_identical(
    result$metrics$weak_calibration$slope_fit$status,
    "not_evaluable"
  )
  expect_identical(
    result$metrics$calibration$error_class,
    "stormbreaker_survival_calibration_boundary_probability"
  )
  expect_identical(
    result$scored_predictions$survival_probability,
    prediction$predictions$survival_probability
  )
})

test_that("external binary metrics match hand calculations with ties", {
  prediction <- make_sb_binary_final_prediction(
    paste0("p", 1:4),
    c(0.9, 0.9, 0.2, 0.1)
  )
  endpoint <- sb_binary_endpoint(
    data.frame(
      id = paste0("p", 1:4),
      y = c(1L, 0L, 1L, 0L)
    ),
    "id",
    "y"
  )
  result <- evaluate_external_original(prediction, endpoint)
  expected_brier <- mean((c(1, 0, 1, 0) - c(0.9, 0.9, 0.2, 0.1))^2)
  expected_log_loss <- mean(c(
    -log(0.9),
    -log1p(-0.9),
    -log(0.2),
    -log1p(-0.1)
  ))

  expect_equal(result$metrics$roc_auc$estimate, 0.625)
  expect_identical(result$metrics$roc_auc$comparable_pairs, 4L)
  expect_identical(result$metrics$roc_auc$correctly_ordered_pairs, 2)
  expect_identical(result$metrics$roc_auc$tied_pairs, 1)
  expect_equal(result$metrics$average_precision$estimate, 7 / 12)
  expect_equal(result$metrics$brier_score$estimate, expected_brier)
  expect_equal(result$metrics$log_loss$estimate, expected_log_loss)
  expect_identical(result$metrics$roc_auc$scope, "pooled_external_cohort")
  expect_identical(
    result$metrics$brier_score$scope,
    "patient_level_external_cohort"
  )
  expect_identical("fold" %in% names(result$metrics$roc_auc), FALSE)
  expect_identical("per_fold" %in% names(result$metrics$roc_auc), FALSE)
  expect_identical("reference" %in% names(result$metrics$brier_score), FALSE)
})

test_that("external calibration is descriptive and never changes probability", {
  probability <- c(0.2, 0.3, 0.7, 0.8)
  prediction <- make_sb_binary_final_prediction(
    paste0("p", 1:4),
    probability
  )
  endpoint <- sb_binary_endpoint(
    data.frame(
      id = paste0("p", 1:4),
      y = c(0L, 1L, 0L, 1L)
    ),
    "id",
    "y"
  )

  result <- evaluate_external_original(prediction, endpoint)
  calibration <- result$metrics$calibration

  expect_identical(calibration$status, "success")
  expect_true(is.finite(calibration$calibration_intercept))
  expect_true(is.finite(calibration$calibration_slope))
  expect_identical(
    result$scored_predictions$probability,
    probability
  )
  expect_null(result$recalibrated_probability)
})

test_that("constant external probabilities retain a typed calibration failure", {
  ids <- paste0("p", seq_len(8L))
  probability <- rep(0.5, length(ids))
  prediction <- make_sb_binary_final_prediction(ids, probability)
  endpoint <- sb_binary_endpoint(
    data.frame(id = ids, y = rep(c(0L, 1L), 4L)),
    "id",
    "y"
  )

  result <- evaluate_external_original(prediction, endpoint)
  calibration <- result$metrics$calibration

  expect_identical(calibration$status, "incomplete")
  expect_identical(calibration$intercept_fit$status, "success")
  expect_identical(calibration$slope_fit$status, "failed")
  expect_identical(
    calibration$slope_fit$error_class,
    "stormbreaker_binary_calibration_nonfinite_coefficients"
  )
  expect_identical(
    calibration$error_class,
    "stormbreaker_binary_calibration_nonfinite_coefficients"
  )
  expect_identical(result$scored_predictions$probability, probability)
})

test_that("external binary moderate calibration uses one aligned locked fit", {
  ids <- paste0("p", seq_len(24L))
  probability <- seq(0.03, 0.97, length.out = length(ids))
  outcome <- rep(c(0L, 1L), length.out = length(ids))
  prediction <- make_sb_binary_final_prediction(
    rev(ids),
    rev(probability)
  )
  endpoint <- sb_binary_endpoint(
    data.frame(id = ids, y = outcome),
    "id",
    "y"
  )
  original_metric <- sb_binary_moderate_calibration_metric
  calls <- 0L
  observed_scored <- NULL
  testthat::local_mocked_bindings(
    sb_binary_moderate_calibration_metric = function(scored, scope) {
      calls <<- calls + 1L
      observed_scored <<- scored
      original_metric(scored, scope)
    }
  )

  result <- evaluate_external_original(prediction, endpoint)
  direct <- original_metric(
    result$scored_predictions,
    "patient_level_external_cohort"
  )

  expect_identical(calls, 1L)
  expect_identical(observed_scored, result$scored_predictions)
  expect_identical(result$metrics$moderate_calibration$status, "success")
  expect_identical(
    result$metrics$moderate_calibration$patient_calibration$patient_id,
    ids
  )
  expect_identical(result$metrics$moderate_calibration$knots, direct$knots)
  expect_identical(result$metrics$moderate_calibration$curve, direct$curve)
  expect_equal(
    result$metrics$moderate_calibration[c("ici", "e50", "e90")],
    direct[c("ici", "e50", "e90")],
    tolerance = 0
  )
  expect_identical(result$scored_predictions$probability, probability)
  expect_match(
    paste(capture.output(print(result)), collapse = "\n"),
    "External binary moderate calibration ICI:",
    fixed = TRUE
  )
})

test_that("external binary moderate calibration is endpoint-order invariant", {
  ids <- paste0("p", seq_len(24L))
  probability <- seq(0.04, 0.96, length.out = length(ids))
  outcome <- rep(c(0L, 1L), length.out = length(ids))
  prediction <- make_sb_binary_final_prediction(ids, probability)
  first <- evaluate_external_original(
    prediction,
    sb_binary_endpoint(data.frame(id = ids, y = outcome), "id", "y")
  )
  order <- c(seq(2L, length(ids), by = 2L), seq(1L, length(ids), by = 2L))
  second <- evaluate_external_original(
    prediction,
    sb_binary_endpoint(
      data.frame(id = ids[order], y = outcome[order]),
      "id",
      "y"
    )
  )

  expect_identical(first$metrics$moderate_calibration$status, "success")
  expect_identical(second$metrics$moderate_calibration$status, "success")
  expect_equal(
    second$metrics$moderate_calibration[c("ici", "e50", "e90")],
    first$metrics$moderate_calibration[c("ici", "e50", "e90")],
    tolerance = 1e-13
  )
  expect_equal(
    second$metrics$moderate_calibration$curve,
    first$metrics$moderate_calibration$curve,
    tolerance = 1e-13
  )
})

test_that("an external binary moderate-calibration error is isolated", {
  ids <- paste0("p", seq_len(24L))
  probability <- seq(0.03, 0.97, length.out = length(ids))
  endpoint <- sb_binary_endpoint(
    data.frame(id = ids, y = rep(c(0L, 1L), 12L)),
    "id",
    "y"
  )
  prediction <- make_sb_binary_final_prediction(ids, probability)
  control_result <- evaluate_external_original(prediction, endpoint)
  calls <- 0L
  testthat::local_mocked_bindings(
    sb_binary_moderate_calibration_metric = function(...) {
      calls <<- calls + 1L
      condition <- simpleError("external binary moderate calibration exploded")
      class(condition) <- c(
        "test_external_binary_moderate_calibration_error",
        class(condition)
      )
      stop(condition)
    }
  )

  observed <- evaluate_external_original(prediction, endpoint)

  expect_identical(calls, 1L)
  expect_identical(observed$status, control_result$status)
  expect_identical(
    observed$scored_predictions,
    control_result$scored_predictions
  )
  expect_identical(
    observed$metrics[names(observed$metrics) != "moderate_calibration"],
    control_result$metrics[
      names(control_result$metrics) != "moderate_calibration"
    ]
  )
  expect_identical(
    observed$metrics$moderate_calibration$status,
    "metric_failed"
  )
  expect_identical(
    observed$metrics$moderate_calibration$error_class,
    "test_external_binary_moderate_calibration_error"
  )
  expect_match(
    paste(capture.output(print(observed)), collapse = "\n"),
    "External binary moderate calibration: unavailable",
    fixed = TRUE
  )
})

test_that("wrong boundary probabilities retain infinite unclipped log loss", {
  prediction <- make_sb_binary_final_prediction(
    c("p1", "p2"),
    c(0, 1)
  )
  endpoint <- sb_binary_endpoint(
    data.frame(id = c("p1", "p2"), y = c(1L, 0L)),
    "id",
    "y"
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_identical(result$metrics$log_loss$status, "success")
  expect_identical(result$metrics$log_loss$estimate, Inf)
  expect_identical(result$metrics$log_loss$n_infinite_loss, 2L)
  expect_identical(result$metrics$log_loss$probability_clipping, "none")
  expect_identical(result$metrics$calibration$status, "incomplete")
  expect_identical(
    result$metrics$calibration$input_status,
    "not_evaluable"
  )
  expect_identical(
    result$metrics$calibration$error_class,
    "stormbreaker_binary_calibration_boundary_probability"
  )
  expect_identical(result$metrics$moderate_calibration$status, "incomplete")
  expect_identical(
    result$metrics$moderate_calibration$n_boundary_probability,
    2L
  )
  expect_identical(
    result$metrics$moderate_calibration$probability_clipping,
    "none"
  )
})

test_that("single-patient and single-class binary evidence stays partial", {
  cases <- list(
    list(id = "p1", y = 1L, p = 0.7),
    list(id = c("p1", "p2"), y = c(1L, 1L), p = c(0.7, 0.8))
  )

  for (case in cases) {
    prediction <- make_sb_binary_final_prediction(case$id, case$p)
    endpoint <- sb_binary_endpoint(
      data.frame(id = case$id, y = case$y),
      "id",
      "y"
    )
    result <- evaluate_external_original(prediction, endpoint)

    expect_identical(result$status, "incomplete")
    expect_identical(result$metrics$roc_auc$status, "not_evaluable")
    expect_identical(
      result$metrics$average_precision$status,
      "not_evaluable"
    )
    expect_identical(result$metrics$brier_score$status, "success")
    expect_identical(result$metrics$log_loss$status, "success")
    expect_identical(result$metrics$calibration$status, "incomplete")
    expect_identical(
      result$metrics$moderate_calibration$status,
      "incomplete"
    )
    expect_identical(
      result$metrics$moderate_calibration$input_status,
      "not_evaluable"
    )
    expect_identical(is.nan(result$metrics$brier_score$estimate), FALSE)
  }
})

test_that("unrelated endpoint input columns cannot change binary scoring", {
  prediction <- make_sb_binary_final_prediction(
    c("p1", "p2", "p3", "p4"),
    c(0.1, 0.4, 0.6, 0.9)
  )
  first_data <- data.frame(
    id = c("p1", "p2", "p3", "p4"),
    y = c(0L, 0L, 1L, 1L),
    gene = c(1, 2, 3, 4)
  )
  second_data <- first_data[c("gene", "y", "id")]
  second_data$gene <- c(Inf, NA, -100, 100)
  first <- evaluate_external_original(
    prediction,
    sb_binary_endpoint(first_data, "id", "y")
  )
  second <- evaluate_external_original(
    prediction,
    sb_binary_endpoint(second_data, "id", "y")
  )

  expect_identical(first$scored_predictions, second$scored_predictions)
  expect_identical(first$metrics, second$metrics)

  changed <- first_data
  changed$y <- 1L - changed$y
  changed_result <- evaluate_external_original(
    prediction,
    sb_binary_endpoint(changed, "id", "y")
  )
  expect_identical(
    changed_result$scored_predictions$probability,
    first$scored_predictions$probability
  )
  expect_false(identical(changed_result$metrics, first$metrics))
})

test_that("external Harrell C uses the complete supplied cohort", {
  ids <- paste0("p", 1:5)
  risk <- c(1.2, -0.5, 0.7, 2.1, 0.1)
  time <- c(1, 4, 3, 2, 5)
  event <- c(1L, 0L, 1L, 1L, 0L)
  prediction <- make_sb_surv_final_prediction(ids, risk)
  endpoint <- sb_surv_endpoint(
    data.frame(id = ids, time = time, event = event),
    "id",
    "time",
    "event",
    "month"
  )
  expected <- survival::concordance(
    survival::Surv(time, event) ~ risk,
    reverse = TRUE,
    timewt = "n"
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_identical(result$status, "success")
  expect_equal(result$metrics$harrell_c$estimate, expected$concordance)
  expect_identical(result$metrics$harrell_c$scope, "external_cohort")
  expect_identical("ipcw_brier" %in% names(result$metrics), FALSE)
  expect_identical("time_auc" %in% names(result$metrics), FALSE)
  expect_identical("uno_c" %in% names(result$metrics), FALSE)
  expect_false("weak_calibration" %in% names(result$metrics))
})

test_that("external Harrell C preserves its count-ratio rounding", {
  scored <- data.frame(
    time = c(1, 7, 5),
    event = c(1L, 1L, 1L),
    risk = c(-0.8, 0.6, -0.8)
  )
  official <- survival::concordance(
    survival::Surv(time, event) ~ risk,
    data = scored,
    reverse = TRUE,
    timewt = "n"
  )
  expected <- unname(
    (official$count[["concordant"]] + 0.5 * official$count[["tied.x"]]) /
      sum(official$count[c("concordant", "discordant", "tied.x")])
  )

  observed <- sb_external_harrell_c_metric(scored)

  expect_identical(observed$estimate, expected)
  expect_false(identical(observed$estimate, official$concordance))
})

test_that("external Harrell C maps all tie counts without swapping", {
  count_fields <- c(
    "concordant", "discordant", "tied_risk", "tied_time",
    "tied_both", "comparable_pairs"
  )
  cases <- list(
    list(risk = c(0.9, 0.8, 0.2), expected = c(2, 0, 0, 1, 0, 2)),
    list(risk = c(0.5, 0.5, 0.2), expected = c(2, 0, 0, 0, 1, 2))
  )

  for (case in cases) {
    observed <- sb_external_harrell_c_metric(data.frame(
      time = c(1, 1, 2),
      event = c(1L, 1L, 0L),
      risk = case$risk
    ))

    expect_equal(
      unname(unlist(observed[count_fields], use.names = FALSE)),
      case$expected,
      tolerance = 0
    )
  }
})

test_that("external IPCW Brier follows the declared patient-level formula", {
  ids <- c("a", "b", "c", "d", "e")
  probability <- c(0.2, 0.4, 0.6, 0.8, 0.3)
  prediction <- make_sb_surv_final_prediction(
    ids,
    risk = c(2, 0, 1, -1, 1.5),
    survival_probability = probability,
    prediction_horizon = 3,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 2, 3, 4, 2),
      event = c(1L, 0L, 0L, 1L, 1L)
    ),
    "id",
    "time",
    "event",
    "month",
    prediction_horizon = 3
  )

  testthat::local_mocked_bindings(
    sb_fit_event_km = function(...) stop("event KM must not run")
  )
  result <- evaluate_external_original(prediction, endpoint)
  brier <- result$metrics$ipcw_brier
  time_auc <- result$metrics$time_auc
  uno_c <- result$metrics$uno_c

  expected_weighted_loss <-
    0.2^2 +
    0.3^2 +
    (1 - 0.6)^2 * 1.5 +
    (1 - 0.8)^2 * 1.5
  expect_identical(brier$status, "success")
  expect_equal(brier$estimate, expected_weighted_loss / 5)
  expect_equal(brier$weighted_loss_sum, expected_weighted_loss)
  expect_equal(brier$weight_sum, 5)
  expect_identical(brier$n_scored, 4L)
  expect_identical(brier$n_censored_at_horizon, 1L)
  expect_identical(brier$censoring_model, "external_reverse_kaplan_meier")
  expect_identical(
    brier$censoring_assumption,
    "marginal_noninformative_censoring"
  )
  expect_identical("reference" %in% names(brier), FALSE)
  expect_identical("skill_score" %in% names(brier), FALSE)
  expect_identical("per_fold" %in% names(brier), FALSE)
  expect_identical(time_auc$status, "success")
  expect_equal(time_auc$estimate, 1, tolerance = 0)
  expect_identical(time_auc$n_cases, 2L)
  expect_identical(time_auc$n_controls, 1L)
  expect_identical(time_auc$n_excluded_censored, 2L)
  expect_equal(time_auc$case_weight_sum, 2, tolerance = 0)
  expect_equal(time_auc$weighted_pair_credit, 2, tolerance = 0)
  expect_equal(time_auc$pair_denominator, 2, tolerance = 0)
  expect_identical(time_auc$marker, "one_minus_survival_probability")
  expect_identical(uno_c$status, "success")
  expect_equal(uno_c$estimate, 1, tolerance = 0)
  expect_identical(uno_c$n_eligible_events, 2L)
  expect_identical(uno_c$n_comparable_pairs, 7L)
  expect_identical(uno_c$n_tied_time_comparable_pairs, 1L)
  expect_equal(uno_c$weighted_pair_credit, 7, tolerance = 0)
  expect_equal(uno_c$weighted_pair_denominator, 7, tolerance = 0)
  expect_identical(uno_c$marker, "risk")
})

test_that("external time AUC reuses one curve without Brier status coupling", {
  ids <- c("case", "control", "censored")
  prediction <- make_sb_surv_final_prediction(
    ids,
    risk = c(2, 0, 1),
    survival_probability = c(0.2, 0.8, 0.5),
    prediction_horizon = 2,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 3, 1.5),
      event = c(1L, 1L, 0L)
    ),
    "id",
    "time",
    "event",
    "month",
    prediction_horizon = 2
  )
  original_fit <- sb_fit_censoring_km
  original_brier <- sb_external_ipcw_brier_metric
  censoring_calls <- 0L
  testthat::local_mocked_bindings(
    sb_fit_censoring_km = function(...) {
      censoring_calls <<- censoring_calls + 1L
      original_fit(...)
    },
    sb_external_ipcw_brier_metric = function(...) {
      metric <- original_brier(...)
      metric$status <- "incomplete"
      metric$estimate <- NA_real_
      metric
    }
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_identical(censoring_calls, 1L)
  expect_identical(result$status, "incomplete")
  expect_identical(result$metrics$ipcw_brier$status, "incomplete")
  expect_identical(result$metrics$time_auc$status, "success")
  expect_equal(result$metrics$time_auc$estimate, 1, tolerance = 0)
  expect_identical(result$metrics$uno_c$status, "success")
  expect_equal(result$metrics$uno_c$estimate, 1, tolerance = 0)
})

test_that("external Uno C uses the aligned external censoring distribution", {
  ids <- paste0("external_g_", 1:5)
  time <- c(1, 2, 3, 4, 5)
  event <- c(1L, 0L, 1L, 0L, 0L)
  risk <- c(0, 1, 3, 2, -1)
  prediction <- make_sb_surv_final_prediction(
    ids,
    risk,
    survival_probability = c(0.8, 0.7, 0.4, 0.3, 0.2),
    prediction_horizon = 5,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(id = ids, time = time, event = event),
    "id",
    "time",
    "event",
    "month",
    prediction_horizon = 5
  )

  result <- evaluate_external_original(prediction, endpoint)
  scored <- result$scored_predictions
  external_curve <- sb_fit_censoring_km(time, event)
  development_curve <- sb_fit_censoring_km(1:5, rep.int(1L, 5L))
  expected <- sb_external_uno_c_metric(
    scored,
    list(time_unit = "month"),
    external_curve
  )
  counterfactual <- sb_external_uno_c_metric(
    scored,
    list(time_unit = "month"),
    development_curve
  )

  expect_identical(result$metrics$uno_c$status, "success")
  expect_equal(result$metrics$uno_c$estimate, expected$estimate, tolerance = 0)
  expect_false(isTRUE(all.equal(expected$estimate, counterfactual$estimate)))
  expect_equal(expected$estimate, 4.55555555555556 / 7.55555555555556)
  expect_equal(counterfactual$estimate, 0.5)
})

test_that("an Uno exception remains inside the external metric", {
  ids <- paste0("uno_error_", 1:6)
  prediction <- make_sb_surv_final_prediction(
    ids,
    risk = c(3, 2, 1, 0, -1, -2),
    survival_probability = c(0.2, 0.3, 0.4, 0.6, 0.7, 0.8),
    prediction_horizon = 3,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = 1:6,
      event = c(1L, 1L, 0L, 1L, 0L, 1L)
    ),
    "id",
    "time",
    "event",
    "month",
    prediction_horizon = 3
  )
  baseline <- evaluate_external_original(prediction, endpoint)
  testthat::local_mocked_bindings(
    sb_uno_c_stratum = function(...) {
      condition <- simpleError("injected Uno failure")
      class(condition) <- c("test_uno_error", class(condition))
      stop(condition)
    }
  )

  observed <- evaluate_external_original(prediction, endpoint)

  expect_identical(observed$status, "incomplete")
  expect_identical(observed$scored_predictions, baseline$scored_predictions)
  expect_identical(
    observed$metrics[names(observed$metrics) != "uno_c"],
    baseline$metrics[names(baseline$metrics) != "uno_c"]
  )
  expect_identical(observed$metrics$uno_c$status, "metric_failed")
  expect_identical(observed$metrics$uno_c$error_class, "test_uno_error")
  expect_identical(
    observed$metrics$uno_c$error_message,
    "injected Uno failure"
  )
})

test_that("all-event IPCW remains valid beyond the largest observed time", {
  ids <- paste0("p", 1:4)
  probability <- c(0.1, 0.2, 0.3, 0.4)
  prediction <- make_sb_surv_final_prediction(
    ids,
    risk = c(4, 3, 2, 1),
    survival_probability = probability,
    prediction_horizon = 5,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(id = ids, time = 1:4, event = 1L),
    "id",
    "time",
    "event",
    "month"
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_identical(result$status, "incomplete")
  expect_identical(result$metrics$ipcw_brier$status, "success")
  expect_identical(result$metrics$time_auc$status, "not_evaluable")
  expect_identical(
    result$metrics$time_auc$error_class,
    "stormbreaker_no_time_auc_control"
  )
  expect_equal(
    result$metrics$ipcw_brier$estimate,
    mean(probability^2)
  )
  expect_true(all(
    result$metrics$ipcw_brier$censoring_curve$survival_after == 1
  ))
})

test_that("survival endpoint information propagates without fabricated metrics", {
  one_prediction <- make_sb_surv_final_prediction("p1", 1)
  one_endpoint <- sb_surv_endpoint(
    data.frame(id = "p1", time = 2, event = 1L),
    "id",
    "time",
    "event",
    "month"
  )
  one <- evaluate_external_original(one_prediction, one_endpoint)
  expect_identical(one$status, "incomplete")
  expect_identical(one$metrics$harrell_c$status, "not_evaluable")
  expect_identical(is.nan(one$metrics$harrell_c$estimate), FALSE)

  after_ids <- c("a", "b", "c")
  after_prediction <- make_sb_surv_final_prediction(
    after_ids,
    risk = c(1, 2, 3),
    survival_probability = c(0.8, 0.7, 0.6),
    prediction_horizon = 3,
    time_unit = "month"
  )
  after_endpoint <- sb_surv_endpoint(
    data.frame(id = after_ids, time = c(4, 5, 6), event = 0L),
    "id",
    "time",
    "event",
    "month"
  )
  after <- evaluate_external_original(after_prediction, after_endpoint)
  expect_identical(after$status, "incomplete")
  expect_identical(after$metrics$harrell_c$status, "not_evaluable")
  expect_identical(after$metrics$ipcw_brier$status, "success")
  expect_equal(
    after$metrics$ipcw_brier$estimate,
    mean((1 - c(0.8, 0.7, 0.6))^2)
  )

  before_prediction <- make_sb_surv_final_prediction(
    after_ids,
    risk = c(1, 2, 3),
    survival_probability = c(0.8, 0.7, 0.6),
    prediction_horizon = 4,
    time_unit = "month"
  )
  before_endpoint <- sb_surv_endpoint(
    data.frame(id = after_ids, time = c(1, 2, 3), event = 0L),
    "id",
    "time",
    "event",
    "month"
  )
  before <- evaluate_external_original(before_prediction, before_endpoint)
  expect_identical(before$status, "incomplete")
  expect_identical(before$metrics$ipcw_brier$status, "incomplete")
  expect_identical(before$metrics$ipcw_brier$n_scored, 0L)
  expect_identical(
    before$metrics$ipcw_brier$error_class,
    "stormbreaker_no_evaluable_horizon_status"
  )
  expect_identical(is.nan(before$metrics$ipcw_brier$estimate), FALSE)
})

test_that("survival horizon and time-unit states are never silently changed", {
  ids <- c("p1", "p2", "p3")
  risk_only <- make_sb_surv_final_prediction(ids, c(3, 2, 1))
  horizon <- make_sb_surv_final_prediction(
    ids,
    c(3, 2, 1),
    survival_probability = c(0.2, 0.5, 0.8),
    prediction_horizon = 3,
    time_unit = "month"
  )
  data <- data.frame(
    id = ids,
    time = c(1, 3, 4),
    event = c(1L, 0L, 1L)
  )
  no_horizon <- sb_surv_endpoint(
    data, "id", "time", "event", "month"
  )
  declared <- sb_surv_endpoint(
    data, "id", "time", "event", "month", 3
  )

  risk_result <- evaluate_external_original(risk_only, no_horizon)
  expect_identical(risk_result$status, "success")
  expect_identical("ipcw_brier" %in% names(risk_result$metrics), FALSE)
  expect_identical("time_auc" %in% names(risk_result$metrics), FALSE)
  expect_identical("uno_c" %in% names(risk_result$metrics), FALSE)
  expect_false("weak_calibration" %in% names(risk_result$metrics))

  risk_mismatch <- evaluate_external_original(risk_only, declared)
  expect_external_failure(
    risk_mismatch,
    "horizon_semantics",
    "stormbreaker_external_horizon_mismatch"
  )

  inherited <- evaluate_external_original(horizon, no_horizon)
  expect_identical(inherited$status, "success")
  expect_identical(
    unique(inherited$scored_predictions$prediction_horizon),
    3
  )

  matched <- evaluate_external_original(horizon, declared)
  expect_identical(matched$status, "success")

  different_horizon <- sb_surv_endpoint(
    data, "id", "time", "event", "month", 2
  )
  horizon_mismatch <- evaluate_external_original(horizon, different_horizon)
  expect_external_failure(
    horizon_mismatch,
    "horizon_semantics",
    "stormbreaker_external_horizon_mismatch"
  )

  different_unit <- sb_surv_endpoint(
    data, "id", "time", "event", "day"
  )
  unit_mismatch <- evaluate_external_original(horizon, different_unit)
  expect_external_failure(
    unit_mismatch,
    "horizon_semantics",
    "stormbreaker_external_time_unit_mismatch"
  )
})

test_that("partial or row-inconsistent survival probability semantics fail", {
  ids <- c("p1", "p2", "p3")
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 4, 5),
      event = c(1L, 0L, 1L)
    ),
    "id",
    "time",
    "event",
    "month"
  )
  base <- make_sb_surv_final_prediction(
    ids,
    c(3, 2, 1),
    survival_probability = c(0.2, 0.5, 0.8),
    prediction_horizon = 3,
    time_unit = "month"
  )
  partial <- base
  partial$predictions$survival_probability <- NULL
  different_horizon <- base
  different_horizon$predictions$prediction_horizon[2L] <- 4
  different_unit <- base
  different_unit$predictions$time_unit[2L] <- "day"

  for (prediction in list(partial, different_horizon, different_unit)) {
    result <- evaluate_external_original(prediction, endpoint)
    expect_external_failure(
      result,
      "prediction_contract",
      "stormbreaker_bad_external_prediction"
    )
  }
})

test_that("a required zero censoring denominator stays incomplete", {
  ids <- c("p1", "p2", "p3")
  prediction <- make_sb_surv_final_prediction(
    ids,
    c(3, 2, 1),
    survival_probability = c(0.2, 0.5, 0.8),
    prediction_horizon = 3,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 4, 5),
      event = c(1L, 0L, 1L)
    ),
    "id",
    "time",
    "event",
    "month"
  )
  testthat::local_mocked_bindings(
    sb_censoring_survival_at = function(...) 0
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_identical(result$status, "incomplete")
  expect_identical(result$metrics$ipcw_brier$status, "incomplete")
  expect_identical(
    result$metrics$ipcw_brier$error_class,
    "stormbreaker_zero_external_censoring_survival"
  )
})

test_that("unrelated survival endpoint input columns never enter scoring", {
  ids <- c("p1", "p2", "p3")
  prediction <- make_sb_surv_final_prediction(
    ids,
    c(3, 2, 1),
    survival_probability = c(0.2, 0.5, 0.8),
    prediction_horizon = 3,
    time_unit = "month"
  )
  first_data <- data.frame(
    id = ids,
    time = c(1, 4, 5),
    event = c(1L, 0L, 1L),
    gene = c(1, 2, 3)
  )
  second_data <- first_data[c("gene", "event", "time", "id")]
  second_data$gene <- c(NA, Inf, -Inf)

  first <- evaluate_external_original(
    prediction,
    sb_surv_endpoint(
      first_data, "id", "time", "event", "month"
    )
  )
  second <- evaluate_external_original(
    prediction,
    sb_surv_endpoint(
      second_data, "id", "time", "event", "month"
    )
  )

  expect_identical(first$scored_predictions, second$scored_predictions)
  expect_identical(first$metrics, second$metrics)
})

test_that("horizon probabilities may vary by patient but remain bounded", {
  ids <- c("p1", "p2", "p3")
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 4, 5),
      event = c(1L, 0L, 1L)
    ),
    "id",
    "time",
    "event",
    "month"
  )
  valid <- make_sb_surv_final_prediction(
    ids,
    c(3, 2, 1),
    survival_probability = c(0.2, 0.5, 0.8),
    prediction_horizon = 3,
    time_unit = "month"
  )
  result <- evaluate_external_original(valid, endpoint)
  expect_identical(result$status, "success")
  expect_identical(
    result$scored_predictions$survival_probability,
    c(0.2, 0.5, 0.8)
  )

  invalid <- valid
  invalid$predictions$survival_probability[2L] <- 2
  failure <- evaluate_external_original(invalid, endpoint)
  expect_external_failure(
    failure,
    "prediction_contract",
    "stormbreaker_bad_external_prediction"
  )
})

test_that("external evaluation survives prediction serialization", {
  binary_prediction <- make_sb_binary_final_prediction(
    c("b1", "b2"),
    c(0.2, 0.8)
  )
  binary_endpoint <- sb_binary_endpoint(
    data.frame(id = c("b1", "b2"), y = c(0L, 1L)),
    "id",
    "y"
  )
  survival_prediction <- make_sb_surv_final_prediction(
    c("s1", "s2", "s3"),
    c(3, 2, 1),
    survival_probability = c(0.2, 0.5, 0.8),
    prediction_horizon = 3,
    time_unit = "month"
  )
  survival_endpoint <- sb_surv_endpoint(
    data.frame(
      id = c("s1", "s2", "s3"),
      time = c(1, 4, 5),
      event = c(1L, 0L, 1L)
    ),
    "id",
    "time",
    "event",
    "month"
  )
  risk_prediction <- make_sb_surv_final_prediction(
    c("s1", "s2", "s3"),
    c(3, 2, 1)
  )

  binary_before <- evaluate_external_original(
    binary_prediction,
    binary_endpoint
  )
  binary_after <- evaluate_external_original(
    unserialize(serialize(binary_prediction, NULL)),
    binary_endpoint
  )
  survival_before <- evaluate_external_original(
    survival_prediction,
    survival_endpoint
  )
  survival_after <- evaluate_external_original(
    unserialize(serialize(survival_prediction, NULL)),
    survival_endpoint
  )
  risk_before <- evaluate_external_original(risk_prediction, survival_endpoint)
  risk_after <- evaluate_external_original(
    unserialize(serialize(risk_prediction, NULL)),
    survival_endpoint
  )

  expect_identical(binary_after, binary_before)
  expect_identical(survival_after, survival_before)
  expect_identical(risk_after, risk_before)
})

test_that("binary external decision curve uses development-declared thresholds", {
  development <- make_sb_binary_fixture(n = 60L, seed = 20260910L)
  external <- make_sb_binary_fixture(n = 30L, seed = 20260911L)
  external$patient_id <- paste0("external_", seq_len(nrow(external)))
  predictors <- c("gene_a", "gene_b", "gene_c")
  thresholds <- c(0.20, 0.40, 0.60)
  task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    predictors,
    decision_thresholds = thresholds
  )
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  prediction <- sb_predict_final(
    fit,
    external[c("patient_id", predictors)],
    "patient_id"
  )
  endpoint <- sb_binary_endpoint(
    external, "patient_id", "outcome"
  )

  result <- evaluate_external_original(prediction, endpoint)
  metric <- result$metrics$decision_curve

  expect_identical(result$status, "success")
  expect_identical(metric$status, "success")
  expect_identical(metric$scope, "patient_level_external_cohort")
  expect_identical(metric$curve$threshold_probability, thresholds)
  for (index in seq_along(thresholds)) {
    threshold <- thresholds[[index]]
    selected <- result$scored_predictions$probability >= threshold
    true_positive <- as.integer(sum(
      selected & result$scored_predictions$outcome == 1L
    ))
    false_positive <- as.integer(sum(
      selected & result$scored_predictions$outcome == 0L
    ))
    expected <- true_positive / nrow(external) -
      false_positive / nrow(external) * threshold / (1 - threshold)

    expect_identical(metric$curve$true_positive[[index]], true_positive)
    expect_identical(metric$curve$false_positive[[index]], false_positive)
    expect_equal(
      metric$curve$model_net_benefit[[index]],
      expected,
      tolerance = 1e-15
    )
  }

  reordered <- external[rev(seq_len(nrow(external))), , drop = FALSE]
  reordered_result <- evaluate_external_original(
    prediction,
    sb_binary_endpoint(reordered, "patient_id", "outcome")
  )
  expect_identical(
    reordered_result$metrics$decision_curve$curve,
    metric$curve
  )
})

test_that("binary external evaluation omits undeclared decision curves", {
  prediction <- make_sb_binary_final_prediction(
    c("p1", "p2", "p3", "p4"),
    c(0.1, 0.4, 0.6, 0.9)
  )
  endpoint <- sb_binary_endpoint(
    data.frame(
      id = c("p1", "p2", "p3", "p4"),
      y = c(0L, 0L, 1L, 1L)
    ),
    "id",
    "y"
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_false("decision_curve" %in% names(result$metrics))
})

test_that("invalid prediction-side decision thresholds fail before scoring", {
  prediction <- make_sb_binary_final_prediction(
    c("p1", "p2"),
    c(0.2, 0.8)
  )
  prediction$decision_thresholds <- c(0.2, 1)
  endpoint <- sb_binary_endpoint(
    data.frame(id = c("p1", "p2"), y = c(0L, 1L)),
    "id",
    "y"
  )

  result <- evaluate_external_original(prediction, endpoint)

  expect_external_failure(
    result,
    "prediction_contract",
    "stormbreaker_bad_decision_thresholds"
  )
})

test_that("external evaluation printing reports metrics without proof claims", {
  prediction <- make_sb_binary_final_prediction(
    c("p1", "p2"),
    c(0.2, 0.8)
  )
  endpoint <- sb_binary_endpoint(
    data.frame(id = c("p1", "p2"), y = c(0L, 1L)),
    "id",
    "y"
  )
  result <- evaluate_external_original(prediction, endpoint)

  output <- paste(capture.output(print(result)), collapse = "\n")

  expect_match(output, "precomputed-prediction external evaluation")
  expect_match(output, "External ROC AUC")
  expect_match(output, "Patients: 2")
  expect_no_match(output, "transportability.*proven", ignore.case = TRUE)
  expect_no_match(output, "clinical utility.*proven", ignore.case = TRUE)
})

test_that("external survival DCA uses fixed predictions after ID alignment", {
  endpoint_ids <- paste0("p", 1:6)
  prediction_order <- c(6L, 2L, 5L, 1L, 4L, 3L)
  time <- 1:6
  event <- c(1L, 0L, 1L, 0L, 1L, 0L)
  event_probability <- c(0.80, 0.70, 0.60, 0.40, 0.55, 0.10)
  prediction <- make_sb_surv_final_prediction(
    endpoint_ids[prediction_order],
    risk = seq_along(prediction_order),
    survival_probability = (1 - event_probability)[prediction_order],
    prediction_horizon = 4,
    time_unit = "month",
    decision_thresholds = 0.50
  )
  endpoint <- sb_surv_endpoint(
    data.frame(id = endpoint_ids, time = time, event = event),
    "id",
    "time",
    "event",
    "month",
    prediction_horizon = 4
  )

  result <- evaluate_external_original(prediction, endpoint)
  metric <- result$metrics$decision_curve
  row <- metric$curve[1L, ]

  expect_identical(result$status, "success")
  expect_identical(metric$status, "success")
  expect_identical(metric$scope, "patient_level_external_cohort")
  expect_identical(result$scored_predictions$patient_id, endpoint_ids)
  expect_equal(
    1 - result$scored_predictions$survival_probability,
    event_probability
  )
  expect_identical(row$n_selected, 4L)
  expect_equal(row$selected_survival_probability, 3 / 8, tolerance = 0)
  expect_equal(row$model_net_benefit, 1 / 6, tolerance = 1e-15)
  expect_equal(row$treat_all_net_benefit, -1 / 4, tolerance = 0)
  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(
    printed,
    "External survival decision curve: 1 declared threshold"
  )
})

test_that("external survival DCA retains endpoint-specific estimability", {
  ids <- c("p1", "p2")
  for (event in list(c(1L, 1L), c(0L, 0L))) {
    prediction <- make_sb_surv_final_prediction(
      ids,
      risk = c(1, 2),
      survival_probability = c(0.40, 0.40),
      prediction_horizon = 2,
      time_unit = "month",
      decision_thresholds = 0.50
    )
    endpoint <- sb_surv_endpoint(
      data.frame(id = ids, time = c(1, 3), event = event),
      "id",
      "time",
      "event",
      "month",
      prediction_horizon = 2
    )
    result <- evaluate_external_original(prediction, endpoint)

    expect_identical(result$metrics$decision_curve$status, "success")
    expect_true(is.finite(
      result$metrics$decision_curve$curve$model_net_benefit
    ))
  }

  prediction <- make_sb_surv_final_prediction(
    c("p1", "p2", "p3"),
    risk = c(1, 2, 3),
    survival_probability = c(0.20, 0.20, 0.90),
    prediction_horizon = 3,
    time_unit = "month",
    decision_thresholds = 0.50
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = c("p1", "p2", "p3"),
      time = c(1, 2, 4),
      event = c(1L, 0L, 0L)
    ),
    "id",
    "time",
    "event",
    "month",
    prediction_horizon = 3
  )
  incomplete <- evaluate_external_original(prediction, endpoint)

  expect_identical(incomplete$status, "incomplete")
  expect_identical(
    incomplete$metrics$decision_curve$curve$horizon_support,
    "insufficient_follow_up"
  )
})

test_that("external survival thresholds are optional and strictly validated", {
  ids <- c("p1", "p2", "p3")
  plain <- make_sb_surv_final_prediction(
    ids,
    risk = c(1, 2, 3),
    survival_probability = c(0.8, 0.5, 0.2),
    prediction_horizon = 2,
    time_unit = "month"
  )
  endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 2, 3),
      event = c(1L, 0L, 1L)
    ),
    "id",
    "time",
    "event",
    "month",
    prediction_horizon = 2
  )
  plain_result <- evaluate_external_original(plain, endpoint)
  expect_false("decision_curve" %in% names(plain_result$metrics))

  invalid <- plain
  invalid$decision_thresholds <- c(0.20, 1)
  invalid_result <- evaluate_external_original(invalid, endpoint)
  expect_external_failure(
    invalid_result,
    "prediction_contract",
    "stormbreaker_bad_decision_thresholds"
  )

  risk_only <- make_sb_surv_final_prediction(
    ids,
    risk = c(1, 2, 3),
    decision_thresholds = 0.50
  )
  risk_endpoint <- sb_surv_endpoint(
    data.frame(
      id = ids,
      time = c(1, 2, 3),
      event = c(1L, 0L, 1L)
    ),
    "id",
    "time",
    "event",
    "month"
  )
  risk_result <- evaluate_external_original(risk_only, risk_endpoint)
  expect_external_failure(
    risk_result,
    "prediction_contract",
    "stormbreaker_bad_external_prediction"
  )
})
