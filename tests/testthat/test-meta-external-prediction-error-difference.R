make_sb_external_prediction_error_difference_bootstrap <- function(
    task_type,
    point,
    standard_error,
    metric = NULL,
    n = 100L,
    reference_pipeline_id = "reference_model",
    candidate_pipeline_id = "candidate_model",
    prediction_target = NULL) {
  if (is.null(metric)) {
    metric <- if (identical(task_type, "binary")) {
      "brier_score"
    } else {
      "ipcw_brier_score"
    }
  }
  if (is.null(prediction_target)) {
    prediction_target <- if (identical(task_type, "binary")) {
      list(
        target_status = "complete",
        has_horizon = FALSE,
        prediction_horizon = NULL,
        time_unit = NULL,
        horizon_source = "not_applicable"
      )
    } else {
      list(
        target_status = "complete",
        has_horizon = TRUE,
        prediction_horizon = 6,
        time_unit = "month",
        horizon_source = "endpoint"
      )
    }
  }
  reference_estimate <- if (identical(metric, "log_loss")) 0.60 else 0.20
  replicate_values <- point + c(-1, 1) * standard_error / sqrt(2)
  candidate_point_estimate <- reference_estimate + point
  candidate_replicate_estimates <- reference_estimate + replicate_values
  point_difference <- candidate_point_estimate - reference_estimate
  replicate_differences <-
    candidate_replicate_estimates - reference_estimate

  result <- list(
    status = "success",
    task_type = task_type,
    reference_pipeline_id = reference_pipeline_id,
    candidate_status = data.frame(
      candidate_pipeline_id = candidate_pipeline_id,
      status = "success",
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    ),
    n = as.integer(n),
    prediction_target = prediction_target,
    reps = 2L,
    metrics = metric,
    point_differences = data.frame(
      reference_pipeline_id = reference_pipeline_id,
      candidate_pipeline_id = candidate_pipeline_id,
      contrast = "candidate_minus_reference",
      metric = metric,
      metric_direction = "lower_is_better",
      reference_estimate = reference_estimate,
      candidate_estimate = candidate_point_estimate,
      estimate = point_difference,
      status = "success",
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    ),
    differences = data.frame(
      replicate = 1:2,
      reference_pipeline_id = reference_pipeline_id,
      candidate_pipeline_id = candidate_pipeline_id,
      contrast = "candidate_minus_reference",
      metric = metric,
      metric_direction = "lower_is_better",
      reference_estimate = rep(reference_estimate, 2L),
      candidate_estimate = candidate_replicate_estimates,
      estimate = replicate_differences,
      status = "success",
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
  )
  class(result) <- "stormbreaker_external_bootstrap_comparison"
  result
}

test_that("paired Brier differences have a fixed numerical estimand", {
  skip_if_not_installed("metafor", minimum_version = "4.0.0")

  points <- c(-0.05, 0.02, 0.08, 0.12)
  standard_errors <- c(0.03, 0.04, 0.025, 0.035)
  comparisons <- Map(
    function(point, standard_error, index) {
      make_sb_external_prediction_error_difference_bootstrap(
        "binary", point, standard_error, n = 90L + index
      )
    },
    points,
    standard_errors,
    seq_along(points)
  )
  names(comparisons) <- paste0("cohort_", letters[1:4])

  observed <- sb_meta_external_prediction_error_difference(
    comparisons,
    candidate_pipeline_id = "candidate_model"
  )

  expect_s3_class(
    observed,
    "stormbreaker_external_prediction_error_difference_meta"
  )
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(observed$reference_pipeline_id, "reference_model")
  expect_identical(observed$candidate_pipeline_id, "candidate_model")
  expect_identical(observed$contrast, "candidate_minus_reference")
  expect_identical(observed$metric, "brier_score")
  expect_identical(observed$cohorts$cohort_id, names(comparisons))
  expect_equal(observed$cohorts$point_difference, points, tolerance = 1e-15)
  expect_equal(
    observed$cohorts$standard_error,
    standard_errors,
    tolerance = 1e-12
  )
  expect_equal(observed$summary$estimate, 0.042549460390, tolerance = 1e-10)
  expect_equal(observed$summary$lower, -0.075960362323, tolerance = 1e-10)
  expect_equal(observed$summary$upper, 0.161059283102, tolerance = 1e-10)
  expect_equal(observed$summary$tau2, 0.004554259189, tolerance = 1e-10)
  expect_identical(
    observed$method$effect_scale,
    "raw_candidate_minus_reference_prediction_error"
  )
  expect_identical(observed$method$tau2_estimator, "REML")
  expect_identical(observed$method$inference, "Knapp-Hartung")
  expect_identical(observed$provenance, list(
    package = "stormbreaker",
    package_version = as.character(utils::packageVersion("stormbreaker")),
    r_version = R.version.string,
    metafor_version = as.character(utils::packageVersion("metafor"))
  ))
  expect_identical(observed$comparisons, comparisons)

  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "candidate - reference", fixed = TRUE)
  expect_match(printed, "paired external prediction error", fixed = TRUE)
  expect_match(printed, "no ranking or model selection", fixed = TRUE)
})

test_that("prediction-error differences keep metric and time semantics", {
  skip_if_not_installed("metafor", minimum_version = "4.0.0")

  log_loss <- list(
    cohort_a = make_sb_external_prediction_error_difference_bootstrap(
      "binary", -0.03, 0.04, metric = "log_loss"
    ),
    cohort_b = make_sb_external_prediction_error_difference_bootstrap(
      "binary", 0.02, 0.05, metric = "log_loss"
    )
  )
  survival <- list(
    cohort_a = make_sb_external_prediction_error_difference_bootstrap(
      "survival", -0.02, 0.03
    ),
    cohort_b = make_sb_external_prediction_error_difference_bootstrap(
      "survival", 0.01, 0.04
    )
  )

  binary_result <- sb_meta_external_prediction_error_difference(
    log_loss,
    candidate_pipeline_id = "candidate_model",
    metric = "log_loss"
  )
  survival_result <- sb_meta_external_prediction_error_difference(
    survival,
    candidate_pipeline_id = "candidate_model"
  )

  expect_identical(binary_result$metric, "log_loss")
  expect_identical(binary_result$metric_target, list(
    has_horizon = FALSE,
    prediction_horizon = NULL,
    time_unit = NULL
  ))
  expect_identical(survival_result$metric, "ipcw_brier_score")
  expect_identical(survival_result$metric_target, list(
    has_horizon = TRUE,
    prediction_horizon = 6,
    time_unit = "month"
  ))
  expect_match(survival_result$estimand, "at 6 month", fixed = TRUE)
  expect_identical(
    survival_result$interpretation$negative,
    "lower_candidate_prediction_error"
  )
})

test_that("prediction-error difference meta rejects invalid target evidence", {
  binary_a <- make_sb_external_prediction_error_difference_bootstrap(
    "binary", -0.03, 0.04
  )
  binary_b <- make_sb_external_prediction_error_difference_bootstrap(
    "binary", 0.02, 0.05
  )

  expect_error(
    sb_meta_external_prediction_error_difference(
      list(only = binary_a),
      candidate_pipeline_id = "candidate_model"
    ),
    class = "stormbreaker_insufficient_external_cohorts"
  )
  expect_error(
    sb_meta_external_prediction_error_difference(
      list(a = binary_a, b = binary_b),
      candidate_pipeline_id = "candidate_model",
      metric = "roc_auc"
    ),
    class = "stormbreaker_bad_external_meta_metric"
  )

  inconsistent_point <- binary_b
  inconsistent_point$point_differences$estimate[[1L]] <- 0.90
  expect_error(
    sb_meta_external_prediction_error_difference(
      list(a = binary_a, b = inconsistent_point),
      candidate_pipeline_id = "candidate_model"
    ),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  inconsistent_replicates <- binary_b
  inconsistent_replicates$differences$estimate <- c(0.20, 0.40)
  expect_error(
    sb_meta_external_prediction_error_difference(
      list(a = binary_a, b = inconsistent_replicates),
      candidate_pipeline_id = "candidate_model"
    ),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  incomplete <- binary_b
  incomplete$differences$status[[1L]] <- "not_evaluable"
  incomplete$differences$estimate[[1L]] <- NA_real_
  expect_error(
    sb_meta_external_prediction_error_difference(
      list(a = binary_a, b = incomplete),
      candidate_pipeline_id = "candidate_model"
    ),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  zero_se <- make_sb_external_prediction_error_difference_bootstrap(
    "binary", 0.02, 0
  )
  expect_error(
    sb_meta_external_prediction_error_difference(
      list(a = binary_a, b = zero_se),
      candidate_pipeline_id = "candidate_model"
    ),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  survival_a <- make_sb_external_prediction_error_difference_bootstrap(
    "survival", -0.02, 0.03
  )
  survival_b <- make_sb_external_prediction_error_difference_bootstrap(
    "survival", 0.01, 0.04
  )
  survival_b$prediction_target$prediction_horizon <- 12
  expect_error(
    sb_meta_external_prediction_error_difference(
      list(a = survival_a, b = survival_b),
      candidate_pipeline_id = "candidate_model"
    ),
    class = "stormbreaker_external_meta_target_mismatch"
  )
})
