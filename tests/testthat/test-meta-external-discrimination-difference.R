make_sb_external_discrimination_difference_bootstrap <- function(
    task_type,
    point,
    standard_error,
    metric = NULL,
    n = 100L,
    reference_pipeline_id = "reference_model",
    candidate_pipeline_id = "candidate_model",
    prediction_target = NULL,
    unrelated_failure = FALSE) {
  if (is.null(metric)) {
    metric <- if (identical(task_type, "binary")) "roc_auc" else "harrell_c"
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
        has_horizon = FALSE,
        prediction_horizon = NULL,
        time_unit = "month",
        horizon_source = "none"
      )
    }
  }
  reference_estimate <- 0.60
  replicate_values <- point + c(-1, 1) * standard_error / sqrt(2)
  candidate_point_estimate <- reference_estimate + point
  point_difference <- candidate_point_estimate - reference_estimate
  candidate_replicate_estimates <- reference_estimate + replicate_values
  replicate_differences <-
    candidate_replicate_estimates - reference_estimate
  candidate_status <- data.frame(
    candidate_pipeline_id = candidate_pipeline_id,
    status = "success",
    error_class = NA_character_,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
  status <- "success"
  if (unrelated_failure) {
    candidate_status <- rbind(
      candidate_status,
      data.frame(
        candidate_pipeline_id = "unrelated_failed_model",
        status = "failed",
        error_class = "stormbreaker_test_failure",
        error_message = "Unrelated candidate failed.",
        stringsAsFactors = FALSE
      )
    )
    status <- "incomplete"
  }
  result <- list(
    status = status,
    task_type = task_type,
    reference_pipeline_id = reference_pipeline_id,
    candidate_status = candidate_status,
    n = as.integer(n),
    prediction_target = prediction_target,
    reps = 2L,
    metrics = metric,
    point_differences = data.frame(
      reference_pipeline_id = reference_pipeline_id,
      candidate_pipeline_id = candidate_pipeline_id,
      contrast = "candidate_minus_reference",
      metric = metric,
      metric_direction = "higher_is_better",
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
      metric_direction = "higher_is_better",
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

make_sb_external_difference_horizon_target <- function(
    prediction_horizon = 6,
    time_unit = "month") {
  list(
    target_status = "complete",
    has_horizon = TRUE,
    prediction_horizon = prediction_horizon,
    time_unit = time_unit,
    horizon_source = "endpoint"
  )
}

test_that("paired discrimination differences have a fixed numerical estimand", {
  skip_if_not_installed("metafor", minimum_version = "4.0.0")

  points <- c(-0.05, 0.02, 0.08, 0.12)
  standard_errors <- c(0.03, 0.04, 0.025, 0.035)
  comparisons <- Map(
    function(point, standard_error, index) {
      make_sb_external_discrimination_difference_bootstrap(
        "binary",
        point,
        standard_error,
        n = 90L + index,
        unrelated_failure = index == 2L
      )
    },
    points,
    standard_errors,
    seq_along(points)
  )
  names(comparisons) <- paste0("cohort_", letters[1:4])

  observed <- sb_meta_external_discrimination_difference(
    comparisons,
    candidate_pipeline_id = "candidate_model"
  )

  expect_s3_class(
    observed,
    "stormbreaker_external_discrimination_difference_meta"
  )
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(observed$reference_pipeline_id, "reference_model")
  expect_identical(observed$candidate_pipeline_id, "candidate_model")
  expect_identical(observed$contrast, "candidate_minus_reference")
  expect_identical(observed$metric, "roc_auc")
  expect_identical(observed$cohorts$cohort_id, names(comparisons))
  expect_equal(observed$cohorts$point_difference, points, tolerance = 1e-15)
  expect_equal(
    observed$cohorts$standard_error,
    standard_errors,
    tolerance = 1e-12
  )
  expect_identical(observed$cohorts$n_bootstrap, rep(2L, 4L))
  expect_equal(observed$summary$estimate, 0.042549460390, tolerance = 1e-10)
  expect_equal(observed$summary$lower, -0.075960362323, tolerance = 1e-10)
  expect_equal(observed$summary$upper, 0.161059283102, tolerance = 1e-10)
  expect_equal(
    observed$summary$prediction_lower,
    -0.202746184516,
    tolerance = 1e-10
  )
  expect_equal(
    observed$summary$prediction_upper,
    0.287845105296,
    tolerance = 1e-10
  )
  expect_equal(observed$summary$tau2, 0.004554259189, tolerance = 1e-10)
  expect_equal(observed$summary$tau, 0.067485251644, tolerance = 1e-10)
  expect_equal(
    observed$summary$i2_percent,
    81.917056396541,
    tolerance = 1e-10
  )
  expect_identical(
    observed$method$effect_scale,
    "raw_candidate_minus_reference_discrimination"
  )
  expect_identical(observed$method$tau2_estimator, "REML")
  expect_identical(observed$method$inference, "Knapp-Hartung")
  expect_identical(observed$method$back_transformation, "none")
  expect_identical(observed$provenance, list(
    package = "stormbreaker",
    package_version = as.character(utils::packageVersion("stormbreaker")),
    r_version = R.version.string,
    metafor_version = as.character(utils::packageVersion("metafor"))
  ))
  expect_identical(observed$comparisons, comparisons)

  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "candidate - reference", fixed = TRUE)
  expect_match(printed, "paired external discrimination", fixed = TRUE)
  expect_match(printed, "no ranking or model selection", fixed = TRUE)
})

test_that("survival discrimination differences keep metric time semantics", {
  skip_if_not_installed("metafor", minimum_version = "4.0.0")

  harrell <- list(
    cohort_a = make_sb_external_discrimination_difference_bootstrap(
      "survival", 0.03, 0.04
    ),
    cohort_b = make_sb_external_discrimination_difference_bootstrap(
      "survival", 0.08, 0.05
    )
  )
  horizon <- make_sb_external_difference_horizon_target()
  time_auc <- list(
    cohort_a = make_sb_external_discrimination_difference_bootstrap(
      "survival", 0.04, 0.03, metric = "time_auc",
      prediction_target = horizon
    ),
    cohort_b = make_sb_external_discrimination_difference_bootstrap(
      "survival", 0.09, 0.04, metric = "time_auc",
      prediction_target = horizon
    )
  )

  horizon_free <- sb_meta_external_discrimination_difference(
    harrell,
    candidate_pipeline_id = "candidate_model"
  )
  fixed_horizon <- sb_meta_external_discrimination_difference(
    time_auc,
    candidate_pipeline_id = "candidate_model",
    metric = "time_auc"
  )

  expect_identical(horizon_free$metric, "harrell_c")
  expect_identical(horizon_free$metric_target, list(
    has_horizon = FALSE,
    prediction_horizon = NULL,
    time_unit = NULL
  ))
  expect_identical(fixed_horizon$metric, "time_auc")
  expect_identical(fixed_horizon$metric_target, list(
    has_horizon = TRUE,
    prediction_horizon = 6,
    time_unit = "month"
  ))
  expect_match(fixed_horizon$estimand, "at 6 month", fixed = TRUE)
})

test_that("difference meta-analysis rejects unlike or incomplete inputs", {
  binary_a <- make_sb_external_discrimination_difference_bootstrap(
    "binary", 0.03, 0.04
  )
  binary_b <- make_sb_external_discrimination_difference_bootstrap(
    "binary", 0.08, 0.05
  )
  survival <- make_sb_external_discrimination_difference_bootstrap(
    "survival", 0.06, 0.05
  )

  expect_error(
    sb_meta_external_discrimination_difference(
      list(only = binary_a),
      candidate_pipeline_id = "candidate_model"
    ),
    class = "stormbreaker_insufficient_external_cohorts"
  )
  expect_error(
    sb_meta_external_discrimination_difference(
      list(binary_a, binary_b),
      candidate_pipeline_id = "candidate_model"
    ),
    class = "stormbreaker_bad_external_meta"
  )
  expect_error(
    sb_meta_external_discrimination_difference(
      list(a = binary_a, b = survival),
      candidate_pipeline_id = "candidate_model"
    ),
    class = "stormbreaker_external_meta_task_mismatch"
  )
  expect_error(
    sb_meta_external_discrimination_difference(
      list(a = binary_a, b = binary_b),
      candidate_pipeline_id = "candidate_model",
      metric = "time_auc"
    ),
    class = "stormbreaker_bad_external_meta_metric"
  )

  other_reference <- binary_b
  other_reference$reference_pipeline_id <- "other_reference"
  other_reference$point_differences$reference_pipeline_id <- "other_reference"
  other_reference$differences$reference_pipeline_id <- "other_reference"
  expect_error(
    sb_meta_external_discrimination_difference(
      list(a = binary_a, b = other_reference),
      candidate_pipeline_id = "candidate_model"
    ),
    class = "stormbreaker_external_meta_model_mismatch"
  )

  failed_candidate <- binary_b
  failed_candidate$candidate_status$status[[1L]] <- "failed"
  expect_error(
    sb_meta_external_discrimination_difference(
      list(a = binary_a, b = failed_candidate),
      candidate_pipeline_id = "candidate_model"
    ),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  duplicate_replicate <- binary_b
  duplicate_replicate$differences$replicate <- c(1L, 1L)
  expect_error(
    sb_meta_external_discrimination_difference(
      list(a = binary_a, b = duplicate_replicate),
      candidate_pipeline_id = "candidate_model"
    ),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  incomplete <- binary_b
  incomplete$differences$status[[1L]] <- "not_evaluable"
  incomplete$differences$estimate[[1L]] <- NA_real_
  expect_error(
    sb_meta_external_discrimination_difference(
      list(a = binary_a, b = incomplete),
      candidate_pipeline_id = "candidate_model"
    ),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  inconsistent_point <- binary_b
  inconsistent_point$point_differences$estimate[[1L]] <- 0.90
  expect_error(
    sb_meta_external_discrimination_difference(
      list(a = binary_a, b = inconsistent_point),
      candidate_pipeline_id = "candidate_model"
    ),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  inconsistent_replicates <- binary_b
  inconsistent_replicates$differences$estimate <- c(0.20, 0.40)
  expect_error(
    sb_meta_external_discrimination_difference(
      list(a = binary_a, b = inconsistent_replicates),
      candidate_pipeline_id = "candidate_model"
    ),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  zero_se <- make_sb_external_discrimination_difference_bootstrap(
    "binary", 0.08, 0
  )
  expect_error(
    sb_meta_external_discrimination_difference(
      list(a = binary_a, b = zero_se),
      candidate_pipeline_id = "candidate_model"
    ),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  time_a <- make_sb_external_discrimination_difference_bootstrap(
    "survival", 0.04, 0.03, metric = "time_auc",
    prediction_target = make_sb_external_difference_horizon_target()
  )
  time_b <- make_sb_external_discrimination_difference_bootstrap(
    "survival", 0.09, 0.04, metric = "time_auc",
    prediction_target = make_sb_external_difference_horizon_target(
      prediction_horizon = 12
    )
  )
  expect_error(
    sb_meta_external_discrimination_difference(
      list(a = time_a, b = time_b),
      candidate_pipeline_id = "candidate_model",
      metric = "time_auc"
    ),
    class = "stormbreaker_external_meta_target_mismatch"
  )
})
