make_sb_external_mean_calibration_bootstrap <- function(
    task_type,
    point,
    log_se,
    pipeline_id = "signature__published_model",
    learner = "published_linear_signature",
    n = 100L,
    unrelated_failure = FALSE,
    prediction_target = NULL) {
  transformed <- log(point) + c(-1, 1) * log_se / sqrt(2)
  point_estimates <- data.frame(
    metric = "mean_calibration_ratio",
    estimate = point,
    status = "success",
    error_class = NA_character_,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
  replicates <- data.frame(
    replicate = 1:2,
    metric = "mean_calibration_ratio",
    estimate = exp(transformed),
    status = "success",
    error_class = NA_character_,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
  if (unrelated_failure) {
    point_estimates <- rbind(
      point_estimates,
      data.frame(
        metric = "unrelated_metric",
        estimate = NA_real_,
        status = "metric_failed",
        error_class = "example_failure",
        error_message = "not used by this estimand",
        stringsAsFactors = FALSE
      )
    )
    replicates <- rbind(
      replicates,
      data.frame(
        replicate = 1:2,
        metric = "unrelated_metric",
        estimate = NA_real_,
        status = "metric_failed",
        error_class = "example_failure",
        error_message = "not used by this estimand",
        stringsAsFactors = FALSE
      )
    )
  }
  result <- list(
    status = if (unrelated_failure) "incomplete" else "success",
    task_type = task_type,
    pipeline_id = pipeline_id,
    learner = learner,
    n = as.integer(n),
    reps = 2L,
    point_estimates = point_estimates,
    replicates = replicates
  )
  if (!is.null(prediction_target)) {
    result$prediction_target <- prediction_target
  }
  structure(result, class = "stormbreaker_external_bootstrap")
}

make_sb_external_mean_calibration_target <- function(
    prediction_horizon = 2,
    time_unit = "year",
    horizon_source = "endpoint") {
  list(
    target_status = "complete",
    has_horizon = TRUE,
    prediction_horizon = prediction_horizon,
    time_unit = time_unit,
    horizon_source = horizon_source
  )
}

test_that("binary mean calibration meta-analysis has a fixed log-ratio estimand", {
  skip_if_not_installed("metafor", minimum_version = "4.0.0")

  points <- c(0.55, 0.80, 1.25, 1.80)
  log_ses <- c(0.12, 0.15, 0.18, 0.14)
  bootstraps <- Map(
    function(point, log_se, index) {
      make_sb_external_mean_calibration_bootstrap(
        "binary",
        point,
        log_se,
        n = 80L + index,
        unrelated_failure = index == 2L
      )
    },
    points,
    log_ses,
    seq_along(points)
  )
  names(bootstraps) <- paste0("cohort_", letters[1:4])

  observed <- sb_meta_external_mean_calibration(bootstraps)

  expect_s3_class(
    observed,
    "stormbreaker_external_mean_calibration_meta"
  )
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(observed$metric, "mean_calibration_ratio")
  expect_identical(observed$metric_target, list(
    has_horizon = FALSE,
    prediction_horizon = NULL,
    time_unit = NULL
  ))
  expect_identical(observed$cohorts$cohort_id, names(bootstraps))
  expect_equal(observed$cohorts$point_estimate, points, tolerance = 0)
  expect_equal(observed$cohorts$log_estimate, log(points), tolerance = 1e-12)
  expect_equal(observed$cohorts$log_se, log_ses, tolerance = 1e-12)
  expect_equal(observed$summary$estimate, 0.992628872121, tolerance = 1e-10)
  expect_equal(observed$summary$lower, 0.433378416291, tolerance = 1e-10)
  expect_equal(observed$summary$upper, 2.273560566770, tolerance = 1e-10)
  expect_equal(
    observed$summary$prediction_lower,
    0.163615219434,
    tolerance = 1e-10
  )
  expect_equal(
    observed$summary$prediction_upper,
    6.022129733259,
    tolerance = 1e-10
  )
  expect_equal(observed$summary$tau2_log, 0.253102260013, tolerance = 1e-10)
  expect_equal(observed$summary$tau_log, 0.503092695249, tolerance = 1e-10)
  expect_equal(
    observed$summary$i2_percent,
    92.340168130849,
    tolerance = 1e-10
  )
  expect_identical(observed$method$tau2_estimator, "REML")
  expect_identical(observed$method$inference, "Knapp-Hartung")
  expect_identical(observed$method$back_transformation, "exponential")
  expect_identical(observed$provenance, list(
    package = "stormbreaker",
    package_version = as.character(utils::packageVersion("stormbreaker")),
    r_version = R.version.string,
    metafor_version = as.character(utils::packageVersion("metafor"))
  ))
  expect_identical(
    observed$estimand,
    paste0(
      "Exponentiated random-effects mean of log observed-to-expected event ",
      "ratio for one locked binary model across an assumed distribution of ",
      "independent exchangeable cohorts"
    )
  )
  expect_identical(observed$bootstraps, bootstraps)

  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "fixed-model external mean calibration", fixed = TRUE)
  expect_match(printed, "O:E below 1 means overprediction", fixed = TRUE)
  expect_match(printed, "does not pool calibration slope", fixed = TRUE)
})

test_that("survival mean calibration meta-analysis keeps one time target", {
  skip_if_not_installed("metafor", minimum_version = "4.0.0")

  bootstraps <- list(
    cohort_a = make_sb_external_mean_calibration_bootstrap(
      "survival",
      0.75,
      0.20,
      prediction_target = make_sb_external_mean_calibration_target()
    ),
    cohort_b = make_sb_external_mean_calibration_bootstrap(
      "survival",
      1.10,
      0.16,
      prediction_target = make_sb_external_mean_calibration_target(
        horizon_source = "prediction"
      )
    )
  )

  observed <- sb_meta_external_mean_calibration(
    bootstraps,
    confidence_level = 0.90
  )

  expect_identical(observed$task_type, "survival")
  expect_identical(observed$confidence_level, 0.90)
  expect_identical(observed$metric_target, list(
    has_horizon = TRUE,
    prediction_horizon = 2,
    time_unit = "year"
  ))
  expect_match(observed$estimand, "at 2 year", fixed = TRUE)
  expect_match(
    paste(capture.output(print(observed)), collapse = "\n"),
    "Target: 2 year",
    fixed = TRUE
  )
})

test_that("mean calibration meta-analysis rejects unlike or incomplete inputs", {
  binary_a <- make_sb_external_mean_calibration_bootstrap(
    "binary", 0.8, 0.18
  )
  binary_b <- make_sb_external_mean_calibration_bootstrap(
    "binary", 1.1, 0.15
  )
  survival <- make_sb_external_mean_calibration_bootstrap(
    "survival",
    0.9,
    0.19,
    prediction_target = make_sb_external_mean_calibration_target()
  )

  expect_error(
    sb_meta_external_mean_calibration(list(only = binary_a)),
    class = "stormbreaker_insufficient_external_cohorts"
  )
  expect_error(
    sb_meta_external_mean_calibration(list(binary_a, binary_b)),
    class = "stormbreaker_bad_external_meta"
  )
  expect_error(
    sb_meta_external_mean_calibration(list(a = binary_a, b = survival)),
    class = "stormbreaker_external_meta_task_mismatch"
  )

  different_model <- binary_b
  different_model$learner <- "different_learner"
  expect_error(
    sb_meta_external_mean_calibration(
      list(a = binary_a, b = different_model)
    ),
    class = "stormbreaker_external_meta_model_mismatch"
  )

  incomplete <- binary_b
  incomplete$replicates$status[[1L]] <- "not_evaluable"
  incomplete$replicates$estimate[[1L]] <- NA_real_
  expect_error(
    sb_meta_external_mean_calibration(list(a = binary_a, b = incomplete)),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  duplicate_replicate <- binary_b
  duplicate_replicate$replicates$replicate <- c(1L, 1L)
  expect_error(
    sb_meta_external_mean_calibration(
      list(a = binary_a, b = duplicate_replicate)
    ),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  zero <- binary_b
  zero$point_estimates$estimate <- 0
  expect_error(
    sb_meta_external_mean_calibration(list(a = binary_a, b = zero)),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  zero_se <- make_sb_external_mean_calibration_bootstrap("binary", 1.1, 0)
  expect_error(
    sb_meta_external_mean_calibration(list(a = binary_a, b = zero_se)),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  survival_missing_target <- make_sb_external_mean_calibration_bootstrap(
    "survival", 1.05, 0.17
  )
  expect_error(
    sb_meta_external_mean_calibration(
      list(a = survival, b = survival_missing_target)
    ),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  survival_other_horizon <- make_sb_external_mean_calibration_bootstrap(
    "survival",
    1.05,
    0.17,
    prediction_target = make_sb_external_mean_calibration_target(
      prediction_horizon = 3
    )
  )
  expect_error(
    sb_meta_external_mean_calibration(
      list(a = survival, b = survival_other_horizon)
    ),
    class = "stormbreaker_external_meta_target_mismatch"
  )
})
