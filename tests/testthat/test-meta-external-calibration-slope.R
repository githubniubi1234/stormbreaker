make_sb_external_calibration_slope_bootstrap <- function(
    task_type,
    point,
    slope_se,
    pipeline_id = "signature__published_model",
    learner = "published_linear_signature",
    n = 100L,
    unrelated_failure = FALSE,
    prediction_target = NULL) {
  replicate_estimates <- point + c(-1, 1) * slope_se / sqrt(2)
  point_estimates <- data.frame(
    metric = "calibration_slope",
    estimate = point,
    status = "success",
    error_class = NA_character_,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
  replicates <- data.frame(
    replicate = 1:2,
    metric = "calibration_slope",
    estimate = replicate_estimates,
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

make_sb_external_calibration_slope_target <- function(
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

test_that("binary calibration-slope meta-analysis uses the original scale", {
  skip_if_not_installed("metafor", minimum_version = "4.0.0")

  points <- c(0.55, 0.85, 1.15, 1.45)
  slope_ses <- c(0.12, 0.16, 0.18, 0.14)
  bootstraps <- Map(
    function(point, slope_se, index) {
      make_sb_external_calibration_slope_bootstrap(
        "binary",
        point,
        slope_se,
        n = 80L + index,
        unrelated_failure = index == 2L
      )
    },
    points,
    slope_ses,
    seq_along(points)
  )
  names(bootstraps) <- paste0("cohort_", letters[1:4])

  observed <- sb_meta_external_calibration_slope(bootstraps)

  expect_s3_class(
    observed,
    "stormbreaker_external_calibration_slope_meta"
  )
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(observed$metric, "calibration_slope")
  expect_identical(observed$metric_target, list(
    has_horizon = FALSE,
    prediction_horizon = NULL,
    time_unit = NULL
  ))
  expect_identical(observed$cohorts$cohort_id, names(bootstraps))
  expect_equal(observed$cohorts$point_estimate, points, tolerance = 0)
  expect_equal(observed$cohorts$standard_error, slope_ses, tolerance = 1e-12)
  expect_equal(observed$summary$estimate, 0.994563728454, tolerance = 1e-10)
  expect_equal(observed$summary$lower, 0.369178754096, tolerance = 1e-10)
  expect_equal(observed$summary$upper, 1.619948702813, tolerance = 1e-10)
  expect_equal(
    observed$summary$prediction_lower,
    -0.337282509459,
    tolerance = 1e-10
  )
  expect_equal(
    observed$summary$prediction_upper,
    2.326409966368,
    tolerance = 1e-10
  )
  expect_equal(observed$summary$tau2, 0.136523783945, tolerance = 1e-10)
  expect_equal(observed$summary$tau, 0.369491250160, tolerance = 1e-10)
  expect_equal(
    observed$summary$i2_percent,
    86.295777724312,
    tolerance = 1e-10
  )
  expect_identical(observed$method$effect_scale, "original_calibration_slope")
  expect_identical(observed$method$tau2_estimator, "REML")
  expect_identical(observed$method$inference, "Knapp-Hartung")
  expect_identical(observed$method$back_transformation, "none")
  expect_identical(observed$provenance, list(
    package = "stormbreaker",
    package_version = as.character(utils::packageVersion("stormbreaker")),
    r_version = R.version.string,
    metafor_version = as.character(utils::packageVersion("metafor"))
  ))
  expect_identical(
    observed$estimand,
    paste0(
      "Random-effects mean calibration slope for one locked binary model ",
      "across an assumed distribution of independent exchangeable cohorts"
    )
  )
  expect_identical(observed$bootstraps, bootstraps)

  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "original-scale calibration-slope", fixed = TRUE)
  expect_match(
    printed,
    "does not establish mean or moderate calibration",
    fixed = TRUE
  )
  expect_match(printed, "no model ranking or recalibration", fixed = TRUE)
})

test_that("calibration-slope meta-analysis retains finite negative slopes", {
  skip_if_not_installed("metafor", minimum_version = "4.0.0")

  observed <- sb_meta_external_calibration_slope(list(
    cohort_a = make_sb_external_calibration_slope_bootstrap(
      "binary", -0.20, 0.10
    ),
    cohort_b = make_sb_external_calibration_slope_bootstrap(
      "binary", 0.40, 0.12
    )
  ))

  expect_identical(observed$cohorts$point_estimate, c(-0.20, 0.40))
  expect_identical(is.finite(observed$summary$estimate), TRUE)
})

test_that("survival calibration-slope meta-analysis keeps one time target", {
  skip_if_not_installed("metafor", minimum_version = "4.0.0")

  bootstraps <- list(
    cohort_a = make_sb_external_calibration_slope_bootstrap(
      "survival",
      0.75,
      0.20,
      prediction_target = make_sb_external_calibration_slope_target()
    ),
    cohort_b = make_sb_external_calibration_slope_bootstrap(
      "survival",
      1.10,
      0.16,
      prediction_target = make_sb_external_calibration_slope_target(
        horizon_source = "prediction"
      )
    )
  )

  observed <- sb_meta_external_calibration_slope(
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

test_that("calibration-slope meta-analysis rejects unlike or incomplete inputs", {
  binary_a <- make_sb_external_calibration_slope_bootstrap(
    "binary", 0.8, 0.18
  )
  binary_b <- make_sb_external_calibration_slope_bootstrap(
    "binary", 1.1, 0.15
  )
  survival <- make_sb_external_calibration_slope_bootstrap(
    "survival",
    0.9,
    0.19,
    prediction_target = make_sb_external_calibration_slope_target()
  )

  expect_error(
    sb_meta_external_calibration_slope(list(only = binary_a)),
    class = "stormbreaker_insufficient_external_cohorts"
  )
  expect_error(
    sb_meta_external_calibration_slope(list(binary_a, binary_b)),
    class = "stormbreaker_bad_external_meta"
  )
  expect_error(
    sb_meta_external_calibration_slope(list(a = binary_a, b = survival)),
    class = "stormbreaker_external_meta_task_mismatch"
  )

  different_model <- binary_b
  different_model$pipeline_id <- "signature__different_model"
  expect_error(
    sb_meta_external_calibration_slope(
      list(a = binary_a, b = different_model)
    ),
    class = "stormbreaker_external_meta_model_mismatch"
  )

  incomplete <- binary_b
  incomplete$replicates$status[[1L]] <- "not_evaluable"
  incomplete$replicates$estimate[[1L]] <- NA_real_
  expect_error(
    sb_meta_external_calibration_slope(list(a = binary_a, b = incomplete)),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  duplicate_replicate <- binary_b
  duplicate_replicate$replicates$replicate <- c(1L, 1L)
  expect_error(
    sb_meta_external_calibration_slope(
      list(a = binary_a, b = duplicate_replicate)
    ),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  nonfinite <- binary_b
  nonfinite$point_estimates$estimate <- Inf
  expect_error(
    sb_meta_external_calibration_slope(list(a = binary_a, b = nonfinite)),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  zero_se <- make_sb_external_calibration_slope_bootstrap("binary", 1.1, 0)
  expect_error(
    sb_meta_external_calibration_slope(list(a = binary_a, b = zero_se)),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  survival_missing_target <- make_sb_external_calibration_slope_bootstrap(
    "survival", 1.05, 0.17
  )
  expect_error(
    sb_meta_external_calibration_slope(
      list(a = survival, b = survival_missing_target)
    ),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  survival_other_unit <- make_sb_external_calibration_slope_bootstrap(
    "survival",
    1.05,
    0.17,
    prediction_target = make_sb_external_calibration_slope_target(
      time_unit = "month"
    )
  )
  expect_error(
    sb_meta_external_calibration_slope(
      list(a = survival, b = survival_other_unit)
    ),
    class = "stormbreaker_external_meta_target_mismatch"
  )
})
