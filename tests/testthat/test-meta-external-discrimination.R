make_sb_external_discrimination_bootstrap <- function(
  task_type,
  point,
  logit_se,
  pipeline_id = "signature__published_model",
  learner = "published_linear_signature",
  n = 100L,
  unrelated_failure = FALSE,
  metric = NULL,
  prediction_target = NULL
) {
  if (is.null(metric)) {
    metric <- if (identical(task_type, "binary")) "roc_auc" else "harrell_c"
  }
  transformed <- stats::qlogis(point) +
    c(-1, 1) * logit_se / sqrt(2)
  point_estimates <- data.frame(
    metric = metric,
    estimate = point,
    status = "success",
    error_class = NA_character_,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
  replicates <- data.frame(
    replicate = 1:2,
    metric = metric,
    estimate = stats::plogis(transformed),
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

make_sb_external_meta_horizon_target <- function(
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

test_that("binary discrimination meta-analysis has a fixed numerical estimand", {
  skip_if_not_installed("metafor", minimum_version = "4.0.0")

  points <- c(0.62, 0.71, 0.79, 0.68)
  logit_ses <- c(0.18, 0.15, 0.22, 0.17)
  bootstraps <- Map(
    function(point, logit_se, index) {
      make_sb_external_discrimination_bootstrap(
        "binary",
        point,
        logit_se,
        n = 80L + index,
        unrelated_failure = index == 2L
      )
    },
    points,
    logit_ses,
    seq_along(points)
  )
  names(bootstraps) <- paste0("cohort_", letters[1:4])

  observed <- sb_meta_external_discrimination(bootstraps)

  expect_s3_class(observed, "stormbreaker_external_discrimination_meta")
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(observed$metric, "roc_auc")
  expect_identical(observed$cohorts$cohort_id, names(bootstraps))
  expect_equal(observed$cohorts$point_estimate, points, tolerance = 0)
  expect_equal(observed$cohorts$logit_estimate, qlogis(points), tolerance = 1e-12)
  expect_equal(observed$cohorts$logit_se, logit_ses, tolerance = 1e-12)
  expect_identical(observed$cohorts$n_bootstrap, rep(2L, 4L))
  expect_equal(observed$summary$estimate, 0.7006585, tolerance = 1e-7)
  expect_equal(observed$summary$lower, 0.5786782, tolerance = 1e-7)
  expect_equal(observed$summary$upper, 0.7995560, tolerance = 1e-7)
  expect_equal(observed$summary$prediction_lower, 0.4614810, tolerance = 1e-7)
  expect_equal(observed$summary$prediction_upper, 0.8647427, tolerance = 1e-7)
  expect_equal(observed$summary$tau2_logit, 0.0716311, tolerance = 1e-7)
  expect_equal(observed$summary$tau_logit, 0.2676399, tolerance = 1e-6)
  expect_equal(observed$summary$i2_percent, 69.5609163, tolerance = 1e-7)
  expect_identical(observed$method$tau2_estimator, "REML")
  expect_identical(observed$method$inference, "Knapp-Hartung")
  expect_identical(observed$provenance, list(
    package = "stormbreaker",
    package_version = as.character(utils::packageVersion("stormbreaker")),
    r_version = R.version.string,
    metafor_version = as.character(utils::packageVersion("metafor"))
  ))
  expect_identical(
    observed$estimand,
    paste0(
      "Back-transformed random-effects mean of logit ROC AUC for one ",
      "locked model across an assumed distribution of independent ",
      "exchangeable cohorts"
    )
  )
  expect_identical(observed$bootstraps, bootstraps)

  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "fixed-model external discrimination", fixed = TRUE)
  expect_match(printed, "no model ranking or selection", fixed = TRUE)
  expect_match(printed, "caller responsibilities", fixed = TRUE)
})

test_that("survival discrimination uses horizon-free Harrell C", {
  skip_if_not_installed("metafor", minimum_version = "4.0.0")

  bootstraps <- list(
    cohort_a = make_sb_external_discrimination_bootstrap(
      "survival", 0.61, 0.20
    ),
    cohort_b = make_sb_external_discrimination_bootstrap(
      "survival", 0.69, 0.16
    )
  )

  observed <- sb_meta_external_discrimination(
    bootstraps,
    confidence_level = 0.90
  )

  expect_identical(observed$task_type, "survival")
  expect_identical(observed$metric, "harrell_c")
  expect_identical(observed$confidence_level, 0.90)
  expect_identical(observed$summary$k, 2L)
  expect_true(observed$summary$estimate > 0 && observed$summary$estimate < 1)
  expect_match(observed$prediction_interval_caution, "Only 2 cohorts", fixed = TRUE)
  expect_match(
    observed$estimand,
    "mean of logit horizon-free Harrell concordance",
    fixed = TRUE
  )
})

test_that("fixed-horizon discrimination keeps one explicit time target", {
  skip_if_not_installed("metafor", minimum_version = "4.0.0")

  points <- c(0.62, 0.71, 0.79, 0.68)
  logit_ses <- c(0.18, 0.15, 0.22, 0.17)
  bootstraps <- Map(function(point, logit_se) {
    make_sb_external_discrimination_bootstrap(
      "survival",
      point,
      logit_se,
      metric = "time_auc",
      prediction_target = make_sb_external_meta_horizon_target()
    )
  }, points, logit_ses)
  names(bootstraps) <- paste0("cohort_", letters[1:4])
  bootstraps[[2L]]$prediction_target$horizon_source <- "prediction"

  observed <- sb_meta_external_discrimination(
    bootstraps,
    metric = "time_auc"
  )

  expect_identical(observed$metric, "time_auc")
  expect_identical(observed$metric_target, list(
    has_horizon = TRUE,
    prediction_horizon = 2,
    time_unit = "year"
  ))
  expect_equal(observed$summary$estimate, 0.7006585, tolerance = 1e-7)
  expect_identical(
    observed$estimand,
    paste0(
      "Back-transformed random-effects mean of logit cumulative/dynamic ",
      "time-dependent AUC at 2 year across an assumed distribution of ",
      "independent exchangeable cohorts"
    )
  )
  expect_match(
    paste(capture.output(print(observed)), collapse = "\n"),
    "Target: 2 year",
    fixed = TRUE
  )

  uno_bootstraps <- lapply(bootstraps[1:2], function(bootstrap) {
    bootstrap$point_estimates$metric <- "uno_c"
    bootstrap$replicates$metric <- "uno_c"
    bootstrap
  })
  uno <- sb_meta_external_discrimination(uno_bootstraps, metric = "uno_c")
  expect_identical(uno$metric, "uno_c")
  expect_identical(uno$metric_target, observed$metric_target)
  expect_identical(
    uno$estimand,
    paste0(
      "Back-transformed random-effects mean of logit truncated Uno ",
      "concordance at 2 year across an assumed distribution of independent ",
      "exchangeable cohorts"
    )
  )
})

test_that("discrimination meta-analysis rejects unlike or incomplete inputs", {
  binary_a <- make_sb_external_discrimination_bootstrap("binary", 0.62, 0.18)
  binary_b <- make_sb_external_discrimination_bootstrap("binary", 0.71, 0.15)
  survival <- make_sb_external_discrimination_bootstrap("survival", 0.65, 0.19)

  expect_error(
    sb_meta_external_discrimination(list(only = binary_a)),
    class = "stormbreaker_insufficient_external_cohorts"
  )
  expect_error(
    sb_meta_external_discrimination(list(binary_a, binary_b)),
    class = "stormbreaker_bad_external_meta"
  )
  expect_error(
    sb_meta_external_discrimination(list(a = binary_a, b = survival)),
    class = "stormbreaker_external_meta_task_mismatch"
  )
  expect_error(
    sb_meta_external_discrimination(
      list(a = binary_a, b = binary_b),
      metric = "time_auc"
    ),
    class = "stormbreaker_bad_external_meta_metric"
  )

  different_model <- binary_b
  different_model$pipeline_id <- "signature__different_model"
  expect_error(
    sb_meta_external_discrimination(list(a = binary_a, b = different_model)),
    class = "stormbreaker_external_meta_model_mismatch"
  )

  incomplete <- binary_b
  incomplete$replicates$status[[1L]] <- "not_evaluable"
  incomplete$replicates$estimate[[1L]] <- NA_real_
  expect_error(
    sb_meta_external_discrimination(list(a = binary_a, b = incomplete)),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  duplicate_replicate <- binary_b
  duplicate_replicate$replicates$replicate <- c(1L, 1L)
  expect_error(
    sb_meta_external_discrimination(list(a = binary_a, b = duplicate_replicate)),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  boundary <- binary_b
  boundary$replicates$estimate[[1L]] <- 1
  expect_error(
    sb_meta_external_discrimination(list(a = binary_a, b = boundary)),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  zero_se <- make_sb_external_discrimination_bootstrap("binary", 0.71, 0)
  expect_error(
    sb_meta_external_discrimination(list(a = binary_a, b = zero_se)),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  time_a <- make_sb_external_discrimination_bootstrap(
    "survival",
    0.62,
    0.18,
    metric = "time_auc",
    prediction_target = make_sb_external_meta_horizon_target()
  )
  time_missing_target <- make_sb_external_discrimination_bootstrap(
    "survival", 0.71, 0.15, metric = "time_auc"
  )
  expect_error(
    sb_meta_external_discrimination(
      list(a = time_a, b = time_missing_target),
      metric = "time_auc"
    ),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  time_risk_only <- make_sb_external_discrimination_bootstrap(
    "survival",
    0.71,
    0.15,
    metric = "time_auc",
    prediction_target = list(
      target_status = "complete",
      has_horizon = FALSE,
      prediction_horizon = NULL,
      time_unit = "year",
      horizon_source = "none"
    )
  )
  expect_error(
    sb_meta_external_discrimination(
      list(a = time_a, b = time_risk_only),
      metric = "time_auc"
    ),
    class = "stormbreaker_external_meta_target_unavailable"
  )

  time_other_horizon <- make_sb_external_discrimination_bootstrap(
    "survival",
    0.71,
    0.15,
    metric = "time_auc",
    prediction_target = make_sb_external_meta_horizon_target(
      prediction_horizon = 3
    )
  )
  expect_error(
    sb_meta_external_discrimination(
      list(a = time_a, b = time_other_horizon),
      metric = "time_auc"
    ),
    class = "stormbreaker_external_meta_target_mismatch"
  )

  time_other_unit <- make_sb_external_discrimination_bootstrap(
    "survival",
    0.71,
    0.15,
    metric = "time_auc",
    prediction_target = make_sb_external_meta_horizon_target(
      time_unit = "month"
    )
  )
  expect_error(
    sb_meta_external_discrimination(
      list(a = time_a, b = time_other_unit),
      metric = "time_auc"
    ),
    class = "stormbreaker_external_meta_target_mismatch"
  )
})
