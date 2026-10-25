make_sb_bootstrap_binary_evaluation <- function(
    outcome,
    probability,
    decision_thresholds = NULL) {
  ids <- sprintf("b%02d", seq_along(outcome))
  prediction <- list(
    status = "success",
    task_type = "binary",
    pipeline_id = "none__logistic",
    learner = "logistic",
    id_col = "id",
    predictions = data.frame(
      patient_id = ids,
      pipeline_id = rep("none__logistic", length(ids)),
      learner = rep("logistic", length(ids)),
      positive_class = rep.int(1L, length(ids)),
      probability = probability,
      status = rep("success", length(ids)),
      error_class = rep(NA_character_, length(ids)),
      error_message = rep(NA_character_, length(ids)),
      stringsAsFactors = FALSE
    )
  )
  if (!is.null(decision_thresholds)) {
    prediction$decision_thresholds <- decision_thresholds
  }
  class(prediction) <- "stormbreaker_final_prediction"
  endpoint <- sb_binary_endpoint(
    data.frame(id = ids, outcome = outcome),
    "id",
    "outcome"
  )
  sb_evaluate_external(prediction, endpoint, reverse = FALSE)
}

make_sb_bootstrap_survival_evaluation <- function(
    time,
    event,
    risk,
    survival_probability = NULL,
    prediction_horizon = NULL,
    time_unit = "month",
    decision_thresholds = NULL,
    endpoint_prediction_horizon = prediction_horizon) {
  ids <- sprintf("s%02d", seq_along(time))
  rows <- data.frame(
    patient_id = ids,
    pipeline_id = rep("none__coxph", length(ids)),
    learner = rep("coxph", length(ids)),
    risk = risk,
    stringsAsFactors = FALSE
  )
  if (!is.null(prediction_horizon)) {
    rows$prediction_horizon <- rep(prediction_horizon, length(ids))
    rows$time_unit <- rep(time_unit, length(ids))
    rows$survival_probability <- survival_probability
  }
  rows$status <- rep("success", length(ids))
  rows$error_class <- rep(NA_character_, length(ids))
  rows$error_message <- rep(NA_character_, length(ids))
  prediction <- list(
    status = "success",
    task_type = "survival",
    pipeline_id = "none__coxph",
    learner = "coxph",
    id_col = "id",
    predictions = rows
  )
  if (!is.null(decision_thresholds)) {
    prediction$decision_thresholds <- decision_thresholds
  }
  class(prediction) <- "stormbreaker_final_prediction"
  endpoint <- sb_surv_endpoint(
    data.frame(id = ids, time = time, event = event),
    "id",
    "time",
    "event",
    time_unit,
    prediction_horizon = endpoint_prediction_horizon
  )
  sb_evaluate_external(prediction, endpoint, reverse = FALSE)
}

manual_external_ipcw_brier <- function(
    time,
    event,
    survival_probability,
    horizon) {
  times <- sort(unique(time))
  current <- 1
  survival_after <- numeric(length(times))
  for (index in seq_along(times)) {
    at_risk <- sum(time >= times[[index]])
    observed_events <- sum(time == times[[index]] & event == 1L)
    censoring_events <- sum(time == times[[index]] & event == 0L)
    censoring_risk <- at_risk - observed_events
    if (censoring_events > 0L) {
      current <- current * (1 - censoring_events / censoring_risk)
    }
    survival_after[[index]] <- current
  }
  g_before <- function(value) {
    eligible <- times < value
    if (!any(eligible)) 1 else utils::tail(survival_after[eligible], 1L)
  }
  event_by_horizon <- time <= horizon & event == 1L
  surviving_at_horizon <- time > horizon |
    (time == horizon & event == 0L)
  if (!any(event_by_horizon | surviving_at_horizon)) {
    return(NA_real_)
  }
  event_g <- vapply(time[event_by_horizon], g_before, numeric(1))
  horizon_g <- g_before(horizon)
  if (any(event_g <= 0) ||
      (any(surviving_at_horizon) && horizon_g <= 0)) {
    return(NA_real_)
  }
  weights <- numeric(length(time))
  weights[event_by_horizon] <- 1 / event_g
  weights[surviving_at_horizon] <- 1 / horizon_g
  losses <- numeric(length(time))
  losses[event_by_horizon] <- survival_probability[event_by_horizon]^2
  losses[surviving_at_horizon] <-
    (1 - survival_probability[surviving_at_horizon])^2
  sum(weights * losses) / length(time)
}

manual_external_time_auc <- function(
    time,
    event,
    survival_probability,
    horizon) {
  cases <- event == 1L & time <= horizon
  controls <- time > horizon
  if (!any(cases) || !any(controls)) {
    return(NA_real_)
  }
  times <- sort(unique(time))
  current <- 1
  survival_after <- numeric(length(times))
  for (index in seq_along(times)) {
    at_risk <- sum(time >= times[[index]])
    observed_events <- sum(time == times[[index]] & event == 1L)
    censoring_events <- sum(time == times[[index]] & event == 0L)
    censoring_risk <- at_risk - observed_events
    if (censoring_events > 0L) {
      current <- current * (1 - censoring_events / censoring_risk)
    }
    survival_after[[index]] <- current
  }
  g_before <- function(value) {
    eligible <- times < value
    if (!any(eligible)) 1 else utils::tail(survival_after[eligible], 1L)
  }
  case_g <- vapply(time[cases], g_before, numeric(1))
  if (any(!is.finite(case_g)) || any(case_g <= 0)) {
    return(NA_real_)
  }
  case_weight <- 1 / case_g
  marker <- 1 - survival_probability
  credit <- outer(
    marker[cases],
    marker[controls],
    function(case, control) {
      (case > control) + 0.5 * (case == control)
    }
  )
  sum(case_weight * rowSums(credit)) /
    (sum(case_weight) * sum(controls))
}

manual_external_uno_c <- function(time, event, risk, horizon) {
  times <- sort(unique(time))
  current <- 1
  survival_after <- numeric(length(times))
  for (index in seq_along(times)) {
    at_risk <- sum(time >= times[[index]])
    observed_events <- sum(time == times[[index]] & event == 1L)
    censoring_events <- sum(time == times[[index]] & event == 0L)
    censoring_risk <- at_risk - observed_events
    if (censoring_events > 0L) {
      current <- current * (1 - censoring_events / censoring_risk)
    }
    survival_after[[index]] <- current
  }
  g_before <- function(value) {
    eligible <- times < value
    if (!any(eligible)) 1 else utils::tail(survival_after[eligible], 1L)
  }
  if (!is.finite(g_before(horizon)) || g_before(horizon) <= 0) {
    return(NA_real_)
  }
  eligible <- which(event == 1L & time < horizon)
  if (length(eligible) == 0L) {
    return(NA_real_)
  }
  numerator <- denominator <- 0
  for (index in eligible) {
    event_g <- g_before(time[[index]])
    if (!is.finite(event_g) || event_g <= 0) {
      return(NA_real_)
    }
    weight <- 1 / event_g^2
    comparators <- which(
      time > time[[index]] |
        (time == time[[index]] & event == 0L)
    )
    numerator <- numerator + weight * sum(
      (risk[[index]] > risk[comparators]) +
        0.5 * (risk[[index]] == risk[comparators])
    )
    denominator <- denominator + weight * length(comparators)
  }
  if (!is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }
  numerator / denominator
}

manual_external_mean_calibration_ratio <- function(
    time,
    event,
    survival_probability,
    horizon) {
  event_times <- sort(unique(time[event == 1L & time <= horizon]))
  observed_survival <- prod(vapply(event_times, function(event_time) {
    1 - sum(time == event_time & event == 1L) / sum(time >= event_time)
  }, numeric(1)))
  observed_event_probability <- 1 - observed_survival
  observed_event_probability / mean(-expm1(log(survival_probability)))
}

manual_external_calibration_slope <- function(
    time,
    event,
    survival_probability,
    horizon) {
  transformed_prediction <- log(-log(survival_probability))
  fit <- survival::coxph(
    survival::Surv(
      pmin(time, horizon),
      as.integer(event == 1L & time <= horizon)
    ) ~ transformed_prediction,
    ties = "efron"
  )
  unname(stats::coef(fit)[["transformed_prediction"]])
}

make_sb_bootstrap_calibration_evaluation <- function(
    n = 80L,
    seed = 20260820L,
    horizon = 2) {
  generated <- withr::with_seed(seed, {
    prognostic_index <- stats::rnorm(n)
    event_time <- stats::rexp(
      n,
      rate = 0.15 * exp(0.45 * prognostic_index)
    )
    censoring_time <- stats::rexp(n, rate = 0.05)
    data.frame(
      time = pmin(event_time, censoring_time),
      event = as.integer(event_time <= censoring_time),
      risk = prognostic_index,
      survival_probability = exp(
        -0.15 * horizon * exp(0.45 * prognostic_index)
      )
    )
  })
  make_sb_bootstrap_survival_evaluation(
    generated$time,
    generated$event,
    generated$risk,
    generated$survival_probability,
    horizon,
    "year"
  )
}

manual_binary_auc <- function(outcome, probability) {
  positive <- probability[outcome == 1L]
  negative <- probability[outcome == 0L]
  if (length(positive) == 0L || length(negative) == 0L) {
    return(NA_real_)
  }
  mean(vapply(positive, function(value) {
    mean((value > negative) + 0.5 * (value == negative))
  }, numeric(1)))
}

manual_binary_average_precision <- function(outcome, probability) {
  n_positive <- sum(outcome == 1L)
  n_negative <- sum(outcome == 0L)
  if (n_positive == 0L || n_negative == 0L) {
    return(NA_real_)
  }
  thresholds <- sort(unique(probability), decreasing = TRUE)
  previous_recall <- 0
  total <- 0
  for (threshold in thresholds) {
    predicted_positive <- probability >= threshold
    true_positive <- sum(predicted_positive & outcome == 1L)
    false_positive <- sum(predicted_positive & outcome == 0L)
    recall <- true_positive / n_positive
    precision <- true_positive / (true_positive + false_positive)
    total <- total + (recall - previous_recall) * precision
    previous_recall <- recall
  }
  total
}

expect_external_bootstrap_failure <- function(result, stage, class_name) {
  expect_s3_class(result, "stormbreaker_external_bootstrap")
  expect_identical(result$status, "failed")
  expect_identical(result$error_stage, stage)
  expect_identical(result$error_class, class_name)
  expect_identical(nrow(result$point_estimates), 0L)
  expect_identical(nrow(result$replicates), 0L)
  expect_identical(nrow(result$intervals), 0L)
}

test_that("external bootstrap exposes only evaluation and sampling arguments", {
  expect_identical(
    names(formals(sb_bootstrap_external)),
    c("evaluation", "reps", "confidence_level", "seed")
  )
  expect_identical(formals(sb_bootstrap_external)$reps, 1000L)
  expect_identical(formals(sb_bootstrap_external)$confidence_level, 0.95)
  expect_identical(formals(sb_bootstrap_external)$seed, 1L)
})

test_that("external bootstrap rejects invalid evaluation and arguments explicitly", {
  evaluation <- make_sb_bootstrap_binary_evaluation(
    c(0L, 0L, 1L, 1L),
    c(0.1, 0.3, 0.7, 0.9)
  )
  invalid_reps <- list(1L, c(2L, 3L), NA_real_, Inf, 2.5,
                       .Machine$integer.max + 1)
  for (value in invalid_reps) {
    result <- sb_bootstrap_external(evaluation, reps = value)
    expect_external_bootstrap_failure(
      result,
      "arguments",
      "stormbreaker_bad_integer"
    )
  }
  for (value in list(0, 1, NA_real_, Inf, c(0.9, 0.95))) {
    result <- sb_bootstrap_external(
      evaluation,
      reps = 2L,
      confidence_level = value
    )
    expect_external_bootstrap_failure(
      result,
      "arguments",
      "stormbreaker_bad_confidence_level"
    )
  }
  for (value in list(0L, -1L, 1.5, NA_real_, c(1L, 2L))) {
    result <- sb_bootstrap_external(evaluation, reps = 2L, seed = value)
    expect_external_bootstrap_failure(
      result,
      "arguments",
      "stormbreaker_bad_integer"
    )
  }

  wrong_class <- unclass(evaluation)
  result <- sb_bootstrap_external(wrong_class, reps = 2L)
  expect_external_bootstrap_failure(
    result,
    "evaluation_contract",
    "stormbreaker_bad_external_evaluation"
  )

  failed <- evaluation
  failed$status <- "failed"
  result <- sb_bootstrap_external(failed, reps = 2L)
  expect_external_bootstrap_failure(
    result,
    "evaluation_contract",
    "stormbreaker_failed_external_evaluation"
  )
})

test_that("incomplete external evaluations retain estimable bootstrap targets", {
  evaluation <- make_sb_bootstrap_binary_evaluation(
    rep.int(1L, 4L),
    c(0.6, 0.7, 0.8, 0.9)
  )
  expect_identical(evaluation$status, "incomplete")

  result <- sb_bootstrap_external(evaluation, reps = 8L, seed = 19L)

  expect_identical(result$status, "incomplete")
  expect_identical(
    result$point_estimates$metric,
    c(
      "roc_auc", "average_precision", "brier_score", "log_loss",
      "mean_calibration_ratio",
      "calibration_intercept", "calibration_slope",
      "moderate_calibration_ici", "moderate_calibration_e50",
      "moderate_calibration_e90"
    )
  )
  expect_identical(
    result$intervals$status[result$intervals$metric == "brier_score"],
    "success"
  )
  expect_identical(
    result$intervals$n_finite_successful[
      result$intervals$metric == "brier_score"
    ],
    8L
  )
  expect_true(all(
    result$replicates$status[
      result$replicates$metric == "roc_auc"
    ] != "success"
  ))
})

test_that("ordinary paired resampling matches a manual fixed-seed oracle", {
  outcome <- c(0L, 0L, 1L, 1L, 1L)
  probability <- c(0.05, 0.25, 0.55, 0.75, 0.95)
  evaluation <- make_sb_bootstrap_binary_evaluation(outcome, probability)
  reps <- 12L
  seed <- 20260820L
  expected <- withr::with_seed(seed, {
    lapply(seq_len(reps), function(unused) {
      index <- sample.int(length(outcome), length(outcome), replace = TRUE)
      c(
        brier = mean((outcome[index] - probability[index])^2),
        positive = sum(outcome[index] == 1L),
        negative = sum(outcome[index] == 0L)
      )
    })
  })
  expected_brier <- vapply(expected, `[[`, numeric(1), "brier")
  expected_positive <- as.integer(vapply(
    expected,
    `[[`,
    numeric(1),
    "positive"
  ))
  expected_negative <- as.integer(vapply(
    expected,
    `[[`,
    numeric(1),
    "negative"
  ))

  result <- sb_bootstrap_external(
    evaluation,
    reps = reps,
    confidence_level = 0.9,
    seed = seed
  )
  brier_rows <- result$replicates[
    result$replicates$metric == "brier_score",
    ,
    drop = FALSE
  ]
  mean_calibration_rows <- result$replicates[
    result$replicates$metric == "mean_calibration_ratio",
    ,
    drop = FALSE
  ]
  expected_mean_calibration <- vapply(
    withr::with_seed(seed, {
      lapply(seq_len(reps), function(unused) {
        index <- sample.int(length(outcome), length(outcome), replace = TRUE)
        sum(outcome[index]) / sum(probability[index])
      })
    }),
    identity,
    numeric(1)
  )

  expect_equal(brier_rows$estimate, expected_brier, tolerance = 0)
  expect_identical(mean_calibration_rows$status, rep("success", reps))
  expect_equal(
    mean_calibration_rows$estimate,
    expected_mean_calibration,
    tolerance = 0
  )
  expect_identical(result$resample_audit$n_positive, expected_positive)
  expect_identical(result$resample_audit$n_negative, expected_negative)
  expect_identical(result$resampling$unit, "patient")
  expect_identical(
    result$resampling$scheme,
    "ordinary_nonparametric_with_replacement"
  )
  expect_identical(result$resampling$paired_endpoint_prediction, TRUE)
  expect_identical(result$resampling$prediction, "fixed")
  expect_identical(result$resampling$model_access, "none")
})

test_that("bootstrap retains the effective external prediction target", {
  binary <- sb_bootstrap_external(
    make_sb_bootstrap_binary_evaluation(
      c(0L, 0L, 1L, 1L),
      c(0.1, 0.3, 0.7, 0.9)
    ),
    reps = 2L,
    seed = 1L
  )
  risk_only <- sb_bootstrap_external(
    make_sb_bootstrap_survival_evaluation(
      time = 1:6,
      event = c(1L, 0L, 1L, 0L, 1L, 0L),
      risk = c(0.8, 0.2, 0.5, -0.1, 0.3, -0.4),
      time_unit = "month"
    ),
    reps = 2L,
    seed = 2L
  )
  endpoint_horizon <- sb_bootstrap_external(
    make_sb_bootstrap_survival_evaluation(
      time = 1:6,
      event = c(1L, 0L, 1L, 0L, 1L, 0L),
      risk = c(0.8, 0.2, 0.5, -0.1, 0.3, -0.4),
      survival_probability = c(0.4, 0.7, 0.5, 0.8, 0.6, 0.9),
      prediction_horizon = 2,
      time_unit = "year"
    ),
    reps = 2L,
    seed = 3L
  )
  prediction_horizon <- sb_bootstrap_external(
    make_sb_bootstrap_survival_evaluation(
      time = 1:6,
      event = c(1L, 0L, 1L, 0L, 1L, 0L),
      risk = c(0.8, 0.2, 0.5, -0.1, 0.3, -0.4),
      survival_probability = c(0.4, 0.7, 0.5, 0.8, 0.6, 0.9),
      prediction_horizon = 2,
      time_unit = "year",
      endpoint_prediction_horizon = NULL
    ),
    reps = 2L,
    seed = 4L
  )

  expect_identical(binary$prediction_target, list(
    target_status = "complete",
    has_horizon = FALSE,
    prediction_horizon = NULL,
    time_unit = NULL,
    horizon_source = "not_applicable"
  ))
  expect_identical(risk_only$prediction_target, list(
    target_status = "complete",
    has_horizon = FALSE,
    prediction_horizon = NULL,
    time_unit = "month",
    horizon_source = "none"
  ))
  expect_identical(endpoint_horizon$prediction_target, list(
    target_status = "complete",
    has_horizon = TRUE,
    prediction_horizon = 2,
    time_unit = "year",
    horizon_source = "endpoint"
  ))
  expect_identical(prediction_horizon$prediction_target, list(
    target_status = "complete",
    has_horizon = TRUE,
    prediction_horizon = 2,
    time_unit = "year",
    horizon_source = "prediction"
  ))
  expect_match(
    paste(capture.output(print(risk_only)), collapse = "\n"),
    "Target: horizon-free risk",
    fixed = TRUE
  )
  expect_match(
    paste(capture.output(print(prediction_horizon)), collapse = "\n"),
    "Horizon source: prediction",
    fixed = TRUE
  )
})

test_that("bootstrap preserves RNG state and records the actual RNG kinds", {
  evaluation <- make_sb_bootstrap_binary_evaluation(
    c(0L, 0L, 1L, 1L),
    c(0.1, 0.4, 0.6, 0.9)
  )
  set.seed(481L)
  before_seed <- .Random.seed
  before_kind <- RNGkind()

  first <- sb_bootstrap_external(evaluation, reps = 10L, seed = 91L)

  expect_identical(.Random.seed, before_seed)
  expect_identical(RNGkind(), before_kind)
  expect_identical(first$provenance$rng_kind, before_kind[[1L]])
  expect_identical(first$provenance$normal_kind, before_kind[[2L]])
  expect_identical(first$provenance$sample_kind, before_kind[[3L]])
  second <- sb_bootstrap_external(evaluation, reps = 10L, seed = 91L)
  changed <- sb_bootstrap_external(evaluation, reps = 10L, seed = 92L)
  expect_identical(first$replicates, second$replicates)
  expect_identical(first$resample_audit, second$resample_audit)
  expect_false(identical(first$replicates, changed$replicates))
})

test_that("bootstrap output stores no model, predictors, or resample indices", {
  evaluation <- make_sb_bootstrap_binary_evaluation(
    c(0L, 0L, 1L, 1L),
    c(0.1, 0.3, 0.7, 0.9)
  )
  result <- sb_bootstrap_external(evaluation, reps = 5L, seed = 8L)
  flattened_names <- unlist(lapply(result, names), use.names = FALSE)

  expect_false(any(c(
    "fit", "model", "runtime", "predictors", "indices", "patient_ids"
  ) %in% c(names(result), flattened_names)))
  expect_identical(result$n, 4L)
  expect_identical(result$reps, 5L)
  expect_identical(nrow(result$resample_audit), 5L)
})

test_that("bootstrap canonicalizes away unrelated scored columns", {
  evaluation <- make_sb_bootstrap_binary_evaluation(
    c(0L, 0L, 1L, 1L),
    c(0.1, 0.3, 0.7, 0.9)
  )
  evaluation$scored_predictions$gene_not_for_bootstrap <- c(3, 4, 5, 6)

  result <- sb_bootstrap_external(evaluation, reps = 4L, seed = 11L)

  expect_false(identical(result$status, "failed"))
  expect_false("gene_not_for_bootstrap" %in% names(result))
  expect_false("gene_not_for_bootstrap" %in% names(result$resample_audit))
  expect_false("gene_not_for_bootstrap" %in% names(result$replicates))
})

test_that("mutated patient and horizon rows fail the evaluation contract", {
  binary <- make_sb_bootstrap_binary_evaluation(
    c(0L, 0L, 1L, 1L),
    c(0.1, 0.3, 0.7, 0.9)
  )
  duplicate <- binary
  duplicate$scored_predictions$patient_id[[2L]] <-
    duplicate$scored_predictions$patient_id[[1L]]
  bad_probability <- binary
  bad_probability$scored_predictions$probability[[1L]] <- 2
  bad_n <- binary
  bad_n$n <- 99L
  for (evaluation in list(duplicate, bad_probability, bad_n)) {
    result <- sb_bootstrap_external(evaluation, reps = 2L)
    expect_external_bootstrap_failure(
      result,
      "evaluation_contract",
      "stormbreaker_bad_external_evaluation"
    )
  }

  survival <- make_sb_bootstrap_survival_evaluation(
    c(1, 2, 3),
    c(1L, 0L, 1L),
    c(0.8, 0.5, 0.2),
    c(0.2, 0.5, 0.8),
    2,
    "month"
  )
  changed_horizon <- survival
  changed_horizon$scored_predictions$prediction_horizon[[2L]] <- 3
  changed_unit <- survival
  changed_unit$scored_predictions$time_unit[[2L]] <- "year"
  uniformly_changed_horizon <- survival
  uniformly_changed_horizon$scored_predictions$prediction_horizon <-
    rep(3, nrow(uniformly_changed_horizon$scored_predictions))
  uniformly_changed_unit <- survival
  uniformly_changed_unit$scored_predictions$time_unit <-
    rep("year", nrow(uniformly_changed_unit$scored_predictions))
  bad_survival_probability <- survival
  bad_survival_probability$scored_predictions$survival_probability[[1L]] <- NA
  for (evaluation in list(
    changed_horizon,
    changed_unit,
    uniformly_changed_horizon,
    uniformly_changed_unit,
    bad_survival_probability
  )) {
    result <- sb_bootstrap_external(evaluation, reps = 2L)
    expect_external_bootstrap_failure(
      result,
      "evaluation_contract",
      "stormbreaker_bad_external_evaluation"
    )
  }
})

test_that("failed horizon bootstrap reports validated censoring scope honestly", {
  evaluation <- make_sb_bootstrap_survival_evaluation(
    time = c(1, 2, 3, 4),
    event = c(1L, 0L, 1L, 0L),
    risk = c(0.8, 0.6, 0.4, 0.2),
    survival_probability = c(0.2, 0.4, 0.6, 0.8),
    prediction_horizon = 3,
    time_unit = "month"
  )

  argument_failure <- sb_bootstrap_external(evaluation, reps = 1L)
  expect_external_bootstrap_failure(
    argument_failure,
    "arguments",
    "stormbreaker_bad_integer"
  )
  expect_identical(
    argument_failure$resampling$censoring,
    "unknown_not_validated"
  )

  testthat::local_mocked_bindings(
    sb_external_bootstrap_resamples = function(...) {
      stop("injected resampling failure")
    }
  )
  resampling_failure <- sb_bootstrap_external(
    evaluation,
    reps = 2L,
    seed = 5L
  )
  expect_external_bootstrap_failure(
    resampling_failure,
    "resampling",
    "simpleError"
  )
  expect_identical(
    resampling_failure$resampling$censoring,
    "reestimated_external_reverse_kaplan_meier_each_replicate"
  )
  expect_identical(resampling_failure$prediction_target, list(
    target_status = "complete",
    has_horizon = TRUE,
    prediction_horizon = 3,
    time_unit = "month",
    horizon_source = "endpoint"
  ))
})

test_that("tied binary AUC and average precision match independent oracles", {
  outcome <- c(0L, 0L, 1L, 1L, 0L, 1L)
  probability <- c(0.2, 0.5, 0.5, 0.8, 0.8, 0.8)
  reps <- 14L
  seed <- 202L
  evaluation <- make_sb_bootstrap_binary_evaluation(outcome, probability)
  indices <- withr::with_seed(seed, {
    lapply(seq_len(reps), function(unused) {
      sample.int(length(outcome), length(outcome), replace = TRUE)
    })
  })

  result <- sb_bootstrap_external(evaluation, reps = reps, seed = seed)
  auc_rows <- result$replicates[result$replicates$metric == "roc_auc", ]
  ap_rows <- result$replicates[
    result$replicates$metric == "average_precision",
  ]
  for (replicate in seq_len(reps)) {
    index <- indices[[replicate]]
    expected_auc <- manual_binary_auc(outcome[index], probability[index])
    expected_ap <- manual_binary_average_precision(
      outcome[index],
      probability[index]
    )
    if (is.finite(expected_auc)) {
      expect_identical(auc_rows$status[[replicate]], "success")
      expect_equal(auc_rows$estimate[[replicate]], expected_auc)
      expect_identical(ap_rows$status[[replicate]], "success")
      expect_equal(ap_rows$estimate[[replicate]], expected_ap)
    } else {
      expect_false(identical(auc_rows$status[[replicate]], "success"))
      expect_false(identical(ap_rows$status[[replicate]], "success"))
    }
  }
})

test_that("point targets are recomputed instead of copied from metric summaries", {
  outcome <- c(0L, 1L, 0L, 1L, 0L, 1L)
  probability <- c(0.15, 0.35, 0.45, 0.65, 0.55, 0.85)
  evaluation <- make_sb_bootstrap_binary_evaluation(outcome, probability)
  original_metrics <- evaluation$metrics
  evaluation$metrics$roc_auc$estimate <- -99
  evaluation$metrics$brier_score$estimate <- -99
  evaluation$metrics$log_loss$estimate <- -99

  result <- sb_bootstrap_external(evaluation, reps = 8L, seed = 12L)
  point <- setNames(result$point_estimates$estimate, result$point_estimates$metric)

  expect_equal(point[["roc_auc"]], original_metrics$roc_auc$estimate)
  expect_equal(
    point[["average_precision"]],
    original_metrics$average_precision$estimate
  )
  expect_equal(point[["brier_score"]], mean((outcome - probability)^2))
  expect_equal(
    point[["log_loss"]],
    mean(-(outcome * log(probability) +
      (1 - outcome) * log1p(-probability)))
  )
  expect_equal(
    point[["mean_calibration_ratio"]],
    sum(outcome) / sum(probability),
    tolerance = 0
  )
  expect_false(any(result$point_estimates$estimate == -99))
})

test_that("calibration intercept and slope retain component-specific status", {
  evaluation <- make_sb_bootstrap_binary_evaluation(
    c(0L, 0L, 1L, 1L),
    c(0.1, 0.2, 0.8, 0.9)
  )
  result <- sb_bootstrap_external(evaluation, reps = 5L, seed = 3L)
  point <- result$point_estimates

  expect_identical(
    point$status[point$metric == "calibration_intercept"],
    "success"
  )
  expect_identical(
    point$status[point$metric == "calibration_slope"],
    "failed"
  )
  expect_identical(
    point$error_class[point$metric == "calibration_slope"],
    "stormbreaker_calibration_separation"
  )
  expect_true(is.finite(
    point$estimate[point$metric == "calibration_intercept"]
  ))
  expect_true(is.na(point$estimate[point$metric == "calibration_slope"]))
})

test_that("binary bootstrap fits one moderate curve per sampled patient table", {
  n <- 80L
  outcome <- rep(c(0L, 1L), length.out = n)
  probability <- seq(0.03, 0.97, length.out = n)
  evaluation <- make_sb_bootstrap_binary_evaluation(outcome, probability)
  reps <- 9L
  seed <- 810L
  original_metric <- sb_binary_moderate_calibration_metric
  calls <- 0L
  seen_ids <- list()
  testthat::local_mocked_bindings(
    sb_binary_moderate_calibration_metric = function(scored, scope) {
      calls <<- calls + 1L
      seen_ids[[calls]] <<- scored$patient_id
      original_metric(scored, scope)
    }
  )

  result <- sb_bootstrap_external(evaluation, reps = reps, seed = seed)
  target_names <- c(
    "moderate_calibration_ici",
    "moderate_calibration_e50",
    "moderate_calibration_e90"
  )
  direct <- original_metric(
    evaluation$scored_predictions,
    "patient_level_external_bootstrap_sample"
  )
  point <- result$point_estimates[
    result$point_estimates$metric %in% target_names,
    ,
    drop = FALSE
  ]

  expect_identical(calls, reps + 1L)
  expect_identical(result$point_estimates$metric, c(
    "roc_auc", "average_precision", "brier_score", "log_loss",
    "mean_calibration_ratio",
    "calibration_intercept", "calibration_slope", target_names
  ))
  expect_identical(point$metric, target_names)
  expect_identical(point$status, rep("success", 3L))
  expect_equal(
    point$estimate,
    c(direct$ici, direct$e50, direct$e90),
    tolerance = 0
  )
  target_counts <- table(
    factor(
      result$replicates$metric[
        result$replicates$metric %in% target_names
      ],
      levels = target_names
    )
  )
  expect_identical(names(target_counts), target_names)
  expect_identical(as.integer(target_counts), rep.int(reps, 3L))
  expect_identical(seen_ids[[1L]], evaluation$scored_predictions$patient_id)
  expect_identical(
    all(vapply(seen_ids[-1L], length, integer(1)) == n),
    TRUE
  )
  expect_identical(
    any(vapply(seen_ids[-1L], anyDuplicated, integer(1)) > 0L),
    TRUE
  )
})

test_that("binary moderate bootstrap errors affect only their three targets", {
  n <- 80L
  evaluation <- make_sb_bootstrap_binary_evaluation(
    rep(c(0L, 1L), length.out = n),
    seq(0.03, 0.97, length.out = n)
  )
  reps <- 8L
  seed <- 811L
  control_result <- sb_bootstrap_external(evaluation, reps = reps, seed = seed)
  calls <- 0L
  set.seed(20260820L)
  rng_before <- .Random.seed
  testthat::local_mocked_bindings(
    sb_binary_moderate_calibration_metric = function(...) {
      calls <<- calls + 1L
      condition <- simpleError("injected bootstrap moderate calibration error")
      class(condition) <- c(
        "test_bootstrap_binary_moderate_calibration_error",
        class(condition)
      )
      stop(condition)
    },
    sb_fit_final = function(...) stop("must not refit"),
    sb_predict_final = function(...) stop("must not repredict"),
    sb_evaluate_external = function(...) stop("must not reevaluate")
  )

  observed <- sb_bootstrap_external(evaluation, reps = reps, seed = seed)
  rng_after <- .Random.seed
  targets <- c(
    "moderate_calibration_ici",
    "moderate_calibration_e50",
    "moderate_calibration_e90"
  )
  point_other <- !control_result$point_estimates$metric %in% targets
  replicate_other <- !control_result$replicates$metric %in% targets
  interval_other <- !control_result$intervals$metric %in% targets

  expect_identical(calls, reps + 1L)
  expect_identical(rng_after, rng_before)
  expect_identical(observed$resample_audit, control_result$resample_audit)
  expect_identical(observed$resampling, control_result$resampling)
  expect_identical(observed$interval_method, control_result$interval_method)
  expect_identical(observed$provenance, control_result$provenance)
  expect_identical(
    observed$point_estimates[point_other, , drop = FALSE],
    control_result$point_estimates[point_other, , drop = FALSE]
  )
  expect_identical(
    observed$replicates[replicate_other, , drop = FALSE],
    control_result$replicates[replicate_other, , drop = FALSE]
  )
  expect_identical(
    observed$intervals[interval_other, , drop = FALSE],
    control_result$intervals[interval_other, , drop = FALSE]
  )
  expect_identical(
    observed$point_estimates$status[
      observed$point_estimates$metric %in% targets
    ],
    rep("metric_failed", 3L)
  )
  expect_identical(
    unique(observed$replicates$error_class[
      observed$replicates$metric %in% targets
    ]),
    "test_bootstrap_binary_moderate_calibration_error"
  )
})

test_that("bootstrap never calls fitting prediction evaluation or event KM", {
  evaluation <- make_sb_bootstrap_binary_evaluation(
    c(0L, 0L, 1L, 1L, 0L, 1L),
    c(0.1, 0.3, 0.6, 0.8, 0.4, 0.7)
  )
  before <- evaluation
  testthat::local_mocked_bindings(
    sb_fit_final = function(...) stop("must not refit"),
    sb_predict_final = function(...) stop("must not repredict"),
    sb_evaluate_external = function(...) stop("must not reevaluate"),
    sb_fit_event_km = function(...) stop("event KM must not run")
  )

  result <- sb_bootstrap_external(evaluation, reps = 6L, seed = 17L)

  expect_false(identical(result$status, "failed"))
  expect_identical(evaluation, before)
  expect_identical(result$resampling$model_access, "none")
})

test_that("risk-only survival replicates match independent Harrell calculations", {
  time <- c(1, 2, 3, 4, 5, 6)
  event <- c(1L, 1L, 0L, 1L, 0L, 1L)
  risk <- c(0.9, 0.7, 0.6, 0.5, 0.2, 0.1)
  reps <- 15L
  seed <- 62L
  evaluation <- make_sb_bootstrap_survival_evaluation(time, event, risk)
  indices <- withr::with_seed(seed, {
    lapply(seq_len(reps), function(unused) {
      sample.int(length(time), length(time), replace = TRUE)
    })
  })

  result <- sb_bootstrap_external(evaluation, reps = reps, seed = seed)
  rows <- result$replicates[result$replicates$metric == "harrell_c", ]

  for (replicate in seq_len(reps)) {
    index <- indices[[replicate]]
    expected <- survival::concordance(
      survival::Surv(time[index], event[index]) ~ risk[index],
      reverse = TRUE,
      timewt = "n"
    )
    comparable <- unname(sum(expected$count[c(
      "concordant", "discordant", "tied.x"
    )]))
    if (is.finite(comparable) && comparable > 0) {
      estimate <- unname(
        (expected$count[["concordant"]] +
          0.5 * expected$count[["tied.x"]]) / comparable
      )
      expect_identical(rows$status[[replicate]], "success")
      expect_equal(rows$estimate[[replicate]], estimate)
    } else {
      expect_false(identical(rows$status[[replicate]], "success"))
      expect_true(is.na(rows$estimate[[replicate]]))
    }
  }
  expect_identical(result$point_estimates$metric, "harrell_c")
  expect_identical(
    result$resampling$censoring,
    "not_applicable_risk_only"
  )
})

test_that("horizon survival shares each patient draw across all scalar metrics", {
  time <- c(1, 2, 2, 3, 4, 5)
  event <- c(1L, 1L, 0L, 0L, 1L, 0L)
  risk <- c(0.9, 0.7, 0.6, 0.5, 0.4, 0.2)
  probability <- c(0.2, 0.3, 0.45, 0.55, 0.7, 0.8)
  horizon <- 3
  reps <- 9L
  seed <- 81L
  evaluation <- make_sb_bootstrap_survival_evaluation(
    time,
    event,
    risk,
    probability,
    horizon,
    "month"
  )
  indices <- withr::with_seed(seed, {
    lapply(seq_len(reps), function(unused) {
      sample.int(length(time), length(time), replace = TRUE)
    })
  })
  original_censoring_fit <- sb_fit_censoring_km
  original_calibration <- sb_survival_calibration_metric
  original_weak_calibration <- sb_survival_weak_calibration_metric
  original_event_fit <- sb_fit_event_km
  censoring_calls <- 0L
  calibration_calls <- 0L
  weak_calibration_calls <- 0L
  event_fit_calls <- 0L
  testthat::local_mocked_bindings(
    sb_fit_censoring_km = function(...) {
      censoring_calls <<- censoring_calls + 1L
      original_censoring_fit(...)
    },
    sb_survival_calibration_metric = function(...) {
      calibration_calls <<- calibration_calls + 1L
      original_calibration(...)
    },
    sb_survival_weak_calibration_metric = function(...) {
      weak_calibration_calls <<- weak_calibration_calls + 1L
      original_weak_calibration(...)
    },
    sb_fit_event_km = function(...) {
      event_fit_calls <<- event_fit_calls + 1L
      original_event_fit(...)
    }
  )

  result <- sb_bootstrap_external(evaluation, reps = reps, seed = seed)

  expect_identical(censoring_calls, reps + 1L)
  expect_identical(calibration_calls, reps + 1L)
  expect_identical(weak_calibration_calls, reps + 1L)
  expect_identical(event_fit_calls, reps + 1L)
  expect_identical(
    result$point_estimates$metric,
    c(
      "harrell_c", "ipcw_brier_score", "time_auc", "uno_c",
      "calibration_ici", "calibration_e50", "calibration_e90",
      "mean_calibration_ratio", "calibration_slope"
    )
  )
  expect_identical(
    result$resampling$censoring,
    "reestimated_external_reverse_kaplan_meier_each_replicate"
  )
  brier_rows <- result$replicates[
    result$replicates$metric == "ipcw_brier_score",
    ,
    drop = FALSE
  ]
  time_auc_rows <- result$replicates[
    result$replicates$metric == "time_auc",
    ,
    drop = FALSE
  ]
  uno_rows <- result$replicates[
    result$replicates$metric == "uno_c",
    ,
    drop = FALSE
  ]
  expect_equal(
    result$point_estimates$estimate[
      result$point_estimates$metric == "time_auc"
    ],
    manual_external_time_auc(time, event, probability, horizon),
    tolerance = 1e-15
  )
  expect_equal(
    result$point_estimates$estimate[
      result$point_estimates$metric == "uno_c"
    ],
    manual_external_uno_c(time, event, risk, horizon),
    tolerance = 1e-15
  )
  for (replicate in seq_len(reps)) {
    index <- indices[[replicate]]
    expected <- manual_external_ipcw_brier(
      time[index],
      event[index],
      probability[index],
      horizon
    )
    if (is.finite(expected)) {
      expect_identical(brier_rows$status[[replicate]], "success")
      expect_equal(brier_rows$estimate[[replicate]], expected)
    } else {
      expect_false(identical(brier_rows$status[[replicate]], "success"))
      expect_true(is.na(brier_rows$estimate[[replicate]]))
    }
    expected_auc <- manual_external_time_auc(
      time[index],
      event[index],
      probability[index],
      horizon
    )
    if (is.finite(expected_auc)) {
      expect_identical(time_auc_rows$status[[replicate]], "success")
      expect_equal(
        time_auc_rows$estimate[[replicate]],
        expected_auc,
        tolerance = 1e-15
      )
    } else {
      expect_false(identical(time_auc_rows$status[[replicate]], "success"))
      expect_true(is.na(time_auc_rows$estimate[[replicate]]))
    }
    expected_uno <- manual_external_uno_c(
      time[index],
      event[index],
      risk[index],
      horizon
    )
    if (is.finite(expected_uno)) {
      expect_identical(uno_rows$status[[replicate]], "success")
      expect_equal(
        uno_rows$estimate[[replicate]],
        expected_uno,
        tolerance = 1e-15
      )
    } else {
      expect_false(identical(uno_rows$status[[replicate]], "success"))
      expect_true(is.na(uno_rows$estimate[[replicate]]))
    }
  }
})

test_that("horizon bootstrap calibration matches point and manual resamples", {
  evaluation <- make_sb_bootstrap_calibration_evaluation()
  reps <- 8L
  seed <- 91L
  result <- sb_bootstrap_external(evaluation, reps = reps, seed = seed)
  target_names <- c("calibration_ici", "calibration_e50", "calibration_e90")
  point <- result$point_estimates[
    result$point_estimates$metric %in% target_names,
    ,
    drop = FALSE
  ]

  expect_identical(point$metric, target_names)
  expect_identical(point$status, rep("success", 3L))
  expect_identical(
    point$estimate,
    unname(unlist(
      evaluation$metrics$calibration[c("ici", "e50", "e90")],
      use.names = FALSE
    ))
  )

  scored <- sb_external_bootstrap_contract(evaluation)$scored
  indices <- withr::with_seed(seed, {
    lapply(seq_len(reps), function(unused) {
      sample.int(nrow(scored), nrow(scored), replace = TRUE)
    })
  })
  successful_duplicate <- FALSE
  for (replicate in seq_len(5L)) {
    index <- indices[[replicate]]
    sampled <- scored[index, , drop = FALSE]
    rownames(sampled) <- NULL
    expected <- sb_survival_calibration_metric(
      patient_id = sampled$patient_id,
      time = sampled$time,
      event = sampled$event,
      survival_probability = sampled$survival_probability,
      prediction_status = sampled$status,
      prediction_horizon = unique(sampled$prediction_horizon)[[1L]],
      time_unit = unique(sampled$time_unit)[[1L]],
      scope = "patient_level_external_cohort_at_locked_horizon"
    )
    observed <- result$replicates[
      result$replicates$replicate == replicate &
        result$replicates$metric %in% target_names,
      ,
      drop = FALSE
    ]

    expect_identical(observed$metric, target_names)
    expect_identical(observed$status, rep(expected$status, 3L))
    expect_identical(
      observed$error_class,
      rep(expected$error_class, 3L)
    )
    expect_identical(
      observed$error_message,
      rep(expected$error_message, 3L)
    )
    if (identical(expected$status, "success")) {
      expect_equal(
        observed$estimate,
        c(expected$ici, expected$e50, expected$e90),
        tolerance = 0
      )
      expect_identical(
        expected$patient_calibration$row_instance,
        seq_len(nrow(sampled))
      )
      expect_identical(
        expected$patient_calibration$patient_id,
        sampled$patient_id
      )
      if (anyDuplicated(sampled$patient_id) > 0L) {
        successful_duplicate <- TRUE
      }
    } else {
      expect_true(all(is.na(observed$estimate)))
    }
  }
  expect_true(successful_duplicate)
})

test_that("horizon weak calibration matches point and one-draw resample oracles", {
  evaluation <- make_sb_bootstrap_calibration_evaluation()
  reps <- 8L
  seed <- 129L
  result <- sb_bootstrap_external(evaluation, reps = reps, seed = seed)
  target_names <- c("mean_calibration_ratio", "calibration_slope")
  point <- result$point_estimates[
    result$point_estimates$metric %in% target_names,
    ,
    drop = FALSE
  ]

  expect_identical(point$metric, target_names)
  expect_identical(point$status, rep("success", 2L))
  expect_equal(
    point$estimate,
    c(
      evaluation$metrics$weak_calibration$mean_calibration_ratio,
      evaluation$metrics$weak_calibration$calibration_slope
    ),
    tolerance = 0
  )

  scored <- sb_external_bootstrap_contract(evaluation)$scored
  indices <- withr::with_seed(seed, {
    lapply(seq_len(reps), function(unused) {
      sample.int(nrow(scored), nrow(scored), replace = TRUE)
    })
  })
  successful_duplicate <- FALSE
  for (replicate in seq_len(5L)) {
    sampled <- scored[indices[[replicate]], , drop = FALSE]
    rownames(sampled) <- NULL
    expected <- c(
      manual_external_mean_calibration_ratio(
        sampled$time,
        sampled$event,
        sampled$survival_probability,
        unique(sampled$prediction_horizon)[[1L]]
      ),
      manual_external_calibration_slope(
        sampled$time,
        sampled$event,
        sampled$survival_probability,
        unique(sampled$prediction_horizon)[[1L]]
      )
    )
    observed <- result$replicates[
      result$replicates$replicate == replicate &
        result$replicates$metric %in% target_names,
      ,
      drop = FALSE
    ]

    expect_identical(observed$metric, target_names)
    expect_identical(observed$status, rep("success", 2L))
    expect_equal(observed$estimate, expected, tolerance = 1e-12)
    if (anyDuplicated(sampled$patient_id) > 0L) {
      successful_duplicate <- TRUE
    }
  }
  expect_true(successful_duplicate)
})

test_that("weak-calibration replicate failures stay component-specific", {
  time <- c(0.5, 1, 1.5, 4, 5, 6, 7, 8)
  event <- c(1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L)
  survival_probability <- exp(-exp(c(
    -0.2, 0.7, -1.1, 1.0, 0.3, -0.5, 1.3, -1.4
  )))
  evaluation <- make_sb_bootstrap_survival_evaluation(
    time,
    event,
    risk = seq(0.8, -0.6, length.out = 8L),
    survival_probability = survival_probability,
    prediction_horizon = 2,
    time_unit = "year"
  )
  reps <- 40L
  seed <- 37L
  result <- sb_bootstrap_external(evaluation, reps = reps, seed = seed)
  scored <- sb_external_bootstrap_contract(evaluation)$scored
  indices <- withr::with_seed(seed, {
    lapply(seq_len(reps), function(unused) {
      sample.int(nrow(scored), nrow(scored), replace = TRUE)
    })
  })
  expected <- lapply(indices, function(index) {
    sampled <- scored[index, , drop = FALSE]
    sb_survival_weak_calibration_metric(
      patient_id = sampled$patient_id,
      time = sampled$time,
      event = sampled$event,
      survival_probability = sampled$survival_probability,
      prediction_status = sampled$status,
      prediction_horizon = unique(sampled$prediction_horizon)[[1L]],
      time_unit = unique(sampled$time_unit)[[1L]],
      scope = "patient_level_external_cohort_at_locked_horizon"
    )
  })
  expected_mean_status <- vapply(
    expected,
    function(metric) metric$mean_calibration$status,
    character(1)
  )
  expected_slope_status <- vapply(
    expected,
    function(metric) metric$slope_fit$status,
    character(1)
  )
  mean_rows <- result$replicates[
    result$replicates$metric == "mean_calibration_ratio",
    ,
    drop = FALSE
  ]
  slope_rows <- result$replicates[
    result$replicates$metric == "calibration_slope",
    ,
    drop = FALSE
  ]

  expect_identical(mean_rows$status, unname(expected_mean_status))
  expect_identical(slope_rows$status, unname(expected_slope_status))
  expect_gt(sum(slope_rows$status == "success"), 0L)
  expect_gt(sum(slope_rows$status != "success"), 0L)
  expect_true(all(is.na(slope_rows$estimate[slope_rows$status != "success"])))
  slope_interval <- result$intervals[
    result$intervals$metric == "calibration_slope",
    ,
    drop = FALSE
  ]
  expect_identical(
    slope_interval$n_finite_successful,
    as.integer(sum(slope_rows$status == "success" &
      is.finite(slope_rows$estimate)))
  )
  expect_identical(
    slope_interval$n_unavailable,
    as.integer(sum(slope_rows$status != "success"))
  )
})

test_that("weak-calibration component status cannot erase neighboring targets", {
  evaluation <- make_sb_bootstrap_calibration_evaluation()
  reps <- 4L
  seed <- 55L
  baseline <- sb_bootstrap_external(evaluation, reps = reps, seed = seed)

  make_mock_metric <- function(mean_success) {
    metric <- sb_survival_weak_calibration_empty_result(
      prediction_horizon = 2,
      time_unit = "year",
      scope = "patient_level_external_cohort_at_locked_horizon",
      n = evaluation$n
    )
    metric$input_status <- "success"
    if (mean_success) {
      metric$mean_calibration$status <- "success"
      metric$mean_calibration$mean_calibration_ratio <- 0.75
      metric$mean_calibration_ratio <- 0.75
      metric$slope_fit$status <- "not_evaluable"
      metric$slope_fit$error_class <- "test_slope_unavailable"
      metric$slope_fit$error_message <- "slope unavailable"
    } else {
      metric$mean_calibration$status <- "not_evaluable"
      metric$mean_calibration$error_class <- "test_mean_unavailable"
      metric$mean_calibration$error_message <- "mean unavailable"
      metric$slope_fit$status <- "success"
      metric$slope_fit$coefficient <- 1.25
      metric$calibration_slope <- 1.25
    }
    metric
  }

  first <- local({
    calls <- 0L
    testthat::local_mocked_bindings(
      sb_survival_weak_calibration_metric = function(...) {
        calls <<- calls + 1L
        make_mock_metric(TRUE)
      }
    )
    value <- sb_bootstrap_external(evaluation, reps = reps, seed = seed)
    expect_identical(calls, reps + 1L)
    value
  })
  second <- local({
    calls <- 0L
    testthat::local_mocked_bindings(
      sb_survival_weak_calibration_metric = function(...) {
        calls <<- calls + 1L
        make_mock_metric(FALSE)
      }
    )
    value <- sb_bootstrap_external(evaluation, reps = reps, seed = seed)
    expect_identical(calls, reps + 1L)
    value
  })

  first_mean <- first$point_estimates$metric == "mean_calibration_ratio"
  first_slope <- first$point_estimates$metric == "calibration_slope"
  expect_identical(first$point_estimates$status[first_mean], "success")
  expect_identical(first$point_estimates$estimate[first_mean], 0.75)
  expect_identical(first$point_estimates$status[first_slope], "not_evaluable")
  expect_identical(
    first$point_estimates$error_class[first_slope],
    "test_slope_unavailable"
  )
  second_mean <- second$point_estimates$metric == "mean_calibration_ratio"
  second_slope <- second$point_estimates$metric == "calibration_slope"
  expect_identical(second$point_estimates$status[second_mean], "not_evaluable")
  expect_identical(
    second$point_estimates$error_class[second_mean],
    "test_mean_unavailable"
  )
  expect_identical(second$point_estimates$status[second_slope], "success")
  expect_identical(second$point_estimates$estimate[second_slope], 1.25)

  weak_targets <- c("mean_calibration_ratio", "calibration_slope")
  expect_identical(
    first$point_estimates[!first$point_estimates$metric %in% weak_targets, ],
    baseline$point_estimates[
      !baseline$point_estimates$metric %in% weak_targets,
    ]
  )
  expect_identical(
    second$replicates[!second$replicates$metric %in% weak_targets, ],
    baseline$replicates[!baseline$replicates$metric %in% weak_targets, ]
  )
  expect_identical(first$resample_audit, baseline$resample_audit)
  expect_identical(second$resample_audit, baseline$resample_audit)
})

test_that("unavailable calibration resamples remain counted and provisional", {
  generated <- withr::with_seed(1L, {
    x <- stats::rnorm(40L)
    data.frame(
      time = c(3, stats::runif(39L, 0.1, 1.9)),
      event = stats::rbinom(40L, 1L, 0.55),
      risk = x,
      survival_probability = stats::plogis(-0.3 - 0.7 * x)
    )
  })
  evaluation <- make_sb_bootstrap_survival_evaluation(
    generated$time,
    generated$event,
    generated$risk,
    generated$survival_probability,
    prediction_horizon = 2,
    time_unit = "year"
  )

  result <- sb_bootstrap_external(evaluation, reps = 20L, seed = 1L)
  target_names <- c("calibration_ici", "calibration_e50", "calibration_e90")
  intervals <- result$intervals[
    result$intervals$metric %in% target_names,
    ,
    drop = FALSE
  ]
  replicate_rows <- result$replicates[
    result$replicates$metric %in% target_names,
    ,
    drop = FALSE
  ]

  expect_identical(intervals$metric, target_names)
  expect_identical(intervals$status, rep("incomplete", 3L))
  expect_identical(intervals$interval_available, rep(TRUE, 3L))
  expect_identical(intervals$replicate_completeness, rep("partial", 3L))
  expect_identical(intervals$n_finite_successful, rep(17L, 3L))
  expect_identical(intervals$n_unavailable, rep(3L, 3L))
  expect_identical(
    intervals$error_class,
    rep("stormbreaker_incomplete_bootstrap_replicates", 3L)
  )
  failed <- replicate_rows$status != "success"
  expect_identical(sum(failed), 9L)
  expect_identical(
    unique(replicate_rows$error_class[failed]),
    "stormbreaker_survival_calibration_horizon_unsupported"
  )
  expect_identical(nrow(result$resample_audit), 20L)
  expect_identical(result$resample_audit$replicate, seq_len(20L))
})

test_that("a calibration exception leaves bootstrap sampling and metrics intact", {
  evaluation <- make_sb_bootstrap_calibration_evaluation()
  reps <- 6L
  seed <- 44L
  baseline <- sb_bootstrap_external(evaluation, reps = reps, seed = seed)
  scored <- sb_external_bootstrap_contract(evaluation)$scored
  indices <- withr::with_seed(seed, {
    lapply(seq_len(reps), function(unused) {
      sample.int(nrow(scored), nrow(scored), replace = TRUE)
    })
  })
  original_calibration <- sb_survival_calibration_metric
  calls <- 0L
  captured_ids <- list()
  testthat::local_mocked_bindings(
    sb_survival_calibration_metric = function(...) {
      calls <<- calls + 1L
      arguments <- list(...)
      captured_ids[[calls]] <<- arguments$patient_id
      if (calls == 3L) {
        condition <- simpleError("injected survival calibration failure")
        class(condition) <- c(
          "test_bootstrap_calibration_error",
          class(condition)
        )
        stop(condition)
      }
      do.call(original_calibration, arguments)
    }
  )
  set.seed(2468L)
  before_seed <- .Random.seed

  observed <- sb_bootstrap_external(evaluation, reps = reps, seed = seed)

  expect_identical(.Random.seed, before_seed)
  expect_identical(calls, reps + 1L)
  for (replicate in seq_len(reps)) {
    expect_identical(
      captured_ids[[replicate + 1L]],
      scored$patient_id[indices[[replicate]]]
    )
  }
  calibration_targets <- c(
    "calibration_ici", "calibration_e50", "calibration_e90"
  )
  non_calibration_point <-
    !baseline$point_estimates$metric %in% calibration_targets
  non_calibration_replicate <-
    !baseline$replicates$metric %in% calibration_targets
  non_calibration_interval <-
    !baseline$intervals$metric %in% calibration_targets
  expect_identical(
    observed$point_estimates[non_calibration_point, , drop = FALSE],
    baseline$point_estimates[non_calibration_point, , drop = FALSE]
  )
  expect_identical(
    observed$replicates[non_calibration_replicate, , drop = FALSE],
    baseline$replicates[non_calibration_replicate, , drop = FALSE]
  )
  expect_identical(
    observed$intervals[non_calibration_interval, , drop = FALSE],
    baseline$intervals[non_calibration_interval, , drop = FALSE]
  )
  expect_identical(observed$resample_audit, baseline$resample_audit)
  injected <- observed$replicates$replicate == 2L &
    observed$replicates$metric %in% calibration_targets
  expect_identical(sum(injected), 3L)
  expect_identical(
    observed$replicates$status[injected],
    rep("metric_failed", 3L)
  )
  expect_identical(
    observed$replicates$error_class[injected],
    rep("test_bootstrap_calibration_error", 3L)
  )
})

test_that("one time-AUC exception does not erase sampled neighbor metrics", {
  evaluation <- make_sb_bootstrap_survival_evaluation(
    time = c(1, 2, 3, 4, 5, 6),
    event = c(1L, 1L, 0L, 1L, 0L, 1L),
    risk = c(0.9, 0.7, 0.6, 0.5, 0.2, 0.1),
    survival_probability = c(0.2, 0.3, 0.45, 0.6, 0.7, 0.8),
    prediction_horizon = 3
  )
  original_auc <- sb_external_time_auc_metric
  auc_calls <- 0L
  testthat::local_mocked_bindings(
    sb_external_time_auc_metric = function(...) {
      auc_calls <<- auc_calls + 1L
      if (auc_calls == 3L) {
        stop("injected time AUC failure")
      }
      original_auc(...)
    }
  )

  result <- sb_bootstrap_external(evaluation, reps = 5L, seed = 44L)
  injected <- !is.na(result$replicates$error_message) &
    result$replicates$error_message == "injected time AUC failure"

  expect_identical(sum(injected, na.rm = TRUE), 1L)
  expect_identical(result$replicates$metric[injected], "time_auc")
  expect_identical(result$replicates$replicate[injected], 2L)
  same_replicate <- result$replicates$replicate == 2L
  expect_false(any(
    result$replicates$error_message[
      same_replicate & result$replicates$metric != "time_auc"
    ] == "injected time AUC failure",
    na.rm = TRUE
  ))
  expect_true(any(result$replicates$replicate > 2L))
})

test_that("one Uno exception does not erase sampled neighbor metrics", {
  evaluation <- make_sb_bootstrap_survival_evaluation(
    time = c(1, 2, 3, 4, 5, 6),
    event = c(1L, 1L, 0L, 1L, 0L, 1L),
    risk = c(0.9, 0.7, 0.6, 0.5, 0.2, 0.1),
    survival_probability = c(0.2, 0.3, 0.45, 0.6, 0.7, 0.8),
    prediction_horizon = 3
  )
  original_uno <- sb_external_uno_c_metric
  uno_calls <- 0L
  testthat::local_mocked_bindings(
    sb_external_uno_c_metric = function(...) {
      uno_calls <<- uno_calls + 1L
      if (uno_calls == 3L) {
        stop("injected Uno bootstrap failure")
      }
      original_uno(...)
    }
  )

  result <- sb_bootstrap_external(evaluation, reps = 5L, seed = 44L)
  injected <- !is.na(result$replicates$error_message) &
    result$replicates$error_message == "injected Uno bootstrap failure"

  expect_identical(sum(injected, na.rm = TRUE), 1L)
  expect_identical(result$replicates$metric[injected], "uno_c")
  expect_identical(result$replicates$replicate[injected], 2L)
  same_replicate <- result$replicates$replicate == 2L
  expect_false(any(
    result$replicates$error_message[
      same_replicate & result$replicates$metric != "uno_c"
    ] == "injected Uno bootstrap failure",
    na.rm = TRUE
  ))
  expect_true(any(result$replicates$replicate > 2L))
})

test_that("one injected metric exception affects only one target replicate", {
  evaluation <- make_sb_bootstrap_binary_evaluation(
    c(0L, 0L, 0L, 1L, 1L, 1L),
    c(0.1, 0.25, 0.4, 0.6, 0.75, 0.9)
  )
  original_auc <- sb_external_binary_auc_metric
  auc_calls <- 0L
  testthat::local_mocked_bindings(
    sb_external_binary_auc_metric = function(scored) {
      auc_calls <<- auc_calls + 1L
      if (auc_calls == 3L) {
        stop("injected AUC failure")
      }
      original_auc(scored)
    }
  )

  result <- sb_bootstrap_external(evaluation, reps = 5L, seed = 44L)
  injected <- !is.na(result$replicates$error_message) &
    result$replicates$error_message == "injected AUC failure"

  expect_identical(sum(injected, na.rm = TRUE), 1L)
  expect_identical(result$replicates$metric[injected], "roc_auc")
  expect_identical(result$replicates$replicate[injected], 2L)
  same_replicate <- result$replicates$replicate == 2L &
    result$replicates$metric != "roc_auc"
  expect_false(any(
    result$replicates$error_message[same_replicate] ==
      "injected AUC failure",
    na.rm = TRUE
  ))
  expect_true(any(result$replicates$replicate > 2L))
})

test_that("percentile intervals use type 7 and exhaustive replicate counts", {
  point <- data.frame(
    metric = c("complete", "provisional", "single", "none", "bad_point"),
    estimate = c(3, 3, 1, 1, Inf),
    status = rep("success", 5L),
    error_class = rep(NA_character_, 5L),
    error_message = rep(NA_character_, 5L),
    stringsAsFactors = FALSE
  )
  replicate_values <- list(
    complete = c(1, 2, 3, 4, 5),
    provisional = c(1, 2, 3, Inf, NA_real_),
    single = c(1, NA_real_, NA_real_, NA_real_, NA_real_),
    none = c(Inf, -Inf, NaN, NA_real_, NA_real_),
    bad_point = c(1, 2, 3, 4, 5)
  )
  replicate_status <- list(
    complete = rep("success", 5L),
    provisional = c(rep("success", 4L), "not_evaluable"),
    single = c("success", rep("not_evaluable", 4L)),
    none = c(rep("success", 4L), "metric_failed"),
    bad_point = rep("success", 5L)
  )
  replicates <- do.call(rbind, lapply(names(replicate_values), function(metric) {
    data.frame(
      replicate = seq_len(5L),
      metric = metric,
      estimate = replicate_values[[metric]],
      status = replicate_status[[metric]],
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(replicates) <- NULL

  intervals <- sb_external_percentile_intervals(
    point,
    replicates,
    confidence_level = 0.8,
    reps = 5L
  )
  row <- function(metric) intervals[intervals$metric == metric, , drop = FALSE]

  expect_equal(
    c(row("complete")$lower, row("complete")$upper),
    unname(stats::quantile(1:5, c(0.1, 0.9), type = 7))
  )
  expect_identical(row("complete")$status, "success")
  expect_identical(row("complete")$replicate_completeness, "complete")
  expect_equal(
    c(row("provisional")$lower, row("provisional")$upper),
    unname(stats::quantile(1:3, c(0.1, 0.9), type = 7))
  )
  expect_identical(row("provisional")$status, "incomplete")
  expect_identical(row("provisional")$interval_available, TRUE)
  expect_identical(row("provisional")$replicate_completeness, "partial")
  expect_identical(row("provisional")$n_finite_successful, 3L)
  expect_identical(row("provisional")$n_nonfinite, 1L)
  expect_identical(row("provisional")$n_unavailable, 1L)
  expect_identical(row("single")$replicate_completeness, "partial")
  expect_identical(row("single")$interval_available, FALSE)
  expect_identical(row("none")$replicate_completeness, "none")
  expect_identical(row("none")$n_finite_successful, 0L)
  expect_identical(row("none")$n_nonfinite, 4L)
  expect_identical(row("none")$n_unavailable, 1L)
  expect_identical(row("bad_point")$interval_available, FALSE)
  expect_identical(row("bad_point")$status, "incomplete")
  expect_true(is.na(row("bad_point")$lower))
  expect_true(is.na(row("bad_point")$upper))
  expect_identical(
    intervals$n_finite_successful + intervals$n_nonfinite +
      intervals$n_unavailable,
    intervals$n_requested
  )
  expect_equal(intervals$valid_fraction, intervals$n_finite_successful / 5)
})

test_that("unclipped infinite log loss never receives a fabricated interval", {
  evaluation <- make_sb_bootstrap_binary_evaluation(
    c(0L, 1L),
    c(1, 0)
  )
  result <- sb_bootstrap_external(evaluation, reps = 7L, seed = 23L)
  log_loss <- result$intervals[
    result$intervals$metric == "log_loss",
    ,
    drop = FALSE
  ]

  expect_identical(result$status, "incomplete")
  expect_identical(log_loss$point_estimate, Inf)
  expect_identical(log_loss$interval_available, FALSE)
  expect_identical(log_loss$replicate_completeness, "none")
  expect_identical(log_loss$n_finite_successful, 0L)
  expect_identical(log_loss$n_nonfinite, 7L)
  expect_identical(log_loss$n_unavailable, 0L)
  expect_true(is.na(log_loss$lower))
  expect_true(is.na(log_loss$upper))
  expect_true(all(is.infinite(
    result$replicates$estimate[result$replicates$metric == "log_loss"]
  )))
})

test_that("one-patient survival remains explicitly not evaluable", {
  evaluation <- make_sb_bootstrap_survival_evaluation(
    time = 2,
    event = 1L,
    risk = 0.4
  )
  result <- sb_bootstrap_external(evaluation, reps = 5L, seed = 7L)

  expect_identical(result$status, "incomplete")
  expect_identical(result$point_estimates$status, "not_evaluable")
  expect_true(is.na(result$point_estimates$estimate))
  expect_identical(result$intervals$replicate_completeness, "none")
  expect_identical(result$intervals$n_unavailable, 5L)
  expect_identical(result$intervals$interval_available, FALSE)
})

test_that("all-event horizon bootstrap keeps censoring survival equal to one", {
  time <- c(1, 2, 3, 4)
  event <- rep.int(1L, 4L)
  risk <- c(0.8, 0.6, 0.4, 0.2)
  probability <- c(0.2, 0.3, 0.4, 0.5)
  horizon <- 6
  evaluation <- make_sb_bootstrap_survival_evaluation(
    time,
    event,
    risk,
    probability,
    horizon
  )
  result <- sb_bootstrap_external(evaluation, reps = 6L, seed = 13L)
  point <- result$point_estimates[
    result$point_estimates$metric == "ipcw_brier_score",
    ,
    drop = FALSE
  ]
  interval <- result$intervals[
    result$intervals$metric == "ipcw_brier_score",
    ,
    drop = FALSE
  ]
  time_auc <- result$intervals[
    result$intervals$metric == "time_auc",
    ,
    drop = FALSE
  ]
  time_auc_point <- result$point_estimates[
    result$point_estimates$metric == "time_auc",
    ,
    drop = FALSE
  ]

  expect_identical(point$status, "success")
  expect_equal(point$estimate, mean(probability^2))
  expect_identical(interval$replicate_completeness, "complete")
  expect_identical(interval$n_finite_successful, 6L)
  expect_identical(interval$n_unavailable, 0L)
  expect_identical(interval$n_nonfinite, 0L)
  expect_identical(time_auc_point$status, "not_evaluable")
  expect_identical(time_auc$replicate_completeness, "none")
  expect_identical(time_auc$n_unavailable, 6L)
})

test_that("all-censored horizon edges distinguish observed and unknown status", {
  risk <- c(0.8, 0.6, 0.4, 0.2)
  probability <- c(0.2, 0.4, 0.6, 0.8)
  after <- make_sb_bootstrap_survival_evaluation(
    time = c(4, 5, 6, 7),
    event = rep.int(0L, 4L),
    risk = risk,
    survival_probability = probability,
    prediction_horizon = 3
  )
  after_result <- sb_bootstrap_external(after, reps = 6L, seed = 4L)
  after_brier <- after_result$intervals[
    after_result$intervals$metric == "ipcw_brier_score",
    ,
    drop = FALSE
  ]
  after_auc <- after_result$intervals[
    after_result$intervals$metric == "time_auc",
    ,
    drop = FALSE
  ]
  after_auc_point <- after_result$point_estimates[
    after_result$point_estimates$metric == "time_auc",
    ,
    drop = FALSE
  ]

  expect_identical(after_result$status, "incomplete")
  expect_identical(after_brier$status, "success")
  expect_equal(after_brier$point_estimate, mean((1 - probability)^2))
  expect_identical(after_brier$replicate_completeness, "complete")
  expect_identical(after_auc_point$status, "not_evaluable")
  expect_identical(after_auc$replicate_completeness, "none")
  expect_identical(after_auc$n_unavailable, 6L)

  before <- make_sb_bootstrap_survival_evaluation(
    time = c(0.5, 1, 1.5, 2),
    event = rep.int(0L, 4L),
    risk = risk,
    survival_probability = probability,
    prediction_horizon = 3
  )
  before_result <- sb_bootstrap_external(before, reps = 6L, seed = 4L)
  before_brier <- before_result$intervals[
    before_result$intervals$metric == "ipcw_brier_score",
    ,
    drop = FALSE
  ]
  before_auc <- before_result$intervals[
    before_result$intervals$metric == "time_auc",
    ,
    drop = FALSE
  ]
  before_auc_point <- before_result$point_estimates[
    before_result$point_estimates$metric == "time_auc",
    ,
    drop = FALSE
  ]

  expect_identical(before_result$status, "incomplete")
  expect_identical(before_brier$interval_available, FALSE)
  expect_identical(before_brier$replicate_completeness, "none")
  expect_identical(before_brier$n_unavailable, 6L)
  expect_true(is.na(before_brier$point_estimate))
  expect_true(is.na(before_brier$lower))
  expect_true(is.na(before_brier$upper))
  expect_identical(before_auc_point$status, "not_evaluable")
  expect_identical(before_auc$replicate_completeness, "none")
  expect_identical(before_auc$n_unavailable, 6L)
})

test_that("bootstrap results serialize without changing numeric evidence", {
  evaluation <- make_sb_bootstrap_binary_evaluation(
    c(0L, 0L, 1L, 1L, 0L, 1L),
    c(0.1, 0.3, 0.55, 0.75, 0.45, 0.9)
  )
  result <- sb_bootstrap_external(
    evaluation,
    reps = 10L,
    confidence_level = 0.9,
    seed = 31L
  )
  restored <- unserialize(serialize(result, NULL))

  expect_s3_class(restored, "stormbreaker_external_bootstrap")
  expect_identical(restored$point_estimates, result$point_estimates)
  expect_identical(restored$replicates, result$replicates)
  expect_identical(restored$intervals, result$intervals)
  expect_identical(restored$resample_audit, result$resample_audit)
})

test_that("bootstrap evidence contains no fold-derived uncertainty fields", {
  evaluation <- make_sb_bootstrap_binary_evaluation(
    c(0L, 0L, 1L, 1L),
    c(0.1, 0.4, 0.6, 0.9)
  )
  result <- sb_bootstrap_external(evaluation, reps = 4L, seed = 9L)
  field_names <- c(
    names(result),
    unlist(lapply(result, names), use.names = FALSE)
  )

  expect_false(any(grepl("fold|standard_deviation|standard_error", field_names)))
  expect_identical(result$interval_method$name, "percentile")
  expect_identical(result$interval_method$quantile_type, 7L)
})

test_that("binary DCA bootstrap stays outside scalar target tables", {
  evaluation <- make_sb_bootstrap_binary_evaluation(
    c(0L, 0L, 1L, 1L, 0L, 1L),
    c(0.1, 0.3, 0.55, 0.75, 0.45, 0.9),
    decision_thresholds = c(0.2, 0.5, 0.8)
  )
  expect_identical(
    evaluation$metrics$decision_curve$status,
    "success"
  )

  result <- sb_bootstrap_external(
    evaluation,
    reps = 6L,
    confidence_level = 0.9,
    seed = 41L
  )
  all_targets <- unique(c(
    result$point_estimates$metric,
    result$intervals$metric,
    result$replicates$metric
  ))

  expect_false(any(grepl(
    "decision|net_benefit|threshold",
    all_targets,
    ignore.case = TRUE
  )))
  expect_identical(result$decision_curve$status, "success")
  expect_identical(
    result$decision_curve$interval_scope,
    "pointwise_not_simultaneous"
  )
})

test_that("survival DCA bootstrap stays nested and preserves scalar evidence", {
  args <- list(
    time = c(1, 2, 3, 4, 5, 6),
    event = c(1L, 0L, 1L, 0L, 1L, 0L),
    risk = c(1, 2, 3, 4, 5, 6),
    survival_probability = c(0.20, 0.30, 0.40, 0.60, 0.45, 0.90),
    prediction_horizon = 4,
    time_unit = "month"
  )
  plain_evaluation <- do.call(
    make_sb_bootstrap_survival_evaluation,
    args
  )
  threshold_evaluation <- do.call(
    make_sb_bootstrap_survival_evaluation,
    c(args, list(decision_thresholds = 0.50))
  )
  plain <- sb_bootstrap_external(
    plain_evaluation,
    reps = 6L,
    confidence_level = 0.9,
    seed = 83L
  )
  thresholded <- sb_bootstrap_external(
    threshold_evaluation,
    reps = 6L,
    confidence_level = 0.9,
    seed = 83L
  )

  expect_identical(
    threshold_evaluation$metrics$decision_curve$status,
    "success"
  )
  expect_identical("decision_curve" %in% names(plain), FALSE)
  expect_identical("decision_curve" %in% names(thresholded), TRUE)
  expect_identical(thresholded$decision_curve$task_type, "survival")
  expect_identical(
    thresholded$decision_curve$interval_scope,
    "pointwise_not_simultaneous"
  )
  expect_identical(thresholded$point_estimates, plain$point_estimates)
  expect_identical(thresholded$intervals, plain$intervals)
  expect_identical(thresholded$replicates, plain$replicates)
  expect_identical(thresholded$resample_audit, plain$resample_audit)
})
