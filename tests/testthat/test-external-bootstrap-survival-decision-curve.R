make_sb_survival_dca_bootstrap_evaluation <- function(
    time,
    event,
    risk,
    survival_probability = NULL,
    prediction_horizon = NULL,
    time_unit = "month",
    decision_thresholds = NULL) {
  ids <- sprintf("sd%02d", seq_along(time))
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
    prediction_horizon = prediction_horizon
  )
  sb_evaluate_external(prediction, endpoint)
}

test_that("survival DCA bootstrap reuses the scalar patient draws", {
  args <- list(
    time = c(1, 2, 5, 5, 4, 6),
    event = c(1L, 0L, 0L, 1L, 0L, 1L),
    risk = c(0.9, 0.7, 0.6, 0.4, 0.2, 0.1),
    survival_probability = c(0.1, 0.1, 0.1, 0.9, 0.9, 0.9),
    prediction_horizon = 4
  )
  plain_evaluation <- do.call(
    make_sb_survival_dca_bootstrap_evaluation,
    args
  )
  threshold_evaluation <- do.call(
    make_sb_survival_dca_bootstrap_evaluation,
    c(args, list(decision_thresholds = 0.5))
  )
  plain <- sb_bootstrap_external(
    plain_evaluation,
    reps = 12L,
    confidence_level = 0.8,
    seed = 1L
  )
  thresholded <- sb_bootstrap_external(
    threshold_evaluation,
    reps = 12L,
    confidence_level = 0.8,
    seed = 1L
  )
  decision_curve <- thresholded$decision_curve
  targets <- c("model_minus_treat_all", "model_minus_treat_none")

  expect_identical(thresholded$point_estimates, plain$point_estimates)
  expect_identical(thresholded$intervals, plain$intervals)
  expect_identical(thresholded$replicates, plain$replicates)
  expect_identical(thresholded$resample_audit, plain$resample_audit)
  expect_identical("decision_curve" %in% names(plain), FALSE)
  expect_identical("decision_curve" %in% names(thresholded), TRUE)
  expect_identical(decision_curve$status, "incomplete")
  expect_identical(decision_curve$task_type, "survival")
  expect_identical(
    decision_curve$scope,
    "paired_patient_bootstrap_external_survival_decision_curve"
  )
  expect_identical(
    decision_curve$interval_scope,
    "pointwise_not_simultaneous"
  )
  expect_identical(
    decision_curve$event_probability_estimator,
    "kaplan_meier_within_threshold_selected_group_each_resample"
  )
  expect_identical(decision_curve$prediction, "fixed")
  expect_identical(
    decision_curve$threshold_source,
    "study_declared_survival_task"
  )
  expect_identical(decision_curve$comparison_targets, targets)
  expect_identical(decision_curve$prediction_horizon, 4)
  expect_identical(decision_curve$time_unit, "month")
  expect_identical(decision_curve$n_thresholds, 1L)
  expect_identical(
    decision_curve$point_comparisons$comparison,
    targets
  )
  expect_equal(
    decision_curve$point_comparisons$estimate,
    c(0.5, -1 / 6),
    tolerance = 1e-15
  )
  expect_identical(
    decision_curve$point_comparisons$status,
    rep("success", 2L)
  )
  expect_identical(nrow(decision_curve$replicates), 12L)
  expect_identical(decision_curve$replicates$replicate, seq_len(12L))
  expect_identical(
    any(c("patient_id", "indices") %in% names(decision_curve$replicates)),
    FALSE
  )
  zero_before_horizon <-
    decision_curve$replicates$horizon_support ==
      "survival_zero_before_horizon"
  expect_identical(any(zero_before_horizon), TRUE)
  expect_identical(
    decision_curve$replicates$model_minus_treat_all_status[
      zero_before_horizon
    ],
    rep("success", sum(zero_before_horizon))
  )
  all_event_replicate <- thresholded$resample_audit$replicate[
    thresholded$resample_audit$n_events == thresholded$n
  ][[1L]]
  all_event_curve <- decision_curve$replicates[
    decision_curve$replicates$replicate == all_event_replicate,
    ,
    drop = FALSE
  ]
  expect_identical(all_event_curve$status, "success")

  indices <- withr::with_seed(1L, {
    lapply(seq_len(12L), function(unused) {
      sample.int(length(args$time), length(args$time), replace = TRUE)
    })
  })
  first <- indices[[1L]]
  sampled_time <- args$time[first]
  sampled_event <- args$event[first]
  sampled_survival <- args$survival_probability[first]
  selected <- 1 - sampled_survival >= 0.5
  manual_survival <- function(time, event, horizon) {
    event_times <- sort(unique(time[event == 1L & time <= horizon]))
    if (length(event_times) == 0L) {
      return(1)
    }
    prod(vapply(event_times, function(value) {
      1 - sum(time == value & event == 1L) / sum(time >= value)
    }, numeric(1)))
  }
  selected_fraction <- sum(selected) / length(selected)
  selected_survival <- manual_survival(
    sampled_time[selected],
    sampled_event[selected],
    4
  )
  overall_survival <- manual_survival(sampled_time, sampled_event, 4)
  true_positive_rate <- selected_fraction * (1 - selected_survival)
  false_positive_rate <- selected_fraction * selected_survival
  model_net_benefit <- true_positive_rate - false_positive_rate
  treat_all_net_benefit <- (1 - overall_survival) - overall_survival
  first_row <- decision_curve$replicates[
    decision_curve$replicates$replicate == 1L,
    ,
    drop = FALSE
  ]

  expect_identical(first_row$n_selected, as.integer(sum(selected)))
  expect_equal(first_row$selected_fraction, selected_fraction, tolerance = 0)
  expect_equal(
    first_row$selected_survival_probability,
    selected_survival,
    tolerance = 1e-15
  )
  expect_equal(
    first_row$overall_survival_probability,
    overall_survival,
    tolerance = 1e-15
  )
  expect_equal(
    first_row$estimated_true_positive_rate,
    true_positive_rate,
    tolerance = 1e-15
  )
  expect_equal(
    first_row$estimated_false_positive_rate,
    false_positive_rate,
    tolerance = 1e-15
  )
  expect_equal(first_row$model_net_benefit, model_net_benefit, tolerance = 1e-15)
  expect_equal(
    first_row$treat_all_net_benefit,
    treat_all_net_benefit,
    tolerance = 1e-15
  )
  expect_equal(
    first_row$model_minus_treat_all,
    model_net_benefit - treat_all_net_benefit,
    tolerance = 1e-15
  )
  expect_equal(
    first_row$model_minus_treat_none,
    model_net_benefit,
    tolerance = 1e-15
  )
  expect_identical(
    thresholded$resample_audit$n_unique_patients[[1L]],
    as.integer(length(unique(first)))
  )

  probabilities <- thresholded$interval_method$probabilities
  for (comparison in targets) {
    interval <- decision_curve$intervals[
      decision_curve$intervals$comparison == comparison,
      ,
      drop = FALSE
    ]
    replicate_rows <- decision_curve$replicates[
      decision_curve$replicates[[paste0(comparison, "_status")]] == "success" &
        is.finite(decision_curve$replicates[[comparison]]),
      ,
      drop = FALSE
    ]
    expected_bounds <- unname(stats::quantile(
      replicate_rows[[comparison]],
      probs = probabilities,
      type = 7
    ))
    expect_identical(interval$status, "incomplete")
    expect_identical(interval$interval_available, TRUE)
    expect_identical(interval$replicate_completeness, "partial")
    expect_identical(interval$n_requested, 12L)
    expect_identical(interval$n_finite_successful, 8L)
    expect_identical(interval$n_unavailable, 4L)
    expect_identical(interval$n_nonfinite, 0L)
    expect_equal(interval$valid_fraction, 8 / 12, tolerance = 0)
    expect_equal(
      c(interval$lower, interval$upper),
      expected_bounds,
      tolerance = 0
    )
    expect_identical(
      interval$error_class,
      "stormbreaker_incomplete_bootstrap_replicates"
    )
  }
  public_messages <- na.omit(decision_curve$intervals$error_message)
  expect_identical(any(grepl("survival_dca_", public_messages)), FALSE)
  expect_match(
    paste(capture.output(print(thresholded)), collapse = "\n"),
    paste0(
      "External survival DCA paired pointwise intervals: ",
      "1 declared thresholds \\| incomplete"
    )
  )
})

test_that("survival DCA bootstrap recomputes curves and rejects bad thresholds", {
  evaluation <- make_sb_survival_dca_bootstrap_evaluation(
    time = c(1, 2, 5, 5, 4, 6),
    event = c(1L, 0L, 0L, 1L, 0L, 1L),
    risk = c(0.9, 0.7, 0.6, 0.4, 0.2, 0.1),
    survival_probability = c(0.1, 0.1, 0.1, 0.9, 0.9, 0.9),
    prediction_horizon = 4,
    decision_thresholds = c(0.2, 0.5)
  )
  changed <- evaluation
  changed$metrics$decision_curve$curve$model_net_benefit <- 999
  changed$metrics$decision_curve$curve$treat_all_net_benefit <- 998
  changed$metrics$decision_curve$curve$treat_none_net_benefit <- 997
  changed$metrics$decision_curve$curve$status <- "not_evaluable"

  expected <- sb_bootstrap_external(
    evaluation,
    reps = 3L,
    confidence_level = 0.8,
    seed = 17L
  )
  observed <- sb_bootstrap_external(
    changed,
    reps = 3L,
    confidence_level = 0.8,
    seed = 17L
  )

  expect_identical(observed$decision_curve, expected$decision_curve)
  expect_identical(observed$point_estimates, expected$point_estimates)
  expect_identical(observed$replicates, expected$replicates)

  invalid <- evaluation
  invalid$metrics$decision_curve$curve$threshold_probability <- c(0.2, 1)
  failed <- sb_bootstrap_external(invalid, reps = 2L, seed = 17L)
  expect_identical(failed$status, "failed")
  expect_identical(failed$error_stage, "evaluation_contract")
  expect_identical(failed$error_class, "stormbreaker_bad_decision_thresholds")
  expect_identical("decision_curve" %in% names(failed), FALSE)
})

test_that("risk-only and threshold-free survival keep the old result shape", {
  threshold_free <- make_sb_survival_dca_bootstrap_evaluation(
    time = c(1, 2, 5, 5, 4, 6),
    event = c(1L, 0L, 0L, 1L, 0L, 1L),
    risk = c(0.9, 0.7, 0.6, 0.4, 0.2, 0.1),
    survival_probability = c(0.1, 0.1, 0.1, 0.9, 0.9, 0.9),
    prediction_horizon = 4
  )
  risk_only <- make_sb_survival_dca_bootstrap_evaluation(
    time = c(1, 2, 3, 4, 5, 6),
    event = c(1L, 0L, 1L, 0L, 1L, 0L),
    risk = c(0.9, 0.7, 0.6, 0.4, 0.2, 0.1)
  )
  source_dca <- make_sb_survival_dca_bootstrap_evaluation(
    time = c(1, 2, 5, 5, 4, 6),
    event = c(1L, 0L, 0L, 1L, 0L, 1L),
    risk = c(0.9, 0.7, 0.6, 0.4, 0.2, 0.1),
    survival_probability = c(0.1, 0.1, 0.1, 0.9, 0.9, 0.9),
    prediction_horizon = 4,
    decision_thresholds = 0.5
  )
  risk_only$metrics$decision_curve <- source_dca$metrics$decision_curve

  threshold_free_result <- sb_bootstrap_external(
    threshold_free,
    reps = 2L,
    seed = 5L
  )
  risk_only_result <- sb_bootstrap_external(risk_only, reps = 2L, seed = 5L)

  expect_identical("decision_curve" %in% names(threshold_free_result), FALSE)
  expect_identical("decision_curve" %in% names(risk_only_result), FALSE)
})

test_that("paired survival DCA contrast avoids near-one cancellation", {
  evaluation <- make_sb_survival_dca_bootstrap_evaluation(
    time = c(rep(2, 99), 1),
    event = c(rep(0L, 99), 1L),
    risk = seq_len(100),
    survival_probability = c(rep(0, 99), 1),
    prediction_horizon = 1.5,
    decision_thresholds = 1 - .Machine$double.eps
  )

  result <- sb_bootstrap_external(
    evaluation,
    reps = 2L,
    confidence_level = 0.8,
    seed = 7L
  )
  row <- result$decision_curve$point_curve
  naive <- row$model_net_benefit - row$treat_all_net_benefit
  direct <-
    (row$estimated_true_positive_rate - row$overall_event_probability) -
    (row$estimated_false_positive_rate -
       row$overall_survival_probability) * row$threshold_odds

  expect_equal(row$selected_fraction, 0.99, tolerance = 0)
  expect_equal(row$selected_survival_probability, 1, tolerance = 0)
  expect_equal(row$overall_survival_probability, 0.99, tolerance = 1e-15)
  expect_equal(naive, 0, tolerance = 0)
  expect_equal(direct, -0.01, tolerance = 1e-15)
  expect_equal(row$model_minus_treat_all, -0.01, tolerance = 1e-15)
})

test_that("survival DCA support is specific to each paired contrast", {
  evaluation <- make_sb_survival_dca_bootstrap_evaluation(
    time = c(1, 2, 4),
    event = c(1L, 0L, 0L),
    risk = c(0.9, 0.5, 0.1),
    survival_probability = c(0, 1, 1),
    prediction_horizon = 3,
    decision_thresholds = 0.5
  )
  scored <- sb_external_bootstrap_contract(evaluation)$scored
  row <- sb_external_bootstrap_survival_decision_curve_replicate(
    scored[c(1L, 2L, 2L), , drop = FALSE],
    0.5,
    1L
  )

  expect_identical(row$status, "not_evaluable")
  expect_equal(row$model_net_benefit, 1 / 3, tolerance = 1e-15)
  expect_identical(is.na(row$treat_all_net_benefit), TRUE)
  expect_identical(row$model_minus_treat_none_status, "success")
  expect_equal(row$model_minus_treat_none, 1 / 3, tolerance = 1e-15)
  expect_identical(row$model_minus_treat_all_status, "not_evaluable")
  expect_identical(is.na(row$model_minus_treat_all), TRUE)
  expect_identical(
    row$model_minus_treat_all_error_class,
    "stormbreaker_insufficient_survival_dca_follow_up"
  )
})

test_that("identical survival DCA strategies have exact zero contrasts", {
  common <- list(
    time = c(1, 2, 2),
    event = c(0L, 0L, 0L),
    risk = c(0.9, 0.5, 0.1),
    prediction_horizon = 3,
    decision_thresholds = 0.5
  )
  all_selected <- do.call(
    make_sb_survival_dca_bootstrap_evaluation,
    c(common, list(survival_probability = c(0, 0, 0)))
  )
  none_selected <- do.call(
    make_sb_survival_dca_bootstrap_evaluation,
    c(common, list(survival_probability = c(1, 1, 1)))
  )
  all_result <- sb_bootstrap_external(
    all_selected,
    reps = 4L,
    confidence_level = 0.8,
    seed = 9L
  )
  none_result <- sb_bootstrap_external(
    none_selected,
    reps = 4L,
    confidence_level = 0.8,
    seed = 9L
  )

  all_point <- all_result$decision_curve$point_curve
  all_interval <- all_result$decision_curve$intervals[
    all_result$decision_curve$intervals$comparison ==
      "model_minus_treat_all",
    ,
    drop = FALSE
  ]
  expect_identical(is.na(all_point$model_net_benefit), TRUE)
  expect_identical(is.na(all_point$treat_all_net_benefit), TRUE)
  expect_identical(all_point$model_minus_treat_all, 0)
  expect_identical(all_point$model_minus_treat_all_status, "success")
  expect_identical(
    all_result$decision_curve$replicates$model_minus_treat_all,
    rep(0, 4L)
  )
  expect_identical(all_interval$status, "success")
  expect_identical(all_interval$n_finite_successful, 4L)
  expect_identical(c(all_interval$lower, all_interval$upper), c(0, 0))

  none_point <- none_result$decision_curve$point_curve
  none_interval <- none_result$decision_curve$intervals[
    none_result$decision_curve$intervals$comparison ==
      "model_minus_treat_none",
    ,
    drop = FALSE
  ]
  expect_identical(none_point$model_net_benefit, 0)
  expect_identical(is.na(none_point$treat_all_net_benefit), TRUE)
  expect_identical(none_point$model_minus_treat_none, 0)
  expect_identical(none_point$model_minus_treat_none_status, "success")
  expect_identical(
    none_result$decision_curve$replicates$model_minus_treat_none,
    rep(0, 4L)
  )
  expect_identical(none_interval$status, "success")
  expect_identical(none_interval$n_finite_successful, 4L)
  expect_identical(c(none_interval$lower, none_interval$upper), c(0, 0))
  expect_identical(all_result$decision_curve$status, "incomplete")
  expect_identical(none_result$decision_curve$status, "incomplete")
  expect_identical(all_result$resample_audit$n_events, rep(0L, 4L))
  expect_identical(none_result$resample_audit$n_events, rep(0L, 4L))
})

test_that("survival DCA intervals expose finite-replicate availability", {
  evaluation <- make_sb_survival_dca_bootstrap_evaluation(
    time = c(1, 2, 5, 5, 4, 6),
    event = c(1L, 0L, 0L, 1L, 0L, 1L),
    risk = c(0.9, 0.7, 0.6, 0.4, 0.2, 0.1),
    survival_probability = c(0.1, 0.1, 0.1, 0.9, 0.9, 0.9),
    prediction_horizon = 4,
    decision_thresholds = 0.5
  )
  scored <- sb_external_bootstrap_contract(evaluation)$scored
  successful <- sb_external_bootstrap_survival_decision_curve_replicate(
    scored,
    0.5,
    1L
  )
  unavailable <- sb_external_bootstrap_survival_decision_curve_replicate(
    scored[c(1L, 2L, 4L, 4L, 4L, 4L), , drop = FALSE],
    0.5,
    2L
  )
  result <- sb_external_bootstrap_survival_decision_curve_result(
    scored,
    0.5,
    rbind(successful, unavailable),
    0.8,
    2L
  )

  expect_identical(result$intervals$n_finite_successful, rep(1L, 2L))
  expect_identical(result$intervals$n_unavailable, rep(1L, 2L))
  expect_identical(result$intervals$interval_available, rep(FALSE, 2L))
  expect_identical(
    result$intervals$error_class,
    rep("stormbreaker_insufficient_bootstrap_replicates", 2L)
  )
})

test_that("survival DCA thresholds retain declared order and serialize", {
  thresholds <- c(0.2, 0.5, 0.8)
  result <- sb_bootstrap_external(
    make_sb_survival_dca_bootstrap_evaluation(
      time = c(1, 2, 5, 5, 4, 6),
      event = c(1L, 0L, 0L, 1L, 0L, 1L),
      risk = c(0.9, 0.7, 0.6, 0.4, 0.2, 0.1),
      survival_probability = c(0.1, 0.1, 0.1, 0.9, 0.9, 0.9),
      prediction_horizon = 4,
      decision_thresholds = thresholds
    ),
    reps = 2L,
    confidence_level = 0.8,
    seed = 2L
  )
  curve <- result$decision_curve

  expect_identical(curve$point_curve$threshold_index, seq_along(thresholds))
  expect_identical(curve$point_curve$threshold_probability, thresholds)
  expect_identical(
    curve$point_comparisons$threshold_index,
    rep(seq_along(thresholds), each = 2L)
  )
  expect_identical(
    curve$point_comparisons$comparison,
    rep(curve$comparison_targets, times = length(thresholds))
  )
  expect_identical(
    curve$replicates$replicate,
    rep(seq_len(2L), each = length(thresholds))
  )
  expect_identical(
    curve$replicates$threshold_index,
    rep(seq_along(thresholds), times = 2L)
  )
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(result, path)
  expect_identical(readRDS(path), result)
})

test_that("survival DCA bootstrap keeps fixed probabilities and avoids refitting", {
  evaluation <- make_sb_survival_dca_bootstrap_evaluation(
    time = c(1, 2, 5, 5, 4, 6),
    event = c(1L, 0L, 0L, 1L, 0L, 1L),
    risk = c(0.9, 0.7, 0.6, 0.4, 0.2, 0.1),
    survival_probability = c(0.1, 0.1, 0.1, 0.9, 0.9, 0.9),
    prediction_horizon = 4,
    decision_thresholds = 0.5
  )
  reps <- 4L
  seed <- 1L
  original_dca <- sb_survival_decision_curve_metric
  original_event_km <- sb_fit_event_km
  dca_calls <- 0L
  event_km_calls <- 0L
  seen_probabilities <- list()
  testthat::local_mocked_bindings(
    sb_survival_decision_curve_metric = function(predictions, ...) {
      dca_calls <<- dca_calls + 1L
      seen_probabilities[[dca_calls]] <<- predictions$survival_probability
      original_dca(predictions, ...)
    },
    sb_fit_event_km = function(time, event) {
      event_km_calls <<- event_km_calls + 1L
      original_event_km(time, event)
    },
    sb_fit_final = function(...) stop("must not refit"),
    sb_predict_final = function(...) stop("must not repredict"),
    sb_evaluate_external = function(...) stop("must not reevaluate")
  )

  result <- sb_bootstrap_external(evaluation, reps = reps, seed = seed)
  probabilities <- evaluation$scored_predictions$survival_probability
  indices <- withr::with_seed(seed, {
    lapply(seq_len(reps), function(unused) {
      sample.int(length(probabilities), length(probabilities), replace = TRUE)
    })
  })

  expect_identical(dca_calls, reps + 1L)
  expect_identical(event_km_calls, 3L * (reps + 1L))
  for (replicate in seq_len(reps)) {
    expect_identical(
      seen_probabilities[[replicate]],
      probabilities[indices[[replicate]]]
    )
  }
  expect_identical(seen_probabilities[[reps + 1L]], probabilities)
  expect_identical(result$decision_curve$prediction, "fixed")
})
