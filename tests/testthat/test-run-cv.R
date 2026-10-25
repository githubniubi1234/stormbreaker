test_that("five-fold CV gives every patient one finite OOF risk", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  result <- sb_run_cv(task, sb_pipeline(), folds = 5, seed = 5201314)

  expect_s3_class(result, "stormbreaker_cv_result")
  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_equal(nrow(result$predictions), nrow(d))
  expect_identical(anyDuplicated(result$predictions$patient_id), 0L)
  expect_setequal(result$predictions$patient_id, d$patient_id)
  expect_equal(result$predictions$time, d$time)
  expect_identical(result$predictions$event, d$event)
  expect_identical(unique(result$predictions$time_unit), "month")
  expect_identical(all(is.finite(result$predictions$risk)), TRUE)
  expect_identical(unique(result$predictions$status), "success")
  expect_identical(nrow(result$runs), 5L)
  expect_identical(unique(result$runs$learner_used), "coxph")
  expect_length(result$coefficients, 5L)
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("success", 5L)
  )
  expect_named(result$metrics, "harrell_c")
  expect_identical(result$metrics$harrell_c$status, "success")
  expect_identical(result$metrics$harrell_c$n_folds, 5L)
  expect_identical(result$metrics$harrell_c$n_evaluable_folds, 5L)
  expect_identical(
    result$metrics$harrell_c$per_fold$fold,
    seq_len(5L)
  )
  expect_identical(
    result$metrics$harrell_c$per_fold$status,
    rep("success", 5L)
  )
  official <- survival::concordance(
    survival::Surv(time, event) ~ risk + survival::strata(fold),
    data = result$predictions,
    reverse = TRUE,
    timewt = "n"
  )
  expect_equal(
    result$metrics$harrell_c$estimate,
    unname(official$concordance),
    tolerance = 1e-15
  )
  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(printed, "Fold-stratified Harrell C:")
  expect_match(printed, "Comparable pairs:")
  expect_identical(
    result$runs$n_analysis_events + result$runs$n_assessment_events,
    rep.int(sum(d$event), 5L)
  )
})

test_that("the public result schema includes audited metrics", {
  d <- make_sb_fixture(n = 20)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  result <- sb_run_cv(task, folds = 2, seed = 2)

  expect_named(
    result,
    c(
      "status", "complete_oof", "predictions", "metrics", "runs",
      "selection_stability", "coefficients", "folds", "task", "pipeline",
      "provenance", "glmnet_active_set_stability"
    ),
    ignore.order = FALSE
  )
  expect_identical(result$selection_stability$status, "not_applicable")
  expect_identical(result$selection_stability$n_folds, 2L)
  expect_false("weak_calibration" %in% names(result$metrics))
  expect_identical(result$selection_stability$n_successful_folds, 0L)
  expect_true(is.na(result$selection_stability$mean_pairwise_jaccard))
  expect_identical(
    result$glmnet_active_set_stability$status,
    "not_applicable"
  )
  expect_true(is.na(result$glmnet_active_set_stability$alpha))
  expect_named(
    result$predictions,
    c(
      "patient_id", "fold", "time", "event", "time_unit", "risk",
      "status", "pipeline_id", "learner"
    ),
    ignore.order = FALSE
  )
  expect_named(
    result$runs,
    c(
      "fold", "seed", "pipeline_id", "selector_requested",
      "learner_requested", "learner_used", "n_analysis",
      "n_analysis_events", "n_assessment", "n_assessment_events",
      "status", "error_class", "error_message"
    ),
    ignore.order = FALSE
  )
})

test_that("a declared horizon adds reconstructible Cox OOF survival probability", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    )
  )

  for (pipeline in pipelines) {
    result <- sb_run_cv(task, pipeline, folds = 5L, seed = 20260819L)

    expect_identical(result$status, "success")
    expect_identical(result$complete_oof, TRUE)
    expect_identical(result$complete_survival_probability_oof, TRUE)
    expect_named(
      result$metrics,
      c(
        "harrell_c", "ipcw_brier", "time_auc", "uno_c", "calibration",
        "weak_calibration"
      )
    )
    expect_identical(result$metrics$ipcw_brier$status, "success")
    expect_identical(result$metrics$ipcw_brier$prediction_horizon, 7)
    expect_true(is.finite(result$metrics$ipcw_brier$estimate))
    expect_identical(
      result$metrics$ipcw_brier$reference$status,
      "success"
    )
    expect_identical(
      result$metrics$ipcw_brier$skill_score$status,
      "success"
    )
    expect_true(is.finite(
      result$metrics$ipcw_brier$skill_score$estimate
    ))
    expect_identical(result$metrics$time_auc$status, "success")
    expect_identical(result$metrics$time_auc$prediction_horizon, 7)
    expect_true(is.finite(result$metrics$time_auc$estimate))
    expect_gt(result$metrics$time_auc$pair_denominator, 0)
    expect_identical(result$metrics$time_auc$n_evaluable_folds, 5L)
    expect_identical(result$metrics$uno_c$status, "success")
    expect_identical(result$metrics$uno_c$prediction_horizon, 7)
    expect_true(is.finite(result$metrics$uno_c$estimate))
    expect_gt(result$metrics$uno_c$weighted_pair_denominator, 0)
    expect_identical(result$metrics$uno_c$n_evaluable_folds, 5L)
    calibration <- result$metrics$calibration
    expect_identical(calibration$status, "success")
    expect_identical(
      calibration$scope,
      "pooled_patient_level_out_of_fold"
    )
    expect_identical(calibration$purpose, "descriptive_evaluation_only")
    expect_identical(calibration$prediction_horizon, 7)
    expect_identical(calibration$time_unit, "month")
    expect_identical(
      calibration$patient_calibration$patient_id,
      result$predictions$patient_id
    )
    direct_calibration <- sb_survival_calibration_metric(
      patient_id = result$predictions$patient_id,
      time = result$predictions$time,
      event = result$predictions$event,
      survival_probability = result$predictions$survival_probability,
      prediction_status = ifelse(
        result$predictions$status == "success" &
          result$predictions$survival_probability_status == "success",
        "success",
        "failed"
      ),
      prediction_horizon = result$task$prediction_horizon,
      time_unit = result$task$time_unit,
      scope = "pooled_patient_level_out_of_fold",
      fold = result$predictions$fold
    )
    expect_identical(calibration$knots, direct_calibration$knots)
    expect_identical(calibration$curve, direct_calibration$curve)
    expect_identical(
      calibration$patient_calibration$absolute_error,
      direct_calibration$patient_calibration$absolute_error
    )
    expect_identical(
      calibration[c("ici", "e50", "e90")],
      direct_calibration[c("ici", "e50", "e90")]
    )
    weak_calibration <- result$metrics$weak_calibration
    direct_weak_calibration <- sb_survival_weak_calibration_metric(
      patient_id = result$predictions$patient_id,
      time = result$predictions$time,
      event = result$predictions$event,
      survival_probability = result$predictions$survival_probability,
      prediction_status = ifelse(
        result$predictions$status == "success" &
          result$predictions$survival_probability_status == "success",
        "success",
        "failed"
      ),
      prediction_horizon = result$task$prediction_horizon,
      time_unit = result$task$time_unit,
      scope = "pooled_patient_level_out_of_fold",
      fold = result$predictions$fold
    )
    expect_identical(weak_calibration$status, "success")
    expect_identical(
      weak_calibration$patient_predictions$patient_id,
      result$predictions$patient_id
    )
    expect_identical(weak_calibration$per_fold$fold, seq_len(5L))
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
    expect_named(
      result$predictions,
      c(
        "patient_id", "fold", "time", "event", "time_unit", "risk",
        "prediction_horizon", "survival_probability",
        "survival_probability_status", "status", "pipeline_id", "learner"
      ),
      ignore.order = FALSE
    )
    expect_identical(
      unique(result$predictions$prediction_horizon),
      7
    )
    expect_identical(
      unique(result$predictions$survival_probability_status),
      "success"
    )
    expect_true(all(
      result$predictions$survival_probability >= 0 &
        result$predictions$survival_probability <= 1
    ))
    expect_length(result$baseline_hazards, 5L)

    for (fold in seq_along(result$folds)) {
      rows <- result$predictions$fold == fold
      coefficient_audit <- result$coefficients[[fold]]
      baseline_audit <- result$baseline_hazards[[fold]]
      x <- d[
        match(
          result$predictions$patient_id[rows],
          d$patient_id
        ),
        coefficient_audit$values$feature,
        drop = FALSE
      ]
      reconstructed_risk <- unname(as.vector(
        as.matrix(x) %*% coefficient_audit$values$coefficient
      ))
      reconstructed_probability <- exp(
        -baseline_audit$cumulative_baseline_hazard *
          exp(reconstructed_risk)
      )

      expect_identical(baseline_audit$status, "success")
      expect_identical(baseline_audit$ties, "efron")
      expect_identical(
        baseline_audit$source,
        "coxph_basehaz_zero_reference"
      )
      expect_equal(
        result$predictions$risk[rows],
        reconstructed_risk,
        tolerance = 1e-12
      )
      expect_equal(
        result$predictions$survival_probability[rows],
        reconstructed_probability,
        tolerance = 1e-12
      )
      reference_curve <-
        result$metrics$ipcw_brier$reference$survival_curves[[fold]]
      reference_per_fold <- result$metrics$ipcw_brier$reference$per_fold
      expect_equal(
        reference_per_fold$reference_survival_horizon[[fold]],
        sb_event_survival_at(reference_curve, 7),
        tolerance = 0
      )
    }
    printed <- paste(capture.output(print(result)), collapse = "\n")
    expect_match(printed, "IPCW Brier reference:", fixed = TRUE)
    expect_match(printed, "IPCW Brier skill score:", fixed = TRUE)
    expect_match(printed, "Time-dependent AUC:", fixed = TRUE)
    expect_match(printed, "Fold-stratified Uno C:", fixed = TRUE)
    expect_match(printed, "Survival calibration ICI:", fixed = TRUE)
    expect_match(printed, "Mean calibration O:E:", fixed = TRUE)
    expect_match(printed, "Survival calibration slope:", fixed = TRUE)
    expect_match(printed, "E50:", fixed = TRUE)
    expect_match(printed, "E90:", fixed = TRUE)
  }
})

test_that("survival weak calibration exceptions stay inside the optional metric", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  baseline <- sb_run_cv(task, sb_pipeline(), folds = 5L, seed = 20260819L)
  testthat::local_mocked_bindings(
    sb_survival_weak_calibration_metric = function(...) {
      condition <- simpleError("weak calibration exploded")
      class(condition) <- c("test_weak_calibration_error", class(condition))
      stop(condition)
    }
  )

  observed <- sb_run_cv(task, sb_pipeline(), folds = 5L, seed = 20260819L)

  expect_identical(observed$status, baseline$status)
  expect_identical(observed$predictions, baseline$predictions)
  expect_identical(observed$folds, baseline$folds)
  expect_identical(observed$coefficients, baseline$coefficients)
  expect_identical(observed$baseline_hazards, baseline$baseline_hazards)
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
    "test_weak_calibration_error"
  )
  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "Mean calibration O:E: unavailable", fixed = TRUE)
  expect_match(
    printed,
    "Survival calibration slope: unavailable",
    fixed = TRUE
  )
})

test_that("survival calibration incompleteness cannot alter the CV run", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  baseline <- sb_run_cv(task, sb_pipeline(), folds = 5L, seed = 20260819L)
  original_fit <- sb_fit_coxph
  learner_calls <- 0L
  testthat::local_mocked_bindings(
    sb_fit_coxph = function(...) {
      learner_calls <<- learner_calls + 1L
      original_fit(...)
    },
    sb_survival_calibration_metric = function(...) {
      metric <- sb_survival_calibration_empty_result(
        prediction_horizon = 7,
        time_unit = "month",
        scope = "pooled_patient_level_out_of_fold",
        n = nrow(d)
      )
      metric$input_status <- "not_evaluable"
      metric$error_class <- "test_calibration_unavailable"
      metric$error_message <- "deliberately unavailable"
      metric
    }
  )

  observed <- sb_run_cv(task, sb_pipeline(), folds = 5L, seed = 20260819L)

  expect_identical(learner_calls, 5L)
  expect_identical(observed$status, baseline$status)
  expect_identical(observed$predictions, baseline$predictions)
  expect_identical(observed$folds, baseline$folds)
  expect_identical(observed$coefficients, baseline$coefficients)
  expect_identical(
    observed$metrics[c("harrell_c", "ipcw_brier", "time_auc", "uno_c")],
    baseline$metrics[c("harrell_c", "ipcw_brier", "time_auc", "uno_c")]
  )
  expect_identical(observed$metrics$calibration$status, "incomplete")
})

test_that("survival calibration exceptions stay inside the optional metric", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  baseline <- sb_run_cv(task, sb_pipeline(), folds = 5L, seed = 20260819L)
  testthat::local_mocked_bindings(
    sb_survival_calibration_metric = function(...) {
      condition <- simpleError("calibration exploded")
      class(condition) <- c("test_calibration_error", class(condition))
      stop(condition)
    }
  )

  observed <- sb_run_cv(task, sb_pipeline(), folds = 5L, seed = 20260819L)

  expect_identical(observed$status, baseline$status)
  expect_identical(observed$predictions, baseline$predictions)
  expect_identical(observed$folds, baseline$folds)
  expect_identical(observed$coefficients, baseline$coefficients)
  expect_identical(
    observed$metrics[c("harrell_c", "ipcw_brier", "time_auc", "uno_c")],
    baseline$metrics[c("harrell_c", "ipcw_brier", "time_auc", "uno_c")]
  )
  expect_identical(observed$metrics$calibration$status, "metric_failed")
  expect_identical(
    observed$metrics$calibration$error_class,
    "test_calibration_error"
  )
  expect_identical(
    observed$metrics$calibration$error_message,
    "calibration exploded"
  )
  expect_match(
    paste(capture.output(print(observed)), collapse = "\n"),
    "Survival calibration: unavailable",
    fixed = TRUE
  )
})

test_that("horizon CV shares one censoring fit without Brier status coupling", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  original_fit <- sb_fit_censoring_km
  original_brier <- sb_survival_ipcw_brier_metric
  censoring_calls <- 0L
  testthat::local_mocked_bindings(
    sb_fit_censoring_km = function(...) {
      censoring_calls <<- censoring_calls + 1L
      original_fit(...)
    },
    sb_survival_ipcw_brier_metric = function(...) {
      metric <- original_brier(...)
      metric$status <- "incomplete"
      metric$estimate <- NA_real_
      metric
    }
  )

  result <- sb_run_cv(task, sb_pipeline(), folds = 5L, seed = 20260819L)

  expect_identical(censoring_calls, 5L)
  expect_identical(result$status, "failed")
  expect_identical(result$metrics$ipcw_brier$status, "incomplete")
  expect_identical(result$metrics$time_auc$status, "success")
  expect_true(is.finite(result$metrics$time_auc$estimate))
  expect_identical(result$metrics$uno_c$status, "success")
  expect_true(is.finite(result$metrics$uno_c$estimate))
})

test_that("glmnet Cox OOF survival probability uses saved Breslow baselines", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 101L
  )
  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 20260819L)

  expect_identical(result$status, "success")
  expect_identical(result$complete_survival_probability_oof, TRUE)
  for (fold in seq_along(result$folds)) {
    rows <- result$predictions$fold == fold
    coefficient_audit <- result$coefficients[[fold]]
    baseline_audit <- result$baseline_hazards[[fold]]
    x <- d[
      match(result$predictions$patient_id[rows], d$patient_id),
      coefficient_audit$values$feature,
      drop = FALSE
    ]
    reconstructed_risk <- unname(as.vector(
      as.matrix(x) %*% coefficient_audit$values$coefficient
    ))
    reconstructed_probability <- exp(
      -baseline_audit$cumulative_baseline_hazard *
        exp(reconstructed_risk)
    )

    expect_identical(baseline_audit$status, "success")
    expect_identical(baseline_audit$ties, "breslow")
    expect_identical(baseline_audit$source, "glmnet_fixed_beta_breslow")
    expect_equal(
      result$predictions$survival_probability[rows],
      reconstructed_probability,
      tolerance = 1e-12
    )
  }
})

test_that("assessment outcomes cannot change fitted survival probabilities", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  predictors <- c("gene_a", "gene_b")
  make_task <- function(data) {
    sb_surv_task(
      data,
      "patient_id",
      "time",
      "event",
      predictors,
      "month",
      prediction_horizon = 5
    )
  }
  split <- list(
    fold = 1L,
    analysis_rows = seq_len(48L),
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[seq_len(48L)],
    assessment_ids = d$patient_id[49:60]
  )
  pipeline <- sb_pipeline()
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  original <- sb_run_one_fold(
    make_task(d), pipeline, runtime, split, seed = 101L
  )

  changed_outcomes <- d
  changed_outcomes$time[split$assessment_rows] <-
    changed_outcomes$time[split$assessment_rows] + 1000
  changed_outcomes$event[split$assessment_rows] <-
    1L - changed_outcomes$event[split$assessment_rows]
  outcome_perturbed <- sb_run_one_fold(
    make_task(changed_outcomes), pipeline, runtime, split, seed = 101L
  )

  expect_identical(
    original$baseline_hazard,
    outcome_perturbed$baseline_hazard
  )
  expect_identical(original$coefficients, outcome_perturbed$coefficients)
  expect_equal(
    original$predictions$survival_probability,
    outcome_perturbed$predictions$survival_probability,
    tolerance = 0
  )
  analysis_inputs <- list(list(
    fold = 1L,
    time = d$time[split$analysis_rows],
    event = d$event[split$analysis_rows]
  ))
  original_metric <- sb_survival_ipcw_brier_metric(
    original$predictions,
    analysis_inputs,
    prediction_horizon = 5
  )
  perturbed_metric <- sb_survival_ipcw_brier_metric(
    outcome_perturbed$predictions,
    analysis_inputs,
    prediction_horizon = 5
  )
  expect_identical(
    original_metric$censoring_curves,
    perturbed_metric$censoring_curves
  )
  expect_identical(
    original_metric$reference$survival_curves,
    perturbed_metric$reference$survival_curves
  )
  expect_identical(
    original_metric$reference$per_fold$reference_survival_horizon,
    perturbed_metric$reference$per_fold$reference_survival_horizon
  )
  expect_gt(
    abs(
      original_metric$weighted_loss_sum -
        perturbed_metric$weighted_loss_sum
    ),
    0
  )
  expect_gt(
    abs(
      original_metric$reference$weighted_loss_sum -
        perturbed_metric$reference$weighted_loss_sum
    ),
    0
  )

  changed_predictors <- d
  changed_predictors$gene_a[split$assessment_rows] <-
    changed_predictors$gene_a[split$assessment_rows] + 5
  predictor_perturbed <- sb_run_one_fold(
    make_task(changed_predictors), pipeline, runtime, split, seed = 101L
  )

  expect_identical(original$baseline_hazard, predictor_perturbed$baseline_hazard)
  expect_identical(original$coefficients, predictor_perturbed$coefficients)
  expect_gt(
    max(abs(
      original$predictions$survival_probability -
        predictor_perturbed$predictions$survival_probability
    )),
    0
  )
  predictor_metric <- sb_survival_ipcw_brier_metric(
    predictor_perturbed$predictions,
    analysis_inputs,
    prediction_horizon = 5
  )
  expect_identical(original_metric$censoring_curves, predictor_metric$censoring_curves)
  expect_identical(original_metric$reference, predictor_metric$reference)
  expect_gt(
    abs(
      original_metric$weighted_loss_sum -
        predictor_metric$weighted_loss_sum
    ),
    0
  )
})

test_that("survival-probability failure preserves successful risk evidence", {
  d <- make_sb_fixture(n = 30L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 5
  )
  split <- list(
    fold = 1L,
    analysis_rows = 1:20,
    assessment_rows = 21:30,
    analysis_ids = d$patient_id[1:20],
    assessment_ids = d$patient_id[21:30]
  )
  pipeline <- sb_pipeline()
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  runtime$learner$predict_survival <- function(...) {
    sb_abort(
      "Synthetic survival-probability failure.",
      "stormbreaker_test_survival_probability_failure"
    )
  }

  result <- sb_run_one_fold(task, pipeline, runtime, split, seed = 1L)

  expect_identical(result$run$status, "success")
  expect_true(all(is.finite(result$predictions$risk)))
  expect_identical(unique(result$predictions$status), "success")
  expect_true(all(is.na(result$predictions$survival_probability)))
  expect_identical(
    unique(result$predictions$survival_probability_status),
    "failed"
  )
  expect_identical(result$coefficients$status, "success")
  expect_identical(result$baseline_hazard$status, "failed")
  expect_identical(
    result$baseline_hazard$error_class,
    "stormbreaker_test_survival_probability_failure"
  )
  metric <- sb_survival_ipcw_brier_metric(
    result$predictions,
    list(list(
      fold = 1L,
      time = d$time[split$analysis_rows],
      event = d$event[split$analysis_rows]
    )),
    prediction_horizon = 5
  )
  expect_identical(metric$status, "incomplete")
  expect_identical(metric$reference$status, "success")
  expect_identical(metric$reference$n_evaluable_folds, 1L)
  expect_identical(metric$skill_score$status, "incomplete")
  expect_identical(is.na(metric$skill_score$estimate), TRUE)
  uno <- sb_survival_uno_c_metric(
    result$predictions,
    list(list(
      fold = 1L,
      time = d$time[split$analysis_rows],
      event = d$event[split$analysis_rows]
    )),
    prediction_horizon = 5,
    censoring_curves = metric$censoring_curves
  )
  expect_identical(uno$status, "success")
  expect_true(is.finite(uno$estimate))
})

test_that("an unavailable requested Uno C preserves all other CV evidence", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  original_uno <- sb_survival_uno_c_metric
  testthat::local_mocked_bindings(
    sb_survival_uno_c_metric = function(...) {
      metric <- original_uno(...)
      metric$status <- "incomplete"
      metric$estimate <- NA_real_
      metric
    }
  )

  result <- sb_run_cv(task, sb_pipeline(), folds = 5L, seed = 20260819L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(result$complete_survival_probability_oof, TRUE)
  expect_identical(result$metrics$harrell_c$status, "success")
  expect_identical(result$metrics$ipcw_brier$status, "success")
  expect_identical(result$metrics$time_auc$status, "success")
  expect_identical(result$metrics$uno_c$status, "incomplete")
  expect_true(is.na(result$metrics$uno_c$estimate))
})

test_that("analysis and assessment patients never overlap", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month"
  )
  result <- sb_run_cv(task, folds = 5, seed = 17)

  for (split in result$folds) {
    expect_length(
      intersect(split$analysis_ids, split$assessment_ids),
      0L
    )
    expect_setequal(
      c(split$analysis_ids, split$assessment_ids),
      d$patient_id
    )
  }
})

test_that("OOF risk is reproduced by fitting only the recorded analysis rows", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month"
  )
  result <- sb_run_cv(task, folds = 5, seed = 88)

  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    audit <- result$coefficients[[fold]]
    analysis <- d[d$patient_id %in% split$analysis_ids, , drop = FALSE]
    assessment <- d[d$patient_id %in% split$assessment_ids, , drop = FALSE]
    manual <- survival::coxph(
      survival::Surv(time, event) ~ gene_a + gene_b,
      data = analysis,
      ties = "efron",
      singular.ok = FALSE,
      x = TRUE,
      y = TRUE
    )
    expected <- unname(stats::predict(
      manual,
      newdata = assessment,
      type = "lp",
      reference = "zero"
    ))
    reconstructed <- unname(as.vector(
      as.matrix(assessment[audit$values$feature]) %*%
        audit$values$coefficient
    ))
    observed <- result$predictions$risk[
      match(assessment$patient_id, result$predictions$patient_id)
    ]
    expect_identical(audit$fold, as.integer(fold))
    expect_identical(audit$learner, "coxph")
    expect_identical(audit$status, "success")
    expect_true(is.na(audit$lambda))
    expect_true(is.na(audit$alpha))
    expect_identical(audit$values$feature, c("gene_a", "gene_b"))
    expect_equal(
      audit$values$coefficient,
      unname(stats::coef(manual)[audit$values$feature]),
      tolerance = 1e-15
    )
    expect_identical(audit$n_nonzero, as.integer(sum(audit$values$nonzero)))
    expect_identical(
      audit$nonzero_features,
      audit$values$feature[audit$values$nonzero]
    )
    expect_equal(observed, expected, tolerance = 1e-12)
    expect_equal(observed, reconstructed, tolerance = 1e-12)
  }
})

test_that("CV declares the common zero-reference linear risk scale", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  result <- sb_run_cv(task, folds = 5, seed = 88)

  expect_identical(result$pipeline$risk_type, "linear_predictor")
  expect_identical(result$pipeline$risk_reference, "zero")
  expect_identical(result$pipeline$risk_formula, "x_beta")
  expect_identical(
    result$pipeline$risk_direction,
    "higher_is_higher_hazard"
  )
})

test_that("the same seed reproduces folds and OOF risks", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month"
  )
  first <- sb_run_cv(task, folds = 5, seed = 101)
  second <- sb_run_cv(task, folds = 5, seed = 101)

  expect_identical(first$predictions, second$predictions)
  expect_identical(first$metrics, second$metrics)
  expect_identical(first$runs, second$runs)
  expect_identical(first$folds, second$folds)
})

test_that("fold assignment does not change the caller's random stream", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month"
  )

  set.seed(404)
  expected <- stats::runif(3)
  set.seed(404)
  sb_run_cv(task, folds = 5, seed = 101)
  observed <- stats::runif(3)

  expect_identical(observed, expected)
})

test_that("failed Cox fits remain visible and never use a fallback", {
  d <- make_sb_fixture(n = 60)
  d$gene_b <- d$gene_a
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month"
  )
  result <- sb_run_cv(task, folds = 5, seed = 9)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$status), "failed")
  expect_identical(unique(result$runs$learner_requested), "coxph")
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(all(nzchar(result$runs$error_message)), TRUE)
  expect_identical(all(is.na(result$predictions$risk)), TRUE)
  expect_identical(result$metrics$harrell_c$status, "incomplete")
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_fitted", 5L)
  )
  expect_identical(
    vapply(result$coefficients, function(audit) nrow(audit$values), integer(1)),
    rep(0L, 5L)
  )
  expect_identical(is.na(result$metrics$harrell_c$estimate), TRUE)
  expect_identical(
    result$metrics$harrell_c$per_fold$status,
    rep("prediction_failed", 5L)
  )
  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(
    printed,
    "Harrell C: unavailable (incomplete)",
    fixed = TRUE
  )
  expect_equal(result$predictions$time, d$time)
  expect_identical(result$predictions$event, d$event)
  expect_identical(unique(result$predictions$status), "failed")
})

test_that("folds require enough events and censored patients", {
  d <- make_sb_fixture(n = 20)
  d$event <- c(rep(1L, 4), rep(0L, 16))
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  error <- capture_sb_error(sb_run_cv(task, folds = 5, seed = 1))

  expect_s3_class(error, "stormbreaker_infeasible_folds")
})

test_that("logical event coding runs through stratified CV", {
  d <- make_sb_fixture(n = 20)
  d$event <- as.logical(d$event)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  result <- sb_run_cv(task, folds = 2, seed = 1)

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_setequal(unique(result$predictions$event), c(0L, 1L))
})

test_that("run arguments must use STORMBREAKER contracts", {
  d <- make_sb_fixture(n = 20)
  task_error <- capture_sb_error(sb_run_cv(d, folds = 2, seed = 1))
  expect_s3_class(task_error, "stormbreaker_bad_task")

  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  pipeline_error <- capture_sb_error(
    sb_run_cv(task, pipeline = list(learner = "coxph"), folds = 2, seed = 1)
  )
  expect_s3_class(pipeline_error, "stormbreaker_bad_pipeline")
})

test_that("mutated pipeline metadata cannot masquerade as the core pipeline", {
  d <- make_sb_fixture(n = 20)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  pipeline <- sb_pipeline()
  pipeline$id <- "forged_pipeline"
  error <- capture_sb_error(
    sb_run_cv(task, pipeline = pipeline, folds = 2, seed = 1)
  )

  expect_s3_class(error, "stormbreaker_unresolved_pipeline")
})

test_that("mutated risk scale metadata cannot masquerade as the core pipeline", {
  d <- make_sb_fixture(n = 20)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  pipeline <- sb_pipeline()
  pipeline$risk_reference <- "sample"
  error <- capture_sb_error(
    sb_run_cv(task, pipeline = pipeline, folds = 2, seed = 1)
  )

  expect_s3_class(error, "stormbreaker_unresolved_pipeline")
})

test_that("mutated task metadata is rejected before resampling", {
  d <- make_sb_fixture(n = 20)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  task$n <- task$n - 1L
  error <- capture_sb_error(sb_run_cv(task, folds = 2, seed = 1))

  expect_s3_class(error, "stormbreaker_mutated_task")
})

test_that("Top-K CV stores one complete selection audit per fold", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    selector_n_features = 2L
  )
  result <- sb_run_cv(task, pipeline, folds = 5, seed = 101)

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_length(result$selections, 5L)
  expect_length(result$coefficients, 5L)
  expect_identical(result$metrics$harrell_c$status, "success")
  for (selection in result$selections) {
    expect_identical(selection$status, "success")
    expect_length(selection$selected_features, 2L)
    expect_identical(nrow(selection$scores), 3L)
  }
  for (fold in seq_along(result$folds)) {
    selection <- result$selections[[fold]]
    audit <- result$coefficients[[fold]]
    assessment_rows <- result$folds[[fold]]$assessment_rows
    reconstructed <- unname(as.vector(
      as.matrix(d[
        assessment_rows,
        audit$values$feature,
        drop = FALSE
      ]) %*% audit$values$coefficient
    ))
    observed <- result$predictions$risk[
      match(
        result$folds[[fold]]$assessment_ids,
        result$predictions$patient_id
      )
    ]

    expect_identical(audit$status, "success")
    expect_identical(audit$values$feature, selection$selected_features)
    expect_equal(reconstructed, observed, tolerance = 1e-12)
  }
  stability <- result$selection_stability
  expected_count <- vapply(task$predictors, function(feature) {
    as.integer(sum(vapply(
      result$selections,
      function(selection) feature %in% selection$selected_features,
      logical(1)
    )))
  }, integer(1))
  expected_jaccard <- vapply(
    seq_len(nrow(stability$pairwise)),
    function(row) {
      selected_i <- result$selections[[
        stability$pairwise$fold_i[[row]]
      ]]$selected_features
      selected_j <- result$selections[[
        stability$pairwise$fold_j[[row]]
      ]]$selected_features
      length(intersect(selected_i, selected_j)) /
        length(union(selected_i, selected_j))
    },
    numeric(1)
  )
  expect_identical(stability$status, "success")
  expect_identical(stability$n_successful_folds, 5L)
  expect_identical(stability$feature_frequency$feature, task$predictors)
  expect_identical(
    stability$feature_frequency$selected_count,
    unname(expected_count)
  )
  expect_equal(stability$pairwise$jaccard, expected_jaccard, tolerance = 0)
  expect_equal(
    stability$mean_pairwise_jaccard,
    mean(expected_jaccard),
    tolerance = 0
  )
  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(printed, "Selection stability", fixed = TRUE)
})

test_that("failed survival selector folds keep stability incomplete", {
  d <- make_sb_fixture(n = 40L, seed = 20260819L)
  d$constant_a <- 1
  d$constant_b <- 2
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("constant_a", "constant_b"),
    "month"
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    selector_n_features = 1L
  )

  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 20260819L)
  stability <- result$selection_stability

  expect_identical(result$status, "failed")
  expect_identical(
    vapply(result$selections, `[[`, "", "status"),
    rep("failed", 4L)
  )
  expect_identical(stability$status, "incomplete")
  expect_identical(stability$n_successful_folds, 0L)
  expect_identical(stability$n_evaluable_pairs, 0L)
  expect_true(is.na(stability$mean_pairwise_jaccard))
  expect_identical(stability$feature_frequency$selected_count, c(0L, 0L))
  expect_true(all(is.na(stability$feature_frequency$selection_frequency)))
  expect_true(all(stability$per_fold$status == "failed"))
  expect_true(all(stability$pairwise$status == "incomplete"))
})

test_that("nested glmnet Cox tunes once inside every outer analysis fold", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260813)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipeline <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 1701L
  )
  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 101L)

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(result$metrics$harrell_c$status, "success")
  expect_identical(unique(result$runs$learner_used), "glmnet_cox")
  expect_length(result$tuning, 5L)
  expect_identical(
    vapply(result$tuning, function(x) x$status, ""),
    rep("success", 5L)
  )
  expect_identical(
    vapply(result$tuning, function(x) x$selection_rule, ""),
    rep("lambda.min", 5L)
  )
  expect_identical(
    vapply(result$tuning, function(x) x$inner_folds, integer(1)),
    rep(3L, 5L)
  )
  expect_identical(
    vapply(result$tuning, function(x) is.finite(x$lambda), logical(1)),
    rep(TRUE, 5L)
  )
  for (fold in seq_along(result$tuning)) {
    expect_setequal(
      names(result$tuning[[fold]]$fold_id),
      result$folds[[fold]]$analysis_ids
    )
  }
  expect_identical(
    result$provenance$glmnet_version,
    as.character(utils::packageVersion("glmnet"))
  )
})

test_that("Elastic Net settings are preserved in every outer-fold audit", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260817)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 1701L,
    glmnet_standardize = FALSE,
    glmnet_lambda_rule = "lambda.1se",
    glmnet_tuning_measure = "C"
  )

  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 101L)

  expect_identical(result$status, "success")
  expect_identical(result$pipeline$learner_params$alpha, 0.5)
  expect_identical(result$pipeline$learner_params$standardize, FALSE)
  expect_identical(
    result$pipeline$learner_params$selection_rule,
    "lambda.1se"
  )
  expect_identical(
    vapply(result$tuning, function(x) x$alpha, numeric(1)),
    rep(0.5, 5L)
  )
  expect_identical(
    vapply(result$coefficients, function(x) x$alpha, numeric(1)),
    rep(0.5, 5L)
  )
  expect_identical(
    vapply(result$tuning, function(x) x$standardize, logical(1)),
    rep(FALSE, 5L)
  )
  expect_identical(
    vapply(result$tuning, function(x) x$selection_rule, character(1)),
    rep("lambda.1se", 5L)
  )
  expect_identical(
    vapply(result$tuning, function(x) x$measure, character(1)),
    rep("C", 5L)
  )
})

test_that("Ridge alpha zero reaches tuning and final coefficient audits", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260818)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0,
    inner_folds = 3L,
    inner_seed = 1701L
  )

  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 101L)

  expect_identical(result$status, "success")
  expect_identical(
    vapply(result$tuning, function(x) x$alpha, numeric(1)),
    rep(0, 5L)
  )
  expect_identical(
    vapply(result$coefficients, function(x) x$alpha, numeric(1)),
    rep(0, 5L)
  )
})

test_that("nested glmnet Cox records outer-fold coefficient audits", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260813)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month"
  )
  pipeline <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 1701L
  )
  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 101L)

  expect_length(result$coefficients, 5L)
  for (fold in seq_along(result$coefficients)) {
    audit <- result$coefficients[[fold]]
    expect_identical(audit$fold, as.integer(fold))
    expect_identical(audit$learner, "glmnet_cox")
    expect_identical(audit$status, "success")
    expect_equal(audit$lambda, result$tuning[[fold]]$lambda, tolerance = 0)
    expect_identical(audit$values$feature, predictors)
    expect_identical(nrow(audit$values), length(predictors))
    expect_identical(
      audit$n_nonzero,
      as.integer(sum(audit$values$nonzero))
    )
    expect_identical(
      audit$nonzero_features,
      audit$values$feature[audit$values$nonzero]
    )
    assessment <- d[
      result$folds[[fold]]$assessment_rows,
      predictors,
      drop = FALSE
    ]
    reconstructed_risk <- unname(as.vector(
      as.matrix(assessment) %*% audit$values$coefficient
    ))
    observed_risk <- result$predictions$risk[
      match(
        result$folds[[fold]]$assessment_ids,
        result$predictions$patient_id
      )
    ]
    expect_equal(reconstructed_risk, observed_risk, tolerance = 1e-12)
  }
})

test_that("outer assessment patients cannot change fitted coefficients", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260813)
  predictors <- c("gene_a", "gene_b", "gene_c")
  pipeline <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 1701L
  )
  make_task <- function(data) {
    sb_surv_task(
      data, "patient_id", "time", "event", predictors, "month"
    )
  }
  original_task <- make_task(d)
  split <- sb_make_stratified_folds(
    original_task,
    folds = 5L,
    seed = 101L
  )$splits[[1]]
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  original <- sb_run_one_fold(
    original_task, pipeline, runtime, split, seed = 101L
  )

  changed <- d
  rows <- split$assessment_rows
  changed[rows, predictors] <- changed[rows, predictors] * 10000
  changed$time[rows] <- changed$time[rows] + 10000
  changed$event[rows] <- 1L - changed$event[rows]
  perturbed <- sb_run_one_fold(
    make_task(changed), pipeline, runtime, split, seed = 101L
  )

  expect_identical(original$coefficients, perturbed$coefficients)
})

test_that("a final glmnet fit failure cannot invent coefficients", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260813)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipeline <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 1701L
  )
  original_adapter <- sb_glmnet_cox_learner_adapter
  local_mocked_bindings(
    sb_glmnet_cox_learner_adapter = function(...) {
      adapter <- original_adapter(...)
      adapter$fit <- function(x, time, event) {
        sb_abort(
          "Synthetic final glmnet fit failure.",
          "stormbreaker_test_glmnet_final_fit_failure"
        )
      }
      adapter
    }
  )

  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 101L)

  expect_identical(
    vapply(result$tuning, function(x) x$status, ""),
    rep("success", 5L)
  )
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_test_glmnet_final_fit_failure"
  )
  expect_identical(
    vapply(result$coefficients, function(x) x$status, ""),
    rep("not_fitted", 5L)
  )
  expect_identical(all(is.na(result$predictions$risk)), TRUE)
})

test_that("each nested lambda ignores its outer assessment patients", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260813)
  predictors <- c("gene_a", "gene_b", "gene_c")
  pipeline <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 1701L
  )
  make_task <- function(data) {
    sb_surv_task(
      data, "patient_id", "time", "event", predictors, "month"
    )
  }
  original_task <- make_task(d)
  split <- sb_make_stratified_folds(
    original_task,
    folds = 5L,
    seed = 101L
  )$splits[[1]]
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  original <- sb_run_one_fold(
    original_task,
    pipeline,
    runtime,
    split,
    seed = 101L
  )

  changed <- d
  rows <- split$assessment_rows
  changed[rows, predictors] <- changed[rows, predictors] * 10000
  changed$time[rows] <- changed$time[rows] + 10000
  changed$event[rows] <- 1L - changed$event[rows]
  perturbed <- sb_run_one_fold(
    make_task(changed),
    pipeline,
    runtime,
    split,
    seed = 101L
  )

  expect_equal(
    original$tuning$lambda,
    perturbed$tuning$lambda,
    tolerance = 0
  )
  expect_identical(
    original$tuning$fold_id,
    perturbed$tuning$fold_id
  )
  expect_identical(
    original$tuning$curve,
    perturbed$tuning$curve
  )
})

test_that("nested tuning failures remain visible without fallback", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260813)
  d$constant_a <- 1
  d$constant_b <- 2
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("constant_a", "constant_b"), "month"
  )
  pipeline <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 1701L
  )
  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 101L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$status), "failed")
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_glmnet_tuning_failure"
  )
  expect_identical(all(is.na(result$predictions$risk)), TRUE)
  expect_identical(
    vapply(result$tuning, function(x) x$status, ""),
    rep("failed", 5L)
  )
  expect_identical(
    vapply(result$tuning, function(x) x$error_class, ""),
    rep("stormbreaker_glmnet_tuning_failure", 5L)
  )
  expect_identical(
    vapply(result$tuning, function(x) nzchar(x$error_message), logical(1)),
    rep(TRUE, 5L)
  )
  expect_identical(
    vapply(result$coefficients, function(x) x$status, ""),
    rep("not_fitted", 5L)
  )
  expect_identical(
    vapply(result$coefficients, function(x) nrow(x$values), integer(1)),
    rep(0L, 5L)
  )
})

test_that("failed tail-warning validation retains the original warning", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260816)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipeline <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 1701L
  )
  tail_warning <- paste0(
    "from glmnet C++ code (error code -27); Convergence for 27th ",
    "lambda value not reached after maxit=100000 iterations; ",
    "solutions for larger lambdas returned"
  )
  local_mocked_bindings(
    cv.glmnet = function(...) {
      warning(tail_warning, call. = FALSE)
      list(
        lambda = c(0.2, 0.1),
        cvm = c(2, 1),
        cvsd = c(0.1, 0.1),
        nzero = c(0L, 1L),
        lambda.min = 0.1
      )
    },
    glmnet = function(...) {
      structure(list(lambda = 0.5), class = "glmnet")
    },
    .package = "glmnet"
  )

  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 101L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_glmnet_selected_lambda_not_converged"
  )
  expect_identical(
    vapply(
      result$tuning,
      function(x) identical(x$warnings, tail_warning),
      logical(1)
    ),
    rep(TRUE, 3L)
  )
  expect_identical(
    vapply(result$tuning, function(x) is.null(x$warning_validation),
           logical(1)),
    rep(TRUE, 3L)
  )
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(all(is.na(result$predictions$risk)), TRUE)
})

test_that("survival CV DCA is reconstructed from fixed OOF probabilities", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  thresholds <- c(0.20, 0.40, 0.60)
  plain_task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month",
    prediction_horizon = 7
  )
  threshold_task <- sb_surv_task(
    d, "patient_id", "time", "event", predictors, "month",
    prediction_horizon = 7,
    decision_thresholds = thresholds
  )
  pipeline <- sb_pipeline()
  plain <- sb_run_cv(plain_task, pipeline, folds = 5L, seed = 20260819L)
  observed <- sb_run_cv(
    threshold_task,
    pipeline,
    folds = 5L,
    seed = 20260819L
  )
  metric <- observed$metrics$decision_curve

  expect_identical(observed$status, "success")
  expect_identical(metric$status, "success")
  expect_identical(metric$scope, "pooled_patient_level_out_of_fold")
  expect_identical(metric$curve$threshold_probability, thresholds)
  expect_identical(observed$task$decision_thresholds, thresholds)
  expect_identical(observed$predictions, plain$predictions)
  expect_identical(observed$folds, plain$folds)
  expect_false("decision_curve" %in% names(plain$metrics))
  expect_false(any(c(
    "decision", "action", "treat", "selected"
  ) %in% names(observed$predictions)))

  threshold <- thresholds[[1L]]
  selected <- 1 - observed$predictions$survival_probability >= threshold
  manual_survival <- function(time, event, horizon) {
    event_times <- sort(unique(time[event == 1L & time <= horizon]))
    if (length(event_times) == 0L) {
      return(1)
    }
    prod(vapply(event_times, function(value) {
      1 - sum(time == value & event == 1L) / sum(time >= value)
    }, numeric(1)))
  }
  selected_survival <- manual_survival(
    observed$predictions$time[selected],
    observed$predictions$event[selected],
    7
  )
  overall_survival <- manual_survival(
    observed$predictions$time,
    observed$predictions$event,
    7
  )
  selected_fraction <- sum(selected) / nrow(observed$predictions)
  odds <- threshold / (1 - threshold)
  expected_model <- selected_fraction * (1 - selected_survival) -
    selected_fraction * selected_survival * odds
  expected_all <- (1 - overall_survival) - overall_survival * odds
  row <- metric$curve[1L, ]

  expect_equal(row$selected_survival_probability, selected_survival)
  expect_equal(row$model_net_benefit, expected_model)
  expect_equal(row$treat_all_net_benefit, expected_all)
  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "OOF survival decision curve: 3 declared thresholds")
})

test_that("survival CV reuses an exact patient-keyed outer assignment", {
  d <- make_sb_fixture(n = 60L, seed = 6101L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 6
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "coxph",
    selector_n_features = 2L
  )
  generated <- sb_run_cv(task, pipeline, folds = 5L, seed = 6102L)
  assignment <- setNames(
    generated$predictions$fold,
    generated$predictions$patient_id
  )
  supplied <- assignment[rev(seq_along(assignment))]

  reused <- sb_run_cv(task, pipeline, folds = supplied, seed = 999L)
  task_ids <- as.character(task$data[[task$id_col]])
  expected_assignment <- data.frame(
    patient_id = task_ids,
    fold = as.integer(supplied[match(task_ids, names(supplied))]),
    stringsAsFactors = FALSE
  )

  expect_identical(reused$folds, generated$folds)
  expect_identical(reused$predictions, generated$predictions)
  expect_identical(reused$metrics, generated$metrics)
  expect_identical(reused$coefficients, generated$coefficients)
  expect_identical(reused$baseline_hazards, generated$baseline_hazards)
  expect_identical(reused$selections, generated$selections)
  expect_identical(
    reused$selection_stability,
    generated$selection_stability
  )
  expect_identical(
    reused$runs[setdiff(names(reused$runs), "seed")],
    generated$runs[setdiff(names(generated$runs), "seed")]
  )
  expect_identical(generated$runs$seed, rep(6102L, 5L))
  expect_identical(reused$runs$seed, rep(NA_integer_, 5L))
  expect_identical(generated$provenance$fold_source, "generated_stratified")
  expect_identical(reused$provenance$fold_source, "preassigned_patient_id")
  expect_identical(generated$provenance$seed, 6102L)
  expect_identical(reused$provenance$seed, NA_integer_)
  expect_identical(generated$provenance$folds, 5L)
  expect_identical(reused$provenance$folds, 5L)
  expect_identical(reused$provenance$fold_assignment, expected_assignment)
  expect_identical(
    generated$provenance$fold_assignment,
    expected_assignment
  )
  printed <- paste(capture.output(print(reused)), collapse = "\n")
  expect_match(printed, "preassigned patient-ID folds", fixed = TRUE)
})

test_that("preassigned survival folds keep glmnet inner tuning in analysis", {
  skip_if_not_installed("glmnet")
  d <- make_sb_fixture(n = 60L, seed = 6111L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 6112L
  )
  generated_folds <- sb_make_stratified_folds(task, 3L, 6113L)
  ids <- as.character(task$data[[task$id_col]])
  supplied <- setNames(generated_folds$assignment, ids)

  result <- sb_run_cv(task, pipeline, folds = supplied, seed = 999L)

  expect_identical(result$status, "success")
  expect_identical(result$provenance$fold_source, "preassigned_patient_id")
  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    inner_ids <- names(result$tuning[[fold]]$fold_id)
    expect_setequal(inner_ids, split$analysis_ids)
    expect_length(intersect(inner_ids, split$assessment_ids), 0L)
  }
})
