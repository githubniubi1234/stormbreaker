make_sb_paired_binary_prediction <- function(patient_ids, probability,
                                             pipeline_id) {
  n <- length(patient_ids)
  result <- list(
    status = "success",
    task_type = "binary",
    pipeline_id = pipeline_id,
    learner = "logistic",
    id_col = "patient_id",
    predictions = data.frame(
      patient_id = patient_ids,
      pipeline_id = rep(pipeline_id, n),
      learner = rep("logistic", n),
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

make_sb_paired_survival_prediction <- function(
    patient_ids,
    risk,
    survival_probability,
    pipeline_id) {
  n <- length(patient_ids)
  result <- list(
    status = "success",
    task_type = "survival",
    pipeline_id = pipeline_id,
    learner = "coxph",
    id_col = "patient_id",
    predictions = data.frame(
      patient_id = patient_ids,
      pipeline_id = rep(pipeline_id, n),
      learner = rep("coxph", n),
      risk = risk,
      prediction_horizon = rep(30, n),
      time_unit = rep("month", n),
      survival_probability = survival_probability,
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

test_that("binary model differences use the same patient bootstrap draws", {
  ids <- sprintf("p%02d", seq_len(40L))
  outcome <- rep(c(0L, 1L), 20L)
  reference_probability <- seq(0.08, 0.92, length.out = length(ids))
  candidate_probability <- pmin(
    0.97,
    pmax(0.03, reference_probability + rep(c(-0.04, 0.06), 20L))
  )
  predictions <- list(
    make_sb_paired_binary_prediction(
      ids, reference_probability, "reference_binary"
    ),
    make_sb_paired_binary_prediction(
      ids, candidate_probability, "candidate_binary"
    ),
    make_sb_paired_binary_prediction(
      ids, rev(candidate_probability), "candidate_binary_2"
    )
  )
  endpoint <- sb_binary_endpoint(
    data.frame(patient_id = ids, outcome = outcome),
    "patient_id",
    "outcome"
  )
  comparison <- sb_compare_external(predictions, endpoint)

  set.seed(9400L)
  caller_state <- .Random.seed
  observed <- sb_bootstrap_external_comparison(
    comparison,
    reference_pipeline_id = "reference_binary",
    reps = 30L,
    confidence_level = 0.90,
    seed = 9401L
  )
  expect_identical(.Random.seed, caller_state)

  reference <- sb_bootstrap_external(
    comparison$evaluations$reference_binary,
    reps = 30L,
    confidence_level = 0.90,
    seed = 9401L
  )
  candidate <- sb_bootstrap_external(
    comparison$evaluations$candidate_binary,
    reps = 30L,
    confidence_level = 0.90,
    seed = 9401L
  )
  reference_rows <- reference$replicates[
    reference$replicates$metric %in% observed$metrics,
    ,
    drop = FALSE
  ]
  candidate_rows <- candidate$replicates[
    candidate$replicates$metric %in% observed$metrics,
    ,
    drop = FALSE
  ]
  expected_key <- paste(reference_rows$replicate, reference_rows$metric)
  candidate_rows <- candidate_rows[match(
    expected_key,
    paste(candidate_rows$replicate, candidate_rows$metric)
  ), , drop = FALSE]
  difference_rows <- observed$differences[
    observed$differences$candidate_pipeline_id == "candidate_binary",
    ,
    drop = FALSE
  ]
  difference_rows <- difference_rows[match(
    expected_key,
    paste(difference_rows$replicate, difference_rows$metric)
  ), , drop = FALSE]

  expect_s3_class(
    observed,
    "stormbreaker_external_bootstrap_comparison"
  )
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(observed$reference_pipeline_id, "reference_binary")
  expect_identical(observed$prediction_target, list(
    target_status = "complete",
    has_horizon = FALSE,
    prediction_horizon = NULL,
    time_unit = NULL,
    horizon_source = "not_applicable"
  ))
  expect_identical(
    observed$metrics,
    c("roc_auc", "average_precision", "brier_score", "log_loss")
  )
  expect_identical(observed$resample_audit, reference$resample_audit)
  expect_equal(
    difference_rows$reference_estimate,
    reference_rows$estimate,
    tolerance = 0
  )
  expect_equal(
    difference_rows$candidate_estimate,
    candidate_rows$estimate,
    tolerance = 0
  )
  expect_equal(
    difference_rows$estimate,
    candidate_rows$estimate - reference_rows$estimate,
    tolerance = 0
  )
  brier_differences <- difference_rows$estimate[
    difference_rows$metric == "brier_score"
  ]
  brier_interval <- observed$intervals[
    observed$intervals$candidate_pipeline_id == "candidate_binary" &
      observed$intervals$metric == "brier_score",
    ,
    drop = FALSE
  ]
  expect_equal(
    unname(c(brier_interval$lower, brier_interval$upper)),
    unname(stats::quantile(
      brier_differences,
      probs = c(0.05, 0.95),
      names = FALSE,
      type = 7
    )),
    tolerance = 0
  )
  expect_identical(
    brier_interval$metric_direction,
    "lower_is_better"
  )
  second_candidate <- observed$differences[
    observed$differences$candidate_pipeline_id == "candidate_binary_2",
    ,
    drop = FALSE
  ]
  second_candidate <- second_candidate[match(
    expected_key,
    paste(second_candidate$replicate, second_candidate$metric)
  ), , drop = FALSE]
  expect_equal(
    second_candidate$reference_estimate,
    difference_rows$reference_estimate,
    tolerance = 0
  )
  expect_identical(
    observed$candidate_status$status,
    c("success", "success")
  )
  expect_identical(
    any(grepl("calibration|decision|rank|winner", observed$metrics)),
    FALSE
  )
})

test_that("survival horizon differences preserve paired censoring draws", {
  ids <- sprintf("s%02d", seq_len(60L))
  time <- as.numeric(seq_len(60L))
  event <- rep.int(1L, length(ids))
  reference_risk <- seq(-1.2, 1.2, length.out = length(ids))
  candidate_risk <- reference_risk + sin(seq_along(ids)) / 8
  predictions <- list(
    make_sb_paired_survival_prediction(
      ids,
      reference_risk,
      stats::plogis(-reference_risk),
      "reference_survival"
    ),
    make_sb_paired_survival_prediction(
      ids,
      candidate_risk,
      stats::plogis(-candidate_risk),
      "candidate_survival"
    )
  )
  endpoint <- sb_surv_endpoint(
    data.frame(patient_id = ids, time = time, event = event),
    "patient_id",
    "time",
    "event",
    "month",
    prediction_horizon = 30
  )
  comparison <- sb_compare_external(predictions, endpoint)

  observed <- sb_bootstrap_external_comparison(
    comparison,
    reference_pipeline_id = "reference_survival",
    reps = 20L,
    confidence_level = 0.90,
    seed = 9501L
  )
  reference <- sb_bootstrap_external(
    comparison$evaluations$reference_survival,
    reps = 20L,
    confidence_level = 0.90,
    seed = 9501L
  )
  candidate <- sb_bootstrap_external(
    comparison$evaluations$candidate_survival,
    reps = 20L,
    confidence_level = 0.90,
    seed = 9501L
  )
  reference_rows <- reference$replicates[
    reference$replicates$metric %in% observed$metrics,
    ,
    drop = FALSE
  ]
  candidate_rows <- candidate$replicates[
    candidate$replicates$metric %in% observed$metrics,
    ,
    drop = FALSE
  ]
  key <- paste(reference_rows$replicate, reference_rows$metric)
  candidate_rows <- candidate_rows[match(
    key,
    paste(candidate_rows$replicate, candidate_rows$metric)
  ), , drop = FALSE]
  difference_rows <- observed$differences[match(
    key,
    paste(observed$differences$replicate, observed$differences$metric)
  ), , drop = FALSE]

  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "survival")
  expect_identical(observed$prediction_target, list(
    target_status = "complete",
    has_horizon = TRUE,
    prediction_horizon = 30,
    time_unit = "month",
    horizon_source = "endpoint"
  ))
  expect_identical(
    observed$metrics,
    c("harrell_c", "ipcw_brier_score", "time_auc", "uno_c")
  )
  expect_identical(observed$resample_audit, reference$resample_audit)
  expect_equal(
    difference_rows$estimate,
    candidate_rows$estimate - reference_rows$estimate,
    tolerance = 0
  )
  expect_identical(
    observed$intervals$metric_direction,
    c(
      "higher_is_better", "lower_is_better",
      "higher_is_better", "higher_is_better"
    )
  )
})

test_that("a failed candidate cannot erase valid paired comparisons", {
  ids <- sprintf("f%02d", seq_len(40L))
  outcome <- rep(c(0L, 1L), 20L)
  base_probability <- seq(0.10, 0.90, length.out = length(ids))
  predictions <- list(
    make_sb_paired_binary_prediction(ids, base_probability, "reference"),
    make_sb_paired_binary_prediction(
      ids,
      rev(base_probability),
      "valid_candidate"
    ),
    make_sb_paired_binary_prediction(
      ids,
      pmin(base_probability + 0.02, 0.99),
      "failed_candidate"
    )
  )
  predictions[[3L]]$predictions$patient_id[[1L]] <- "wrong_patient"
  endpoint <- sb_binary_endpoint(
    data.frame(patient_id = ids, outcome = outcome),
    "patient_id",
    "outcome"
  )
  comparison <- sb_compare_external(predictions, endpoint)

  observed <- sb_bootstrap_external_comparison(
    comparison,
    reference_pipeline_id = "reference",
    reps = 12L,
    confidence_level = 0.90,
    seed = 9601L
  )

  expect_identical(comparison$status, "incomplete")
  expect_identical(observed$status, "incomplete")
  expect_identical(
    observed$candidate_status$candidate_pipeline_id,
    c("valid_candidate", "failed_candidate")
  )
  expect_identical(
    observed$candidate_status$status,
    c("success", "failed")
  )
  expect_identical(
    observed$candidate_status$error_class,
    c(NA_character_, "stormbreaker_failed_external_evaluation")
  )
  expect_identical(
    unique(observed$intervals$candidate_pipeline_id),
    "valid_candidate"
  )
  expect_identical(
    unique(observed$differences$candidate_pipeline_id),
    "valid_candidate"
  )
  expect_identical(
    "failed_candidate" %in% observed$model_replicates$pipeline_id,
    FALSE
  )
  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "No ranking, winner, or p value", fixed = TRUE)
})
