make_sb_binary_external_comparison <- function(
  prefix,
  reverse = FALSE,
  second_signature_id = "published_b",
  fail_second = FALSE
) {
  index <- seq_len(24L)
  data <- data.frame(
    patient_id = paste0(prefix, "_", index),
    outcome = rep(c(0L, 1L), 12L),
    gene_a = seq(-1, 1, length.out = 24L),
    gene_b = cos(index),
    stringsAsFactors = FALSE
  )
  signatures <- list(
    sb_linear_signature(
      "published_a", "binary", c(gene_a = 0.8, gene_b = -0.2),
      intercept = -0.1
    ),
    sb_linear_signature(
      second_signature_id, "binary", c(gene_a = -0.3, gene_b = 0.6),
      intercept = 0.2
    )
  )
  if (reverse) {
    signatures <- rev(signatures)
  }
  predictions <- lapply(signatures, function(signature) {
    sb_predict_signature(
      signature,
      data[c("patient_id", "gene_a", "gene_b")],
      "patient_id"
    )
  })
  if (fail_second) {
    predictions[[2L]]$predictions$patient_id[[1L]] <- "wrong_patient"
  }
  sb_compare_external(
    predictions,
    sb_binary_endpoint(data, "patient_id", "outcome")
  )
}

make_sb_survival_external_comparison <- function(
  prefix,
  time_unit,
  horizon = NULL,
  endpoint_horizon = TRUE,
  fail_prediction = FALSE
) {
  index <- seq_len(24L)
  data <- data.frame(
    patient_id = paste0(prefix, "_", index),
    time = index,
    event = rep(c(1L, 0L, 1L), 8L),
    gene_a = seq(-1, 1, length.out = 24L),
    stringsAsFactors = FALSE
  )
  signature_args <- list(
    signature_id = "published_survival",
    task_type = "survival",
    coefficients = c(gene_a = 0.4),
    intercept = 0
  )
  if (!is.null(horizon)) {
    signature_args$prediction_horizon <- horizon
    signature_args$time_unit <- time_unit
    signature_args$baseline_survival <- 0.8
  }
  signature <- do.call(sb_linear_signature, signature_args)
  predictor_data <- if (fail_prediction) {
    data["patient_id"]
  } else {
    data[c("patient_id", "gene_a")]
  }
  prediction <- sb_predict_signature(
    signature,
    predictor_data,
    "patient_id"
  )
  endpoint_args <- list(
    data = data,
    id_col = "patient_id",
    time_col = "time",
    event_col = "event",
    time_unit = time_unit
  )
  if (!is.null(horizon) && endpoint_horizon) {
    endpoint_args$prediction_horizon <- horizon
  }
  sb_compare_external(
    list(prediction),
    do.call(
      sb_surv_endpoint,
      endpoint_args
    )
  )
}

test_that("survival collections retain the common evaluation target", {
  horizon_collection <- sb_collect_external_comparisons(list(
    cohort_a = make_sb_survival_external_comparison("horizon_a", "month", 6),
    cohort_b = make_sb_survival_external_comparison("horizon_b", "month", 6)
  ))
  risk_collection <- sb_collect_external_comparisons(list(
    cohort_a = make_sb_survival_external_comparison("risk_a", "month"),
    cohort_b = make_sb_survival_external_comparison("risk_b", "month")
  ))
  inherited_collection <- sb_collect_external_comparisons(list(
    cohort_a = make_sb_survival_external_comparison(
      "inherited_a", "month", 6, endpoint_horizon = FALSE
    ),
    cohort_b = make_sb_survival_external_comparison(
      "inherited_b", "month", 6, endpoint_horizon = FALSE
    )
  ))
  failed_inherited <- make_sb_survival_external_comparison(
    "failed_inherited",
    "month",
    6,
    endpoint_horizon = FALSE,
    fail_prediction = TRUE
  )
  mixed_status_collection <- sb_collect_external_comparisons(list(
    successful = inherited_collection$comparisons$cohort_a,
    failed = failed_inherited
  ))

  expect_identical(horizon_collection$status, "success")
  expect_identical(horizon_collection$task_type, "survival")
  expect_identical(horizon_collection$cohorts$time_unit, c("month", "month"))
  expect_identical(horizon_collection$cohorts$prediction_horizon, c(6, 6))
  expect_identical(
    horizon_collection$results$cohort_id,
    c("cohort_a", "cohort_b")
  )
  expect_identical(risk_collection$status, "success")
  expect_identical(
    "prediction_horizon" %in% names(risk_collection$cohorts),
    FALSE
  )
  expect_identical(
    inherited_collection$cohorts$prediction_horizon,
    c(6, 6)
  )
  expect_identical(
    c("ipcw_brier", "time_auc", "uno_c") %in%
      names(inherited_collection$results),
    rep(TRUE, 3L)
  )
  expect_identical(mixed_status_collection$status, "incomplete")
  expect_identical(
    mixed_status_collection$cohorts$prediction_horizon,
    c(6, 6)
  )
})

test_that("external comparisons remain cohort-specific and ordered", {
  cohort_a <- make_sb_binary_external_comparison("a")
  cohort_b <- make_sb_binary_external_comparison("b", reverse = TRUE)

  observed <- sb_collect_external_comparisons(list(
    cohort_a = cohort_a,
    cohort_b = cohort_b
  ))

  expect_s3_class(observed, "stormbreaker_external_collection")
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expected_ids <- c("signature__published_a", "signature__published_b")
  expect_identical(observed$pipeline_ids, expected_ids)
  expect_identical(observed$cohorts$cohort_id, c("cohort_a", "cohort_b"))
  expect_identical(observed$cohorts$n_patients, c(24L, 24L))
  expect_identical(
    observed$results$cohort_id,
    rep(c("cohort_a", "cohort_b"), each = 2L)
  )
  expect_identical(
    observed$results$pipeline_id,
    rep(expected_ids, 2L)
  )
  expect_equal(
    observed$results$roc_auc[1:2],
    cohort_a$comparison$roc_auc,
    tolerance = 0
  )
  expect_equal(
    observed$results$roc_auc[3:4],
    rev(cohort_b$comparison$roc_auc),
    tolerance = 0
  )
  expect_identical(observed$comparisons, list(
    cohort_a = cohort_a,
    cohort_b = cohort_b
  ))

  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "No metric is pooled", fixed = TRUE)
  expect_match(printed, "no pipeline is ranked or selected", fixed = TRUE)
})

test_that("incomplete cohort comparisons remain visible", {
  cohort_a <- make_sb_binary_external_comparison("complete")
  cohort_b <- make_sb_binary_external_comparison(
    "incomplete",
    fail_second = TRUE
  )

  observed <- sb_collect_external_comparisons(list(
    complete = cohort_a,
    incomplete = cohort_b
  ))

  expect_identical(observed$status, "incomplete")
  expect_identical(observed$cohorts$status, c("success", "incomplete"))
  expect_identical(
    observed$results$status,
    c("success", "success", "success", "failed")
  )
  expect_identical(
    observed$results$error_class[[4L]],
    "stormbreaker_external_patient_mismatch"
  )
})

test_that("external collection rejects ambiguous or incompatible inputs", {
  binary <- make_sb_binary_external_comparison("binary")
  other_pipelines <- make_sb_binary_external_comparison(
    "other",
    second_signature_id = "different_signature"
  )
  month_6 <- make_sb_survival_external_comparison("month_6", "month", 6)
  month_12 <- make_sb_survival_external_comparison("month_12", "month", 12)
  year_6 <- make_sb_survival_external_comparison("year_6", "year", 6)
  risk_only <- make_sb_survival_external_comparison("risk_only", "month")
  inherited_6 <- make_sb_survival_external_comparison(
    "inherited_6", "month", 6, endpoint_horizon = FALSE
  )
  inherited_12 <- make_sb_survival_external_comparison(
    "inherited_12", "month", 12, endpoint_horizon = FALSE
  )
  failed_inherited_6 <- make_sb_survival_external_comparison(
    "failed_inherited_6",
    "month",
    6,
    endpoint_horizon = FALSE,
    fail_prediction = TRUE
  )
  failed_inherited_12 <- make_sb_survival_external_comparison(
    "failed_inherited_12",
    "month",
    12,
    endpoint_horizon = FALSE,
    fail_prediction = TRUE
  )
  unknown_target <- failed_inherited_6
  unknown_target$prediction_target <- list(
    target_status = "unknown",
    has_horizon = NA,
    prediction_horizon = NULL,
    time_unit = "month",
    horizon_source = "unknown"
  )

  empty <- capture_sb_error(sb_collect_external_comparisons(list()))
  unnamed <- capture_sb_error(sb_collect_external_comparisons(list(binary)))
  direct <- capture_sb_error(sb_collect_external_comparisons(binary))
  duplicate_names <- capture_sb_error(sb_collect_external_comparisons(
    structure(list(binary, binary), names = c("same", "same"))
  ))
  task_mismatch <- capture_sb_error(sb_collect_external_comparisons(list(
    cohort_a = binary,
    cohort_b = month_6
  )))
  pipeline_mismatch <- capture_sb_error(sb_collect_external_comparisons(list(
    cohort_a = binary,
    cohort_b = other_pipelines
  )))
  horizon_mismatch <- capture_sb_error(sb_collect_external_comparisons(list(
    cohort_a = month_6,
    cohort_b = month_12
  )))
  horizon_presence_mismatch <- capture_sb_error(
    sb_collect_external_comparisons(list(
      cohort_a = month_6,
      cohort_b = risk_only
    ))
  )
  unit_mismatch <- capture_sb_error(sb_collect_external_comparisons(list(
    cohort_a = month_6,
    cohort_b = year_6
  )))
  inherited_horizon_mismatch <- capture_sb_error(
    sb_collect_external_comparisons(list(
      cohort_a = inherited_6,
      cohort_b = inherited_12
    ))
  )
  failed_horizon_mismatch <- capture_sb_error(
    sb_collect_external_comparisons(list(
      cohort_a = failed_inherited_6,
      cohort_b = failed_inherited_12
    ))
  )
  unknown_target_error <- capture_sb_error(
    sb_collect_external_comparisons(list(cohort_a = unknown_target))
  )

  expect_s3_class(empty, "stormbreaker_bad_external_collection")
  expect_s3_class(unnamed, "stormbreaker_bad_external_collection")
  expect_s3_class(direct, "stormbreaker_bad_external_collection")
  expect_s3_class(duplicate_names, "stormbreaker_bad_external_collection")
  expect_s3_class(task_mismatch, "stormbreaker_external_task_mismatch")
  expect_s3_class(
    pipeline_mismatch,
    "stormbreaker_external_pipeline_mismatch"
  )
  expect_s3_class(
    horizon_mismatch,
    "stormbreaker_external_target_mismatch"
  )
  expect_s3_class(
    horizon_presence_mismatch,
    "stormbreaker_external_target_mismatch"
  )
  expect_s3_class(unit_mismatch, "stormbreaker_external_target_mismatch")
  expect_s3_class(
    inherited_horizon_mismatch,
    "stormbreaker_external_target_mismatch"
  )
  expect_s3_class(
    failed_horizon_mismatch,
    "stormbreaker_external_target_mismatch"
  )
  expect_s3_class(
    unknown_target_error,
    "stormbreaker_external_target_mismatch"
  )
})
