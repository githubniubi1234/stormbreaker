test_that("a locked binary signature computes its published logistic formula", {
  signature <- sb_linear_signature(
    signature_id = "published_binary_demo",
    task_type = "binary",
    coefficients = c(gene_b = -0.4, gene_a = 0.7),
    intercept = -0.2,
    source = "Example publication formula",
    input_description = "Preprocessed values on the publication scale"
  )
  new_data <- data.frame(
    patient_id = paste0("p", seq_len(8L)),
    gene_a = seq(-1, 1, length.out = 8L),
    unrelated_outcome = rep(c(0L, 1L), 4L),
    gene_b = seq(0.5, -0.5, length.out = 8L),
    stringsAsFactors = FALSE
  )

  observed <- sb_predict_signature(signature, new_data, "patient_id")
  changed_outcome <- new_data
  changed_outcome$unrelated_outcome <- 1L - changed_outcome$unrelated_outcome
  repeated <- sb_predict_signature(
    signature,
    changed_outcome[c("gene_b", "patient_id", "unrelated_outcome", "gene_a")],
    "patient_id"
  )
  expected_score <- -0.2 - 0.4 * new_data$gene_b + 0.7 * new_data$gene_a

  expect_s3_class(signature, "stormbreaker_linear_signature")
  expect_identical(signature$features, c("gene_b", "gene_a"))
  expect_identical(signature$source, "Example publication formula")
  expect_s3_class(observed, "stormbreaker_final_prediction")
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(observed$pipeline_id, "signature__published_binary_demo")
  expect_identical(observed$learner, "locked_linear_signature")
  expect_equal(
    observed$predictions$probability,
    stats::plogis(expected_score),
    tolerance = 0
  )
  expect_identical(
    repeated$predictions$probability,
    observed$predictions$probability
  )
  expect_identical(
    observed$signature$coefficients,
    c(gene_b = -0.4, gene_a = 0.7)
  )

  endpoint <- sb_binary_endpoint(
    data.frame(
      patient_id = new_data$patient_id,
      outcome = rep(c(0L, 1L), 4L)
    ),
    "patient_id",
    "outcome"
  )
  evaluation <- sb_evaluate_external(observed, endpoint)
  expect_identical(evaluation$status, "success")
  expect_identical(
    evaluation$scored_predictions$probability,
    observed$predictions$probability
  )
})

test_that("binary probabilities require a stated intercept", {
  omitted <- capture_sb_error(sb_linear_signature(
    "omitted_binary_intercept",
    "binary",
    c(gene_a = 1)
  ))
  explicit_zero <- sb_linear_signature(
    "explicit_zero_binary_intercept",
    "binary",
    c(gene_a = 1),
    intercept = 0
  )
  prediction <- sb_predict_signature(
    explicit_zero,
    data.frame(patient_id = "p1", gene_a = 0),
    "patient_id"
  )
  survival_default <- sb_linear_signature(
    "survival_zero_reference",
    "survival",
    c(gene_a = 1)
  )

  expect_s3_class(omitted, "stormbreaker_missing_signature_intercept")
  expect_identical(prediction$predictions$probability, 0.5)
  expect_identical(survival_default$intercept, 0)
})

test_that("a survival signature can remain risk-only or use a stated baseline", {
  coefficients <- c(gene_a = 0.5, gene_b = -0.25)
  risk_signature <- sb_linear_signature(
    "published_survival_risk",
    "survival",
    coefficients,
    intercept = 0.1
  )
  horizon_signature <- sb_linear_signature(
    "published_survival_horizon",
    "survival",
    coefficients,
    intercept = 0.1,
    prediction_horizon = 6,
    time_unit = "month",
    baseline_survival = 0.82
  )
  new_data <- data.frame(
    patient_id = paste0("s", seq_len(16L)),
    gene_b = seq(0.8, -0.8, length.out = 16L),
    event = rep(c(1L, 0L), 8L),
    gene_a = seq(-1.2, 1.2, length.out = 16L),
    stringsAsFactors = FALSE
  )
  expected_risk <- 0.1 + 0.5 * new_data$gene_a - 0.25 * new_data$gene_b

  risk_prediction <- sb_predict_signature(
    risk_signature,
    new_data,
    "patient_id"
  )
  horizon_prediction <- sb_predict_signature(
    horizon_signature,
    new_data,
    "patient_id"
  )

  expect_identical(risk_prediction$status, "success")
  expect_equal(risk_prediction$predictions$risk, expected_risk, tolerance = 0)
  expect_identical(
    "survival_probability" %in% names(risk_prediction$predictions),
    FALSE
  )
  expect_equal(
    horizon_prediction$predictions$risk,
    expected_risk,
    tolerance = 0
  )
  expect_equal(
    horizon_prediction$predictions$survival_probability,
    0.82^exp(expected_risk),
    tolerance = 0
  )
  expect_identical(
    unique(horizon_prediction$predictions$prediction_horizon),
    6
  )
  expect_identical(
    unique(horizon_prediction$predictions$time_unit),
    "month"
  )

  endpoint_data <- data.frame(
    patient_id = new_data$patient_id,
    time = as.numeric(seq_len(16L)),
    event = rep(c(1L, 0L), 8L)
  )
  evaluation <- sb_evaluate_external(
    horizon_prediction,
    sb_surv_endpoint(
      endpoint_data,
      "patient_id",
      "time",
      "event",
      "month",
      prediction_horizon = 6
    )
  )
  expect_identical(evaluation$status, "success")
  expect_identical(evaluation$metrics$ipcw_brier$status, "success")
  expect_identical(evaluation$metrics$time_auc$status, "success")
  expect_identical(evaluation$metrics$uno_c$status, "success")

  zero_baseline <- sb_linear_signature(
    "zero_baseline",
    "survival",
    c(gene_a = 1),
    prediction_horizon = 6,
    time_unit = "month",
    baseline_survival = 0
  )
  zero_prediction <- sb_predict_signature(
    zero_baseline,
    data.frame(patient_id = "extreme", gene_a = -1000),
    "patient_id"
  )
  expect_identical(
    zero_prediction$predictions$survival_probability,
    0
  )
})

test_that("signature prediction failures stay explicit without substitution", {
  signature <- sb_linear_signature(
    "missing_feature_demo",
    "binary",
    c(gene_a = 0.5, gene_b = -0.2),
    intercept = 0
  )
  missing <- data.frame(
    patient_id = c("p1", "p2"),
    gene_a = c(0.1, 0.2)
  )
  nonnumeric <- data.frame(
    patient_id = c("p1", "p2"),
    gene_a = c(0.1, 0.2),
    gene_b = c("low", "high")
  )

  missing_result <- sb_predict_signature(signature, missing, "patient_id")
  nonnumeric_result <- sb_predict_signature(
    signature,
    nonnumeric,
    "patient_id"
  )
  partial_horizon <- capture_sb_error(sb_linear_signature(
    "partial_horizon",
    "survival",
    c(gene_a = 0.5),
    prediction_horizon = 12,
    time_unit = "month"
  ))
  unnamed <- capture_sb_error(sb_linear_signature(
    "unnamed",
    "binary",
    coefficients = c(0.5, -0.2)
  ))

  expect_identical(missing_result$status, "failed")
  expect_identical(missing_result$error_stage, "features")
  expect_identical(
    missing_result$error_class,
    "stormbreaker_missing_signature_features"
  )
  expect_identical(
    missing_result$predictions$status,
    rep("failed", 2L)
  )
  expect_identical(missing_result$learner, "locked_linear_signature")
  expect_identical(nonnumeric_result$status, "failed")
  expect_identical(nonnumeric_result$error_stage, "features")
  expect_s3_class(partial_horizon, "stormbreaker_bad_signature_horizon")
  expect_s3_class(unnamed, "stormbreaker_bad_signature_coefficients")
})

test_that("locked signatures reuse external comparison and paired bootstrap", {
  ids <- sprintf("e%02d", seq_len(40L))
  data <- data.frame(
    patient_id = ids,
    outcome = rep(c(0L, 1L), 20L),
    gene_a = seq(-1, 1, length.out = 40L),
    gene_b = cos(seq_len(40L) / 5),
    stringsAsFactors = FALSE
  )
  signatures <- list(
    sb_linear_signature(
      "published_reference",
      "binary",
      c(gene_a = 0.6, gene_b = -0.2),
      intercept = -0.1
    ),
    sb_linear_signature(
      "published_candidate",
      "binary",
      c(gene_a = 0.4, gene_b = 0.3),
      intercept = 0.2
    )
  )
  predictions <- lapply(signatures, function(signature) {
    sb_predict_signature(
      signature,
      data[c("patient_id", "gene_a", "gene_b")],
      "patient_id"
    )
  })
  comparison <- sb_compare_external(
    predictions,
    sb_binary_endpoint(data, "patient_id", "outcome")
  )
  paired <- sb_bootstrap_external_comparison(
    comparison,
    reference_pipeline_id = predictions[[1L]]$pipeline_id,
    reps = 12L,
    confidence_level = 0.90,
    seed = 9901L
  )

  expect_identical(comparison$status, "success")
  expect_identical(
    comparison$comparison$learner,
    rep("locked_linear_signature", 2L)
  )
  expect_identical(paired$status, "success")
  expect_identical(
    unique(paired$intervals$candidate_pipeline_id),
    predictions[[2L]]$pipeline_id
  )
  expect_identical(
    paired$pairing$draw_sharing,
    "same_patient_indices_for_all_models"
  )
})

test_that("signature predictions print their formula-only origin", {
  signature <- sb_linear_signature(
    "printed_signature",
    "binary",
    c(gene_a = 0.5),
    intercept = 0
  )
  prediction <- sb_predict_signature(
    signature,
    data.frame(patient_id = "patient_1", gene_a = 0.2),
    "patient_id"
  )

  printed <- paste(capture.output(print(prediction)), collapse = "\n")

  expect_match(printed, "locked formula-only prediction", fixed = TRUE)
  expect_match(printed, "printed_signature", fixed = TRUE)
  expect_identical(
    grepl("final state-only prediction", printed, fixed = TRUE),
    FALSE
  )
})
