test_that("binary locked predictions share one external endpoint", {
  development <- make_sb_binary_fixture(n = 60L, seed = 8101L)
  development_task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  )
  fits <- lapply(pipelines, function(pipeline) {
    sb_fit_final(development_task, pipeline)
  })

  external <- make_sb_binary_fixture(n = 36L, seed = 8102L)
  external$patient_id <- paste0("external_", external$patient_id)
  predictors <- external[c(
    "patient_id", "gene_a", "gene_b", "gene_c"
  )]
  predictions <- lapply(fits, function(fit) {
    sb_predict_final(fit, predictors, "patient_id")
  })
  endpoint <- sb_binary_endpoint(
    external,
    "patient_id",
    "outcome"
  )

  observed <- sb_compare_external(predictions, endpoint)
  direct <- lapply(predictions, sb_evaluate_external, endpoint = endpoint, reverse = FALSE)
  pipeline_ids <- vapply(predictions, `[[`, character(1L), "pipeline_id")
  names(direct) <- pipeline_ids

  expect_s3_class(observed, "stormbreaker_external_comparison")
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(observed$n_patients, 36L)
  expect_identical(observed$evaluations, direct)
  expect_identical(observed$comparison$pipeline_id, pipeline_ids)
  expect_equal(
    observed$comparison$roc_auc,
    unname(vapply(
      direct, function(x) x$metrics$roc_auc$estimate, numeric(1L)
    )),
    tolerance = 0
  )
  expect_equal(
    observed$comparison$brier_score,
    unname(vapply(
      direct, function(x) x$metrics$brier_score$estimate, numeric(1L)
    )),
    tolerance = 0
  )
  expect_equal(
    observed$comparison$mean_calibration_ratio,
    unname(vapply(
      direct,
      function(x) x$metrics$mean_calibration$mean_calibration_ratio,
      numeric(1L)
    )),
    tolerance = 0
  )
  expect_identical(
    unname(vapply(
      observed$evaluations,
      function(x) x$patient_match$output_order,
      character(1L)
    )),
    rep("external_endpoint", 2L)
  )
  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "No pipeline is ranked or selected", fixed = TRUE)
})

test_that("survival horizon comparisons copy existing external metrics", {
  development <- make_sb_fixture(n = 72L, seed = 8201L)
  development_task <- sb_surv_task(
    development,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 6
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    )
  )
  fits <- lapply(pipelines, function(pipeline) {
    sb_fit_final(development_task, pipeline)
  })

  external <- make_sb_fixture(n = 60L, seed = 8202L)
  external$patient_id <- paste0("external_", external$patient_id)
  predictors <- external[c(
    "patient_id", "gene_a", "gene_b", "gene_c"
  )]
  predictions <- lapply(fits, function(fit) {
    sb_predict_final(fit, predictors, "patient_id")
  })
  endpoint <- sb_surv_endpoint(
    external,
    "patient_id",
    "time",
    "event",
    "month",
    prediction_horizon = 6
  )

  observed <- sb_compare_external(predictions, endpoint)
  direct <- lapply(predictions, sb_evaluate_external, endpoint = endpoint, reverse = FALSE)
  pipeline_ids <- vapply(predictions, `[[`, character(1L), "pipeline_id")
  names(direct) <- pipeline_ids

  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "survival")
  expect_identical(observed$endpoint$prediction_horizon, 6)
  expect_identical(observed$evaluations, direct)
  expect_equal(
    observed$comparison$harrell_c,
    unname(vapply(
      direct, function(x) x$metrics$harrell_c$estimate, numeric(1L)
    )),
    tolerance = 0
  )
  expect_equal(
    observed$comparison$ipcw_brier,
    unname(vapply(
      direct, function(x) x$metrics$ipcw_brier$estimate, numeric(1L)
    )),
    tolerance = 0
  )
  expect_equal(
    observed$comparison$time_auc,
    unname(vapply(
      direct, function(x) x$metrics$time_auc$estimate, numeric(1L)
    )),
    tolerance = 0
  )
  expect_equal(
    observed$comparison$uno_c,
    unname(vapply(
      direct, function(x) x$metrics$uno_c$estimate, numeric(1L)
    )),
    tolerance = 0
  )
})

test_that("external comparison retains a horizon inherited from predictions", {
  index <- seq_len(30L)
  data <- data.frame(
    patient_id = paste0("inherited_", index),
    time = index,
    event = rep(c(1L, 0L, 1L), 10L),
    gene_a = seq(-1, 1, length.out = 30L),
    stringsAsFactors = FALSE
  )
  signatures <- list(
    sb_linear_signature(
      "inherited_a", "survival", c(gene_a = 0.3), intercept = 0,
      prediction_horizon = 6, time_unit = "month",
      baseline_survival = 0.8
    ),
    sb_linear_signature(
      "inherited_b", "survival", c(gene_a = -0.2), intercept = 0,
      prediction_horizon = 6, time_unit = "month",
      baseline_survival = 0.75
    )
  )
  predictions <- lapply(signatures, function(signature) {
    sb_predict_signature(
      signature,
      data[c("patient_id", "gene_a")],
      "patient_id"
    )
  })
  endpoint <- sb_surv_endpoint(
    data, "patient_id", "time", "event", "month"
  )

  observed <- sb_compare_external(predictions, endpoint)
  direct <- lapply(predictions, sb_evaluate_external, endpoint = endpoint, reverse = FALSE)

  expect_identical(observed$prediction_target$has_horizon, TRUE)
  expect_identical(observed$prediction_target$prediction_horizon, 6)
  expect_identical(observed$prediction_target$time_unit, "month")
  expect_identical(observed$prediction_target$horizon_source, "prediction")
  expect_identical("prediction_horizon" %in% names(observed$endpoint), FALSE)
  expect_identical(
    c("ipcw_brier", "time_auc", "uno_c") %in% names(observed$comparison),
    rep(TRUE, 3L)
  )
  expect_equal(
    observed$comparison$ipcw_brier,
    vapply(direct, function(x) x$metrics$ipcw_brier$estimate, numeric(1L)),
    tolerance = 0
  )
})

test_that("comparison rejects unlike inherited survival targets", {
  index <- seq_len(24L)
  data <- data.frame(
    patient_id = paste0("mixed_", index),
    time = index,
    event = rep(c(1L, 0L, 1L), 8L),
    gene_a = seq(-1, 1, length.out = 24L),
    stringsAsFactors = FALSE
  )
  make_prediction <- function(signature_id, horizon = NULL) {
    arguments <- list(
      signature_id = signature_id,
      task_type = "survival",
      coefficients = c(gene_a = 0.3),
      intercept = 0
    )
    if (!is.null(horizon)) {
      arguments$prediction_horizon <- horizon
      arguments$time_unit <- "month"
      arguments$baseline_survival <- 0.8
    }
    sb_predict_signature(
      do.call(sb_linear_signature, arguments),
      data[c("patient_id", "gene_a")],
      "patient_id"
    )
  }
  endpoint <- sb_surv_endpoint(
    data, "patient_id", "time", "event", "month"
  )
  different_horizons <- capture_sb_error(sb_compare_external(
    list(make_prediction("h6", 6), make_prediction("h12", 12)),
    endpoint
  ))
  mixed_presence <- capture_sb_error(sb_compare_external(
    list(make_prediction("risk"), make_prediction("horizon", 6)),
    endpoint
  ))

  expect_s3_class(
    different_horizons,
    "stormbreaker_external_target_mismatch"
  )
  expect_s3_class(mixed_presence, "stormbreaker_external_target_mismatch")
})

test_that("failed predictions retain recoverable survival targets", {
  index <- seq_len(20L)
  data <- data.frame(
    patient_id = paste0("failed_", index),
    time = index,
    event = rep(c(1L, 0L), 10L),
    gene_a = seq(-1, 1, length.out = 20L),
    stringsAsFactors = FALSE
  )
  signature <- sb_linear_signature(
    "failed_horizon", "survival", c(gene_a = 0.4), intercept = 0,
    prediction_horizon = 6, time_unit = "month",
    baseline_survival = 0.8
  )
  failed <- sb_predict_signature(
    signature,
    data["patient_id"],
    "patient_id"
  )
  endpoint <- sb_surv_endpoint(
    data, "patient_id", "time", "event", "month"
  )

  observed <- sb_compare_external(list(failed), endpoint)

  expect_identical(failed$status, "failed")
  expect_identical(observed$status, "incomplete")
  expect_identical(observed$prediction_target$target_status, "complete")
  expect_identical(observed$prediction_target$has_horizon, TRUE)
  expect_identical(observed$prediction_target$prediction_horizon, 6)
  expect_identical(observed$prediction_target$time_unit, "month")

  unknown <- failed
  unknown$signature <- NULL
  unknown$predictions <- unknown$predictions[0, , drop = FALSE]
  unknown_comparison <- sb_compare_external(list(unknown), endpoint)
  expect_identical(
    unknown_comparison$prediction_target$target_status,
    "unknown"
  )
  expect_identical(
    is.na(unknown_comparison$prediction_target$has_horizon),
    TRUE
  )
})

test_that("a failed prediction with unknown task makes the target partial", {
  index <- seq_len(20L)
  data <- data.frame(
    patient_id = paste0("partial_", index),
    time = index,
    event = rep(c(1L, 0L), 10L),
    gene_a = seq(-1, 1, length.out = 20L),
    stringsAsFactors = FALSE
  )
  known_signature <- sb_linear_signature(
    "known_horizon", "survival", c(gene_a = 0.4), intercept = 0,
    prediction_horizon = 6, time_unit = "month",
    baseline_survival = 0.8
  )
  known <- sb_predict_signature(
    known_signature,
    data[c("patient_id", "gene_a")],
    "patient_id"
  )
  unknown <- sb_predict_signature(
    list(
      signature_id = "unknown_task",
      task_type = "not_a_task",
      coefficients = c(gene_a = 0.4),
      intercept = 0
    ),
    data[c("patient_id", "gene_a")],
    "patient_id"
  )
  endpoint <- sb_surv_endpoint(
    data, "patient_id", "time", "event", "month"
  )

  observed <- sb_compare_external(list(known, unknown), endpoint)
  collection_error <- capture_sb_error(
    sb_collect_external_comparisons(list(cohort_a = observed))
  )

  expect_identical(unknown$status, "failed")
  expect_identical(is.na(unknown$task_type), TRUE)
  expect_identical(observed$status, "incomplete")
  expect_identical(observed$prediction_target$target_status, "partial")
  expect_s3_class(
    collection_error,
    "stormbreaker_external_target_mismatch"
  )
})

test_that("one invalid external prediction remains visible without fallback", {
  development <- make_sb_binary_fixture(n = 60L, seed = 8301L)
  development_task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  )
  fits <- lapply(pipelines, function(pipeline) {
    sb_fit_final(development_task, pipeline)
  })

  external <- make_sb_binary_fixture(n = 30L, seed = 8302L)
  external$patient_id <- paste0("external_", external$patient_id)
  predictors <- external[c(
    "patient_id", "gene_a", "gene_b", "gene_c"
  )]
  predictions <- lapply(fits, function(fit) {
    sb_predict_final(fit, predictors, "patient_id")
  })
  predictions[[2L]]$predictions$patient_id[[1L]] <- "wrong_patient"
  endpoint <- sb_binary_endpoint(
    external,
    "patient_id",
    "outcome"
  )

  observed <- sb_compare_external(predictions, endpoint)

  expect_identical(observed$status, "incomplete")
  expect_identical(observed$comparison$status, c("success", "failed"))
  expect_identical(
    observed$comparison$error_stage,
    c(NA_character_, "patient_match")
  )
  expect_identical(
    observed$comparison$error_class,
    c(NA_character_, "stormbreaker_external_patient_mismatch")
  )
  expect_identical(is.na(observed$comparison$roc_auc), c(FALSE, TRUE))
  expect_identical(observed$evaluations[[1L]]$status, "success")
  expect_identical(observed$evaluations[[2L]]$status, "failed")
  expect_identical(
    observed$evaluations[[2L]]$error_stage,
    "patient_match"
  )
})

test_that("external comparison requires distinct final predictions", {
  development <- make_sb_binary_fixture(n = 30L, seed = 8401L)
  task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  prediction <- sb_predict_final(
    fit,
    development[c("patient_id", "gene_a", "gene_b", "gene_c")],
    "patient_id"
  )
  endpoint <- sb_binary_endpoint(development, "patient_id", "outcome")

  empty <- capture_sb_error(sb_compare_external(list(), endpoint))
  duplicate <- capture_sb_error(sb_compare_external(
    list(prediction, prediction),
    endpoint
  ))

  expect_s3_class(empty, "stormbreaker_bad_prediction_collection")
  expect_s3_class(duplicate, "stormbreaker_duplicate_pipeline_id")
})
