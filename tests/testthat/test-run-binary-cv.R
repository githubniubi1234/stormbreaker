test_that("binary CV gives every patient one finite OOF probability", {
  d <- make_sb_binary_fixture(n = 60)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
  )
  result <- sb_run_cv(task, folds = 5, seed = 20260818)

  expect_s3_class(result, "stormbreaker_binary_cv_result")
  expect_s3_class(result, "stormbreaker_cv_result")
  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(result$pipeline$learner, "logistic")
  expect_identical(result$pipeline$task_type, "binary")
  expect_identical(nrow(result$predictions), nrow(d))
  expect_identical(anyDuplicated(result$predictions$patient_id), 0L)
  expect_equal(result$predictions$outcome, d$outcome)
  expect_identical(
    all(result$predictions$probability >= 0 &
          result$predictions$probability <= 1),
    TRUE
  )
  expect_named(
    result$metrics,
    c(
      "roc_auc", "average_precision", "brier_score", "log_loss",
      "calibration", "moderate_calibration"
    )
  )
  expect_identical(result$metrics$roc_auc$status, "success")
  expect_identical(result$metrics$average_precision$status, "success")
  expect_equal(
    result$metrics$average_precision$estimate,
    sum(
      result$metrics$average_precision$per_fold$n_positive *
        result$metrics$average_precision$per_fold$estimate
    ) / result$metrics$average_precision$n_positive,
    tolerance = 1e-15
  )
  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(
    printed,
    "Fold-stratified average precision:",
    fixed = TRUE
  )
  expect_identical(result$metrics$brier_score$status, "success")
  expect_equal(
    result$metrics$brier_score$estimate,
    mean(
      (result$predictions$outcome - result$predictions$probability)^2
    ),
    tolerance = 0
  )
  case_loss <- ifelse(
    result$predictions$outcome == 1L,
    -log(result$predictions$probability),
    -log1p(-result$predictions$probability)
  )
  expect_identical(result$metrics$log_loss$status, "success")
  expect_equal(
    result$metrics$log_loss$estimate,
    mean(case_loss),
    tolerance = 0
  )
  analysis_prevalence <- result$runs$n_analysis_positive /
    result$runs$n_analysis
  reference_probability <- analysis_prevalence[result$predictions$fold]
  reference_case_loss <- ifelse(
    result$predictions$outcome == 1L,
    -log(reference_probability),
    -log1p(-reference_probability)
  )
  expect_equal(
    result$metrics$log_loss$reference$estimate,
    mean(reference_case_loss),
    tolerance = 0
  )
  expect_equal(
    result$metrics$log_loss$skill_score$estimate,
    1 - sum(case_loss) / sum(reference_case_loss),
    tolerance = 0
  )
  expect_match(printed, "OOF log loss:", fixed = TRUE)
  expect_match(
    printed,
    "Log-loss skill vs analysis prevalence:",
    fixed = TRUE
  )
  logit_probability <- stats::qlogis(result$predictions$probability)
  direct_calibration_intercept <- stats::glm(
    result$predictions$outcome ~ 1,
    family = stats::binomial(),
    offset = logit_probability
  )
  direct_calibration_slope <- stats::glm(
    result$predictions$outcome ~ logit_probability,
    family = stats::binomial()
  )
  expect_identical(result$metrics$calibration$status, "success")
  expect_equal(
    result$metrics$calibration$calibration_intercept,
    unname(stats::coef(direct_calibration_intercept)[[1L]]),
    tolerance = 1e-12
  )
  expect_equal(
    result$metrics$calibration$calibration_slope,
    unname(stats::coef(direct_calibration_slope)[[2L]]),
    tolerance = 1e-12
  )
  expect_match(printed, "OOF calibration intercept:", fixed = TRUE)
  direct_moderate_calibration <- sb_binary_moderate_calibration_metric(
    result$predictions,
    "pooled_patient_level_out_of_fold"
  )
  expect_identical(result$metrics$moderate_calibration$status, "success")
  expect_identical(
    result$metrics$moderate_calibration$patient_calibration$patient_id,
    result$predictions$patient_id
  )
  expect_equal(
    result$metrics$moderate_calibration[c("ici", "e50", "e90")],
    direct_moderate_calibration[c("ici", "e50", "e90")],
    tolerance = 0
  )
  expect_match(
    printed,
    "OOF binary moderate calibration ICI:",
    fixed = TRUE
  )
  expect_length(result$coefficients, 5L)
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("success", 5L)
  )
})

test_that("binary OOF probabilities reproduce analysis-only glm fits", {
  d <- make_sb_binary_fixture(n = 60)
  predictors <- c("gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  result <- sb_run_cv(task, folds = 5, seed = 101)

  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    analysis_rows <- match(split$analysis_ids, d$patient_id)
    assessment_rows <- match(split$assessment_ids, d$patient_id)
    manual <- stats::glm(
      outcome ~ gene_a + gene_b,
      data = d[analysis_rows, , drop = FALSE],
      family = stats::binomial(link = "logit")
    )
    expected <- unname(stats::predict(
      manual,
      newdata = d[assessment_rows, , drop = FALSE],
      type = "response"
    ))
    observed <- result$predictions$probability[
      match(split$assessment_ids, result$predictions$patient_id)
    ]
    expect_equal(observed, expected, tolerance = 1e-12)

    audit <- result$coefficients[[fold]]$values
    design <- stats::model.matrix(
      ~ gene_a + gene_b,
      data = d[assessment_rows, , drop = FALSE]
    )
    reconstructed <- unname(drop(stats::plogis(
      design %*% audit$coefficient[match(colnames(design), audit$term)]
    )))
    expect_equal(observed, reconstructed, tolerance = 1e-12)
  }
})

test_that("binary CV is reproducible and folds never overlap", {
  d <- make_sb_binary_fixture(n = 60)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b")
  )
  first <- sb_run_cv(task, folds = 5, seed = 77)
  second <- sb_run_cv(task, folds = 5, seed = 77)

  expect_identical(first$predictions, second$predictions)
  expect_identical(first$metrics, second$metrics)
  expect_identical(first$runs, second$runs)
  expect_identical(first$folds, second$folds)
  for (split in first$folds) {
    expect_length(intersect(split$analysis_ids, split$assessment_ids), 0L)
    expect_setequal(
      c(split$analysis_ids, split$assessment_ids),
      d$patient_id
    )
  }
})

test_that("binary fit failures remain visible without fallback", {
  d <- make_sb_binary_fixture(n = 60)
  d$gene_b <- d$gene_a
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b")
  )
  result <- sb_run_cv(task, folds = 5, seed = 9)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$status), "failed")
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(all(is.na(result$predictions$probability)), TRUE)
  expect_identical(result$metrics$roc_auc$status, "incomplete")
  expect_identical(is.na(result$metrics$roc_auc$estimate), TRUE)
  expect_identical(
    result$metrics$average_precision$status,
    "incomplete"
  )
  expect_identical(
    is.na(result$metrics$average_precision$estimate),
    TRUE
  )
  expect_identical(result$metrics$brier_score$status, "incomplete")
  expect_identical(is.na(result$metrics$brier_score$estimate), TRUE)
  expect_identical(result$metrics$log_loss$status, "incomplete")
  expect_identical(is.na(result$metrics$log_loss$estimate), TRUE)
  expect_identical(result$metrics$log_loss$reference$status, "success")
  expect_identical(result$metrics$log_loss$skill_score$status, "incomplete")
  expect_identical(result$metrics$calibration$status, "incomplete")
  expect_identical(
    result$metrics$calibration$input_status,
    "prediction_failed"
  )
  expect_identical(
    result$metrics$moderate_calibration$status,
    "incomplete"
  )
  expect_identical(
    result$metrics$moderate_calibration$input_status,
    "prediction_failed"
  )
  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(
    printed,
    "Fold-stratified average precision: unavailable (incomplete)",
    fixed = TRUE
  )
  expect_match(
    printed,
    "OOF log loss: unavailable (incomplete)",
    fixed = TRUE
  )
  expect_match(
    printed,
    "Log-loss skill: unavailable (incomplete)",
    fixed = TRUE
  )
  expect_match(
    printed,
    "OOF calibration: unavailable (incomplete)",
    fixed = TRUE
  )
  expect_match(
    printed,
    "OOF binary moderate calibration: unavailable",
    fixed = TRUE
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_fitted", 5L)
  )
})

test_that("binary CV fits one optional moderate curve over pooled OOF rows", {
  d <- make_sb_binary_fixture(n = 60L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  original_metric <- sb_binary_moderate_calibration_metric
  calls <- 0L
  observed_scored <- NULL
  testthat::local_mocked_bindings(
    sb_binary_moderate_calibration_metric = function(scored, scope) {
      calls <<- calls + 1L
      observed_scored <<- scored
      original_metric(scored, scope)
    }
  )

  result <- sb_run_cv(task, folds = 5L, seed = 20260820L)

  expect_identical(calls, 1L)
  expect_identical(observed_scored, result$predictions)
  expect_identical(result$metrics$moderate_calibration$status, "success")
  expect_identical(
    result$metrics$moderate_calibration$per_fold$status,
    rep("success", 5L)
  )
  expect_identical(
    result$metrics$moderate_calibration$n_folds,
    5L
  )
})

test_that("a binary moderate-calibration error does not erase OOF evidence", {
  d <- make_sb_binary_fixture(n = 60L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  control_result <- sb_run_cv(task, folds = 5L, seed = 114L)
  testthat::local_mocked_bindings(
    sb_binary_moderate_calibration_metric = function(...) {
      condition <- simpleError("injected binary moderate calibration failure")
      class(condition) <- c(
        "test_binary_moderate_calibration_error",
        class(condition)
      )
      stop(condition)
    }
  )

  observed <- sb_run_cv(task, folds = 5L, seed = 114L)

  expect_identical(observed$status, control_result$status)
  expect_identical(observed$complete_oof, control_result$complete_oof)
  expect_identical(observed$predictions, control_result$predictions)
  expect_identical(observed$runs, control_result$runs)
  expect_identical(
    observed$metrics[names(observed$metrics) != "moderate_calibration"],
    control_result$metrics[
      names(control_result$metrics) != "moderate_calibration"
    ]
  )
  expect_identical(
    observed$metrics$moderate_calibration$status,
    "metric_failed"
  )
  expect_identical(
    observed$metrics$moderate_calibration$error_class,
    "test_binary_moderate_calibration_error"
  )
})

test_that("binary folds require enough patients in both classes", {
  d <- make_sb_binary_fixture(n = 20)
  d$outcome <- c(rep(1L, 4), rep(0L, 16))
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")

  error <- capture_sb_error(sb_run_cv(task, folds = 5, seed = 1))

  expect_s3_class(error, "stormbreaker_infeasible_folds")
})

test_that("task and pipeline outcome types cannot be mixed", {
  binary_data <- make_sb_binary_fixture(n = 20)
  binary_task <- sb_binary_task(
    binary_data, "patient_id", "outcome", "gene_a"
  )
  survival_data <- make_sb_fixture(n = 20)
  survival_task <- sb_surv_task(
    survival_data, "patient_id", "time", "event", "gene_a", "month"
  )

  binary_error <- capture_sb_error(
    sb_run_cv(binary_task, sb_pipeline(), folds = 2, seed = 1)
  )
  survival_error <- capture_sb_error(
    sb_run_cv(
      survival_task,
      sb_pipeline(learner = "logistic"),
      folds = 2,
      seed = 1
    )
  )

  expect_s3_class(binary_error, "stormbreaker_pipeline_task_mismatch")
  expect_s3_class(survival_error, "stormbreaker_pipeline_task_mismatch")
})

test_that("assessment outcomes cannot change fixed-split probabilities", {
  d <- make_sb_binary_fixture(n = 60)
  split <- list(
    fold = 1L,
    analysis_rows = 1:48,
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[1:48],
    assessment_ids = d$patient_id[49:60]
  )
  pipeline <- sb_pipeline(learner = "logistic")
  make_task <- function(data) {
    sb_binary_task(
      data, "patient_id", "outcome", c("gene_a", "gene_b")
    )
  }
  original <- sb_run_one_binary_fold(
    make_task(d), pipeline, sb_resolve_pipeline_runtime(pipeline), split,
    seed = 1L
  )
  changed <- d
  changed$outcome[49:60] <- 1L - changed$outcome[49:60]
  perturbed <- sb_run_one_binary_fold(
    make_task(changed), pipeline, sb_resolve_pipeline_runtime(pipeline),
    split, seed = 1L
  )

  expect_equal(
    original$predictions$probability,
    perturbed$predictions$probability,
    tolerance = 0
  )

  brier_from_fold <- function(fold_result) {
    reference <- data.frame(
      fold = fold_result$run$fold,
      analysis_prevalence =
        fold_result$run$n_analysis_positive / fold_result$run$n_analysis
    )
    sb_binary_brier_metric(fold_result$predictions, reference)
  }
  original_brier <- brier_from_fold(original)
  perturbed_brier <- brier_from_fold(perturbed)
  expect_identical(
    original_brier$reference$per_fold$analysis_prevalence,
    perturbed_brier$reference$per_fold$analysis_prevalence
  )
  reference_probability <-
    original_brier$reference$per_fold$analysis_prevalence
  expect_equal(
    original_brier$reference$estimate,
    mean((original$predictions$outcome - reference_probability)^2),
    tolerance = 0
  )
  expect_equal(
    perturbed_brier$reference$estimate,
    mean((perturbed$predictions$outcome - reference_probability)^2),
    tolerance = 0
  )

  log_loss_from_fold <- function(fold_result) {
    reference <- data.frame(
      fold = fold_result$run$fold,
      analysis_prevalence =
        fold_result$run$n_analysis_positive / fold_result$run$n_analysis
    )
    sb_binary_log_loss_metric(fold_result$predictions, reference)
  }
  original_log_loss <- log_loss_from_fold(original)
  perturbed_log_loss <- log_loss_from_fold(perturbed)
  expect_identical(
    original_log_loss$reference$per_fold$analysis_prevalence,
    perturbed_log_loss$reference$per_fold$analysis_prevalence
  )
  original_reference_loss <- ifelse(
    original$predictions$outcome == 1L,
    -log(reference_probability),
    -log1p(-reference_probability)
  )
  perturbed_reference_loss <- ifelse(
    perturbed$predictions$outcome == 1L,
    -log(reference_probability),
    -log1p(-reference_probability)
  )
  expect_equal(
    original_log_loss$reference$estimate,
    mean(original_reference_loss),
    tolerance = 0
  )
  expect_equal(
    perturbed_log_loss$reference$estimate,
    mean(perturbed_reference_loss),
    tolerance = 0
  )
})

test_that("binary probability references use each outer analysis prevalence", {
  d <- make_sb_binary_fixture(n = 60L, seed = 2101L)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b")
  )
  result <- sb_run_cv(task, folds = 5L, seed = 2201L)
  reference <- result$metrics$brier_score$reference
  skill <- result$metrics$brier_score$skill_score
  analysis_prevalence <- vapply(result$folds, function(split) {
    mean(d$outcome[match(split$analysis_ids, d$patient_id)])
  }, numeric(1))
  patient_reference <- analysis_prevalence[result$predictions$fold]
  reference_loss <- (
    result$predictions$outcome - patient_reference
  )^2
  model_loss <- (
    result$predictions$outcome - result$predictions$probability
  )^2

  expect_identical(reference$status, "success")
  expect_equal(
    reference$per_fold$analysis_prevalence,
    analysis_prevalence,
    tolerance = 0
  )
  expect_equal(reference$estimate, mean(reference_loss), tolerance = 0)
  expect_identical(skill$status, "success")
  expect_equal(
    skill$estimate,
    1 - sum(model_loss) / sum(reference_loss),
    tolerance = 0
  )

  log_loss <- result$metrics$log_loss
  reference_log_loss <- ifelse(
    result$predictions$outcome == 1L,
    -log(patient_reference),
    -log1p(-patient_reference)
  )
  model_log_loss <- ifelse(
    result$predictions$outcome == 1L,
    -log(result$predictions$probability),
    -log1p(-result$predictions$probability)
  )
  expect_identical(log_loss$reference$status, "success")
  expect_equal(
    log_loss$reference$per_fold$analysis_prevalence,
    analysis_prevalence,
    tolerance = 0
  )
  expect_equal(
    log_loss$reference$estimate,
    mean(reference_log_loss),
    tolerance = 0
  )
  expect_identical(log_loss$skill_score$status, "success")
  expect_equal(
    log_loss$skill_score$estimate,
    1 - sum(model_log_loss) / sum(reference_log_loss),
    tolerance = 1e-14
  )
})

test_that("the binary result schema is explicit", {
  d <- make_sb_binary_fixture(n = 20)
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
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
      "patient_id", "fold", "outcome", "positive_class", "probability",
      "status", "pipeline_id", "learner"
    ),
    ignore.order = FALSE
  )
})

test_that("mutated binary task metadata is rejected", {
  d <- make_sb_binary_fixture(n = 20)
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  task$n_positive <- task$n_positive + 1L

  error <- capture_sb_error(sb_run_cv(task, folds = 2, seed = 1))

  expect_s3_class(error, "stormbreaker_mutated_task")
})

test_that("nested glmnet logistic supports predictors exceeding patients", {
  skip_if_not_installed("glmnet")

  d <- make_sb_high_dim_binary_fixture(n = 60, p = 80)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 1701L
  )
  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 101L)

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(nrow(result$predictions), 60L)
  expect_identical(
    all(is.finite(result$predictions$probability)),
    TRUE
  )
  expect_length(result$tuning, 5L)
  expect_identical(
    vapply(result$tuning, `[[`, "", "status"),
    rep("success", 5L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("success", 5L)
  )
  expect_identical(
    result$provenance$glmnet_version,
    as.character(utils::packageVersion("glmnet"))
  )

  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    audit <- result$tuning[[fold]]
    expect_setequal(names(audit$fold_id), split$analysis_ids)
    expect_length(intersect(names(audit$fold_id), split$assessment_ids), 0L)
    expect_identical(all(is.finite(audit$curve$mean_deviance)), TRUE)
  }
})

test_that("nested binary OOF probabilities match direct glmnet execution", {
  skip_if_not_installed("glmnet")

  d <- make_sb_high_dim_binary_fixture(n = 60, p = 20)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 0.4,
    inner_folds = 3L,
    inner_seed = 301L,
    glmnet_standardize = FALSE,
    glmnet_lambda_rule = "lambda.1se",
    glmnet_tuning_measure = "auc"
  )
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 401L)
  split <- result$folds[[1L]]
  analysis_rows <- match(split$analysis_ids, d$patient_id)
  assessment_rows <- match(split$assessment_ids, d$patient_id)
  tuning <- result$tuning[[1L]]
  manual_cv <- glmnet::cv.glmnet(
    x = as.matrix(d[analysis_rows, predictors, drop = FALSE]),
    y = d$outcome[analysis_rows],
    family = "binomial",
    alpha = 0.4,
    foldid = unname(tuning$fold_id[split$analysis_ids]),
    type.measure = "auc",
    standardize = FALSE
  )
  manual_fit <- glmnet::glmnet(
    x = as.matrix(d[analysis_rows, predictors, drop = FALSE]),
    y = d$outcome[analysis_rows],
    family = "binomial",
    alpha = 0.4,
    standardize = FALSE
  )
  expected <- unname(drop(stats::predict(
    manual_fit,
    newx = as.matrix(d[assessment_rows, predictors, drop = FALSE]),
    s = manual_cv$lambda.1se,
    type = "response"
  )))
  observed <- result$predictions$probability[
    match(split$assessment_ids, result$predictions$patient_id)
  ]

  expect_equal(tuning$lambda, unname(manual_cv$lambda.1se), tolerance = 0)
  expect_identical(tuning$selection_rule, "lambda.1se")
  expect_identical(tuning$standardize, FALSE)
  expect_identical(tuning$measure, "auc")
  expect_equal(observed, expected, tolerance = 1e-12)

  coefficient_audit <- result$coefficients[[1L]]$values
  design <- cbind(
    "(Intercept)" = 1,
    as.matrix(d[assessment_rows, predictors, drop = FALSE])
  )
  reconstructed <- unname(drop(stats::plogis(
    design %*% coefficient_audit$coefficient[
      match(colnames(design), coefficient_audit$term)
    ]
  )))
  expect_equal(observed, reconstructed, tolerance = 1e-12)
})

test_that("outer assessment changes cannot alter nested binary fitting", {
  skip_if_not_installed("glmnet")

  d <- make_sb_high_dim_binary_fixture(n = 60, p = 20)
  predictors <- grep("^gene_", names(d), value = TRUE)
  split <- list(
    fold = 1L,
    analysis_rows = 1:48,
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[1:48],
    assessment_ids = d$patient_id[49:60]
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 1,
    inner_folds = 3L,
    inner_seed = 501L
  )
  make_task <- function(data) {
    sb_binary_task(data, "patient_id", "outcome", predictors)
  }
  original <- sb_run_one_binary_fold(
    make_task(d), pipeline, sb_resolve_pipeline_runtime(pipeline), split,
    seed = 1L
  )
  changed <- d
  changed$outcome[49:60] <- 1L - changed$outcome[49:60]
  perturbed <- sb_run_one_binary_fold(
    make_task(changed), pipeline, sb_resolve_pipeline_runtime(pipeline),
    split, seed = 1L
  )

  expect_identical(original$tuning, perturbed$tuning)
  expect_equal(
    original$predictions$probability,
    perturbed$predictions$probability,
    tolerance = 0
  )
})

test_that("nested binary tuning failures remain visible without fallback", {
  skip_if_not_installed("glmnet")

  d <- make_sb_binary_fixture(n = 30)
  d$constant_a <- 1
  d$constant_b <- 2
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("constant_a", "constant_b")
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 1L
  )
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 1L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$status), "failed")
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_glmnet_logistic_tuning_failure"
  )
  expect_identical(
    vapply(result$tuning, `[[`, "", "status"),
    rep("failed", 3L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_fitted", 3L)
  )
  expect_identical(all(is.na(result$predictions$probability)), TRUE)
  expect_identical(result$metrics$roc_auc$status, "incomplete")
  expect_identical(result$metrics$brier_score$status, "incomplete")
})

test_that("binary AUC Top-K selects and audits features inside every outer fold", {
  d <- make_sb_high_dim_binary_fixture(n = 80L, p = 30L, seed = 1201L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 3L
  )
  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 1301L)

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_length(result$selections, 4L)
  expect_length(result$coefficients, 4L)
  expect_identical(
    vapply(result$selections, `[[`, "", "status"),
    rep("success", 4L)
  )
  stability <- result$selection_stability
  expected_count <- vapply(predictors, function(feature) {
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
  expect_identical(stability$n_successful_folds, 4L)
  expect_identical(stability$feature_frequency$feature, predictors)
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

  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    audit <- result$selections[[fold]]
    analysis_rows <- match(split$analysis_ids, d$patient_id)
    assessment_rows <- match(split$assessment_ids, d$patient_id)
    direct <- sb_fit_univ_binary_auc_top_k(
      d[analysis_rows, predictors, drop = FALSE],
      d$outcome[analysis_rows],
      n_features = 3L
    )

    expect_identical(audit$fold, fold)
    expect_identical(audit$selector, "univ_binary_auc_top_k")
    expect_identical(audit$selected_features, direct$selected_features)
    expect_identical(audit$scores, direct$scores)
    expect_length(audit$selected_features, 3L)
    expect_identical(sum(audit$scores$selected), 3L)
    expect_identical(nrow(audit$scores), length(predictors))

    coefficients <- result$coefficients[[fold]]$values
    expect_identical(
      coefficients$term,
      c("(Intercept)", audit$selected_features)
    )
    design <- cbind(
      "(Intercept)" = 1,
      as.matrix(d[assessment_rows, audit$selected_features, drop = FALSE])
    )
    rebuilt <- unname(drop(stats::plogis(
      design %*% coefficients$coefficient[
        match(colnames(design), coefficients$term)
      ]
    )))
    observed <- result$predictions$probability[
      match(split$assessment_ids, result$predictions$patient_id)
    ]
    expect_equal(observed, rebuilt, tolerance = 1e-12)
  }
})

test_that("outer assessment data cannot alter binary AUC selection or fitting", {
  d <- make_sb_high_dim_binary_fixture(n = 60L, p = 12L, seed = 1401L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  split <- list(
    fold = 1L,
    analysis_rows = 1:48,
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[1:48],
    assessment_ids = d$patient_id[49:60]
  )
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 3L
  )
  run_fixed <- function(data) {
    task <- sb_binary_task(data, "patient_id", "outcome", predictors)
    sb_run_one_binary_fold(
      task,
      pipeline,
      sb_resolve_pipeline_runtime(pipeline),
      split,
      seed = 1L
    )
  }
  original <- run_fixed(d)

  changed_outcome <- d
  changed_outcome$outcome[split$assessment_rows] <-
    1L - changed_outcome$outcome[split$assessment_rows]
  outcome_run <- run_fixed(changed_outcome)

  expect_identical(original$selection, outcome_run$selection)
  expect_identical(original$coefficients, outcome_run$coefficients)
  expect_equal(
    original$predictions$probability,
    outcome_run$predictions$probability,
    tolerance = 0
  )

  changed_predictors <- d
  selected <- original$selection$selected_features[[1L]]
  changed_predictors[split$assessment_rows, selected] <-
    changed_predictors[split$assessment_rows, selected] + 100
  predictor_run <- run_fixed(changed_predictors)

  expect_identical(original$selection, predictor_run$selection)
  expect_identical(original$coefficients, predictor_run$coefficients)
  expect_false(isTRUE(all.equal(
    original$predictions$probability,
    predictor_run$predictions$probability,
    tolerance = 0
  )))
})

test_that("learned binary selector results expose selection audit explicitly", {
  d <- make_sb_binary_fixture(n = 40L, seed = 1501L)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
  )
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 1L
  )
  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 1601L)

  expect_named(
    result,
    c(
      "status", "complete_oof", "predictions", "metrics", "runs",
      "selections", "selection_stability", "coefficients", "folds", "task",
      "pipeline", "provenance", "glmnet_active_set_stability"
    ),
    ignore.order = FALSE
  )
})

test_that("binary selector failures remain visible without fallback", {
  d <- make_sb_binary_fixture(n = 40L, seed = 1901L)
  d$constant_a <- 1
  d$constant_b <- 2
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("constant_a", "constant_b")
  )
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 1L
  )
  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 2001L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$status), "failed")
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_insufficient_valid_features"
  )
  expect_identical(
    vapply(result$selections, `[[`, "", "status"),
    rep("failed", 4L)
  )
  expect_identical(
    vapply(result$selections, `[[`, "", "error_class"),
    rep("stormbreaker_insufficient_valid_features", 4L)
  )
  expect_identical(result$selection_stability$status, "incomplete")
  expect_identical(result$selection_stability$n_successful_folds, 0L)
  expect_identical(result$selection_stability$n_evaluable_pairs, 0L)
  expect_true(is.na(result$selection_stability$mean_pairwise_jaccard))
  expect_identical(
    result$selection_stability$per_fold$error_class,
    rep("stormbreaker_insufficient_valid_features", 4L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_fitted", 4L)
  )
  expect_identical(all(is.na(result$predictions$probability)), TRUE)
  expect_identical(result$metrics$roc_auc$status, "incomplete")
  expect_identical(result$metrics$brier_score$status, "incomplete")
})

test_that("Gaussian Naive Bayes supports high-dimensional numeric candidates", {
  d <- make_sb_high_dim_binary_fixture(n = 60L, p = 80L, seed = 2301L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipeline <- sb_pipeline(
    learner = "gaussian_naive_bayes",
    naive_bayes_var_smoothing = 1e-9
  )
  result <- sb_run_cv(task, pipeline, folds = 5L, seed = 2401L)

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(all(is.finite(result$predictions$probability)), TRUE)
  expect_length(result$distribution_parameters, 5L)
  expect_identical(
    vapply(result$distribution_parameters, `[[`, "", "status"),
    rep("success", 5L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_applicable", 5L)
  )
  expect_identical(result$metrics$roc_auc$status, "success")
  expect_identical(result$metrics$brier_score$status, "success")
  expect_identical(
    result$metrics$brier_score$reference$status,
    "success"
  )

  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    audit <- result$distribution_parameters[[fold]]
    assessment_rows <- match(split$assessment_ids, d$patient_id)
    means <- do.call(rbind, lapply(c(0L, 1L), function(class) {
      values <- audit$values[audit$values$class == class, , drop = FALSE]
      stats::setNames(values$mean, values$feature)
    }))
    rownames(means) <- c("0", "1")
    variances <- do.call(rbind, lapply(c(0L, 1L), function(class) {
      values <- audit$values[audit$values$class == class, , drop = FALSE]
      stats::setNames(values$variance, values$feature)
    }))
    rownames(variances) <- c("0", "1")
    priors <- stats::setNames(audit$priors$prior, audit$priors$class)
    rebuilt <- manual_gaussian_nb_probability(
      as.matrix(d[assessment_rows, predictors, drop = FALSE]),
      priors,
      means,
      variances
    )
    observed <- result$predictions$probability[
      match(split$assessment_ids, result$predictions$patient_id)
    ]
    expect_equal(observed, rebuilt, tolerance = 1e-14)
    expect_identical(nrow(audit$values), 2L * length(predictors))
  }
})

test_that("assessment data cannot alter Gaussian Naive Bayes fitting", {
  d <- make_sb_high_dim_binary_fixture(n = 60L, p = 12L, seed = 2501L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  split <- list(
    fold = 1L,
    analysis_rows = 1:48,
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[1:48],
    assessment_ids = d$patient_id[49:60]
  )
  pipeline <- sb_pipeline(learner = "gaussian_naive_bayes")
  run_fixed <- function(data) {
    task <- sb_binary_task(data, "patient_id", "outcome", predictors)
    sb_run_one_binary_fold(
      task,
      pipeline,
      sb_resolve_pipeline_runtime(pipeline),
      split,
      seed = 1L
    )
  }
  original <- run_fixed(d)

  changed_outcome <- d
  changed_outcome$outcome[split$assessment_rows] <-
    1L - changed_outcome$outcome[split$assessment_rows]
  outcome_run <- run_fixed(changed_outcome)
  expect_identical(
    original$distribution_parameters,
    outcome_run$distribution_parameters
  )
  expect_equal(
    original$predictions$probability,
    outcome_run$predictions$probability,
    tolerance = 0
  )

  changed_predictors <- d
  changed_predictors[split$assessment_rows, predictors[[1L]]] <-
    changed_predictors[split$assessment_rows, predictors[[1L]]] + 100
  predictor_run <- run_fixed(changed_predictors)
  expect_identical(
    original$distribution_parameters,
    predictor_run$distribution_parameters
  )
  expect_false(isTRUE(all.equal(
    original$predictions$probability,
    predictor_run$predictions$probability,
    tolerance = 0
  )))
})

test_that("Gaussian Naive Bayes variance failures have no fallback", {
  d <- make_sb_binary_fixture(n = 40L, seed = 2601L)
  d$constant_a <- 1
  d$constant_b <- 2
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("constant_a", "constant_b")
  )
  pipeline <- sb_pipeline(learner = "gaussian_naive_bayes")
  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 2701L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$status), "failed")
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_invalid_naive_bayes_variance"
  )
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(
    vapply(result$distribution_parameters, `[[`, "", "status"),
    rep("not_fitted", 4L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_applicable", 4L)
  )
  expect_identical(all(is.na(result$predictions$probability)), TRUE)
  expect_identical(result$metrics$brier_score$status, "incomplete")
  expect_identical(
    result$metrics$brier_score$reference$status,
    "success"
  )
  expect_identical(
    result$metrics$brier_score$skill_score$status,
    "incomplete"
  )
})

test_that("fixed LightGBM supports high-dimensional numeric candidates", {
  skip_if_not_installed("lightgbm", minimum_version = "4.0.0")

  d <- make_sb_high_dim_binary_fixture(n = 48L, p = 64L, seed = 3201L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipeline <- sb_pipeline(
    learner = "lightgbm_binary",
    lightgbm_nrounds = 8L,
    lightgbm_learning_rate = 0.1,
    lightgbm_num_leaves = 7L,
    lightgbm_min_data_in_leaf = 3L,
    lightgbm_seed = 3301L
  )
  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 3401L)

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(all(is.finite(result$predictions$probability)), TRUE)
  expect_length(result$tree_models, 4L)
  expect_identical(
    vapply(result$tree_models, `[[`, "", "status"),
    rep("success", 4L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_applicable", 4L)
  )
  expect_identical(result$metrics$roc_auc$status, "success")
  expect_identical(result$metrics$brier_score$status, "success")
  expect_identical(
    result$metrics$brier_score$reference$status,
    "success"
  )

  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    audit <- result$tree_models[[fold]]
    assessment_rows <- match(split$assessment_ids, d$patient_id)
    restored <- lightgbm::lgb.load(model_str = audit$model_string)
    rebuilt <- as.numeric(stats::predict(
      restored,
      as.matrix(d[assessment_rows, audit$features, drop = FALSE]),
      type = "response",
      num_iteration = audit$nrounds,
      params = list(num_threads = 1L)
    ))
    observed <- result$predictions$probability[
      match(split$assessment_ids, result$predictions$patient_id)
    ]
    expect_equal(observed, rebuilt, tolerance = 0)
    expect_identical(audit$features, predictors)
    expect_identical(audit$feature_importance$feature, predictors)
  }
})

test_that("assessment data cannot alter fixed LightGBM fitting", {
  skip_if_not_installed("lightgbm", minimum_version = "4.0.0")

  d <- make_sb_high_dim_binary_fixture(n = 60L, p = 12L, seed = 3501L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  split <- list(
    fold = 1L,
    analysis_rows = 1:48,
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[1:48],
    assessment_ids = d$patient_id[49:60]
  )
  pipeline <- sb_pipeline(
    learner = "lightgbm_binary",
    lightgbm_nrounds = 8L,
    lightgbm_num_leaves = 7L,
    lightgbm_min_data_in_leaf = 3L,
    lightgbm_seed = 3601L
  )
  run_fixed <- function(data) {
    task <- sb_binary_task(data, "patient_id", "outcome", predictors)
    sb_run_one_binary_fold(
      task,
      pipeline,
      sb_resolve_pipeline_runtime(pipeline),
      split,
      seed = 1L
    )
  }
  original <- run_fixed(d)

  changed_outcome <- d
  changed_outcome$outcome[split$assessment_rows] <-
    1L - changed_outcome$outcome[split$assessment_rows]
  outcome_run <- run_fixed(changed_outcome)
  expect_identical(original$tree_model, outcome_run$tree_model)
  expect_equal(
    original$predictions$probability,
    outcome_run$predictions$probability,
    tolerance = 0
  )

  changed_predictors <- d
  changed_predictors[split$assessment_rows, predictors[[1L]]] <-
    changed_predictors[split$assessment_rows, predictors[[1L]]] + 100
  predictor_run <- run_fixed(changed_predictors)
  expect_identical(original$tree_model, predictor_run$tree_model)
  expect_false(isTRUE(all.equal(
    original$predictions$probability,
    predictor_run$predictions$probability,
    tolerance = 0
  )))
})

test_that("missing fixed LightGBM dependency has no fallback", {
  d <- make_sb_binary_fixture(n = 40L, seed = 3701L)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
  )
  pipeline <- sb_pipeline(
    learner = "lightgbm_binary",
    lightgbm_nrounds = 5L,
    lightgbm_min_data_in_leaf = 2L
  )
  local_mocked_bindings(sb_lightgbm_available = function() FALSE)

  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 3801L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$status), "failed")
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_missing_lightgbm"
  )
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(
    vapply(result$tree_models, `[[`, "", "status"),
    rep("not_fitted", 4L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_applicable", 4L)
  )
  expect_identical(all(is.na(result$predictions$probability)), TRUE)
  expect_identical(result$metrics$roc_auc$status, "incomplete")
  expect_identical(result$metrics$brier_score$status, "incomplete")
})

test_that("fixed CatBoost supports high-dimensional numeric candidates", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_high_dim_binary_fixture(n = 48L, p = 64L, seed = 4601L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipeline <- sb_pipeline(
    learner = "catboost_binary",
    catboost_iterations = 8L,
    catboost_learning_rate = 0.1,
    catboost_depth = 3L,
    catboost_l2_leaf_reg = 3,
    catboost_seed = 4701L
  )
  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 4801L)

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(all(is.finite(result$predictions$probability)), TRUE)
  expect_length(result$tree_models, 4L)
  expect_identical(
    vapply(result$tree_models, `[[`, "", "status"),
    rep("success", 4L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_applicable", 4L)
  )
  expect_identical(result$metrics$roc_auc$status, "success")
  expect_identical(result$metrics$brier_score$status, "success")
  expect_identical(
    result$metrics$brier_score$reference$status,
    "success"
  )
  expect_identical(
    result$provenance$catboost_version,
    as.character(utils::packageVersion("catboost"))
  )

  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    audit <- result$tree_models[[fold]]
    assessment_rows <- match(split$assessment_ids, d$patient_id)
    restored <- unserialize(audit$model_raw)
    restored_pool <- catboost::catboost.load_pool(as.matrix(
      d[assessment_rows, audit$features, drop = FALSE]
    ))
    rebuilt <- as.numeric(catboost::catboost.predict(
      restored,
      restored_pool,
      prediction_type = "Probability",
      ntree_start = 0L,
      ntree_end = audit$iterations,
      thread_count = 1L
    ))
    observed <- result$predictions$probability[
      match(split$assessment_ids, result$predictions$patient_id)
    ]
    expect_equal(observed, rebuilt, tolerance = 0)
    expect_identical(audit$tree_count, 8L)
    expect_identical(audit$features, predictors)
    expect_identical(audit$feature_importance$feature, predictors)
  }
})

test_that("assessment data cannot alter fixed CatBoost fitting", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_high_dim_binary_fixture(n = 60L, p = 12L, seed = 4901L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  split <- list(
    fold = 1L,
    analysis_rows = 1:48,
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[1:48],
    assessment_ids = d$patient_id[49:60]
  )
  pipeline <- sb_pipeline(
    learner = "catboost_binary",
    catboost_iterations = 8L,
    catboost_depth = 3L,
    catboost_seed = 5001L
  )
  run_fixed <- function(data) {
    task <- sb_binary_task(data, "patient_id", "outcome", predictors)
    sb_run_one_binary_fold(
      task,
      pipeline,
      sb_resolve_pipeline_runtime(pipeline),
      split,
      seed = 1L
    )
  }
  predict_saved <- function(audit, data) {
    model <- unserialize(audit$model_raw)
    pool <- catboost::catboost.load_pool(as.matrix(
      data[split$assessment_rows, audit$features, drop = FALSE]
    ))
    as.numeric(catboost::catboost.predict(
      model,
      pool,
      prediction_type = "Probability",
      ntree_start = 0L,
      ntree_end = audit$iterations,
      thread_count = 1L
    ))
  }
  original <- run_fixed(d)

  changed_outcome <- d
  changed_outcome$outcome[split$assessment_rows] <-
    1L - changed_outcome$outcome[split$assessment_rows]
  outcome_run <- run_fixed(changed_outcome)
  expect_equal(
    original$predictions$probability,
    outcome_run$predictions$probability,
    tolerance = 0
  )
  expect_equal(
    predict_saved(original$tree_model, d),
    predict_saved(outcome_run$tree_model, d),
    tolerance = 0
  )
  expect_identical(
    original$tree_model$feature_importance,
    outcome_run$tree_model$feature_importance
  )

  changed_predictors <- d
  changed_predictors[split$assessment_rows, predictors] <-
    -changed_predictors[split$assessment_rows, predictors] + 100
  predictor_run <- run_fixed(changed_predictors)
  expect_equal(
    predict_saved(original$tree_model, d),
    predict_saved(predictor_run$tree_model, d),
    tolerance = 0
  )
  expect_identical(
    original$tree_model$feature_importance,
    predictor_run$tree_model$feature_importance
  )
  expect_false(isTRUE(all.equal(
    original$predictions$probability,
    predictor_run$predictions$probability,
    tolerance = 0
  )))
})

test_that("missing fixed CatBoost dependency has no fallback", {
  d <- make_sb_binary_fixture(n = 40L, seed = 5101L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c"),
    decision_thresholds = c(0.2, 0.5)
  )
  pipeline <- sb_pipeline(
    learner = "catboost_binary",
    catboost_iterations = 5L,
    catboost_depth = 2L
  )
  local_mocked_bindings(sb_catboost_available = function() FALSE)

  result <- sb_run_cv(task, pipeline, folds = 4L, seed = 5201L)

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(unique(result$runs$status), "failed")
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_missing_catboost"
  )
  expect_identical(all(is.na(result$runs$learner_used)), TRUE)
  expect_identical(
    vapply(result$tree_models, `[[`, "", "status"),
    rep("not_fitted", 4L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_applicable", 4L)
  )
  expect_identical(all(is.na(result$predictions$probability)), TRUE)
  expect_identical(result$metrics$roc_auc$status, "incomplete")
  expect_identical(result$metrics$brier_score$status, "incomplete")
  expect_identical(result$metrics$decision_curve$status, "incomplete")
  expect_identical(result$metrics$decision_curve$n_scored, 0L)
})

test_that("binary CV decision curve is reconstructed from OOF patients", {
  d <- make_sb_binary_fixture(n = 60L, seed = 5601L)
  thresholds <- c(0.25, 0.50, 0.75)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c"),
    decision_thresholds = thresholds
  )

  result <- sb_run_cv(
    task,
    sb_pipeline(learner = "logistic"),
    folds = 3L,
    seed = 5701L
  )
  metric <- result$metrics$decision_curve

  expect_identical(result$status, "success")
  expect_identical(metric$status, "success")
  expect_identical(metric$scope, "pooled_patient_level_out_of_fold")
  expect_identical(metric$n_scored, task$n)
  expect_identical(result$task$decision_thresholds, thresholds)
  expect_identical(metric$curve$threshold_probability, thresholds)

  for (index in seq_along(thresholds)) {
    threshold <- thresholds[[index]]
    selected <- result$predictions$probability >= threshold
    true_positive <- as.integer(sum(
      selected & result$predictions$outcome == 1L
    ))
    false_positive <- as.integer(sum(
      selected & result$predictions$outcome == 0L
    ))
    true_negative <- as.integer(sum(
      !selected & result$predictions$outcome == 0L
    ))
    false_negative <- as.integer(sum(
      !selected & result$predictions$outcome == 1L
    ))
    expected_net_benefit <- true_positive / task$n -
      false_positive / task$n * threshold / (1 - threshold)

    expect_identical(metric$curve$true_positive[[index]], true_positive)
    expect_identical(metric$curve$false_positive[[index]], false_positive)
    expect_identical(metric$curve$true_negative[[index]], true_negative)
    expect_identical(metric$curve$false_negative[[index]], false_negative)
    expect_equal(
      metric$curve$model_net_benefit[[index]],
      expected_net_benefit,
      tolerance = 1e-15
    )
  }
})

test_that("binary CV omits decision curve without a threshold declaration", {
  d <- make_sb_binary_fixture(n = 60L, seed = 5801L)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
  )

  result <- sb_run_cv(
    task,
    sb_pipeline(learner = "logistic"),
    folds = 3L,
    seed = 5901L
  )

  expect_false("decision_curve" %in% names(result$metrics))
  expect_false("decision_thresholds" %in% names(result$task))
})

test_that("binary CV reuses an exact patient-keyed outer assignment", {
  d <- make_sb_binary_fixture(n = 60L, seed = 6201L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 2L
  )
  generated <- sb_run_cv(task, pipeline, folds = 5L, seed = 6202L)
  assignment <- setNames(
    generated$predictions$fold,
    generated$predictions$patient_id
  )
  supplied <- assignment[c(
    seq(2L, length(assignment), by = 2L),
    seq(1L, length(assignment), by = 2L)
  )]

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
  expect_identical(reused$selections, generated$selections)
  expect_identical(
    reused$selection_stability,
    generated$selection_stability
  )
  expect_identical(
    reused$runs[setdiff(names(reused$runs), "seed")],
    generated$runs[setdiff(names(generated$runs), "seed")]
  )
  expect_identical(generated$runs$seed, rep(6202L, 5L))
  expect_identical(reused$runs$seed, rep(NA_integer_, 5L))
  expect_identical(generated$provenance$fold_source, "generated_stratified")
  expect_identical(reused$provenance$fold_source, "preassigned_patient_id")
  expect_identical(generated$provenance$seed, 6202L)
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

test_that("preassigned binary folds are not silently repaired", {
  d <- make_sb_binary_fixture(n = 20L, seed = 6211L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  supplied <- setNames(
    ifelse(d$outcome == 1L, 1L, 2L),
    d$patient_id
  )

  result <- sb_run_cv(
    task,
    sb_pipeline(learner = "logistic"),
    folds = supplied,
    seed = 6212L
  )

  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(result$provenance$fold_source, "preassigned_patient_id")
  expect_identical(result$provenance$seed, NA_integer_)
  expect_identical(
    result$provenance$fold_assignment$fold,
    as.integer(supplied)
  )
  expect_identical(
    result$runs$error_class,
    rep("stormbreaker_logistic_single_class", 2L)
  )
  expect_identical(result$runs$status, rep("failed", 2L))
})
