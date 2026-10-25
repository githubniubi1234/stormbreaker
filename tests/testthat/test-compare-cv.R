test_that("comparison uses identical outer folds and copied OOF metrics", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260816)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    ),
    sb_pipeline(
      learner = "glmnet_cox",
      inner_folds = 3L,
      inner_seed = 1701L
    )
  )

  comparison <- sb_compare_cv(
    task,
    pipelines,
    folds = 5L,
    seed = 101L
  )

  expect_s3_class(comparison, "stormbreaker_cv_comparison")
  expect_identical(comparison$status, "success")
  expect_identical(comparison$n_patients, 60L)
  expect_identical(comparison$folds, 5L)
  expect_identical(comparison$seed, 101L)
  required_columns <- c(
    "pipeline_id", "selector", "learner", "status", "complete_oof",
    "harrell_c_status", "harrell_c", "comparable_pairs",
    "selection_stability_status", "mean_pairwise_jaccard",
    "n_successful_selection_folds",
    "n_folds", "n_successful_folds",
    "glmnet_active_set_stability_status",
    "glmnet_active_set_mean_pairwise_jaccard",
    "n_successful_glmnet_active_set_folds"
  )
  expect_length(
    setdiff(required_columns, names(comparison$comparison)),
    0L
  )
  expected_ids <- vapply(pipelines, function(x) x$id, character(1))
  expect_identical(comparison$comparison$pipeline_id, expected_ids)
  expect_identical(names(comparison$results), expected_ids)
  expect_identical(comparison$comparison$complete_oof, rep(TRUE, 3L))
  expect_identical(
    comparison$comparison$harrell_c_status,
    rep("success", 3L)
  )
  expect_identical(
    comparison$comparison$selection_stability_status,
    c("not_applicable", "success", "not_applicable")
  )
  expect_identical(
    is.finite(comparison$comparison$mean_pairwise_jaccard),
    c(FALSE, TRUE, FALSE)
  )
  expect_identical(
    comparison$comparison$n_successful_selection_folds,
    c(0L, 5L, 0L)
  )
  expect_equal(
    comparison$comparison$mean_pairwise_jaccard,
    unname(vapply(
      comparison$results,
      function(result) result$selection_stability$mean_pairwise_jaccard,
      numeric(1)
    )),
    tolerance = 0
  )
  expected_active_status <- unname(vapply(
    comparison$results,
    function(result) result$glmnet_active_set_stability$status,
    character(1)
  ))
  expect_identical(
    comparison$comparison$glmnet_active_set_stability_status,
    expected_active_status
  )
  expect_identical(
    expected_active_status[1:2],
    rep("not_applicable", 2L)
  )
  expect_match(expected_active_status[[3L]], "^(success|not_evaluable)$")
  expect_equal(
    comparison$comparison$glmnet_active_set_mean_pairwise_jaccard,
    unname(vapply(
      comparison$results,
      function(result) {
        result$glmnet_active_set_stability$mean_pairwise_jaccard
      },
      numeric(1)
    )),
    tolerance = 0
  )
  expect_identical(
    comparison$comparison$n_successful_glmnet_active_set_folds,
    c(0L, 0L, 5L)
  )
  expect_identical(
    vapply(comparison$results[[1L]]$coefficients, `[[`, "", "status"),
    rep("success", 5L)
  )
  expect_identical(
    vapply(comparison$results[[2L]]$coefficients, `[[`, "", "status"),
    rep("success", 5L)
  )

  assignments <- lapply(
    comparison$results,
    function(x) x$predictions[c("patient_id", "fold")]
  )
  for (assignment in assignments[-1L]) {
    expect_identical(assignment, assignments[[1L]])
  }
  expected_c <- vapply(
    comparison$results,
    function(x) x$metrics$harrell_c$estimate,
    numeric(1)
  )
  expected_pairs <- vapply(
    comparison$results,
    function(x) x$metrics$harrell_c$comparable_pairs,
    numeric(1)
  )
  expect_equal(
    comparison$comparison$harrell_c,
    unname(expected_c),
    tolerance = 0
  )
  expect_equal(
    comparison$comparison$comparable_pairs,
    unname(expected_pairs),
    tolerance = 0
  )

  printed <- paste(capture.output(print(comparison)), collapse = "\n")
  expect_match(
    printed,
    "No pipeline is selected by this comparison.",
    fixed = TRUE
  )
})

test_that("survival comparison copies declared-horizon IPCW Brier", {
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
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    )
  )

  compared <- sb_compare_cv(
    task,
    pipelines,
    folds = 5L,
    seed = 20260819L
  )

  expect_identical(compared$status, "success")
  expect_identical(
    compared$comparison$complete_survival_probability_oof,
    c(TRUE, TRUE)
  )
  expect_identical(compared$comparison$prediction_horizon, c(7, 7))
  expect_identical(
    compared$comparison$ipcw_brier_status,
    c("success", "success")
  )
  expect_identical(
    compared$comparison$reference_ipcw_brier_status,
    c("success", "success")
  )
  expect_identical(
    compared$comparison$ipcw_brier_skill_status,
    c("success", "success")
  )
  expect_identical(
    compared$comparison$time_auc_status,
    c("success", "success")
  )
  expect_identical(
    compared$comparison$uno_c_status,
    c("success", "success")
  )
  expect_identical(
    compared$comparison$calibration_status,
    c("success", "success")
  )
  expect_identical(
    compared$comparison$mean_calibration_status,
    c("success", "success")
  )
  expect_identical(
    compared$comparison$calibration_slope_status,
    c("success", "success")
  )
  expected <- unname(vapply(
    compared$results,
    function(result) result$metrics$ipcw_brier$estimate,
    numeric(1)
  ))
  expect_equal(compared$comparison$ipcw_brier, expected, tolerance = 0)
  expected_reference <- unname(vapply(
    compared$results,
    function(result) result$metrics$ipcw_brier$reference$estimate,
    numeric(1)
  ))
  expected_skill <- unname(vapply(
    compared$results,
    function(result) result$metrics$ipcw_brier$skill_score$estimate,
    numeric(1)
  ))
  expect_equal(
    compared$comparison$reference_ipcw_brier,
    expected_reference,
    tolerance = 0
  )
  expect_equal(
    compared$comparison$ipcw_brier_skill_score,
    expected_skill,
    tolerance = 0
  )
  expected_time_auc <- unname(vapply(
    compared$results,
    function(result) result$metrics$time_auc$estimate,
    numeric(1)
  ))
  expected_time_auc_denominator <- unname(vapply(
    compared$results,
    function(result) result$metrics$time_auc$pair_denominator,
    numeric(1)
  ))
  expect_equal(
    compared$comparison$time_auc,
    expected_time_auc,
    tolerance = 0
  )
  expect_equal(
    compared$comparison$time_auc_pair_denominator,
    expected_time_auc_denominator,
    tolerance = 0
  )
  expected_uno_c <- unname(vapply(
    compared$results,
    function(result) result$metrics$uno_c$estimate,
    numeric(1)
  ))
  expected_uno_c_denominator <- unname(vapply(
    compared$results,
    function(result) result$metrics$uno_c$weighted_pair_denominator,
    numeric(1)
  ))
  expect_equal(compared$comparison$uno_c, expected_uno_c, tolerance = 0)
  expect_equal(
    compared$comparison$uno_c_weighted_pair_denominator,
    expected_uno_c_denominator,
    tolerance = 0
  )
  expected_calibration_ici <- unname(vapply(
    compared$results,
    function(result) result$metrics$calibration$ici,
    numeric(1)
  ))
  expected_calibration_e50 <- unname(vapply(
    compared$results,
    function(result) result$metrics$calibration$e50,
    numeric(1)
  ))
  expected_calibration_e90 <- unname(vapply(
    compared$results,
    function(result) result$metrics$calibration$e90,
    numeric(1)
  ))
  expect_identical(
    compared$comparison$calibration_ici,
    expected_calibration_ici
  )
  expect_identical(
    compared$comparison$calibration_e50,
    expected_calibration_e50
  )
  expect_identical(
    compared$comparison$calibration_e90,
    expected_calibration_e90
  )
  expected_mean_calibration <- unname(vapply(
    compared$results,
    function(result) result$metrics$weak_calibration$mean_calibration_ratio,
    numeric(1)
  ))
  expected_calibration_slope <- unname(vapply(
    compared$results,
    function(result) result$metrics$weak_calibration$calibration_slope,
    numeric(1)
  ))
  expect_identical(
    compared$comparison$mean_calibration_ratio,
    expected_mean_calibration
  )
  expect_identical(
    compared$comparison$calibration_slope,
    expected_calibration_slope
  )
  first_references <- compared$results[[1L]]$metrics$ipcw_brier$reference
  for (result in compared$results[-1L]) {
    observed_reference <- result$metrics$ipcw_brier$reference
    expect_identical(
      observed_reference$per_fold$reference_survival_horizon,
      first_references$per_fold$reference_survival_horizon
    )
    expect_identical(
      observed_reference$per_fold$weight_sum,
      first_references$per_fold$weight_sum
    )
    expect_identical(
      observed_reference$per_fold$weighted_loss_sum,
      first_references$per_fold$weighted_loss_sum
    )
    expect_identical(
      observed_reference$survival_curves,
      first_references$survival_curves
    )
  }
  printed <- paste(capture.output(print(compared)), collapse = "\n")
  expect_match(printed, "ipcw_brier", fixed = TRUE)
  expect_match(printed, "ipcw_brier_skill_score", fixed = TRUE)
  expect_match(printed, "time_auc", fixed = TRUE)
  expect_match(printed, "calibration_ici", fixed = TRUE)
  expect_match(printed, "calibration_e50", fixed = TRUE)
  expect_match(printed, "calibration_e90", fixed = TRUE)
  expect_match(printed, "mean_calibration_ratio", fixed = TRUE)
  expect_match(printed, "calibration_slope", fixed = TRUE)
})

test_that("survival comparison copies saved calibration without recomputing it", {
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
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    )
  )
  saved_results <- lapply(
    pipelines,
    function(pipeline) {
      sb_run_cv(task, pipeline, folds = 5L, seed = 20260819L)
    }
  )
  next_result <- 0L
  testthat::local_mocked_bindings(
    sb_run_cv = function(...) {
      next_result <<- next_result + 1L
      saved_results[[next_result]]
    },
    sb_survival_calibration_metric = function(...) {
      stop("comparison must not refit survival calibration")
    },
    sb_survival_weak_calibration_metric = function(...) {
      stop("comparison must not refit survival weak calibration")
    }
  )

  compared <- sb_compare_cv(
    task,
    pipelines,
    folds = 5L,
    seed = 20260819L
  )

  expect_identical(compared$status, "success")
  expect_identical(next_result, 2L)
  expect_identical(
    compared$comparison$pipeline_id,
    vapply(pipelines, `[[`, "", "id")
  )
  expect_false("winner" %in% names(compared$comparison))
  expect_identical(
    compared$comparison$calibration_ici,
    unname(vapply(saved_results, function(x) x$metrics$calibration$ici, 0))
  )
  expect_identical(
    compared$comparison$mean_calibration_status,
    unname(vapply(
      saved_results,
      function(x) x$metrics$weak_calibration$mean_calibration$status,
      ""
    ))
  )
  expect_identical(
    compared$comparison$mean_calibration_ratio,
    unname(vapply(
      saved_results,
      function(x) x$metrics$weak_calibration$mean_calibration_ratio,
      0
    ))
  )
  expect_identical(
    compared$comparison$calibration_slope_status,
    unname(vapply(
      saved_results,
      function(x) x$metrics$weak_calibration$slope_fit$status,
      ""
    ))
  )
  expect_identical(
    compared$comparison$calibration_slope,
    unname(vapply(
      saved_results,
      function(x) x$metrics$weak_calibration$calibration_slope,
      0
    ))
  )
})

test_that("all survival learners share the fold-local KM reference", {
  skip_if_not_installed("glmnet")

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
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    ),
    sb_pipeline(
      learner = "glmnet_cox",
      inner_folds = 3L,
      inner_seed = 1701L
    )
  )

  compared <- sb_compare_cv(
    task,
    pipelines,
    folds = 5L,
    seed = 20260819L
  )
  reference <- compared$results[[1L]]$metrics$ipcw_brier$reference

  expect_identical(compared$status, "success")
  for (result in compared$results[-1L]) {
    observed <- result$metrics$ipcw_brier$reference
    expect_identical(observed$per_fold, reference$per_fold)
    expect_identical(observed$survival_curves, reference$survival_curves)
    expect_identical(observed$weighted_loss_sum, reference$weighted_loss_sum)
    expect_identical(observed$estimate, reference$estimate)
  }
})

test_that("unsupported survival horizon stays incomplete without metric fallback", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 1e6
  )

  compared <- sb_compare_cv(
    task,
    list(sb_pipeline()),
    folds = 5L,
    seed = 20260819L
  )

  expect_identical(compared$status, "incomplete")
  expect_identical(compared$comparison$status, "failed")
  expect_identical(compared$comparison$complete_oof, TRUE)
  expect_identical(
    compared$comparison$complete_survival_probability_oof,
    FALSE
  )
  expect_identical(compared$comparison$harrell_c_status, "success")
  expect_identical(compared$comparison$ipcw_brier_status, "incomplete")
  expect_identical(compared$comparison$time_auc_status, "incomplete")
  expect_identical(compared$comparison$calibration_status, "incomplete")
  expect_identical(
    compared$comparison$mean_calibration_status,
    "not_evaluated"
  )
  expect_identical(
    compared$comparison$calibration_slope_status,
    "not_fitted"
  )
  expect_identical(
    compared$comparison$reference_ipcw_brier_status,
    "incomplete"
  )
  expect_identical(
    compared$comparison$ipcw_brier_skill_status,
    "incomplete"
  )
  expect_true(is.na(compared$comparison$ipcw_brier))
  expect_true(is.na(compared$comparison$reference_ipcw_brier))
  expect_true(is.na(compared$comparison$ipcw_brier_skill_score))
  expect_true(is.na(compared$comparison$time_auc))
  expect_true(is.na(compared$comparison$calibration_ici))
  expect_true(is.na(compared$comparison$calibration_e50))
  expect_true(is.na(compared$comparison$calibration_e90))
  expect_true(is.na(compared$comparison$mean_calibration_ratio))
  expect_true(is.na(compared$comparison$calibration_slope))
  expect_true(is.finite(compared$comparison$harrell_c))
})

test_that("comparison copies unavailable Uno evidence without ranking", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
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

  compared <- sb_compare_cv(
    task,
    list(sb_pipeline()),
    folds = 5L,
    seed = 20260819L
  )

  expect_identical(compared$status, "incomplete")
  expect_identical(compared$comparison$status, "failed")
  expect_identical(compared$comparison$harrell_c_status, "success")
  expect_identical(compared$comparison$ipcw_brier_status, "success")
  expect_identical(compared$comparison$time_auc_status, "success")
  expect_identical(compared$comparison$uno_c_status, "incomplete")
  expect_true(is.na(compared$comparison$uno_c))
  expect_false("winner" %in% names(compared$comparison))
})

test_that("comparison keeps a failed pipeline visible without fallback", {
  d <- make_sb_fixture(n = 60, seed = 20260816)
  d$gene_b <- d$gene_a
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month"
  )
  pipelines <- list(
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 1L
    ),
    sb_pipeline()
  )

  comparison <- sb_compare_cv(task, pipelines, folds = 5L, seed = 101L)

  expect_identical(comparison$status, "incomplete")
  expect_identical(
    comparison$comparison$pipeline_id,
    c("univ_cox_top_k__coxph__k1", "none__coxph")
  )
  expect_identical(comparison$comparison$status, c("success", "failed"))
  expect_identical(comparison$comparison$complete_oof, c(TRUE, FALSE))
  expect_identical(
    comparison$comparison$harrell_c_status,
    c("success", "incomplete")
  )
  expect_identical(is.finite(comparison$comparison$harrell_c), c(TRUE, FALSE))
  failed_result <- comparison$results[["none__coxph"]]
  expect_identical(all(is.na(failed_result$runs$learner_used)), TRUE)
  expect_identical(all(is.na(failed_result$predictions$risk)), TRUE)
})

test_that("comparison admits different alphas on identical outer folds", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260816)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipelines <- list(
    sb_pipeline(
      learner = "glmnet_cox",
      glmnet_alpha = 0.3,
      inner_folds = 3L,
      inner_seed = 101L
    ),
    sb_pipeline(
      learner = "glmnet_cox",
      glmnet_alpha = 0.7,
      inner_folds = 3L,
      inner_seed = 101L
    )
  )

  comparison <- sb_compare_cv(task, pipelines, folds = 5L, seed = 101L)

  expect_identical(comparison$status, "success")
  expect_identical(
    comparison$comparison$pipeline_id,
    c(
      "none__glmnet_cox__alpha0.3__inner3__seed101",
      "none__glmnet_cox__alpha0.7__inner3__seed101"
    )
  )
  expect_length(comparison$results, 2L)
  assignments <- lapply(
    comparison$results,
    function(x) x$predictions[c("patient_id", "fold")]
  )
  expect_identical(assignments[[1L]], assignments[[2L]])
})

test_that("comparison requires a usable collection of unique pipelines", {
  d <- make_sb_fixture(n = 20)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )

  empty_error <- capture_sb_error(
    sb_compare_cv(task, list(), folds = 2L, seed = 1L)
  )
  expect_s3_class(empty_error, "stormbreaker_bad_pipeline_collection")

  member_error <- capture_sb_error(
    sb_compare_cv(task, list("not a pipeline"), folds = 2L, seed = 1L)
  )
  expect_s3_class(member_error, "stormbreaker_bad_pipeline_collection")

  duplicate_error <- capture_sb_error(
    sb_compare_cv(
      task,
      list(sb_pipeline(), sb_pipeline()),
      folds = 2L,
      seed = 1L
    )
  )
  expect_s3_class(duplicate_error, "stormbreaker_duplicate_pipeline_id")
})

test_that("binary comparison uses identical outer folds and saved AUC", {
  skip_if_not_installed("glmnet")

  d <- make_sb_binary_fixture(n = 60)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes"),
    sb_pipeline(
      learner = "glmnet_logistic",
      glmnet_alpha = 0,
      inner_folds = 3L,
      inner_seed = 101L
    ),
    sb_pipeline(
      learner = "glmnet_logistic",
      glmnet_alpha = 1,
      inner_folds = 3L,
      inner_seed = 101L
    )
  )
  comparison <- sb_compare_cv(
    task, pipelines, folds = 3L, seed = 20260819L
  )

  expect_identical(comparison$status, "success")
  expect_identical(comparison$task_type, "binary")
  expect_identical(
    comparison$comparison$pipeline_id,
    vapply(pipelines, `[[`, "", "id")
  )
  expect_named(
    comparison$comparison,
    c(
      "pipeline_id", "selector", "learner", "status", "complete_oof",
      "roc_auc_status", "roc_auc", "comparable_pairs",
      "average_precision_status", "average_precision", "brier_status",
      "brier_score", "n_scored", "reference_brier_status",
      "reference_brier", "brier_skill_status", "brier_skill_score",
      "log_loss_status", "log_loss", "reference_log_loss_status",
      "reference_log_loss", "log_loss_skill_status", "log_loss_skill_score",
      "calibration_status", "calibration_intercept", "calibration_slope",
      "moderate_calibration_status", "moderate_calibration_ici",
      "moderate_calibration_e50", "moderate_calibration_e90",
      "selection_stability_status", "mean_pairwise_jaccard",
      "n_successful_selection_folds",
      "n_folds", "n_successful_folds",
      "glmnet_active_set_stability_status",
      "glmnet_active_set_mean_pairwise_jaccard",
      "n_successful_glmnet_active_set_folds", "prediction_variation"
    ),
    ignore.order = FALSE
  )
  assignments <- lapply(
    comparison$results,
    function(x) x$predictions[c("patient_id", "fold")]
  )
  expect_identical(assignments[[1L]], assignments[[2L]])
  expect_identical(assignments[[1L]], assignments[[3L]])
  expect_identical(assignments[[1L]], assignments[[4L]])
  expected_auc <- vapply(
    comparison$results,
    function(x) x$metrics$roc_auc$estimate,
    numeric(1)
  )
  expect_equal(
    comparison$comparison$roc_auc,
    unname(expected_auc),
    tolerance = 0
  )
  expected_average_precision <- vapply(
    comparison$results,
    function(x) x$metrics$average_precision$estimate,
    numeric(1)
  )
  expect_equal(
    comparison$comparison$average_precision,
    unname(expected_average_precision),
    tolerance = 0
  )
  expected_brier <- vapply(
    comparison$results,
    function(x) x$metrics$brier_score$estimate,
    numeric(1)
  )
  manual_brier <- vapply(
    comparison$results,
    function(x) {
      mean((x$predictions$outcome - x$predictions$probability)^2)
    },
    numeric(1)
  )
  expect_equal(expected_brier, manual_brier, tolerance = 0)
  expect_equal(
    comparison$comparison$brier_score,
    unname(expected_brier),
    tolerance = 0
  )
  expected_reference <- vapply(
    comparison$results,
    function(x) x$metrics$brier_score$reference$estimate,
    numeric(1)
  )
  expected_skill <- vapply(
    comparison$results,
    function(x) x$metrics$brier_score$skill_score$estimate,
    numeric(1)
  )
  expect_equal(
    comparison$comparison$reference_brier,
    unname(expected_reference),
    tolerance = 0
  )
  expect_length(unique(comparison$comparison$reference_brier), 1L)
  expect_equal(
    comparison$comparison$brier_skill_score,
    unname(expected_skill),
    tolerance = 0
  )
  expected_log_loss <- vapply(
    comparison$results,
    function(x) x$metrics$log_loss$estimate,
    numeric(1)
  )
  expected_reference_log_loss <- vapply(
    comparison$results,
    function(x) x$metrics$log_loss$reference$estimate,
    numeric(1)
  )
  expected_log_loss_skill <- vapply(
    comparison$results,
    function(x) x$metrics$log_loss$skill_score$estimate,
    numeric(1)
  )
  expect_equal(
    comparison$comparison$log_loss,
    unname(expected_log_loss),
    tolerance = 0
  )
  expect_equal(
    comparison$comparison$reference_log_loss,
    unname(expected_reference_log_loss),
    tolerance = 0
  )
  expect_equal(
    comparison$comparison$log_loss_skill_score,
    unname(expected_log_loss_skill),
    tolerance = 0
  )
  expected_calibration_intercept <- vapply(
    comparison$results,
    function(x) x$metrics$calibration$calibration_intercept,
    numeric(1)
  )
  expected_calibration_slope <- vapply(
    comparison$results,
    function(x) x$metrics$calibration$calibration_slope,
    numeric(1)
  )
  expect_equal(
    comparison$comparison$calibration_intercept,
    unname(expected_calibration_intercept),
    tolerance = 0
  )
  expect_equal(
    comparison$comparison$calibration_slope,
    unname(expected_calibration_slope),
    tolerance = 0
  )
  expected_moderate_status <- vapply(
    comparison$results,
    function(x) x$metrics$moderate_calibration$status,
    character(1)
  )
  expect_identical(
    comparison$comparison$moderate_calibration_status,
    unname(expected_moderate_status)
  )
  for (metric in c("ici", "e50", "e90")) {
    expected <- unname(vapply(
      comparison$results,
      function(x) x$metrics$moderate_calibration[[metric]],
      numeric(1)
    ))
    expect_equal(
      comparison$comparison[[paste0("moderate_calibration_", metric)]],
      expected,
      tolerance = 0
    )
  }
  expect_identical(
    comparison$comparison$selection_stability_status,
    rep("not_applicable", length(pipelines))
  )
  expect_true(all(is.na(comparison$comparison$mean_pairwise_jaccard)))
  expect_identical(
    comparison$comparison$n_successful_selection_folds,
    rep(0L, length(pipelines))
  )
  expected_active_status <- unname(vapply(
    comparison$results,
    function(result) result$glmnet_active_set_stability$status,
    character(1)
  ))
  expect_identical(
    comparison$comparison$glmnet_active_set_stability_status,
    expected_active_status
  )
  expect_identical(
    expected_active_status[1:3],
    rep("not_applicable", 3L)
  )
  expect_match(expected_active_status[[4L]], "^(success|not_evaluable)$")
  expect_equal(
    comparison$comparison$glmnet_active_set_mean_pairwise_jaccard,
    unname(vapply(
      comparison$results,
      function(result) {
        result$glmnet_active_set_stability$mean_pairwise_jaccard
      },
      numeric(1)
    )),
    tolerance = 0
  )
  expect_identical(
    comparison$comparison$n_successful_glmnet_active_set_folds,
    c(0L, 0L, 0L, 3L)
  )
  printed <- paste(capture.output(print(comparison)), collapse = "\n")
  expect_match(printed, "average_precision", fixed = TRUE)
  expect_match(printed, "brier_score", fixed = TRUE)
  expect_match(printed, "brier_skill_score", fixed = TRUE)
  expect_match(printed, "log_loss", fixed = TRUE)
  expect_match(printed, "log_loss_skill_score", fixed = TRUE)
  expect_match(printed, "calibration_intercept", fixed = TRUE)
  expect_match(printed, "calibration_slope", fixed = TRUE)
  expect_match(printed, "moderate_calibration_ici", fixed = TRUE)
  expect_match(printed, "moderate_calibration_e50", fixed = TRUE)
  expect_match(printed, "moderate_calibration_e90", fixed = TRUE)
  expect_match(printed, "No pipeline is selected", fixed = TRUE)
})

test_that("binary comparison copies saved moderate calibration without refitting", {
  d <- make_sb_binary_fixture(n = 60L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  )
  saved_results <- lapply(
    pipelines,
    function(pipeline) {
      sb_run_cv(task, pipeline, folds = 5L, seed = 20260820L)
    }
  )
  saved_results[[1L]]$metrics$moderate_calibration$ici <- 0.123456789
  next_result <- 0L
  testthat::local_mocked_bindings(
    sb_run_cv = function(...) {
      next_result <<- next_result + 1L
      saved_results[[next_result]]
    },
    sb_binary_moderate_calibration_metric = function(...) {
      stop("comparison must not refit binary moderate calibration")
    }
  )

  compared <- sb_compare_cv(
    task,
    pipelines,
    folds = 5L,
    seed = 20260820L
  )

  expect_identical(next_result, 2L)
  expect_identical(
    compared$comparison$moderate_calibration_ici[[1L]],
    0.123456789
  )
  expect_identical(
    compared$comparison$moderate_calibration_e50,
    unname(vapply(
      saved_results,
      function(x) x$metrics$moderate_calibration$e50,
      numeric(1)
    ))
  )
  expect_false("winner" %in% names(compared$comparison))
})

test_that("binary comparison continues after an earlier failed logistic", {
  skip_if_not_installed("glmnet")

  d <- make_sb_binary_fixture(n = 60)
  d$gene_b <- d$gene_a
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b")
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(
      learner = "glmnet_logistic",
      inner_folds = 3L,
      inner_seed = 1L
    )
  )
  comparison <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1L)

  expect_identical(comparison$status, "incomplete")
  expect_identical(
    names(comparison$results),
    vapply(pipelines, `[[`, character(1L), "id")
  )
  expect_identical(comparison$comparison$status, c("failed", "success"))
  expect_identical(
    comparison$comparison$roc_auc_status,
    c("incomplete", "success")
  )
  expect_identical(
    comparison$comparison$average_precision_status,
    c("incomplete", "success")
  )
  expect_identical(
    comparison$comparison$brier_status,
    c("incomplete", "success")
  )
  expect_identical(
    comparison$comparison$reference_brier_status,
    c("success", "success")
  )
  expect_identical(
    comparison$comparison$brier_skill_status,
    c("incomplete", "success")
  )
  expect_identical(
    comparison$comparison$log_loss_status,
    c("incomplete", "success")
  )
  expect_identical(
    comparison$comparison$reference_log_loss_status,
    c("success", "success")
  )
  expect_identical(
    comparison$comparison$log_loss_skill_status,
    c("incomplete", "success")
  )
  expect_identical(
    comparison$comparison$calibration_status,
    c("incomplete", "success")
  )
  expect_identical(
    comparison$comparison$moderate_calibration_status,
    c("incomplete", "success")
  )
  failed <- comparison$results[[pipelines[[1L]]$id]]
  succeeded <- comparison$results[[pipelines[[2L]]$id]]
  expect_identical(all(is.na(failed$runs$learner_used)), TRUE)
  expect_identical(all(is.na(failed$predictions$probability)), TRUE)
  expect_identical(succeeded$status, "success")
  expect_true(all(is.finite(succeeded$predictions$probability)))
})

test_that("binary comparison admits fold-local AUC selection without ranking", {
  d <- make_sb_high_dim_binary_fixture(n = 60L, p = 8L, seed = 1701L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "logistic",
      selector_n_features = 2L
    ),
    sb_pipeline(learner = "gaussian_naive_bayes"),
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "gaussian_naive_bayes",
      selector_n_features = 2L
    )
  )

  comparison <- sb_compare_cv(
    task, pipelines, folds = 3L, seed = 1801L
  )

  expect_identical(comparison$status, "success")
  expect_identical(
    comparison$comparison$pipeline_id,
    vapply(pipelines, `[[`, "", "id")
  )
  fold_assignments <- lapply(
    comparison$results,
    function(result) result$predictions[c("patient_id", "fold")]
  )
  for (index in 2:length(fold_assignments)) {
    expect_identical(fold_assignments[[1L]], fold_assignments[[index]])
  }
  expect_equal(
    comparison$comparison$roc_auc,
    unname(vapply(
      comparison$results,
      function(result) result$metrics$roc_auc$estimate,
      numeric(1)
    )),
    tolerance = 0
  )
  expect_equal(
    comparison$comparison$brier_score,
    unname(vapply(
      comparison$results,
      function(result) result$metrics$brier_score$estimate,
      numeric(1)
    )),
    tolerance = 0
  )
  expect_length(comparison$results[[2L]]$selections, 3L)
  expect_length(comparison$results[[4L]]$selections, 3L)
  expect_identical(
    comparison$comparison$selection_stability_status,
    c("not_applicable", "success", "not_applicable", "success")
  )
  expect_identical(
    is.finite(comparison$comparison$mean_pairwise_jaccard),
    c(FALSE, TRUE, FALSE, TRUE)
  )
  expect_identical(
    comparison$comparison$n_successful_selection_folds,
    c(0L, 3L, 0L, 3L)
  )
  expect_equal(
    comparison$comparison$mean_pairwise_jaccard,
    unname(vapply(
      comparison$results,
      function(result) result$selection_stability$mean_pairwise_jaccard,
      numeric(1)
    )),
    tolerance = 0
  )
  selected_audit_matches <- Map(
    function(audit, selection) {
      length(unique(audit$values$feature)) == 2L &&
        setequal(unique(audit$values$feature), selection$selected_features)
    },
    comparison$results[[4L]]$distribution_parameters,
    comparison$results[[4L]]$selections
  )
  expect_true(all(unlist(selected_audit_matches, use.names = FALSE)))
  printed <- paste(capture.output(print(comparison)), collapse = "\n")
  expect_match(printed, "No pipeline is selected", fixed = TRUE)
})

test_that("fixed LightGBM comparison reuses folds with and without Top-K", {
  skip_if_not_installed("lightgbm", minimum_version = "4.0.0")

  d <- make_sb_high_dim_binary_fixture(n = 48L, p = 10L, seed = 3901L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  fixed_args <- list(
    lightgbm_nrounds = 6L,
    lightgbm_learning_rate = 0.1,
    lightgbm_num_leaves = 7L,
    lightgbm_min_data_in_leaf = 3L,
    lightgbm_seed = 4001L
  )
  pipelines <- list(
    do.call(
      sb_pipeline,
      c(list(learner = "lightgbm_binary"), fixed_args)
    ),
    do.call(
      sb_pipeline,
      c(
        list(
          selector = "univ_binary_auc_top_k",
          learner = "lightgbm_binary",
          selector_n_features = 3L
        ),
        fixed_args
      )
    )
  )

  comparison <- sb_compare_cv(
    task, pipelines, folds = 3L, seed = 4101L
  )

  expect_identical(comparison$status, "success")
  expect_identical(
    comparison$comparison$pipeline_id,
    vapply(pipelines, `[[`, "", "id")
  )
  fold_assignments <- lapply(
    comparison$results,
    function(result) result$predictions[c("patient_id", "fold")]
  )
  expect_identical(fold_assignments[[1L]], fold_assignments[[2L]])
  expect_length(comparison$results[[1L]]$tree_models, 3L)
  expect_length(comparison$results[[2L]]$tree_models, 3L)
  expect_length(comparison$results[[2L]]$selections, 3L)
  for (fold in seq_len(3L)) {
    expect_setequal(
      comparison$results[[2L]]$tree_models[[fold]]$features,
      comparison$results[[2L]]$selections[[fold]]$selected_features
    )
  }
  printed <- paste(capture.output(print(comparison)), collapse = "\n")
  expect_match(printed, "No pipeline is selected", fixed = TRUE)
})

test_that("fixed CatBoost comparison reuses folds with and without Top-K", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_high_dim_binary_fixture(n = 48L, p = 10L, seed = 5301L)
  predictors <- grep("^gene_", names(d), value = TRUE)
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  fixed_args <- list(
    catboost_iterations = 6L,
    catboost_learning_rate = 0.1,
    catboost_depth = 3L,
    catboost_l2_leaf_reg = 3,
    catboost_seed = 5401L
  )
  pipelines <- list(
    do.call(
      sb_pipeline,
      c(list(learner = "catboost_binary"), fixed_args)
    ),
    do.call(
      sb_pipeline,
      c(
        list(
          selector = "univ_binary_auc_top_k",
          learner = "catboost_binary",
          selector_n_features = 3L
        ),
        fixed_args
      )
    )
  )

  comparison <- sb_compare_cv(
    task, pipelines, folds = 3L, seed = 5501L
  )

  expect_identical(comparison$status, "success")
  expect_identical(
    comparison$comparison$pipeline_id,
    vapply(pipelines, `[[`, "", "id")
  )
  fold_assignments <- lapply(
    comparison$results,
    function(result) result$predictions[c("patient_id", "fold")]
  )
  expect_identical(fold_assignments[[1L]], fold_assignments[[2L]])
  expect_length(comparison$results[[1L]]$tree_models, 3L)
  expect_length(comparison$results[[2L]]$tree_models, 3L)
  expect_length(comparison$results[[2L]]$selections, 3L)
  for (fold in seq_len(3L)) {
    expect_setequal(
      comparison$results[[2L]]$tree_models[[fold]]$features,
      comparison$results[[2L]]$selections[[fold]]$selected_features
    )
  }
  printed <- paste(capture.output(print(comparison)), collapse = "\n")
  expect_match(printed, "No pipeline is selected", fixed = TRUE)
})

test_that("binary comparison copies declared decision curves without selection", {
  d <- make_sb_binary_fixture(n = 60L, seed = 5601L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c"),
    decision_thresholds = c(0.2, 0.5, 0.8)
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  )

  compared <- sb_compare_cv(
    task,
    pipelines,
    folds = 3L,
    seed = 5701L
  )

  expect_identical(
    names(compared$decision_curves),
    vapply(pipelines, `[[`, character(1), "id")
  )
  for (index in seq_along(pipelines)) {
    expect_identical(
      compared$decision_curves[[index]],
      compared$results[[index]]$metrics$decision_curve
    )
  }
  expect_false(any(grepl(
    "best|winner|selected_threshold",
    names(compared),
    ignore.case = TRUE
  )))
  expect_false(any(grepl(
    "decision",
    names(compared$comparison),
    ignore.case = TRUE
  )))
})

test_that("binary comparison omits decision curves when undeclared", {
  d <- make_sb_binary_fixture(n = 60L, seed = 5801L)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
  )
  compared <- sb_compare_cv(
    task,
    list(
      sb_pipeline(learner = "logistic"),
      sb_pipeline(learner = "gaussian_naive_bayes")
    ),
    folds = 3L,
    seed = 5901L
  )

  expect_false("decision_curves" %in% names(compared))
})

test_that("survival comparison copies declared curves without selection", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7,
    decision_thresholds = c(0.20, 0.40, 0.60)
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    )
  )
  compared <- sb_compare_cv(
    task,
    pipelines,
    folds = 5L,
    seed = 20260819L
  )

  expect_identical(compared$status, "success")
  expect_named(compared$decision_curves, vapply(pipelines, `[[`, "", "id"))
  for (index in seq_along(pipelines)) {
    expect_identical(
      compared$decision_curves[[index]],
      compared$results[[index]]$metrics$decision_curve
    )
  }
  expect_false(any(grepl(
    "decision|winner|best|threshold",
    names(compared$comparison),
    ignore.case = TRUE
  )))
})

test_that("survival comparison reuses one patient-keyed outer assignment", {
  d <- make_sb_fixture(n = 30L, seed = 6301L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    )
  )
  generated <- sb_make_stratified_folds(task, 3L, 6302L)
  supplied <- setNames(generated$assignment, d$patient_id)
  supplied <- supplied[c(30L, seq_len(29L))]

  first <- sb_compare_cv(task, pipelines, folds = supplied, seed = 1L)
  second <- sb_compare_cv(task, pipelines, folds = supplied, seed = 999L)

  expect_identical(first, second)
  expect_identical(first$fold_source, "preassigned_patient_id")
  expect_identical(first$seed, NA_integer_)
  expect_identical(first$folds, 3L)
  expect_identical(
    first$fold_assignment,
    data.frame(
      patient_id = d$patient_id,
      fold = as.integer(supplied[match(d$patient_id, names(supplied))]),
      stringsAsFactors = FALSE
    )
  )
  for (index in seq_along(pipelines)) {
    direct <- sb_run_cv(
      task,
      pipelines[[index]],
      folds = supplied,
      seed = 123L
    )
    expect_identical(first$results[[index]], direct)
    expect_identical(
      first$results[[index]]$provenance$fold_assignment,
      first$fold_assignment
    )
  }
  printed <- paste(capture.output(print(first)), collapse = "\n")
  expect_match(printed, "preassigned patient-ID folds", fixed = TRUE)
})

test_that("binary comparison reuses one patient-keyed outer assignment", {
  d <- make_sb_binary_fixture(n = 30L, seed = 6311L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  )
  generated <- sb_make_stratified_folds(task, 3L, 6312L)
  supplied <- setNames(generated$assignment, d$patient_id)
  supplied <- supplied[rev(seq_along(supplied))]

  first <- sb_compare_cv(task, pipelines, folds = supplied, seed = 1L)
  second <- sb_compare_cv(task, pipelines, folds = supplied, seed = 999L)

  expect_identical(first, second)
  expect_identical(first$fold_source, "preassigned_patient_id")
  expect_identical(first$seed, NA_integer_)
  expect_identical(first$folds, 3L)
  expect_identical(
    first$fold_assignment,
    data.frame(
      patient_id = d$patient_id,
      fold = as.integer(supplied[match(d$patient_id, names(supplied))]),
      stringsAsFactors = FALSE
    )
  )
  for (index in seq_along(pipelines)) {
    direct <- sb_run_cv(
      task,
      pipelines[[index]],
      folds = supplied,
      seed = 123L
    )
    expect_identical(first$results[[index]], direct)
    expect_identical(
      first$results[[index]]$provenance$fold_assignment,
      first$fold_assignment
    )
  }
  printed <- paste(capture.output(print(first)), collapse = "\n")
  expect_match(printed, "preassigned patient-ID folds", fixed = TRUE)
})
