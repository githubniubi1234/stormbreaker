test_that("IECV preserves group roles and rejects patients shared by cohorts", {
  tasks <- lapply(c("a", "b"), function(cohort) {
    d <- sb_example_binary(40L, 3754L)
    d$patient_id <- paste0(cohort, d$patient_id)
    group_col <- paste0(cohort, "_group")
    d[[group_col]] <- paste0(cohort, rep(1:20, each = 2))
    sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"),
                   group_col = group_col)
  })
  names(tasks) <- c("a", "b")
  validated <- sb_validate_iecv_cohorts(tasks)
  pooled <- sb_iecv_development_task(validated$tasks, tasks[[1]], "binary")
  expect_identical(pooled$group_col, "a_group")
  expect_identical(pooled$data$a_group,
    c(tasks$a$data$a_group, tasks$b$data$b_group))
  tasks$b$data$b_group[1:2] <- tasks$a$data$a_group[1]
  error <- capture_sb_error(sb_validate_iecv_cohorts(tasks))
  expect_s3_class(error, "stormbreaker_internal_external_patient_overlap")
})

make_sb_iecv_binary_task <- function(cohort_id, seed, predictors,
                                     id_col = "patient_id",
                                     outcome_col = "outcome") {
  data <- make_sb_binary_fixture(n = 36L, seed = seed)
  data$patient_id <- paste0(cohort_id, "_", data$patient_id)
  names(data)[names(data) == "patient_id"] <- id_col
  names(data)[names(data) == "outcome"] <- outcome_col
  sb_binary_task(
    data,
    id_col,
    outcome_col,
    predictors,
    decision_thresholds = c(0.2, 0.4)
  )
}

make_sb_iecv_survival_task <- function(cohort_id, seed, predictors,
                                       horizon = 6) {
  data <- make_sb_fixture(n = 36L, seed = seed)
  data$patient_id <- paste0(cohort_id, "_", data$patient_id)
  sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = horizon,
    decision_thresholds = c(0.2, 0.4)
  )
}

test_that("binary IECV holds out each cohort from fitting", {
  cohorts <- list(
    cohort_a = make_sb_iecv_binary_task(
      "a", 9101L, c("gene_a", "gene_b", "gene_c")
    ),
    cohort_b = make_sb_iecv_binary_task(
      "b", 9102L, c("gene_c", "gene_a", "gene_b"),
      id_col = "subject_id", outcome_col = "label"
    ),
    cohort_c = make_sb_iecv_binary_task(
      "c", 9103L, c("gene_b", "gene_c", "gene_a")
    )
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  )

  observed <- sb_run_internal_external_cv(cohorts, pipelines)
  pipeline_ids <- vapply(pipelines, `[[`, character(1L), "id")

  expect_s3_class(observed, "stormbreaker_internal_external_cv")
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(
    observed$validation_scheme,
    "leave_one_cohort_out_internal_external_cv"
  )
  expect_identical(observed$pipeline_ids, pipeline_ids)
  expect_identical(names(observed$folds), names(cohorts))
  expect_identical(nrow(observed$results), 6L)
  expect_identical(
    observed$results$heldout_cohort_id,
    rep(names(cohorts), each = 2L)
  )
  expect_identical(observed$results$pipeline_id, rep(pipeline_ids, 3L))
  expect_identical(
    observed$results$n_development_patients,
    rep(72L, 6L)
  )

  all_ids <- lapply(cohorts, function(task) {
    as.character(task$data[[task$id_col]])
  })
  for (heldout_id in names(cohorts)) {
    fold <- observed$folds[[heldout_id]]
    development_ids <- unlist(
      all_ids[fold$development_cohort_ids],
      use.names = FALSE
    )
    heldout_ids <- all_ids[[heldout_id]]

    expect_identical(fold$heldout_cohort_id, heldout_id)
    expect_setequal(fold$development_patient_ids, development_ids)
    expect_length(intersect(fold$development_patient_ids, heldout_ids), 0L)
    expect_identical(fold$heldout_patient_ids, heldout_ids)
    expect_identical(fold$comparison$status, "success")
    expect_identical(names(fold$fits), pipeline_ids)
    expect_identical(names(fold$predictions), pipeline_ids)
    for (pipeline_id in pipeline_ids) {
      expect_identical(
        fold$fits[[pipeline_id]]$development$patient_ids,
        fold$development_patient_ids
      )
      expect_identical(
        fold$predictions[[pipeline_id]]$predictions$patient_id,
        heldout_ids
      )
    }
  }

  expect_false(any(c(
    "pooled_metric", "ranking", "winner", "final_fit"
  ) %in% names(observed)))
  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "different fitted model", fixed = TRUE)
  expect_match(
    printed,
    "not final independent external validation",
    fixed = TRUE
  )
})

test_that("held-out binary outcomes cannot change their fold predictions", {
  cohorts <- list(
    cohort_a = make_sb_iecv_binary_task(
      "a", 9201L, c("gene_a", "gene_b")
    ),
    cohort_b = make_sb_iecv_binary_task(
      "b", 9202L, c("gene_b", "gene_a")
    ),
    cohort_c = make_sb_iecv_binary_task(
      "c", 9203L, c("gene_a", "gene_b")
    )
  )
  pipeline <- sb_pipeline(learner = "logistic")
  original <- sb_run_internal_external_cv(cohorts, list(pipeline))

  changed_data <- cohorts$cohort_c$data
  changed_data$outcome <- 1L - changed_data$outcome
  changed <- cohorts
  changed$cohort_c <- sb_binary_task(
    changed_data,
    "patient_id",
    "outcome",
    cohorts$cohort_c$predictors,
    decision_thresholds = cohorts$cohort_c$decision_thresholds
  )
  altered <- sb_run_internal_external_cv(changed, list(pipeline))

  original_fold <- original$folds$cohort_c
  altered_fold <- altered$folds$cohort_c
  pipeline_id <- pipeline$id
  expect_identical(
    original_fold$development_patient_ids,
    altered_fold$development_patient_ids
  )
  expect_identical(
    original_fold$fits[[pipeline_id]]$learner_state,
    altered_fold$fits[[pipeline_id]]$learner_state
  )
  expect_equal(
    original_fold$predictions[[pipeline_id]]$predictions$probability,
    altered_fold$predictions[[pipeline_id]]$predictions$probability,
    tolerance = 0
  )
})

test_that("binary stack IECV fits only the non-held-out cohorts", {
  skip_if_not_installed("glmnet")

  cohorts <- list(
    cohort_a = make_sb_iecv_binary_task(
      "a", 9251L, c("gene_a", "gene_b")
    ),
    cohort_b = make_sb_iecv_binary_task(
      "b", 9252L, c("gene_b", "gene_a")
    ),
    cohort_c = make_sb_iecv_binary_task(
      "c", 9253L, c("gene_a", "gene_b")
    )
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 9254L)

  observed <- sb_run_internal_external_cv(cohorts, stack)

  expect_identical(observed$status, "success")
  expect_identical(observed$strategy_type, "stack")
  expect_identical(observed$pipeline_ids, stack$id)
  expect_identical(observed$results$pipeline_id, rep(stack$id, 3L))
  for (heldout_id in names(cohorts)) {
    fold <- observed$folds[[heldout_id]]
    heldout_ids <- as.character(
      cohorts[[heldout_id]]$data[[cohorts[[heldout_id]]$id_col]]
    )
    development_ids <- unlist(lapply(
      cohorts[names(cohorts) != heldout_id],
      function(task) as.character(task$data[[task$id_col]])
    ), use.names = FALSE)
    fit <- fold$fits[[stack$id]]
    prediction <- fold$predictions[[stack$id]]

    expect_identical(names(fold$fits), stack$id)
    expect_identical(names(fold$predictions), stack$id)
    expect_setequal(fit$development$patient_ids, development_ids)
    expect_length(intersect(fit$development$patient_ids, heldout_ids), 0L)
    expect_identical(names(fit$member_fits), stack$pipeline_ids)
    expect_identical(
      prediction$predictions$patient_id,
      heldout_ids
    )
    expect_identical(fold$comparison$comparison$pipeline_id, stack$id)
  }

  changed_data <- cohorts$cohort_c$data
  changed_data$outcome <- 1L - changed_data$outcome
  changed <- cohorts
  changed$cohort_c <- sb_binary_task(
    changed_data,
    cohorts$cohort_c$id_col,
    cohorts$cohort_c$outcome_col,
    cohorts$cohort_c$predictors,
    decision_thresholds = cohorts$cohort_c$decision_thresholds
  )
  altered <- sb_run_internal_external_cv(changed, stack)
  original_fold <- observed$folds$cohort_c
  altered_fold <- altered$folds$cohort_c

  expect_identical(
    original_fold$fits[[stack$id]]$meta$state,
    altered_fold$fits[[stack$id]]$meta$state
  )
  expect_identical(
    lapply(
      original_fold$fits[[stack$id]]$member_fits,
      `[[`,
      "learner_state"
    ),
    lapply(
      altered_fold$fits[[stack$id]]$member_fits,
      `[[`,
      "learner_state"
    )
  )
  expect_equal(
    original_fold$predictions[[stack$id]]$predictions$probability,
    altered_fold$predictions[[stack$id]]$predictions$probability,
    tolerance = 0
  )
})

test_that("survival IECV retains one declared time target", {
  cohorts <- list(
    cohort_a = make_sb_iecv_survival_task(
      "a", 9301L, c("gene_a", "gene_b", "gene_c")
    ),
    cohort_b = make_sb_iecv_survival_task(
      "b", 9302L, c("gene_c", "gene_b", "gene_a")
    ),
    cohort_c = make_sb_iecv_survival_task(
      "c", 9303L, c("gene_b", "gene_a", "gene_c")
    )
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "coxph",
    selector_n_features = 2L
  )

  observed <- sb_run_internal_external_cv(cohorts, list(pipeline))

  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "survival")
  expect_identical(observed$cohorts$time_unit, rep("month", 3L))
  expect_identical(observed$cohorts$prediction_horizon, rep(6, 3L))
  expect_identical(
    c("ipcw_brier", "time_auc", "uno_c") %in% names(observed$results),
    rep(TRUE, 3L)
  )
  for (fold in observed$folds) {
    expect_identical(fold$comparison$prediction_target$has_horizon, TRUE)
    expect_identical(
      fold$comparison$prediction_target$prediction_horizon,
      6
    )
    expect_identical(
      fold$comparison$prediction_target$time_unit,
      "month"
    )
    expect_identical(
      fold$predictions[[pipeline$id]]$decision_thresholds,
      c(0.2, 0.4)
    )
  }
})

test_that("survival stack IECV keeps its baseline inside development cohorts", {
  skip_if_not_installed("glmnet")

  cohorts <- list(
    cohort_a = make_sb_iecv_survival_task(
      "sa", 9321L, c("gene_a", "gene_b", "gene_c")
    ),
    cohort_b = make_sb_iecv_survival_task(
      "sb", 9322L, c("gene_c", "gene_a", "gene_b")
    ),
    cohort_c = make_sb_iecv_survival_task(
      "sc", 9323L, c("gene_b", "gene_c", "gene_a")
    )
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "coxph"),
    sb_pipeline(learner = "survreg_weibull")
  ), inner_folds = 3L, inner_seed = 9324L)

  observed <- sb_run_internal_external_cv(cohorts, stack)

  expect_identical(observed$status, "success")
  expect_identical(observed$strategy_type, "stack")
  for (heldout_id in names(cohorts)) {
    fold <- observed$folds[[heldout_id]]
    fit <- fold$fits[[stack$id]]
    prediction <- fold$predictions[[stack$id]]
    heldout_ids <- as.character(
      cohorts[[heldout_id]]$data[[cohorts[[heldout_id]]$id_col]]
    )

    expect_length(intersect(fit$development$patient_ids, heldout_ids), 0L)
    expect_identical(
      fit$meta$state$baseline_hazard$prediction_horizon,
      6
    )
    expect_identical(
      prediction$predictions$prediction_horizon,
      rep(6, length(heldout_ids))
    )
    expect_false(any(!is.finite(prediction$predictions$risk)))
    expect_false(any(!is.finite(
      prediction$predictions$survival_probability
    )))
  }
})

test_that("risk-only survival IECV does not invent a horizon", {
  make_risk_task <- function(cohort_id, seed) {
    data <- make_sb_fixture(n = 30L, seed = seed)
    data$patient_id <- paste0(cohort_id, "_", data$patient_id)
    sb_surv_task(
      data,
      "patient_id",
      "time",
      "event",
      c("gene_a", "gene_b"),
      "month"
    )
  }
  cohorts <- list(
    cohort_a = make_risk_task("a", 9351L),
    cohort_b = make_risk_task("b", 9352L)
  )

  observed <- sb_run_internal_external_cv(
    cohorts,
    list(sb_pipeline())
  )

  expect_identical(observed$status, "success")
  expect_identical(
    "prediction_horizon" %in% names(observed$cohorts),
    FALSE
  )
  expect_identical(
    c("ipcw_brier", "time_auc", "uno_c") %in% names(observed$results),
    rep(FALSE, 3L)
  )
  for (fold in observed$folds) {
    expect_identical(fold$comparison$prediction_target$has_horizon, FALSE)
    expect_identical(
      fold$comparison$prediction_target$prediction_horizon,
      NULL
    )
  }
})

test_that("IECV validates cohort isolation and common targets", {
  binary_a <- make_sb_iecv_binary_task(
    "a", 9401L, c("gene_a", "gene_b")
  )
  binary_b <- make_sb_iecv_binary_task(
    "b", 9402L, c("gene_a", "gene_b")
  )
  different_features <- make_sb_iecv_binary_task(
    "c", 9403L, c("gene_a", "gene_c")
  )
  different_thresholds <- binary_b
  different_thresholds$decision_thresholds <- c(0.3, 0.5)
  different_thresholds <- sb_binary_task(
    different_thresholds$data,
    different_thresholds$id_col,
    different_thresholds$outcome_col,
    different_thresholds$predictors,
    different_thresholds$decision_thresholds
  )
  overlapping <- binary_b
  overlapping$data[[overlapping$id_col]] <- binary_a$data[[binary_a$id_col]]
  overlapping <- sb_binary_task(
    overlapping$data,
    overlapping$id_col,
    overlapping$outcome_col,
    overlapping$predictors,
    overlapping$decision_thresholds
  )
  survival <- make_sb_iecv_survival_task(
    "s", 9404L, c("gene_a", "gene_b")
  )
  pipeline <- sb_pipeline(learner = "logistic")

  expect_s3_class(
    capture_sb_error(sb_run_internal_external_cv(
      list(binary_a), list(pipeline)
    )),
    "stormbreaker_bad_internal_external_cohorts"
  )
  expect_s3_class(
    capture_sb_error(sb_run_internal_external_cv(
      list(a = binary_a, s = survival), list(pipeline)
    )),
    "stormbreaker_internal_external_task_mismatch"
  )
  expect_s3_class(
    capture_sb_error(sb_run_internal_external_cv(
      list(a = binary_a, c = different_features), list(pipeline)
    )),
    "stormbreaker_internal_external_predictor_mismatch"
  )
  expect_s3_class(
    capture_sb_error(sb_run_internal_external_cv(
      list(a = binary_a, b = different_thresholds), list(pipeline)
    )),
    "stormbreaker_internal_external_target_mismatch"
  )
  expect_s3_class(
    capture_sb_error(sb_run_internal_external_cv(
      list(a = binary_a, b = overlapping), list(pipeline)
    )),
    "stormbreaker_internal_external_patient_overlap"
  )
  expect_s3_class(
    capture_sb_error(sb_run_internal_external_cv(
      list(a = binary_a, b = binary_b), list(pipeline, pipeline)
    )),
    "stormbreaker_duplicate_pipeline_id"
  )
})

test_that("IECV continues later pipelines after each earlier fold-fit failure", {
  make_high_dimensional_task <- function(cohort_id, seed) {
    data <- make_sb_fixture(n = 14L, seed = seed)
    data$patient_id <- paste0(cohort_id, "_", data$patient_id)
    withr::with_seed(seed + 1000L, {
      for (index in 4:31) {
        data[[sprintf("gene_%02d", index)]] <- stats::rnorm(nrow(data))
      }
    })
    predictors <- c(
      "gene_a", "gene_b", "gene_c", sprintf("gene_%02d", 4:31)
    )
    sb_surv_task(
      data, "patient_id", "time", "event", predictors, "month"
    )
  }
  cohorts <- list(
    cohort_a = make_high_dimensional_task("a", 9501L),
    cohort_b = make_high_dimensional_task("b", 9502L),
    cohort_c = make_high_dimensional_task("c", 9503L)
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      learner = "coxph",
      selector_n_features = 2L
    )
  )

  observed <- sb_run_internal_external_cv(cohorts, pipelines)
  pipeline_ids <- vapply(pipelines, `[[`, character(1L), "id")

  expect_identical(observed$status, "incomplete")
  expect_identical(
    observed$results$status,
    rep(c("failed", "success"), 3L)
  )
  for (fold in observed$folds) {
    expect_identical(fold$fits[[pipeline_ids[[1L]]]]$status, "failed")
    expect_identical(
      fold$predictions[[pipeline_ids[[1L]]]]$status,
      "failed"
    )
    expect_identical(fold$fits[[pipeline_ids[[2L]]]]$status, "success")
    expect_identical(
      fold$predictions[[pipeline_ids[[2L]]]]$status,
      "success"
    )
  }
})

test_that("IECV keeps a failed stack intact without member fallback", {
  data_a <- make_sb_iecv_binary_task(
    "fa", 9551L, c("gene_a", "gene_b")
  )
  data_b <- make_sb_iecv_binary_task(
    "fb", 9552L, c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 9553L)
  local_mocked_bindings(sb_glmnet_available = function() FALSE)

  observed <- sb_run_internal_external_cv(
    list(cohort_a = data_a, cohort_b = data_b),
    stack
  )

  expect_identical(observed$status, "incomplete")
  expect_identical(observed$results$status, rep("failed", 2L))
  for (fold in observed$folds) {
    expect_identical(fold$fits[[stack$id]]$status, "failed")
    expect_identical(fold$fits[[stack$id]]$error_stage, "meta_fit")
    expect_length(fold$fits[[stack$id]]$member_fits, 0L)
    expect_identical(fold$predictions[[stack$id]]$status, "failed")
    expect_identical(
      fold$predictions[[stack$id]]$error_stage,
      "fit"
    )
  }
})
