registry_test_arguments <- function(row) {
  arguments <- list(
    selector = row$selector,
    learner = row$learner
  )
  if (identical(row$selector_parameter, "selector_n_features")) {
    arguments$selector_n_features <- 2L
  }
  if (identical(row$selector_parameter, "stepcox_direction")) {
    arguments$stepcox_direction <- "forward"
  }
  if (row$selector_parameter %in% c(
    "selector_rsf_nodesize_and_seed",
    "selector_rsf_nodesize_seed_and_stepcox_direction"
  )) {
    arguments$selector_rsf_nodesize <- 5L
    arguments$selector_seed <- 101L
  }
  if (row$selector_parameter %in% c(
    "selector_inner_folds_and_seed",
    "selector_inner_folds_seed_and_stepcox_direction"
  )) {
    arguments$selector_inner_folds <- 3L
    arguments$selector_seed <- 101L
  }
  if (row$selector_parameter %in% c(
    "selector_rsf_nodesize_seed_and_stepcox_direction",
    "selector_inner_folds_seed_and_stepcox_direction"
  )) {
    arguments$stepcox_direction <- "forward"
  }
  if (identical(
    row$selector_parameter,
    "selector_features_and_feature_set_id"
  )) {
    arguments$selector_features <- c("gene_a", "gene_b")
    arguments$feature_set_id <- "registry_pair"
  }
  if (identical(
    row$tuning_semantics,
    "lambda_selected_by_inner_cv"
  )) {
    arguments$inner_folds <- 3L
    arguments$inner_seed <- 101L
  }
  if (identical(row$tuning_semantics, "analysis_fold_internal_cv")) {
    arguments$learner_tuning_folds <- 3L
    arguments$learner_seed <- 101L
    if (identical(row$learner, "gbm_cox_cv")) {
      arguments$gbm_cv_max_trees <- 60L
    }
  }
  arguments
}

removed_python_learner_ids <- function() {
  c(
    "coxtime_fixed", "deephit_fixed", "deepsurv_fixed",
    "logistic_hazard_fixed", "pc_hazard_fixed"
  )
}

test_that("the public registry lists only admitted executable combinations", {
  registry <- sb_pipeline_registry()

  expect_s3_class(registry, "data.frame")
  expect_identical(
    names(registry),
    c(
      "source", "task_type", "selector", "learner",
      "admission_status", "selector_parameter", "tuning_semantics",
      "optional_dependency", "parameter_semantics"
    )
  )
  expect_setequal(
    unique(registry$source),
    c(
      "stormbreaker_builtin", "stormbreaker_abess_runtime",
      "stormbreaker_forest_selector_runtime",
      "stormbreaker_high_dimensional_runtime",
      "stormbreaker_binary_model_selector_runtime",
      "stormbreaker_survival_boosting_selector_runtime",
      "stormbreaker_pamr_selector_runtime",
      "stormbreaker_survival_weight_selector_runtime",
      "stormbreaker_mime117_runtime"
    )
  )
  expect_identical(unique(registry$admission_status), "admitted")
  expect_false(anyDuplicated(registry[c("selector", "learner")]) > 0L)
  expect_identical(
    nrow(sb_pipeline_registry(level = "learner")),
    176L
  )
  expect_false(any(registry$learner %in% removed_python_learner_ids()))
  for (learner in removed_python_learner_ids()) {
    expect_error(
      sb_pipeline(learner = learner),
      class = "stormbreaker_unsupported_learner"
    )
  }

  observed <- paste(
    registry$task_type,
    registry$selector,
    registry$learner,
    sep = "::"
  )
  survival_full_learners <- c(
    "coxph", "glmnet_cox", "catboost_cox",
    "catboost_aft_survival", "gbm_cox",
    sb_mboost_cox_learner_ids(),
    sb_mboost_aft_learner_ids(),
    sb_mboost_ranking_learner_ids(),
    sb_ranger_survival_learner_ids(), sb_rlt_survival_learner_ids(),
    sb_coxphf_learner_ids(), sb_bagged_knn_survival_learner_ids(),
    sb_ahaz_survival_learner_ids(), sb_bujar_linear_learner_ids(), sb_akritas_learner_ids(),
    "coxboost", "plsrcox",
    "superpc", sb_survivalsvm_learner_ids(),
    sb_xgboost_cox_learner_ids(),
    sb_xgboost_aft_learner_ids(),
    sb_survreg_learner_ids(),
    sb_flexsurv_parametric_learner_ids(),
    sb_flexsurvspline_learner_ids(),
    sb_mtlr_survival_learner_ids(),
    sb_pamr_survival_learner_ids(), "rpart_survival",
    "rsf_survival",
    "ctree_survival", "cforest_survival", "bart_survival",
    "aorsf_survival", sb_timereg_survival_learner_ids(),
    sb_ncvreg_survival_learner_ids(),
    "abess_cox_gic",
    sb_analysis_tuned_survival_learner_ids(),
    "bart_aft_survival"
  )
  survival_restricted_learners <- c(
    "pre_survival_rule_ensemble",
    "penaft_gehan_elastic_net"
  )
  survival_learners <- c(
    survival_full_learners,
    survival_restricted_learners
  )
  binary_full_learners <- c(
    "logistic", "glmnet_logistic", "gaussian_naive_bayes",
    sb_omics_specialist_binary_learner_ids(),
    sb_lightgbm_binary_learner_ids(),
    "catboost_binary", sb_ranger_binary_learner_ids(),
    sb_rlt_binary_learner_ids(), sb_wsrf_binary_learner_ids(),
    sb_random_glm_binary_learner_ids(),
    "adaboost_binary", "logitboost_binary",
    sb_classical_binary_learner_ids(),
    sb_specialized_binary_linear_learner_ids(),
    sb_xgboost_binary_learner_ids(), sb_ncvreg_binary_learner_ids(),
    "abess_binary_gic",
    sb_mclust_da_learner_ids(), sb_mixomics_da_learner_ids(),
    sb_hdda_binary_learner_ids()
  )
  binary_restricted_learners <- c(
    "pre_binary_rule_ensemble",
    sb_slope_binary_learner_ids(),
    "cancerclass_binary"
  )
  binary_learners <- c(binary_full_learners, binary_restricted_learners)
  paired_keys <- function(task, selectors, learners) {
    unlist(lapply(learners, function(learner) {
      paste(task, selectors, learner, sep = "::")
    }), use.names = FALSE)
  }
  expected <- c(
    paired_keys(
      "survival",
      c("none", "univ_cox_top_k"),
      survival_full_learners
    ),
    paste(
      "survival",
      "none",
      survival_restricted_learners,
      sep = "::"
    ),
    paired_keys(
      "binary",
      c("none", "univ_binary_auc_top_k"),
      binary_full_learners
    ),
    paste("binary", "none", binary_restricted_learners, sep = "::"),
    paste(
      "survival",
      "stepcox",
      survival_full_learners,
      sep = "::"
    ),
    paste(
      "survival",
      "rsf_minimal_depth",
      c(
        "coxboost_cv", "glmnet_cox", "gbm_cox_cv",
        "plsrcox_cv", "superpc_cv", sb_survivalsvm_learner_ids()
      ),
      sep = "::"
    ),
    paste(
      "survival",
      "coxboost_nonzero",
      c(
        "glmnet_cox", "gbm_cox_cv", "plsrcox_cv",
        "superpc_cv", sb_survivalsvm_learner_ids()
      ),
      sep = "::"
    ),
    paste(
      "survival",
      "lasso_nonzero",
      c(
        "coxboost_cv", "gbm_cox_cv", "plsrcox_cv",
        "rsf_survival", "superpc_cv", sb_survivalsvm_learner_ids()
      ),
      sep = "::"
    ),
    paste(
      "survival",
      c(
        "rsf_then_stepcox", "coxboost_then_stepcox",
        "lasso_then_stepcox"
      ),
      "coxph",
      sep = "::"
    ),
    paste(
      "binary",
      "abess_binary_gic_nonzero",
      setdiff(binary_full_learners, "abess_binary_gic"),
      sep = "::"
    ),
    paste(
      "survival",
      "abess_cox_gic_nonzero",
      setdiff(
        survival_full_learners,
        c("abess_cox_gic", sb_analysis_tuned_survival_learner_ids())
      ),
      sep = "::"
    ),
    paired_keys(
      "binary",
      c(
        "boruta_binary_roughfix", "varselrf_binary_oob_1se",
        "lasso_binary_nonzero", "stepglm_binary_aic",
        "glmboost_binary_nonzero", "rf_binary_importance",
        "enet_binary_nonzero", "gbm_binary_importance", "svm_rfe_binary", "gamboost_binary_nonzero", "xgboost_binary_importance", "plsrglm_binary_weights",
        setdiff(
          sb_high_dimensional_binary_selector_ids(),
          "stabsel_binary_lasso"
        )
      ),
      binary_full_learners
    ),
    paste(
      "binary",
      "stabsel_binary_lasso",
      setdiff(binary_full_learners, "glmnet_logistic"),
      sep = "::"
    ),
    paste(
      "survival",
      "boruta_survival_roughfix",
      setdiff(
        survival_full_learners,
        sb_analysis_tuned_survival_learner_ids()
      ),
      sep = "::"
    ),
    paste(
      "survival",
      "lpc_cox_nlogn",
      setdiff(
        survival_full_learners,
        sb_analysis_tuned_survival_learner_ids()
      ),
      sep = "::"
    ),
    paste(
      "survival",
      "stabsel_cox_lasso",
      setdiff(
        survival_full_learners,
        c(sb_analysis_tuned_survival_learner_ids(), "glmnet_cox")
      ),
      sep = "::"
    ),
    paste("survival", "fixed_features", survival_learners, sep = "::"),
    paste("binary", "fixed_features", binary_learners, sep = "::")
  )
  expected <- c(expected,
    paired_keys("binary", "pamr_binary_nonzero", binary_full_learners),
    paired_keys("survival", "pamr_survival_nonzero",
      setdiff(survival_full_learners, sb_analysis_tuned_survival_learner_ids())),
    paired_keys("survival", c("lasso_nonzero", "enet_cox_nonzero", "coxboost_nonzero", "rsf_minimal_depth"),
      setdiff(survival_full_learners, sb_analysis_tuned_survival_learner_ids())),
    paired_keys("survival", c("ridge_cox_weights", "plsrcox_weights", "svm_survival_weights"),
      setdiff(survival_full_learners, sb_analysis_tuned_survival_learner_ids())),
    paired_keys("survival", c("gbm_cox_importance", "glmboost_cox_nonzero", "superpc_survival_scores", "gamboost_cox_nonzero", "xgboost_cox_importance"),
      setdiff(survival_full_learners, sb_analysis_tuned_survival_learner_ids())),
    paired_keys("survival", "enet_cox_nonzero",
      unique(c("coxph", "glmnet_cox", "ranger_survival", "coxboost_cv",
        "gbm_cox_cv", "plsrcox_cv", "rsf_survival", "superpc_cv", sb_survivalsvm_learner_ids()))))
  expect_setequal(observed, expected)
})

test_that("every registry row resolves through the existing runtime", {
  registry <- sb_pipeline_registry()
  ids <- character(nrow(registry))

  for (index in seq_len(nrow(registry))) {
    row <- registry[index, , drop = FALSE]
    pipeline <- do.call(sb_pipeline, registry_test_arguments(row))
    runtime <- sb_resolve_pipeline_runtime(pipeline)

    expect_identical(pipeline$task_type, row$task_type)
    expect_identical(pipeline$selector, row$selector)
    expect_identical(pipeline$learner, row$learner)
    expect_identical(
      sb_validate_pipeline_contract(pipeline, row$task_type),
      pipeline
    )
    expect_identical(
      sb_validate_pipeline_runtime(runtime, pipeline),
      runtime
    )
    ids[[index]] <- pipeline$id
  }

  expect_identical(anyDuplicated(ids), 0L)
})

test_that("the registry can summarize one row per learner", {
  routes <- sb_pipeline_registry()
  learners <- sb_pipeline_registry(level = "learner")
  route_counts <- table(routes$learner)
  expected_selectors <- vapply(
    learners$learner,
    function(learner) {
      paste(routes$selector[routes$learner == learner], collapse = ", ")
    },
    character(1)
  )

  expect_identical(
    names(learners),
    c(
      "task_type", "learner", "pipeline_count", "selectors",
      "tuning_semantics", "optional_dependency", "parameter_semantics"
    )
  )
  expect_setequal(learners$learner, unique(routes$learner))
  expect_identical(
    learners$pipeline_count,
    unname(as.integer(route_counts[learners$learner]))
  )
  expect_identical(learners$selectors, unname(expected_selectors))

  binary <- sb_pipeline_registry(task_type = "binary", level = "learner")
  expect_identical(unique(binary$task_type), "binary")
  expect_setequal(
    binary$learner,
    unique(routes$learner[routes$task_type == "binary"])
  )
})

test_that("learner summaries can report installed R dependencies", {
  declared <- sb_pipeline_registry(level = "learner")
  checked <- sb_pipeline_registry(
    level = "learner",
    check_installed = TRUE
  )

  expect_identical(checked[names(declared)], declared)
  expect_identical(
    tail(names(checked), 4L),
    c(
      "r_packages", "python_modules", "missing_r_packages",
      "dependency_status"
    )
  )
  expect_identical(
    checked$dependency_status[checked$learner == "coxph"],
    "no_optional_dependency"
  )
  expect_identical(
    checked$dependency_status[
      checked$learner == "survreg_weibull"
    ],
    "r_packages_installed"
  )
  expect_identical(
    any(checked$learner %in% removed_python_learner_ids()),
    FALSE
  )
  expect_true(
    all(is.na(checked$python_modules) | checked$python_modules == "")
  )
  expect_identical(
    any(checked$dependency_status == "python_not_checked"),
    FALSE
  )
  expect_identical(
    setdiff(
      unique(checked$dependency_status),
      c(
        "no_optional_dependency", "r_packages_installed",
        "missing_r_packages", "python_not_checked"
      )
    ),
    character()
  )
})

test_that("fixed nonlinear learners are not printed as elastic nets", {
  binary <- utils::capture.output(print(sb_pipeline(
    learner = "bart_binary"
  )))
  survival <- utils::capture.output(print(sb_pipeline(
    learner = "bart_survival"
  )))

  expect_true(any(grepl("fixed before resampling", binary, fixed = TRUE)))
  expect_true(any(grepl("fixed before resampling", survival, fixed = TRUE)))
  expect_false(any(grepl("Elastic Net", c(binary, survival), fixed = TRUE)))
})

test_that("the linter admits a raw specification or constructed pipeline", {
  raw_report <- sb_lint_pipeline(list(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 2L
  ))
  object_report <- sb_lint_pipeline(
    sb_pipeline(learner = "logistic"),
    task_type = "binary"
  )

  expect_identical(
    names(raw_report),
    c(
      "status", "pipeline_id", "task_type", "selector", "learner",
      "issue_code", "message"
    )
  )
  expect_identical(raw_report$status, "admitted")
  expect_identical(
    raw_report$pipeline_id,
    "univ_binary_auc_top_k__logistic__k2"
  )
  expect_identical(raw_report$task_type, "binary")
  expect_true(is.na(raw_report$issue_code))
  expect_true(is.na(raw_report$message))

  expect_identical(object_report$status, "admitted")
  expect_identical(object_report$pipeline_id, "none__logistic")
})

test_that("the linter reports expected specification and contract errors", {
  missing_inner <- sb_lint_pipeline(list(learner = "glmnet_cox"))
  mismatch <- sb_lint_pipeline(
    list(learner = "logistic"),
    task_type = "survival"
  )
  unknown_argument <- sb_lint_pipeline(list(
    learner = "logistic",
    mystery = 1
  ))
  bad_container <- sb_lint_pipeline("none__coxph")
  bad_learner <- sb_lint_pipeline(list(learner = 123))
  bad_selector <- sb_lint_pipeline(list(
    selector = c("none", "bad")
  ))
  mutated <- sb_pipeline()
  mutated$risk_formula <- "changed"
  mutation <- sb_lint_pipeline(mutated)

  expect_identical(missing_inner$status, "rejected")
  expect_identical(
    missing_inner$issue_code,
    "stormbreaker_missing_inner_parameter"
  )
  expect_identical(
    mismatch$issue_code,
    "stormbreaker_pipeline_task_mismatch"
  )
  expect_identical(
    unknown_argument$issue_code,
    "stormbreaker_unknown_pipeline_argument"
  )
  expect_identical(
    bad_container$issue_code,
    "stormbreaker_bad_pipeline_spec"
  )
  expect_true(is.na(bad_container$selector))
  expect_true(is.na(bad_container$learner))
  expect_identical(bad_learner$selector, "none")
  expect_true(is.na(bad_learner$learner))
  expect_true(is.na(bad_selector$selector))
  expect_identical(bad_selector$learner, "coxph")
  expect_identical(
    mutation$issue_code,
    "stormbreaker_unresolved_pipeline"
  )
  expect_true(all(nzchar(c(
    missing_inner$message,
    mismatch$message,
    unknown_argument$message,
    bad_container$message,
    mutation$message
  ))))
})

test_that("runtime registry and adapter drift remains an internal error", {
  local_mocked_bindings(
    sb_resolve_pipeline_runtime = function(pipeline) {
      sb_abort(
        "injected registry drift",
        "stormbreaker_unresolved_learner_adapter"
      )
    }
  )

  error <- tryCatch(
    sb_lint_pipeline(list(learner = "logistic")),
    error = function(error) error
  )

  expect_s3_class(error, "stormbreaker_unresolved_learner_adapter")
  expect_identical(conditionMessage(error), "injected registry drift")
})

test_that("unexpected contract stormbreaker errors remain visible", {
  local_mocked_bindings(
    sb_validate_pipeline_contract = function(pipeline, task_type) {
      sb_abort(
        "injected internal contract failure",
        "stormbreaker_internal_contract_failure"
      )
    }
  )

  error <- tryCatch(
    sb_lint_pipeline(list(learner = "logistic")),
    error = function(error) error
  )

  expect_s3_class(error, "stormbreaker_internal_contract_failure")
  expect_identical(
    conditionMessage(error),
    "injected internal contract failure"
  )
})

test_that("unexpected constructor stormbreaker errors remain visible", {
  pipeline <- sb_pipeline(learner = "logistic")
  local_mocked_bindings(
    sb_pipeline = function(...) {
      sb_abort(
        "injected internal constructor failure",
        "stormbreaker_internal_constructor_failure"
      )
    }
  )

  error <- tryCatch(
    sb_lint_pipeline(pipeline),
    error = function(error) error
  )

  expect_s3_class(error, "stormbreaker_internal_constructor_failure")
  expect_identical(
    conditionMessage(error),
    "injected internal constructor failure"
  )
})

test_that("mutated pipeline reports use exact fields without defaults", {
  missing <- sb_pipeline(learner = "logistic")
  missing[["learner"]] <- NULL
  missing_report <- sb_lint_pipeline(missing)

  partial <- sb_pipeline(learner = "logistic")
  partial[["learner"]] <- NULL
  partial[["learners"]] <- "spoofed"
  partial_report <- sb_lint_pipeline(partial)

  expect_identical(missing_report$status, "rejected")
  expect_identical(
    missing_report$issue_code,
    "stormbreaker_unresolved_pipeline"
  )
  expect_identical(missing_report$selector, "none")
  expect_true(is.na(missing_report$learner))
  expect_identical(partial_report$status, "rejected")
  expect_identical(
    partial_report$issue_code,
    "stormbreaker_unresolved_pipeline"
  )
  expect_identical(partial_report$selector, "none")
  expect_true(is.na(partial_report$learner))
})

test_that("admission does not claim that optional packages are installed", {
  local_mocked_bindings(sb_lightgbm_available = function() FALSE)

  report <- sb_lint_pipeline(list(learner = "lightgbm_binary"))
  registry_row <- sb_pipeline_registry()
  registry_row <- registry_row[
    registry_row$learner == "lightgbm_binary" &
      registry_row$selector == "none", ,
    drop = FALSE
  ]

  expect_identical(report$status, "admitted")
  expect_identical(registry_row$optional_dependency, "lightgbm")
  expect_identical(registry_row$admission_status, "admitted")
})

test_that("the linter does not disguise unexpected internal errors", {
  local_mocked_bindings(
    sb_resolve_pipeline_runtime = function(pipeline) {
      stop("unexpected internal failure")
    }
  )

  error <- tryCatch(
    sb_lint_pipeline(list(learner = "logistic")),
    error = function(error) error
  )

  expect_s3_class(error, "simpleError")
  expect_false(inherits(error, "stormbreaker_error"))
  expect_identical(conditionMessage(error), "unexpected internal failure")

  pipeline <- sb_pipeline()
  local_mocked_bindings(
    sb_pipeline = function(...) {
      stop("unexpected reconstruction failure")
    }
  )
  reconstruction_error <- tryCatch(
    sb_lint_pipeline(pipeline),
    error = function(error) error
  )

  expect_s3_class(reconstruction_error, "simpleError")
  expect_false(inherits(reconstruction_error, "stormbreaker_error"))
  expect_identical(
    conditionMessage(reconstruction_error),
    "unexpected reconstruction failure"
  )
})
