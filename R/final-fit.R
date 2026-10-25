sb_final_fit_attempt <- function(code) {
  tryCatch(
    list(value = force(code), error = NULL),
    error = function(error) list(value = NULL, error = error)
  )
}

sb_as_final_fit <- function(result) {
  class(result) <- "stormbreaker_final_fit"
  result
}

sb_final_fit_failure <- function(result, stage, error) {
  result$status <- "failed"
  result$error_stage <- stage
  result$error_class <- class(error)[1L]
  result$error_message <- conditionMessage(error)
  if (identical(stage, "selection")) {
    result$selection$status <- "failed"
    result$selection$error_class <- class(error)[1L]
    result$selection$error_message <- conditionMessage(error)
  }
  if (identical(stage, "tuning")) {
    result$tuning$status <- "failed"
    result$tuning$error_class <- class(error)[1L]
    result$tuning$error_message <- conditionMessage(error)
    if (is.character(error$glmnet_warnings)) {
      result$tuning$warnings <- error$glmnet_warnings
    }
  }
  if (identical(stage, "learner_fit")) {
    result$learner$status <- "failed"
    result$learner$error_class <- class(error)[1L]
    result$learner$error_message <- conditionMessage(error)
  }
  sb_as_final_fit(result)
}

sb_empty_final_tuning <- function(pipeline) {
  is_glmnet <- pipeline$learner %in% c(
    "glmnet_cox", "glmnet_logistic"
  )
  params <- if (is_glmnet) pipeline$learner_params else list()
  result <- list(
    learner = pipeline$learner,
    status = if (is_glmnet) "not_fitted" else "not_required",
    lambda = NA_real_,
    alpha = if (is_glmnet) params$alpha else NA_real_,
    selection_rule = if (is_glmnet) {
      params$selection_rule
    } else {
      NA_character_
    },
    measure = if (is_glmnet) params$measure else NA_character_,
    inner_folds = if (is_glmnet) params$inner_folds else NA_integer_,
    inner_seed = if (is_glmnet) params$inner_seed else NA_integer_,
    fold_id = integer(),
    alignment = if (is_glmnet) "lambda" else NA_character_,
    inner_selection_method = if (is_glmnet) {
      sb_pipeline_inner_selection_method(pipeline)
    } else {
      NA_character_
    },
    inner_selections = list(),
    curve = NULL,
    warnings = character(),
    warning_validation = NULL,
    error_class = NA_character_,
    error_message = NA_character_,
    standardize = if (is_glmnet) params$standardize else NA
  )
  if (pipeline$learner %in% sb_svm_probability_learner_ids()) {
    result$status <- "not_fitted"
    result$inner_folds <- 5L
    result$inner_seed <- pipeline$learner_params$seed
    result$selection_rule <- "libsvm_platt"
    result$measure <- "smoothed_binomial_likelihood"
    result$alignment <- "decision_value"
    result$inner_selection_method <- sb_pipeline_inner_selection_method(pipeline)
    result$standardize <- pipeline$learner_params$standardize
  }
  result
}

sb_successful_final_tuning <- function(pipeline, tuned) {
  result <- list(
    learner = pipeline$learner,
    status = "success",
    lambda = tuned$lambda,
    alpha = tuned$alpha,
    selection_rule = tuned$selection_rule,
    measure = tuned$measure,
    inner_folds = tuned$inner_folds,
    inner_seed = tuned$seed,
    fold_id = tuned$fold_id,
    alignment = tuned$alignment,
    inner_selection_method = tuned$inner_selection_method,
    inner_selections = tuned$inner_selections,
    curve = tuned$curve,
    warnings = tuned$warnings,
    warning_validation = if ("warning_validation" %in% names(tuned)) {
      tuned$warning_validation
    } else {
      NULL
    },
    error_class = NA_character_,
    error_message = NA_character_,
    standardize = tuned$standardize
  )
  if (!is.null(tuned$calibration)) result$calibration <- tuned$calibration
  result
}

sb_final_fit_provenance <- function(pipeline) {
  provenance <- list(
    package = "stormbreaker",
    package_version = as.character(utils::packageVersion("stormbreaker")),
    r_version = R.version.string,
    survival_version = as.character(utils::packageVersion("survival"))
  )
  if (pipeline$learner %in% c("glmnet_cox", "glmnet_logistic")) {
    provenance$glmnet_version <- if (sb_glmnet_available()) {
      as.character(utils::packageVersion("glmnet"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_lightgbm_binary_learner_ids()) {
    provenance$lightgbm_version <- if (sb_lightgbm_available()) {
      as.character(utils::packageVersion("lightgbm"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% c(
    "catboost_binary",
    "catboost_cox",
    "catboost_aft_survival"
  )) {
    provenance$catboost_version <- if (sb_catboost_available()) {
      as.character(utils::packageVersion("catboost"))
    } else {
      NA_character_
    }
  }
  if (identical(pipeline$learner, "cancerclass_binary")) {
    provenance$cancerclass_version <- if (sb_cancerclass_available()) {
      as.character(utils::packageVersion("cancerclass"))
    } else {
      NA_character_
    }
    provenance$biobase_version <- if (sb_cancerclass_available()) {
      as.character(utils::packageVersion("Biobase"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% c("gbm_cox", "gbm_cox_cv")) {
    provenance$gbm_version <- if (sb_gbm_available()) {
      as.character(utils::packageVersion("gbm"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_ranger_learner_ids()) {
    provenance$ranger_version <- if (sb_ranger_available()) {
      as.character(utils::packageVersion("ranger"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in%
    sb_classical_binary_boosting_learner_ids()) {
    dependency <- sb_classical_binary_boosting_dependency(pipeline$learner)
    provenance[[dependency$field]] <- if (dependency$available) {
      as.character(utils::packageVersion(dependency$package))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_ncvreg_learner_ids()) {
    provenance$ncvreg_version <- if (sb_ncvreg_available()) {
      as.character(utils::packageVersion("ncvreg"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_rlt_learner_ids()) {
    provenance$rlt_version <- if (sb_rlt_available()) {
      as.character(utils::packageVersion("RLT"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_coxphf_learner_ids()) {
    provenance$coxphf_version <- if (sb_coxphf_available()) {
      as.character(utils::packageVersion("coxphf"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_akritas_learner_ids()) {
    for (package in c("survivalmodels", "distr6"))
      provenance[[paste0(package, "_version")]] <- if (requireNamespace(package, quietly = TRUE))
        as.character(utils::packageVersion(package)) else NA_character_
  }
  if (pipeline$learner %in% sb_bagged_knn_survival_learner_ids()) {
    provenance$bnn_survival_version <- if (
      sb_bagged_knn_survival_available()
    ) {
      as.character(utils::packageVersion("bnnSurvival"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_ahaz_survival_learner_ids()) {
    provenance$ahaz_version <- if (sb_ahaz_available()) {
      as.character(utils::packageVersion("ahaz"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_bujar_linear_learner_ids()) {
    provenance$bujar_version <- if (sb_bujar_available()) {
      as.character(utils::packageVersion("bujar"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_rule_aft_learner_ids()) {
    dependency <- sb_rule_aft_dependency(pipeline$learner)
    for (package in dependency$packages) {
      field <- paste0(tolower(package), "_version")
      provenance[[field]] <- if (requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        NA_character_
      }
    }
  }
  if (pipeline$learner %in% sb_wsrf_binary_learner_ids()) {
    provenance$wsrf_version <- if (sb_wsrf_available()) {
      as.character(utils::packageVersion("wsrf"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_slope_randomglm_binary_learner_ids()) {
    dependency <- sb_slope_randomglm_binary_dependency(pipeline$learner)
    field <- paste0(tolower(dependency$package), "_version")
    provenance[[field]] <- if (dependency$available) {
      as.character(utils::packageVersion(dependency$package))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_mboost_cox_learner_ids()) {
    dependency <- sb_mboost_cox_dependency(pipeline$learner)
    for (package in dependency$packages) {
      field <- paste0(tolower(package), "_version")
      provenance[[field]] <- if (requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        NA_character_
      }
    }
  }
  if (pipeline$learner %in% sb_mboost_aft_learner_ids()) {
    dependency <- sb_mboost_aft_dependency(pipeline$learner)
    for (package in dependency$packages) {
      field <- paste0(tolower(package), "_version")
      provenance[[field]] <- if (requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        NA_character_
      }
    }
  }
  if (pipeline$learner %in% sb_mboost_ranking_learner_ids()) {
    dependency <- sb_mboost_ranking_dependency(pipeline$learner)
    for (package in dependency$packages) {
      field <- paste0(tolower(package), "_version")
      provenance[[field]] <- if (requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        NA_character_
      }
    }
  }
  if (pipeline$learner %in% sb_xgboost_learner_ids() || pipeline$selector %in% c("xgboost_binary_importance", "xgboost_cox_importance")) {
    provenance$xgboost_version <- if (sb_xgboost_available()) {
      as.character(utils::packageVersion("xgboost"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_classical_binary_learner_ids()) {
    dependency <- sb_classical_binary_dependency(pipeline$learner)
    if (!is.na(dependency$package)) {
      field <- paste0(tolower(dependency$package), "_version")
      provenance[[field]] <- if (dependency$available) {
        as.character(utils::packageVersion(dependency$package))
      } else {
        NA_character_
      }
    }
  }
  if (pipeline$learner %in% sb_specialized_binary_linear_learner_ids()) {
    dependency <- sb_specialized_binary_linear_dependency(pipeline$learner)
    if (!is.na(dependency$package)) {
      field <- paste0(tolower(dependency$package), "_version")
      provenance[[field]] <- if (dependency$available) {
        as.character(utils::packageVersion(dependency$package))
      } else {
        NA_character_
      }
    }
  }
  if (pipeline$learner %in% sb_omics_specialist_binary_learner_ids()) {
    dependency <- sb_omics_specialist_binary_dependency(pipeline$learner)
    field <- paste0(tolower(dependency$package), "_version")
    provenance[[field]] <- if (dependency$available) {
      as.character(utils::packageVersion(dependency$package))
    } else {
      NA_character_
    }
  }
  classical_dependency <- if (
    pipeline$learner %in% sb_survivalsvm_learner_ids()) {
    list(
      field = "survivalsvm_version",
      package = "survivalsvm",
      available = sb_survivalsvm_available()
    )
  } else {
    switch(pipeline$learner,
      coxboost = list(
        field = "coxboost_version",
        package = "CoxBoost",
        available = sb_coxboost_available()
      ),
      plsrcox = list(
        field = "plsrcox_version",
        package = "plsRcox",
        available = sb_plsrcox_available()
      ),
      superpc = list(
        field = "superpc_version",
        package = "superpc",
        available = sb_superpc_available()
      ),
      NULL
    )
  }
  if (!is.null(classical_dependency)) {
    provenance[[classical_dependency$field]] <- if (
      classical_dependency$available) {
      as.character(utils::packageVersion(classical_dependency$package))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_additional_survival_learner_ids()) {
    dependency <- sb_additional_survival_dependency(pipeline$learner)
    for (package in dependency$packages) {
      field <- paste0(tolower(package), "_version")
      if (!field %in% names(provenance)) {
        provenance[[field]] <- if (requireNamespace(package, quietly = TRUE)) {
          as.character(utils::packageVersion(package))
        } else {
          NA_character_
        }
      }
    }
  }
  if (pipeline$learner %in% c(
    sb_flexsurv_parametric_learner_ids(),
    sb_flexsurvspline_learner_ids()
  )) {
    provenance$flexsurv_version <- if (requireNamespace(
      "flexsurv",
      quietly = TRUE
    )) {
      as.character(utils::packageVersion("flexsurv"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_mtlr_survival_learner_ids()) {
    provenance$mtlr_version <- if (sb_mtlr_survival_available()) {
      as.character(utils::packageVersion("MTLR"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_pamr_survival_learner_ids()) {
    provenance$pamr_version <- if (sb_pamr_survival_available()) {
      as.character(utils::packageVersion("pamr"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_timereg_survival_learner_ids()) {
    provenance$timereg_version <- if (sb_timereg_available()) {
      as.character(utils::packageVersion("timereg"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_abess_learner_ids() ||
    pipeline$selector %in% sb_abess_selector_ids()) {
    provenance$abess_version <- if (sb_abess_available()) {
      as.character(utils::packageVersion("abess"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_mgcv_gam_learner_ids()) {
    provenance$mgcv_version <- if (sb_mgcv_available()) {
      as.character(utils::packageVersion("mgcv"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_mclust_da_learner_ids()) {
    provenance$mclust_version <- if (sb_mclust_available()) {
      as.character(utils::packageVersion("mclust"))
    } else {
      NA_character_
    }
  }
  if (pipeline$selector %in% sb_forest_selector_ids()) {
    selector_provenance <- sb_forest_selector_provenance(pipeline$selector)
    for (field in names(selector_provenance)) {
      provenance[[field]] <- selector_provenance[[field]]
    }
  }
  if (pipeline$selector == "plsrglm_binary_weights") {
    provenance$plsrglm_version <- if (requireNamespace("plsRglm", quietly = TRUE))
      as.character(utils::packageVersion("plsRglm")) else NA_character_
  }
  if (pipeline$selector %in% sb_pamr_selector_ids()) {
    provenance$pamr_version <- if (requireNamespace("pamr", quietly = TRUE))
      as.character(utils::packageVersion("pamr")) else NA_character_
  }
  high_dimensional_methods <- intersect(
    c(pipeline$selector, pipeline$learner),
    sb_high_dimensional_method_ids()
  )
  for (method in high_dimensional_methods) {
    method_provenance <- sb_high_dimensional_method_provenance(method)
    for (field in names(method_provenance)) {
      provenance[[field]] <- method_provenance[[field]]
    }
  }
  provenance
}

sb_final_fit_prediction_semantics <- function(task, pipeline) {
  if (inherits(task, "stormbreaker_binary_task")) {
    semantics <- list(
      output = "probability",
      positive_class = pipeline$positive_class,
      direction = pipeline$prediction_direction
    )
    if ("decision_thresholds" %in% names(task)) {
      semantics$decision_thresholds <- task$decision_thresholds
    }
    return(semantics)
  }
  semantics <- list(
    output = "risk",
    risk_type = pipeline$risk_type,
    risk_reference = pipeline$risk_reference,
    risk_formula = pipeline$risk_formula,
    risk_direction = pipeline$risk_direction
  )
  if ("prediction_horizon" %in% names(task)) {
    semantics$output <- c("risk", "survival_probability")
    semantics$prediction_horizon <- task$prediction_horizon
    semantics$time_unit <- task$time_unit
    if ("decision_thresholds" %in% names(task)) {
      semantics$decision_thresholds <- task$decision_thresholds
    }
  }
  semantics
}

sb_final_fit_development <- function(task) {
  if (inherits(task, "stormbreaker_binary_task")) {
    endpoint <- list(
      outcome_col = task$outcome_col,
      positive_class = task$positive_class,
      negative_class = task$negative_class,
      n_positive = task$n_positive,
      n_negative = task$n_negative
    )
    if ("decision_thresholds" %in% names(task)) {
      endpoint$decision_thresholds <- task$decision_thresholds
    }
  } else {
    endpoint <- list(
      time_col = task$time_col,
      event_col = task$event_col,
      time_unit = task$time_unit,
      n_events = task$n_events,
      n_censored = task$n_censored
    )
    if ("prediction_horizon" %in% names(task)) {
      endpoint <- append(
        endpoint,
        list(prediction_horizon = task$prediction_horizon),
        after = match("time_unit", names(endpoint))
      )
      if ("decision_thresholds" %in% names(task)) {
        endpoint <- append(
          endpoint,
          list(decision_thresholds = task$decision_thresholds),
          after = match("prediction_horizon", names(endpoint))
        )
      }
    }
  }
  development <- list(
    patient_ids = as.character(task$data[[task$id_col]]),
    n = task$n,
    independent_unit = task$independent_unit,
    candidate_features = task$predictors,
    endpoint = endpoint
  )
  if (!is.null(task$group_col)) {
    development$group_col <- task$group_col
    development$group_ids <- stats::setNames(
      sb_task_group_ids(task), development$patient_ids
    )
  }
  development
}

#' Fit one explicit pipeline on all development patients
#'
#' `sb_fit_final()` refits one user-supplied pipeline on the complete
#' development task. A learned selector is refitted on all development rows,
#' while a fixed feature block is applied unchanged without reading outcomes.
#' A task's declared `group_col` is retained for glmnet inner tuning while
#' unique row IDs remain separate. Grouping does not change metric weighting.
#' A glmnet learner repeats its declared inner cross-validation there. For a
#' supervised selected glmnet pipeline, the selector is additionally refitted
#' on every inner analysis fold before lambda is chosen; for a fixed feature
#' block, lambda is tuned only on those predeclared columns. Glmnet's declared
#' predictor-standardization setting is applied separately in every inner fit
#' and the final full-development refit. Its predeclared minimum-error or
#' one-standard-error lambda rule is reused without comparing rules on
#' development or external performance. The returned learner state is
#' self-contained for later state-only prediction, and its contents are
#' learner-specific. Coefficient-based learners retain prediction summaries;
#' engine-backed learners retain a pruned prediction object. KNN necessarily
#' retains the standardized development feature matrix and binary outcomes,
#' while SVM retains support vectors, coefficients, and class metadata.
#' LVQ states retain analysis-fold scaling values, patient-derived prototypes,
#' and outcome-derived prototype class labels. They contain no patient IDs or
#' complete development table, but a prototype can equal an initialized
#' standardized training row.
#' randomForestSRC retains event-time metadata inside its prediction forest.
#' CatBoost SurvivalAft states retain the serialized booster plus fixed Normal
#' distribution metadata, while cancerclass states retain its all-feature class
#' centroids and inner-OOF logistic calibration coefficients. Neither state
#' retains patient IDs or row-level outcomes.
#' XGBoost AFT states retain the serialized raw booster and fixed normal-loss
#' prediction metadata, but no patient IDs, follow-up times, or event indicators.
#' Mboost AFT states retain the serialized native engine, including the training
#' data that mboost requires for prediction, and must therefore be handled as
#' sensitive patient-derived objects.
#' The timereg learners retain their native event-time coefficient grid but no
#' patient IDs, censoring-time vector, event-indicator vector, or fitted engine.
#' partykit conditional-inference survival states retain their fitted
#' conditional-distribution data because the engine requires it for prediction.
#' Akritas states retain the analysis predictor basis and paired follow-up/event
#' values required for new-patient conditional product-limit estimation, without
#' patient IDs or row names. Treat these states as sensitive training data.
#' They use the native empirical joint CDF and fixed bandwidth 0.5; scalar risk
#' is negative step-RMST through maximum analysis follow-up, with direct native
#' survival probabilities rather than Cox recalibration.
#' Bagged kNN survival states necessarily retain the development feature matrix,
#' follow-up times, and event indicators required for nearest-neighbor prediction.
#' Ahaz and Bujar linear states retain portable coefficients and their
#' analysis-fold horizon-calibration summaries, not fitted engines, patient IDs,
#' raw follow-up vectors, or event-indicator vectors. The Bujar residual
#' Kaplan-Meier step curve is patient-derived calibration data, so final states
#' should still be handled as sensitive objects.
#' Weighted-subspace forest states remove training outcomes and patient-level OOB
#' indices while retaining the native trees required for prediction.
#' Fixed flexsurv parametric states retain only fitted regression and ancillary
#' distribution parameters, not the native engine or analysis outcomes.
#' Fixed flexsurvspline states retain only spline/regression coefficients and
#' analysis-fold knot locations, not the native engine or analysis outcomes.
#' Final
#' fit objects should therefore be handled as potentially patient-derived
#' sensitive data. A task's optional decision-threshold declaration is retained
#' only as prediction metadata and does not change fitting or produce action
#' labels.
#'
#' This operation does not estimate out-of-fold or external performance, does
#' not compare pipelines, and does not select a winner. The pipeline must be an
#' explicit specification created by [sb_pipeline()]. Execution failures are
#' returned as stage-aware failed fit objects without substituting another
#' selector or learner.
#'
#' @param task A task returned by [sb_surv_task()] or [sb_binary_task()].
#' @param pipeline One explicit task-matched pipeline returned by
#'   [sb_pipeline()].
#' @param capture_diagnostics Keep native learner and selector objects, fitted
#'   predictor data and training endpoints for optional algorithm diagnostics.
#'   Default FALSE preserves the existing compact result. TRUE may produce a
#'   large, sensitive, engine-version-dependent object; do not share it as a
#'   de-identified figure or portable prediction signature. No additional
#'   fitting is performed and saved prediction state is unchanged.
#'
#' @return A `stormbreaker_final_fit` containing the development-cohort audit,
#'   fitted selector state, inner-tuning audit when required, portable learner
#'   state, prediction semantics (including any declared decision
#'   thresholds), provenance, and any execution failure.
#' @export
#'
#' @examples
#' d <- sb_example_survival(n = 30, seed = 1)
#' task <- sb_surv_task(
#'   d, "patient_id", "time", "event",
#'   predictors = c("gene_a", "gene_b"),
#'   time_unit = "month"
#' )
#' fit <- sb_fit_final(task, sb_pipeline())
#' fit
sb_fit_final <- function(task, pipeline, capture_diagnostics = FALSE) {
  pipeline <- sb_upgrade_legacy_names(pipeline)

  if (!is.logical(capture_diagnostics) || length(capture_diagnostics) != 1L ||
      is.na(capture_diagnostics)) {
    sb_abort("capture_diagnostics must be TRUE or FALSE.", "stormbreaker_bad_visualization_argument")
  }
  if (inherits(task, "stormbreaker_binary_task")) {
    task <- sb_validate_binary_task_contract(task)
    task_type <- "binary"
    endpoint <- list(
      task_type = task_type,
      outcome = as.integer(task$data[[task$outcome_col]]),
      group_ids = sb_task_group_ids(task)
    )
  } else if (inherits(task, "stormbreaker_surv_task")) {
    task <- sb_validate_task_contract(task)
    task_type <- "survival"
    endpoint <- list(
      task_type = task_type,
      time = task$data[[task$time_col]],
      event = as.integer(task$data[[task$event_col]]),
      group_ids = sb_task_group_ids(task)
    )
  } else {
    sb_abort(
      "`task` must be created by `sb_surv_task()` or `sb_binary_task()`.",
      "stormbreaker_bad_task"
    )
  }
  pipeline <- sb_validate_pipeline_contract(pipeline, task_type)
  runtime <- sb_validate_pipeline_runtime(
    sb_resolve_pipeline_runtime(pipeline),
    pipeline
  )
  if (capture_diagnostics) {
    runtime$selector <- sb_captured_selector_adapter(pipeline$selector, pipeline$selector_params)
  }
  x <- task$data[task$predictors]
  patient_ids <- as.character(task$data[[task$id_col]])

  result <- list(
    status = "not_fitted",
    task_type = task_type,
    pipeline = pipeline,
    development = sb_final_fit_development(task),
    selection = list(
      selector = pipeline$selector,
      status = "not_fitted",
      selected_features = character(),
      scores = NULL,
      state = NULL,
      error_class = NA_character_,
      error_message = NA_character_
    ),
    tuning = sb_empty_final_tuning(pipeline),
    learner = list(
      status = "not_fitted",
      learner = pipeline$learner,
      error_class = NA_character_,
      error_message = NA_character_
    ),
    learner_state = NULL,
    required_features = character(),
    prediction = sb_final_fit_prediction_semantics(task, pipeline),
    provenance = sb_final_fit_provenance(pipeline),
    error_stage = NA_character_,
    error_class = NA_character_,
    error_message = NA_character_
  )

  selection_attempt <- sb_final_fit_attempt(
    sb_fit_runtime_selector(runtime$selector, x, endpoint)
  )
  if (!is.null(selection_attempt$error)) {
    return(sb_final_fit_failure(
      result,
      "selection",
      selection_attempt$error
    ))
  }
  fitted_selector <- selection_attempt$value
  result$selection <- list(
    selector = pipeline$selector,
    status = "success",
    selected_features = fitted_selector$selected_features,
    scores = if ("scores" %in% names(fitted_selector)) {
      fitted_selector$scores
    } else {
      NULL
    },
    state = fitted_selector,
    error_class = NA_character_,
    error_message = NA_character_
  )
  result$required_features <- fitted_selector$selected_features

  application_attempt <- sb_final_fit_attempt(
    runtime$selector$apply(fitted_selector, x)
  )
  if (!is.null(application_attempt$error)) {
    return(sb_final_fit_failure(
      result,
      "selection_application",
      application_attempt$error
    ))
  }
  selected_x <- application_attempt$value

  preparation_attempt <- sb_final_fit_attempt(
    sb_prepare_runtime_learner(
      runtime$learner,
      pipeline,
      selected_x,
      endpoint,
      patient_ids,
      group_ids = sb_task_group_ids(task),
      selection_context = list(
        selector = runtime$selector,
        fitted_selector = fitted_selector,
        x = x
      )
    )
  )
  if (!is.null(preparation_attempt$error)) {
    stage <- if (pipeline$learner %in% c(
      "glmnet_cox", "glmnet_logistic", sb_svm_probability_learner_ids()
    )) {
      "tuning"
    } else {
      "learner_preparation"
    }
    return(sb_final_fit_failure(
      result,
      stage,
      preparation_attempt$error
    ))
  }
  prepared <- preparation_attempt$value
  learner <- prepared$learner
  if (!is.null(prepared$tuning)) {
    result$tuning <- sb_successful_final_tuning(
      pipeline,
      prepared$tuning
    )
  }

  learner_attempt <- sb_final_fit_attempt(
    sb_fit_final_runtime_learner(learner, selected_x, endpoint,
      patient_ids, pipeline$selector)
  )
  if (!is.null(learner_attempt$error)) {
    return(sb_final_fit_failure(
      result,
      "learner_fit",
      learner_attempt$error
    ))
  }
  model <- learner_attempt$value
  result$learner$status <- "success"
  fit_warnings <- attr(model, "stormbreaker_fit_warnings", exact = TRUE)
  if (length(fit_warnings) > 0L) {
    result$tuning$warnings <- c(
      result$tuning$warnings, paste0("Learner refit: ", fit_warnings)
    )
  }

  context <- if (identical(task_type, "survival")) {
    survival_context <- list(
      x = selected_x,
      time = endpoint$time,
      event = endpoint$event
    )
    if ("prediction_horizon" %in% names(task)) {
      survival_context$prediction_horizon <- task$prediction_horizon
      survival_context$time_unit <- task$time_unit
    }
    survival_context
  } else {
    list()
  }
  state_attempt <- sb_final_fit_attempt(
    learner$extract_state(model, context)
  )
  if (!is.null(state_attempt$error)) {
    return(sb_final_fit_failure(
      result,
      "state_extraction",
      state_attempt$error
    ))
  }

  result$learner_state <- state_attempt$value
  if (capture_diagnostics) {
    result$diagnostic_artifacts <- list(learner = model, selector = fitted_selector,
      predictors = selected_x, endpoint = endpoint,
      source = "same final development fit; no additional fitting",
      sensitive_training_data = TRUE)
  }
  result$status <- "success"
  sb_as_final_fit(result)
}

#' @export
print.stormbreaker_final_fit <- function(x, ...) {
  cat("<stormbreaker_final_fit:", x$status, ">\n")
  cat("final development fit\n")
  cat("Pipeline:", x$pipeline$id, "| Task:", x$task_type, "\n")
  cat("Development patients:", x$development$n, "\n")
  cat("Required features:", length(x$required_features), "\n")
  if (identical(x$status, "failed")) {
    cat("Failed stage:", x$error_stage, "|", x$error_class, "\n")
  }
  invisible(x)
}

sb_predict_final_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (inherits(state, "stormbreaker_constant_sum_state")) {
    x <- sb_validate_learner_state_x(state, x, learner)
    return(sb_resolve_learner_adapter(learner, state$params)$predict_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_linear_survival_learner_ids()) {
    return(sb_predict_cox_learner_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_survivalsvm_kernel_learner_ids()) {
    return(sb_predict_survivalsvm_kernel_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_mtlr_survival_learner_ids()) {
    return(sb_predict_mtlr_survival_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_pamr_survival_learner_ids()) {
    return(sb_predict_pamr_survival_state(state, x))
  }
  if (learner %in% c(
    "logistic", "glmnet_logistic", sb_abess_binary_learner_ids(),
    sb_ncvreg_binary_learner_ids()
  )) {
    return(sb_predict_logistic_learner_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_rlt_learner_ids()) {
    return(sb_predict_rlt_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_bagged_knn_survival_learner_ids()) {
    return(sb_predict_bagged_knn_survival_state(state, x))
  }
  if (learner %in% sb_akritas_learner_ids()) return(sb_predict_akritas_state(state, x))
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_ahaz_survival_learner_ids()) {
    return(sb_predict_ahaz_survival_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_bujar_linear_learner_ids()) {
    return(sb_predict_bujar_linear_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_rule_aft_learner_ids()) {
    return(sb_predict_rule_aft_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_wsrf_binary_learner_ids()) {
    return(sb_predict_wsrf_binary_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_slope_randomglm_binary_learner_ids()) {
    return(sb_predict_slope_randomglm_binary_state(state, x))
  }
  if (identical(learner, "gaussian_naive_bayes")) {
    return(sb_predict_naive_bayes_learner_state(state, x))
  }
  if (learner %in% sb_lightgbm_binary_learner_ids()) {
    return(sb_predict_lightgbm_learner_state(state, x))
  }
  if (identical(learner, "catboost_binary")) {
    return(sb_predict_catboost_learner_state(state, x))
  }
  if (identical(learner, "catboost_cox")) {
    return(sb_predict_catboost_cox_learner_state(state, x))
  }
  if (identical(learner, "catboost_aft_survival")) {
    return(sb_predict_catboost_aft_survival_learner_state(state, x))
  }
  if (identical(learner, "cancerclass_binary")) {
    return(sb_predict_cancerclass_binary_learner_state(state, x))
  }
  if (learner %in% c("gbm_cox", "gbm_cox_cv")) {
    return(sb_predict_gbm_cox_learner_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_mboost_cox_learner_ids()) {
    return(sb_predict_mboost_cox_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_mboost_aft_learner_ids()) {
    return(sb_predict_mboost_aft_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_mboost_ranking_learner_ids()) {
    return(sb_predict_mboost_ranking_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_xgboost_learner_ids()) {
    return(sb_predict_xgboost_learner_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_additional_survival_learner_ids()) {
    return(sb_predict_additional_survival_learner_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_flexsurv_parametric_learner_ids()) {
    return(sb_predict_flexsurv_parametric_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_flexsurvspline_learner_ids()) {
    return(sb_predict_flexsurvspline_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_timereg_survival_learner_ids()) {
    return(sb_predict_timereg_survival_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_ranger_binary_learner_ids()) {
    return(sb_predict_ranger_binary_learner_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_classical_binary_boosting_learner_ids()) {
    return(sb_predict_classical_binary_boosting_learner_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_classical_binary_learner_ids()) {
    return(sb_predict_classical_binary_learner_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_specialized_binary_linear_learner_ids()) {
    return(sb_predict_specialized_binary_linear_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_omics_specialist_binary_learner_ids()) {
    return(sb_predict_omics_specialist_binary_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_mgcv_gam_learner_ids()) {
    return(sb_predict_mgcv_gam_learner_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_mclust_da_learner_ids()) {
    return(sb_predict_mclust_da_learner_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_mixomics_da_learner_ids()) {
    return(sb_predict_mixomics_da_learner_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_hdda_binary_learner_ids()) {
    return(sb_predict_hdda_binary_state(state, x))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_ranger_survival_learner_ids()) {
    return(sb_predict_ranger_survival_learner_state(state, x))
  }
  sb_abort(
    "The final fit has no supported portable learner state.",
    "stormbreaker_bad_final_fit"
  )
}

sb_final_prediction_metadata <- function(fit) {
  metadata <- list(
    task_type = NA_character_,
    pipeline_id = NA_character_,
    learner = NA_character_,
    decision_thresholds = NULL,
    prediction_horizon = NULL,
    time_unit = NULL
  )
  if (!is.list(fit)) {
    return(metadata)
  }
  task_type <- fit[["task_type"]]
  if (is.character(task_type) && length(task_type) == 1L &&
    !is.na(task_type)) {
    metadata$task_type <- task_type
  }
  pipeline <- fit[["pipeline"]]
  if (is.list(pipeline)) {
    if (is.character(pipeline[["id"]]) &&
      length(pipeline[["id"]]) == 1L &&
      !is.na(pipeline[["id"]])) {
      metadata$pipeline_id <- pipeline[["id"]]
    }
    if (is.character(pipeline[["learner"]]) &&
      length(pipeline[["learner"]]) == 1L &&
      !is.na(pipeline[["learner"]])) {
      metadata$learner <- pipeline[["learner"]]
    }
  }
  prediction <- fit[["prediction"]]
  has_declared_horizon <- is.list(prediction) && identical(
    prediction[["output"]],
    c("risk", "survival_probability")
  )
  if (has_declared_horizon) {
    horizon <- prediction[["prediction_horizon"]]
    time_unit <- prediction[["time_unit"]]
    valid_horizon <- is.numeric(horizon) && length(horizon) == 1L &&
      !is.na(horizon) && is.finite(horizon) && horizon > 0
    valid_time_unit <- is.character(time_unit) &&
      length(time_unit) == 1L && !is.na(time_unit) && nzchar(time_unit)
    if (valid_horizon && valid_time_unit) {
      metadata$prediction_horizon <- horizon
      metadata$time_unit <- time_unit
    }
  }
  if (is.list(prediction)) {
    thresholds <- prediction[["decision_thresholds"]]
    if (sb_valid_decision_thresholds(thresholds)) {
      metadata$decision_thresholds <- unname(as.numeric(thresholds))
    }
  }
  metadata
}

sb_final_prediction_rows <- function(fit, patient_ids, status,
                                     error_class, error_message,
                                     values = NULL) {
  metadata <- sb_final_prediction_metadata(fit)
  n <- length(patient_ids)
  rows <- data.frame(
    patient_id = patient_ids,
    pipeline_id = rep(metadata$pipeline_id, n),
    learner = rep(metadata$learner, n),
    stringsAsFactors = FALSE
  )
  if (identical(metadata$task_type, "binary")) {
    rows$positive_class <- rep.int(1L, n)
    rows$probability <- if (is.null(values)) {
      rep(NA_real_, n)
    } else {
      values$probability
    }
  } else if (identical(metadata$task_type, "survival")) {
    rows$risk <- if (is.null(values)) rep(NA_real_, n) else values$risk
    if (!is.null(metadata$prediction_horizon)) {
      rows$prediction_horizon <- rep(metadata$prediction_horizon, n)
      rows$time_unit <- rep(metadata$time_unit, n)
      rows$survival_probability <- if (is.null(values)) {
        rep(NA_real_, n)
      } else {
        values$survival_probability
      }
      stable <- sb_final_cox_log_cumulative_hazard(fit, values)
      if (!is.null(stable)) rows$log_cumulative_hazard <- stable
    }
  }
  rows$status <- rep(status, n)
  rows$error_class <- rep(error_class, n)
  rows$error_message <- rep(error_message, n)
  rows
}

sb_as_final_prediction <- function(result) {
  result$prediction_variation <- sb_prediction_variation(result$predictions, result$task_type)
  class(result) <- "stormbreaker_final_prediction"
  result
}

sb_final_prediction_failure <- function(fit, id_col, patient_ids,
                                        stage, error) {
  metadata <- sb_final_prediction_metadata(fit)
  error_class <- class(error)[1L]
  error_message <- conditionMessage(error)
  result <- list(
    status = "failed",
    task_type = metadata$task_type,
    pipeline_id = metadata$pipeline_id,
    learner = metadata$learner,
    id_col = if (is.character(id_col) && length(id_col) == 1L &&
      !is.na(id_col)) {
      id_col
    } else {
      NA_character_
    },
    predictions = sb_final_prediction_rows(
      fit,
      patient_ids,
      "failed",
      error_class,
      error_message
    ),
    error_stage = stage,
    error_class = error_class,
    error_message = error_message
  )
  if (!is.null(metadata$decision_thresholds)) {
    result <- append(
      result,
      list(decision_thresholds = metadata$decision_thresholds),
      after = match("id_col", names(result))
    )
  }
  sb_as_final_prediction(result)
}

sb_extract_final_prediction_ids <- function(new_data, id_col) {
  if (!is.data.frame(new_data) || nrow(new_data) < 1L) {
    sb_abort(
      "`new_data` must be a non-empty data frame.",
      "stormbreaker_bad_prediction_data"
    )
  }
  sb_check_name(id_col, "id_col")
  id_positions <- which(names(new_data) == id_col)
  if (length(id_positions) == 0L) {
    sb_abort(
      paste0("Prediction data are missing ID column `", id_col, "`."),
      "stormbreaker_missing_prediction_id"
    )
  }
  if (length(id_positions) > 1L) {
    sb_abort(
      paste0("Prediction ID column `", id_col, "` is duplicated."),
      "stormbreaker_duplicate_prediction_id_column"
    )
  }
  patient_ids <- as.character(new_data[[id_positions]])
  if (anyNA(patient_ids) || any(!nzchar(trimws(patient_ids)))) {
    sb_abort(
      "Prediction patient identifiers contain missing or blank values.",
      "stormbreaker_bad_prediction_id"
    )
  }
  if (anyDuplicated(patient_ids)) {
    sb_abort(
      "Prediction patient identifiers must be unique.",
      "stormbreaker_duplicate_prediction_id"
    )
  }
  patient_ids
}

sb_final_direct_survival_probability_source <- function(learner) {
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_rlt_survival_learner_ids()) {
    return(sb_rlt_direct_survival_source())
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_bagged_knn_survival_learner_ids()) {
    return(sb_bagged_knn_direct_survival_source())
  }
  if (learner %in% sb_akritas_learner_ids()) return(sb_akritas_direct_source())
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_bujar_linear_learner_ids()) {
    return(sb_bujar_direct_survival_source())
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_rule_aft_survival_learner_ids()) {
    return(sb_rule_aft_direct_survival_source(learner))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_ranger_survival_learner_ids()) {
    return("ranger_survival_forest_direct")
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_xgboost_aft_learner_ids()) {
    return("xgboost_aft_normal_native")
  }
  if (identical(learner, "catboost_aft_survival")) {
    return("catboost_aft_normal_native")
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_mboost_aft_learner_ids()) {
    return(sb_mboost_aft_direct_survival_source(learner))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_mboost_ranking_learner_ids()) {
    return(sb_mboost_ranking_direct_survival_source(learner))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_additional_survival_learner_ids()) {
    return(sb_additional_direct_survival_source(learner))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_flexsurv_parametric_learner_ids()) {
    return(sb_flexsurv_parametric_direct_survival_source(learner))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_flexsurvspline_learner_ids()) {
    return(sb_flexsurvspline_direct_survival_source(learner))
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_mtlr_survival_learner_ids()) {
    return(sb_mtlr_direct_survival_source())
  }
  if (is.character(learner) && length(learner) == 1L && !is.na(learner) &&
    learner %in% sb_timereg_survival_learner_ids()) {
    return(sb_timereg_direct_survival_source(learner))
  }
  NULL
}

sb_final_direct_survival_probability_ties <- function(learner) {
  if (learner %in% sb_akritas_learner_ids()) return("native_grouped_event_times")
  "not_applicable"
}

sb_validate_final_survival_state_semantics <- function(fit) {
  prediction <- fit[["prediction"]]
  state <- fit[["learner_state"]]
  baseline <- state[["baseline_hazard"]]
  direct_probability <- state[["survival_probability"]]
  direct_source <- sb_final_direct_survival_probability_source(
    fit[["pipeline"]][["learner"]]
  )
  uses_direct_probability <- !is.null(direct_source)
  output <- if (is.list(prediction)) prediction[["output"]] else NULL
  if (identical(output, "risk")) {
    if (!is.null(baseline) || !is.null(direct_probability)) {
      sb_abort(
        paste0(
          "A risk-only final fit cannot carry a prediction-horizon ",
          "state."
        ),
        "stormbreaker_final_state_semantics_mismatch"
      )
    }
    return(invisible(TRUE))
  }

  horizon_output <- identical(
    output,
    c("risk", "survival_probability")
  )
  declared_horizon <- if (is.list(prediction)) {
    prediction[["prediction_horizon"]]
  } else {
    NULL
  }
  declared_unit <- if (is.list(prediction)) {
    prediction[["time_unit"]]
  } else {
    NULL
  }
  valid_declaration <- horizon_output &&
    is.numeric(declared_horizon) && length(declared_horizon) == 1L &&
    !is.na(declared_horizon) && is.finite(declared_horizon) &&
    declared_horizon > 0 &&
    is.character(declared_unit) && length(declared_unit) == 1L &&
    !is.na(declared_unit) && nzchar(declared_unit)
  if (uses_direct_probability) {
    saved_horizon <- if (is.list(direct_probability)) {
      direct_probability[["prediction_horizon"]]
    } else {
      NULL
    }
    saved_unit <- if (is.list(direct_probability)) {
      direct_probability[["time_unit"]]
    } else {
      NULL
    }
    valid_direct <- isTRUE(valid_declaration) && is.null(baseline) &&
      is.list(direct_probability) &&
      is.numeric(saved_horizon) && length(saved_horizon) == 1L &&
      !is.na(saved_horizon) && is.finite(saved_horizon) &&
      identical(as.numeric(saved_horizon), as.numeric(declared_horizon)) &&
      identical(saved_unit, declared_unit) &&
      identical(direct_probability[["source"]], direct_source) &&
      identical(
        direct_probability[["ties"]],
        sb_final_direct_survival_probability_ties(
          fit[["pipeline"]][["learner"]]
        )
      )
    if (!valid_direct) {
      sb_abort(
        paste0(
          "The declared survival output and saved direct-probability ",
          "state do not match."
        ),
        "stormbreaker_final_state_semantics_mismatch"
      )
    }
    return(invisible(TRUE))
  }
  saved_horizon <- if (is.list(baseline)) {
    baseline[["prediction_horizon"]]
  } else {
    NULL
  }
  saved_unit <- if (is.list(baseline)) baseline[["time_unit"]] else NULL
  saved_hazard <- if (is.list(baseline)) {
    baseline[["cumulative_baseline_hazard"]]
  } else {
    NULL
  }
  valid_baseline <- isTRUE(valid_declaration) &&
    is.null(direct_probability) && is.list(baseline) &&
    is.numeric(saved_horizon) && length(saved_horizon) == 1L &&
    !is.na(saved_horizon) && is.finite(saved_horizon) &&
    identical(as.numeric(saved_horizon), as.numeric(declared_horizon)) &&
    identical(saved_unit, declared_unit) &&
    is.numeric(saved_hazard) && length(saved_hazard) == 1L &&
    !is.na(saved_hazard) && is.finite(saved_hazard) && saved_hazard >= 0
  if (!valid_declaration || !valid_baseline) {
    sb_abort(
      paste0(
        "The declared survival output and saved baseline-hazard state ",
        "do not match."
      ),
      "stormbreaker_final_state_semantics_mismatch"
    )
  }
  invisible(TRUE)
}

sb_validate_final_decision_semantics <- function(fit) {
  prediction <- fit[["prediction"]]
  development <- fit[["development"]]
  endpoint <- if (is.list(development)) {
    development[["endpoint"]]
  } else {
    NULL
  }
  has_prediction_thresholds <- is.list(prediction) &&
    "decision_thresholds" %in% names(prediction)
  has_development_thresholds <- is.list(endpoint) &&
    "decision_thresholds" %in% names(endpoint)
  if (!identical(
    has_prediction_thresholds,
    has_development_thresholds
  )) {
    sb_abort(
      "The final fit decision-threshold declarations do not match.",
      "stormbreaker_final_state_semantics_mismatch"
    )
  }
  if (!has_prediction_thresholds) {
    return(invisible(TRUE))
  }

  prediction_thresholds <- sb_validate_decision_thresholds(
    prediction[["decision_thresholds"]]
  )
  development_thresholds <- sb_validate_decision_thresholds(
    endpoint[["decision_thresholds"]]
  )
  if (!identical(prediction_thresholds, development_thresholds)) {
    sb_abort(
      "The final fit decision-threshold declarations do not match.",
      "stormbreaker_final_state_semantics_mismatch"
    )
  }
  if (identical(fit[["task_type"]], "survival")) {
    prediction_horizon <- prediction[["prediction_horizon"]]
    development_horizon <- endpoint[["prediction_horizon"]]
    valid_horizon <- identical(
      prediction[["output"]],
      c("risk", "survival_probability")
    ) && is.numeric(prediction_horizon) &&
      length(prediction_horizon) == 1L && !is.na(prediction_horizon) &&
      is.finite(prediction_horizon) && prediction_horizon > 0 &&
      is.numeric(development_horizon) &&
      length(development_horizon) == 1L &&
      !is.na(development_horizon) && is.finite(development_horizon) &&
      development_horizon > 0 &&
      identical(
        as.numeric(prediction_horizon),
        as.numeric(development_horizon)
      )
    if (!isTRUE(valid_horizon)) {
      sb_abort(
        "Survival decision thresholds require one matching horizon.",
        "stormbreaker_final_state_semantics_mismatch"
      )
    }
  }
  invisible(TRUE)
}

sb_validate_final_fit_for_prediction <- function(fit) {
  if (!inherits(fit, "stormbreaker_final_fit") || !is.list(fit)) {
    sb_abort(
      "`fit` must be returned by `sb_fit_final()`.",
      "stormbreaker_bad_final_fit"
    )
  }
  if (identical(fit[["status"]], "failed")) {
    sb_abort(
      paste0(
        "The final fit failed at `", fit[["error_stage"]],
        "` and cannot predict."
      ),
      "stormbreaker_failed_final_fit"
    )
  }
  pipeline <- fit[["pipeline"]]
  state <- fit[["learner_state"]]
  required_features <- fit[["required_features"]]
  task_type <- fit[["task_type"]]
  valid <- identical(fit[["status"]], "success") &&
    is.character(task_type) && length(task_type) == 1L &&
    !is.na(task_type) && task_type %in% c("survival", "binary") &&
    inherits(pipeline, "stormbreaker_pipeline") &&
    is.character(pipeline[["id"]]) && length(pipeline[["id"]]) == 1L &&
    is.character(pipeline[["learner"]]) &&
    length(pipeline[["learner"]]) == 1L &&
    is.character(required_features) && length(required_features) > 0L &&
    !anyNA(required_features) && !anyDuplicated(required_features) &&
    is.list(state) && identical(state[["learner"]], pipeline[["learner"]])
  if (!isTRUE(valid)) {
    sb_abort(
      "The final fit is incomplete or inconsistent for prediction.",
      "stormbreaker_bad_final_fit"
    )
  }
  survival_learner <- pipeline[["learner"]] %in% sb_survival_learner_ids()
  if (!identical(
    fit[["task_type"]],
    if (survival_learner) "survival" else "binary"
  )) {
    sb_abort(
      "The final fit task and learner do not match.",
      "stormbreaker_bad_final_fit"
    )
  }
  if (survival_learner) {
    sb_validate_final_survival_state_semantics(fit)
  }
  sb_validate_final_decision_semantics(fit)
  fit
}

sb_select_final_feature_frame <- function(data, required_features,
                                          data_label, argument) {
  if (!is.data.frame(data) || nrow(data) < 1L) {
    sb_abort(
      paste0("`", argument, "` must be a non-empty data frame."),
      "stormbreaker_bad_adapter_predictors"
    )
  }
  counts <- vapply(required_features, function(feature) {
    sum(names(data) == feature)
  }, integer(1))
  if (any(counts == 0L)) {
    sb_abort(
      paste0(
        data_label, " are missing fitted features: ",
        paste(required_features[counts == 0L], collapse = ", "),
        "."
      ),
      "stormbreaker_missing_final_features"
    )
  }
  if (any(counts > 1L)) {
    sb_abort(
      paste0(
        data_label, " duplicate fitted features: ",
        paste(required_features[counts > 1L], collapse = ", "),
        "."
      ),
      "stormbreaker_duplicate_final_feature_columns"
    )
  }
  positions <- vapply(required_features, function(feature) {
    which(names(data) == feature)
  }, integer(1))
  x <- data[positions]
  names(x) <- required_features
  sb_validate_predictor_frame(x, argument)
}

sb_select_final_prediction_features <- function(new_data, id_col,
                                                required_features) {
  if (id_col %in% required_features) {
    sb_abort(
      "The prediction ID column cannot also be a fitted feature.",
      "stormbreaker_overlapping_prediction_roles"
    )
  }
  sb_select_final_feature_frame(
    new_data,
    required_features,
    "Prediction data",
    "new_data fitted features"
  )
}

sb_validate_final_prediction_values <- function(fit, values, n) {
  if (identical(fit[["task_type"]], "binary")) {
    probability <- if (is.list(values)) values[["probability"]] else NULL
    valid <- is.numeric(probability) && length(probability) == n &&
      !anyNA(probability) && all(is.finite(probability)) &&
      all(probability >= 0 & probability <= 1)
    if (!valid) {
      sb_abort(
        "The saved final state produced invalid class-one probabilities.",
        "stormbreaker_invalid_final_probability"
      )
    }
    return(list(probability = unname(probability)))
  }

  risk <- if (is.list(values)) values[["risk"]] else NULL
  valid_risk <- is.numeric(risk) && length(risk) == n &&
    !anyNA(risk) && all(is.finite(risk))
  if (!valid_risk) {
    sb_abort(
      "The saved final state produced invalid survival risks.",
      "stormbreaker_invalid_final_risk"
    )
  }
  result <- list(risk = unname(risk))
  horizon_output <- identical(
    fit[["prediction"]][["output"]],
    c("risk", "survival_probability")
  )
  if (horizon_output) {
    probability <- values[["survival_probability"]]
    valid_probability <- is.numeric(probability) &&
      length(probability) == n && !anyNA(probability) &&
      all(is.finite(probability)) &&
      all(probability >= 0 & probability <= 1)
    if (!valid_probability) {
      sb_abort(
        "The saved final state produced invalid survival probabilities.",
        "stormbreaker_invalid_final_survival_probability"
      )
    }
    result$survival_probability <- unname(probability)
  }
  result
}

#' Predict new patients from a saved final fit without reading outcomes
#'
#' `sb_predict_final()` uses only the declared ID column and the exact numeric
#' features saved by [sb_fit_final()]. Other columns are ignored, including
#' columns named like outcomes, survival times, or events. Required features
#' are restored to their fitted order, and prediction calls only the portable
#' learner state; no selector, learner, or baseline hazard is refitted.
#'
#' Binary fits return the probability of class `1`. Survival fits return their
#' declared risk score: `x_beta` for ordinary Cox, glmnet Cox, fixed CoxBoost,
#' plsRcox, SuperPC, or calibrated linear survival-SVM; the calibrated native
#' kernel score for additive or RBF survival-SVM; the native raw margin for fixed
#' CatBoost Cox; the negative native log-time margin for CatBoost SurvivalAft
#' or XGBoost AFT; the native link-scale log hazard for fixed GBM-Cox or
#' XGBoost-Cox; the
#' negative native distributional log-time link for mboost AFT; the negative
#' native ranking link for mboost Gehan or Uno-C; the native additive-hazard
#' linear predictor for ahaz; negative predicted log time for Bujar; native
#' forest mortality for ranger split-rule forests or randomForestSRC; negative fitted log time for
#' SurvReg; log relative rate for survival Rpart;
#' final conditional cumulative hazard for partykit Ctree and CForest;
#' cumulative hazard at the last analysis event for Aalen; or `x_beta` for
#' proportional-odds survival regression.
#' When the development task declared a prediction horizon, learners with a
#' saved baseline hazard return the corresponding survival probability; ranger,
#' randomForestSRC, CatBoost SurvivalAft, XGBoost AFT, mboost AFT, mboost Gehan,
#' SurvReg, Ctree,
#' CForest, Aalen, and proportional odds read it directly from their saved
#' fitted survival representation. Mboost Uno-C instead uses its saved
#' analysis-fold Cox/Breslow horizon calibration. An invalid native
#' Aalen probability is not clipped or recalibrated. When a binary task, or a horizon-
#' enabled survival task, declared decision thresholds, the endpoint-free
#' prediction carries that unchanged vector as top-level metadata but does not
#' create class or action labels. Input or state failures return an explicit
#' failed prediction audit and never invoke a fallback learner.
#'
#' @param fit A successful `stormbreaker_final_fit` from [sb_fit_final()].
#' @param new_data A non-empty data frame containing one row per new patient,
#'   the declared ID column, and every fitted numeric feature. Other columns
#'   are ignored.
#' @param id_col Name of the unique patient identifier column in `new_data`.
#'
#' @return A `stormbreaker_final_prediction` with patient-order-preserving
#'   state-only predictions and explicit success or failure information.
#'   `prediction_variation` labels exactly constant risk/probability batches
#'   without changing predictions or fit status. A one-patient batch is
#'   insufficient to assess variation, not evidence of a constant model.
#' @details
#' For ordinary Cox final predictions whose probabilities round to zero or one,
#' the optional log_cumulative_hazard column preserves log(H0(h)) + x beta from
#' the saved model when finite. External time AUC and complementary-log-log
#' calibration, including bootstrap, use this monotone scale. Probabilities,
#' Brier scores and decision curves are not clipped or changed. Without this
#' saved scale, boundary-probability calibration keeps its existing unavailable
#' status; arbitrary risk scores are never substituted.
#' @export
#'
#' @examples
#' d <- sb_example_binary(n = 30, seed = 1)
#' task <- sb_binary_task(
#'   d, "patient_id", "outcome",
#'   predictors = c("gene_a", "gene_b")
#' )
#' fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
#' new_data <- d[1:3, c("patient_id", "gene_a", "gene_b")]
#' sb_predict_final(fit, new_data, "patient_id")
sb_predict_final <- function(fit, new_data, id_col) {
  fit <- sb_upgrade_legacy_names(fit)

  id_attempt <- sb_final_fit_attempt(
    sb_extract_final_prediction_ids(new_data, id_col)
  )
  if (!is.null(id_attempt$error)) {
    return(sb_final_prediction_failure(
      fit,
      id_col,
      character(),
      "id",
      id_attempt$error
    ))
  }
  patient_ids <- id_attempt$value

  fit_attempt <- sb_final_fit_attempt(
    sb_validate_final_fit_for_prediction(fit)
  )
  if (!is.null(fit_attempt$error)) {
    return(sb_final_prediction_failure(
      fit,
      id_col,
      patient_ids,
      "fit",
      fit_attempt$error
    ))
  }
  fit <- fit_attempt$value

  feature_attempt <- sb_final_fit_attempt(
    sb_select_final_prediction_features(
      new_data,
      id_col,
      fit$required_features
    )
  )
  if (!is.null(feature_attempt$error)) {
    return(sb_final_prediction_failure(
      fit,
      id_col,
      patient_ids,
      "features",
      feature_attempt$error
    ))
  }
  x <- feature_attempt$value

  prediction_attempt <- sb_final_fit_attempt({
    values <- sb_predict_final_state(fit$learner_state, x)
    sb_validate_final_prediction_values(fit, values, nrow(x))
  })
  if (!is.null(prediction_attempt$error)) {
    return(sb_final_prediction_failure(
      fit,
      id_col,
      patient_ids,
      "prediction",
      prediction_attempt$error
    ))
  }
  metadata <- sb_final_prediction_metadata(fit)
  result <- list(
    status = "success",
    task_type = metadata$task_type,
    pipeline_id = metadata$pipeline_id,
    learner = metadata$learner,
    id_col = id_col,
    predictions = sb_final_prediction_rows(
      fit,
      patient_ids,
      "success",
      NA_character_,
      NA_character_,
      prediction_attempt$value
    ),
    error_stage = NA_character_,
    error_class = NA_character_,
    error_message = NA_character_
  )
  if (!is.null(metadata$decision_thresholds)) {
    result <- append(
      result,
      list(decision_thresholds = metadata$decision_thresholds),
      after = match("id_col", names(result))
    )
  }
  sb_as_final_prediction(result)
}

#' @export
print.stormbreaker_final_prediction <- function(x, ...) {
  cat("<stormbreaker_final_prediction:", x$status, ">\n")
  origin <- if (identical(x$learner, "locked_linear_signature")) {
    "locked formula-only prediction"
  } else {
    "final state-only prediction"
  }
  cat(origin, "\n")
  cat("Pipeline:", x$pipeline_id, "| Task:", x$task_type, "\n")
  cat("Patients:", nrow(x$predictions), "\n")
  if (!is.null(x$prediction_variation)) {
    cat("Prediction variation:", x$prediction_variation$status, "\n")
  }
  if ("decision_thresholds" %in% names(x)) {
    cat(
      "Declared decision thresholds:",
      length(x$decision_thresholds), "\n"
    )
  }
  if (identical(x$status, "failed")) {
    cat("Failed stage:", x$error_stage, "|", x$error_class, "\n")
  }
  invisible(x)
}
