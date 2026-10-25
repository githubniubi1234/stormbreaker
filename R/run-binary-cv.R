sb_run_one_binary_fold <- function(task, pipeline, runtime, split, seed) {
  id <- as.character(task$data[[task$id_col]])
  outcome <- as.integer(task$data[[task$outcome_col]])
  assessment_ids <- id[split$assessment_rows]
  empty_predictions <- data.frame(
    patient_id = assessment_ids,
    fold = split$fold,
    outcome = outcome[split$assessment_rows],
    positive_class = task$positive_class,
    probability = NA_real_,
    status = "failed",
    pipeline_id = pipeline$id,
    learner = pipeline$learner,
    stringsAsFactors = FALSE
  )
  is_glmnet <- identical(pipeline$learner, "glmnet_logistic")
  is_svm <- pipeline$learner %in% sb_svm_probability_learner_ids()
  is_naive_bayes <- identical(
    pipeline$learner,
    "gaussian_naive_bayes"
  )
  is_lightgbm <- pipeline$learner %in% sb_lightgbm_binary_learner_ids()
  is_catboost <- identical(
    pipeline$learner,
    "catboost_binary"
  )
  is_ranger <- pipeline$learner %in% sb_ranger_binary_learner_ids()
  is_classical_boosting <- pipeline$learner %in%
    sb_classical_binary_boosting_learner_ids()
  empty_selection <- list(
    fold = split$fold,
    selector = pipeline$selector,
    status = "failed",
    selected_features = character(),
    scores = NULL,
    error_class = NA_character_,
    error_message = NA_character_
  )
  empty_tuning <- list(
    fold = split$fold,
    learner = pipeline$learner,
    status = if (is_glmnet) "failed" else "not_required",
    lambda = NA_real_,
    alpha = if (is_glmnet) pipeline$learner_params$alpha else NA_real_,
    selection_rule = if (is_glmnet) {
      pipeline$learner_params$selection_rule
    } else {
      NA_character_
    },
    measure = if (is_glmnet) {
      pipeline$learner_params$measure
    } else {
      NA_character_
    },
    inner_folds = if (is_glmnet) {
      pipeline$learner_params$inner_folds
    } else {
      NA_integer_
    },
    inner_seed = if (is_glmnet) {
      pipeline$learner_params$inner_seed
    } else {
      NA_integer_
    },
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
    error_class = NA_character_,
    error_message = NA_character_,
    standardize = if (is_glmnet) {
      pipeline$learner_params$standardize
    } else {
      NA
    }
  )
  if (is_svm) {
    empty_tuning <- c(list(fold = split$fold), sb_empty_final_tuning(pipeline))
    empty_tuning$status <- "failed"
  }
  empty_coefficients <- if (is_glmnet) {
    list(
      fold = split$fold,
      learner = pipeline$learner,
      status = "not_fitted",
      lambda = NA_real_,
      alpha = pipeline$learner_params$alpha,
      n_nonzero = 0L,
      nonzero_features = character(),
      values = data.frame(
        term = character(),
        coefficient = numeric(),
        nonzero = logical(),
        stringsAsFactors = FALSE
      )
    )
  } else if (!is.function(runtime$learner$extract_coefficients)) {
    list(
      fold = split$fold,
      learner = pipeline$learner,
      status = "not_applicable",
      values = data.frame(
        term = character(),
        coefficient = numeric(),
        stringsAsFactors = FALSE
      )
    )
  } else {
    list(
      fold = split$fold,
      learner = pipeline$learner,
      status = "not_fitted",
      values = data.frame(
        term = character(),
        coefficient = numeric(),
        stringsAsFactors = FALSE
      )
    )
  }
  empty_distribution_parameters <- list(
    fold = split$fold,
    learner = pipeline$learner,
    status = "not_fitted",
    var_smoothing = if (is_naive_bayes) {
      pipeline$learner_params$var_smoothing
    } else {
      NA_real_
    },
    epsilon = NA_real_,
    priors = data.frame(
      class = integer(),
      n = integer(),
      prior = numeric(),
      stringsAsFactors = FALSE
    ),
    values = data.frame(
      feature = character(),
      class = integer(),
      mean = numeric(),
      variance = numeric(),
      standard_deviation = numeric(),
      stringsAsFactors = FALSE
    )
  )
  empty_tree_model <- if (is_classical_boosting) {
    list(
      fold = split$fold,
      learner = pipeline$learner,
      status = "not_fitted",
      iterations = pipeline$learner_params$iterations,
      fitted_iterations = NA_integer_,
      features = character(),
      params = pipeline$learner_params,
      importance_type = "split_count",
      model_serialization = NA_character_,
      model_raw = raw(),
      feature_importance = data.frame(
        feature = character(),
        importance = numeric(),
        stringsAsFactors = FALSE
      )
    )
  } else if (is_catboost) {
    list(
      fold = split$fold,
      learner = pipeline$learner,
      status = "not_fitted",
      iterations = pipeline$learner_params$iterations,
      tree_count = NA_integer_,
      features = character(),
      params = pipeline$learner_params,
      model_serialization = NA_character_,
      model_raw = raw(),
      feature_importance = data.frame(
        feature = character(),
        importance = numeric(),
        stringsAsFactors = FALSE
      )
    )
  } else if (is_ranger) {
    list(
      fold = split$fold,
      learner = pipeline$learner,
      status = "not_fitted",
      num_trees = pipeline$learner_params$num_trees,
      effective_mtry = NA_integer_,
      features = character(),
      params = pipeline$learner_params,
      importance_type = pipeline$learner_params$importance,
      model_serialization = NA_character_,
      model_raw = raw(),
      feature_importance = data.frame(
        feature = character(),
        importance = numeric(),
        stringsAsFactors = FALSE
      )
    )
  } else {
    list(
      fold = split$fold,
      learner = pipeline$learner,
      status = "not_fitted",
      nrounds = if (is_lightgbm) {
        pipeline$learner_params$nrounds
      } else {
        NA_integer_
      },
      features = character(),
      params = if (is_lightgbm) pipeline$learner_params else NULL,
      model_string = NA_character_,
      feature_importance = data.frame(
        feature = character(),
        gain = numeric(),
        cover = numeric(),
        frequency = integer(),
        stringsAsFactors = FALSE
      )
    )
  }
  selection <- empty_selection
  tuning <- empty_tuning
  coefficients <- empty_coefficients
  distribution_parameters <- empty_distribution_parameters
  tree_model <- empty_tree_model

  tryCatch(
    {
      runtime <- sb_validate_pipeline_runtime(runtime, pipeline)
      x_analysis <- task$data[
        split$analysis_rows,
        task$predictors,
        drop = FALSE
      ]
      x_assessment <- task$data[
        split$assessment_rows,
        task$predictors,
        drop = FALSE
      ]
      analysis_endpoint <- list(
        task_type = "binary",
        outcome = outcome[split$analysis_rows],
        group_ids = sb_task_group_ids(task)[split$analysis_rows]
      )
      fitted_selector <- sb_fit_runtime_selector(
        runtime$selector,
        x_analysis,
        analysis_endpoint
      )
      selection <- list(
        fold = split$fold,
        selector = runtime$selector$id,
        status = "success",
        selected_features = fitted_selector$selected_features,
        scores = if (!is.null(fitted_selector$scores)) {
          fitted_selector$scores
        } else {
          NULL
        },
        error_class = NA_character_,
        error_message = NA_character_
      )
      selected_analysis <- runtime$selector$apply(
        fitted_selector,
        x_analysis
      )
      selected_assessment <- runtime$selector$apply(
        fitted_selector,
        x_assessment
      )
      prepared <- sb_prepare_runtime_learner(
        runtime$learner,
        pipeline,
        selected_analysis,
        analysis_endpoint,
        id[split$analysis_rows],
        group_ids = sb_task_group_ids(task, split$analysis_rows),
        selection_context = list(
          selector = runtime$selector,
          fitted_selector = fitted_selector,
          x = x_analysis
        )
      )
      learner <- prepared$learner
      if (!is.null(prepared$tuning)) {
        tuned <- prepared$tuning
        tuning <- list(
          fold = split$fold,
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
          error_class = NA_character_,
          error_message = NA_character_,
          standardize = tuned$standardize
        )
        if (!is.null(tuned$calibration)) tuning$calibration <- tuned$calibration
      }
      model <- sb_fit_reusable_runtime_learner(
        learner,
        selected_analysis,
        analysis_endpoint,
        id[split$analysis_rows],
        pipeline$selector
      )
      fit_warnings <- attr(model, "stormbreaker_fit_warnings", exact = TRUE)
      if (length(fit_warnings) > 0L) {
        tuning$warnings <- c(
          tuning$warnings, paste0("Learner refit: ", fit_warnings)
        )
      }
      if (is.function(learner$extract_parameters)) {
        extracted <- learner$extract_parameters(model)
        distribution_parameters <- c(
          list(
            fold = split$fold,
            learner = pipeline$learner,
            status = "success"
          ),
          extracted
        )
      } else if (is.function(learner$extract_tree_model)) {
        extracted <- learner$extract_tree_model(model)
        tree_model <- c(
          list(
            fold = split$fold,
            learner = pipeline$learner,
            status = "success"
          ),
          extracted
        )
      } else if (is.function(learner$extract_coefficients)) {
        coefficient_values <- learner$extract_coefficients(model)
        coefficients <- if (is_glmnet) {
          predictor_nonzero <- coefficient_values$nonzero &
            coefficient_values$term != "(Intercept)"
          list(
            fold = split$fold,
            learner = pipeline$learner,
            status = "success",
            lambda = tuned$lambda,
            alpha = tuned$alpha,
            n_nonzero = as.integer(sum(predictor_nonzero)),
            nonzero_features = coefficient_values$term[predictor_nonzero],
            values = coefficient_values
          )
        } else {
          list(
            fold = split$fold,
            learner = pipeline$learner,
            status = "success",
            values = coefficient_values
          )
        }
      }
      probability <- learner$predict(model, selected_assessment)
      successful_predictions <- empty_predictions
      successful_predictions$probability <- probability
      successful_predictions$status <- "success"

      list(
        predictions = successful_predictions,
        selection = selection,
        tuning = tuning,
        coefficients = coefficients,
        distribution_parameters = distribution_parameters,
        tree_model = tree_model,
        run = data.frame(
          fold = split$fold,
          seed = seed,
          pipeline_id = pipeline$id,
          selector_requested = pipeline$selector,
          learner_requested = pipeline$learner,
          learner_used = pipeline$learner,
          n_analysis = length(split$analysis_rows),
          n_analysis_positive = sum(outcome[split$analysis_rows]),
          n_assessment = length(split$assessment_rows),
          n_assessment_positive = sum(outcome[split$assessment_rows]),
          status = "success",
          error_class = NA_character_,
          error_message = NA_character_,
          stringsAsFactors = FALSE
        )
      )
    },
    error = function(error) {
      failed_selection <- selection
      if (identical(failed_selection$status, "failed")) {
        failed_selection$error_class <- class(error)[1]
        failed_selection$error_message <- conditionMessage(error)
      }
      failed_tuning <- tuning
      if ((is_glmnet || is_svm) && identical(failed_tuning$status, "failed")) {
        failed_tuning$error_class <- class(error)[1]
        failed_tuning$error_message <- conditionMessage(error)
      }
      list(
        predictions = empty_predictions,
        selection = failed_selection,
        tuning = failed_tuning,
        coefficients = coefficients,
        distribution_parameters = distribution_parameters,
        tree_model = tree_model,
        run = data.frame(
          fold = split$fold,
          seed = seed,
          pipeline_id = pipeline$id,
          selector_requested = pipeline$selector,
          learner_requested = pipeline$learner,
          learner_used = NA_character_,
          n_analysis = length(split$analysis_rows),
          n_analysis_positive = sum(outcome[split$analysis_rows]),
          n_assessment = length(split$assessment_rows),
          n_assessment_positive = sum(outcome[split$assessment_rows]),
          status = "failed",
          error_class = class(error)[1],
          error_message = conditionMessage(error),
          stringsAsFactors = FALSE
        )
      )
    }
  )
}

sb_run_binary_cv <- function(task, pipeline, folds, seed, progress = FALSE) {
  task <- sb_validate_binary_task_contract(task)
  if (is.null(pipeline)) {
    pipeline <- sb_pipeline(learner = "logistic")
  }
  pipeline <- sb_validate_pipeline_contract(pipeline, "binary")
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  split_object <- sb_make_outer_folds(task, folds, seed)
  fold_assignment <- data.frame(
    patient_id = as.character(task$data[[task$id_col]]),
    fold = split_object$assignment,
    stringsAsFactors = FALSE
  )
  fold_results <- lapply(split_object$splits, function(split) {
    sb_progress_message(progress, "Fold ", split$fold, "/", length(split_object$splits), " started")
    result <- sb_run_one_binary_fold(
      task,
      pipeline,
      runtime,
      split,
      split_object$seed
    )
    sb_progress_message(progress, "Fold ", split$fold, "/", length(split_object$splits),
      " complete: ", result$run$status)
    result
  })
  predictions <- do.call(
    rbind,
    lapply(fold_results, function(result) result$predictions)
  )
  runs <- do.call(rbind, lapply(fold_results, function(result) result$run))
  selections <- lapply(fold_results, function(result) result$selection)
  selection_stability <- sb_selection_stability(
    pipeline$selector,
    selections,
    task$predictors
  )
  tuning <- lapply(fold_results, function(result) result$tuning)
  coefficients <- lapply(
    fold_results,
    function(result) result$coefficients
  )
  glmnet_active_set_stability <- sb_glmnet_active_set_stability(
    pipeline,
    coefficients,
    task$predictors
  )
  distribution_parameters <- lapply(
    fold_results,
    function(result) result$distribution_parameters
  )
  tree_models <- lapply(
    fold_results,
    function(result) result$tree_model
  )
  rownames(predictions) <- NULL
  rownames(runs) <- NULL

  original_ids <- as.character(task$data[[task$id_col]])
  predictions <- predictions[
    match(original_ids, predictions$patient_id), ,
    drop = FALSE
  ]
  rownames(predictions) <- NULL
  exactly_once <- nrow(predictions) == task$n &&
    !anyNA(predictions$patient_id) &&
    anyDuplicated(predictions$patient_id) == 0L &&
    setequal(predictions$patient_id, original_ids)
  complete_oof <- exactly_once && all(predictions$status == "success") &&
    all(is.finite(predictions$probability)) &&
    all(predictions$probability >= 0 & predictions$probability <= 1)
  brier_reference <- data.frame(
    fold = runs$fold,
    analysis_prevalence = runs$n_analysis_positive / runs$n_analysis
  )
  moderate_calibration <- tryCatch(
    sb_binary_moderate_calibration_metric(
      predictions,
      scope = "pooled_patient_level_out_of_fold"
    ),
    error = function(error) {
      sb_binary_moderate_calibration_error_metric(
        error,
        scope = "pooled_patient_level_out_of_fold",
        n = nrow(predictions)
      )
    }
  )
  metrics <- list(
    roc_auc = sb_binary_auc_metric(predictions),
    average_precision = sb_binary_average_precision_metric(predictions),
    brier_score = sb_binary_brier_metric(predictions, brier_reference),
    log_loss = sb_binary_log_loss_metric(predictions, brier_reference),
    calibration = sb_binary_calibration_metric(predictions),
    moderate_calibration = moderate_calibration
  )
  if ("decision_thresholds" %in% names(task)) {
    metrics$decision_curve <- sb_binary_decision_curve_metric(
      predictions,
      task$decision_thresholds,
      scope = "pooled_patient_level_out_of_fold"
    )
  }

  result <- list(
    status = if (all(runs$status == "success") && complete_oof) {
      "success"
    } else {
      "failed"
    },
    complete_oof = complete_oof,
    predictions = predictions,
    metrics = metrics,
    runs = runs,
    selection_stability = selection_stability,
    coefficients = coefficients,
    folds = split_object$splits,
    task = list(
      n = task$n,
      n_positive = task$n_positive,
      n_negative = task$n_negative,
      id_col = task$id_col,
      outcome_col = task$outcome_col,
      positive_class = task$positive_class,
      negative_class = task$negative_class,
      predictors = task$predictors,
      independent_unit = task$independent_unit
    ),
    pipeline = pipeline,
    provenance = list(
      seed = split_object$seed,
      folds = split_object$n_folds,
      fold_source = split_object$source,
      fold_assignment = fold_assignment,
      package = "stormbreaker",
      package_version = as.character(utils::packageVersion("stormbreaker")),
      r_version = R.version.string
    )
  )
  if (!is.null(task$group_col)) {
    result$task$group_col <- task$group_col
    result$provenance$group_ids <- stats::setNames(sb_task_group_ids(task), original_ids)
  }
  if ("decision_thresholds" %in% names(task)) {
    result$task$decision_thresholds <- task$decision_thresholds
  }
  if (!identical(pipeline$selector, "none")) {
    result <- append(
      result,
      list(selections = selections),
      after = match("runs", names(result))
    )
  }
  if (identical(pipeline$learner, "gaussian_naive_bayes")) {
    result <- append(
      result,
      list(distribution_parameters = distribution_parameters),
      after = match("coefficients", names(result))
    )
  }
  if (pipeline$learner %in% c(
    sb_lightgbm_binary_learner_ids(),
    "catboost_binary", sb_ranger_binary_learner_ids(),
    sb_classical_binary_boosting_learner_ids()
  )) {
    result <- append(
      result,
      list(tree_models = tree_models),
      after = match("coefficients", names(result))
    )
    if (pipeline$learner %in% sb_lightgbm_binary_learner_ids()) {
      result$provenance$lightgbm_version <- if (sb_lightgbm_available()) {
        as.character(utils::packageVersion("lightgbm"))
      } else {
        NA_character_
      }
    } else if (identical(pipeline$learner, "catboost_binary")) {
      result$provenance$catboost_version <- if (sb_catboost_available()) {
        as.character(utils::packageVersion("catboost"))
      } else {
        NA_character_
      }
    } else if (pipeline$learner %in% sb_ranger_binary_learner_ids()) {
      result$provenance$ranger_version <- if (sb_ranger_available()) {
        as.character(utils::packageVersion("ranger"))
      } else {
        NA_character_
      }
    } else {
      dependency <- sb_classical_binary_boosting_dependency(pipeline$learner)
      result$provenance[[dependency$field]] <- if (dependency$available) {
        as.character(utils::packageVersion(dependency$package))
      } else {
        NA_character_
      }
    }
  }
  if (identical(pipeline$learner, "cancerclass_binary")) {
    result$provenance$cancerclass_version <- if (sb_cancerclass_available()) {
      as.character(utils::packageVersion("cancerclass"))
    } else {
      NA_character_
    }
    result$provenance$biobase_version <- if (sb_cancerclass_available()) {
      as.character(utils::packageVersion("Biobase"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% c("glmnet_logistic", sb_svm_probability_learner_ids()) ||
      any(vapply(tuning, function(x) length(x$warnings) > 0L, logical(1)))) {
    result <- append(
      result,
      list(tuning = tuning),
      after = match("runs", names(result))
    )
  }
  if (identical(pipeline$learner, "glmnet_logistic")) {
    result$provenance$glmnet_version <- if (sb_glmnet_available()) {
      as.character(utils::packageVersion("glmnet"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_classical_binary_learner_ids()) {
    dependency <- sb_classical_binary_dependency(pipeline$learner)
    if (!is.na(dependency$package)) {
      field <- paste0(tolower(dependency$package), "_version")
      result$provenance[[field]] <- if (dependency$available) {
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
      result$provenance[[field]] <- if (dependency$available) {
        as.character(utils::packageVersion(dependency$package))
      } else {
        NA_character_
      }
    }
  }
  if (pipeline$learner %in% sb_xgboost_binary_learner_ids() || pipeline$selector == "xgboost_binary_importance") {
    result$provenance$xgboost_version <- if (sb_xgboost_available()) {
      as.character(utils::packageVersion("xgboost"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_abess_learner_ids() ||
    pipeline$selector %in% sb_abess_selector_ids()) {
    result$provenance$abess_version <- if (sb_abess_available()) {
      as.character(utils::packageVersion("abess"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_ncvreg_binary_learner_ids()) {
    result$provenance$ncvreg_version <- if (sb_ncvreg_available()) {
      as.character(utils::packageVersion("ncvreg"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_rlt_binary_learner_ids()) {
    result$provenance$rlt_version <- if (sb_rlt_available()) {
      as.character(utils::packageVersion("RLT"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_wsrf_binary_learner_ids()) {
    result$provenance$wsrf_version <- if (sb_wsrf_available()) {
      as.character(utils::packageVersion("wsrf"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_slope_randomglm_binary_learner_ids()) {
    dependency <- sb_slope_randomglm_binary_dependency(pipeline$learner)
    field <- paste0(tolower(dependency$package), "_version")
    result$provenance[[field]] <- if (dependency$available) {
      as.character(utils::packageVersion(dependency$package))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_rule_aft_learner_ids()) {
    dependency <- sb_rule_aft_dependency(pipeline$learner)
    for (package in dependency$packages) {
      field <- paste0(tolower(package), "_version")
      result$provenance[[field]] <- if (
        requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        NA_character_
      }
    }
  }
  if (pipeline$learner %in% sb_mgcv_gam_learner_ids()) {
    result$provenance$mgcv_version <- if (sb_mgcv_available()) {
      as.character(utils::packageVersion("mgcv"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_mclust_da_learner_ids()) {
    result$provenance$mclust_version <- if (sb_mclust_available()) {
      as.character(utils::packageVersion("mclust"))
    } else {
      NA_character_
    }
  }
  if (pipeline$selector %in% sb_forest_selector_ids()) {
    selector_provenance <- sb_forest_selector_provenance(pipeline$selector)
    for (field in names(selector_provenance)) {
      result$provenance[[field]] <- selector_provenance[[field]]
    }
  }
  if (pipeline$selector == "plsrglm_binary_weights") {
    result$provenance$plsrglm_version <- if (requireNamespace("plsRglm", quietly = TRUE))
      as.character(utils::packageVersion("plsRglm")) else NA_character_
  }
  if (pipeline$selector %in% sb_pamr_selector_ids()) {
    result$provenance$pamr_version <- if (requireNamespace("pamr", quietly = TRUE))
      as.character(utils::packageVersion("pamr")) else NA_character_
  }
  high_dimensional_methods <- intersect(
    c(pipeline$selector, pipeline$learner),
    sb_high_dimensional_method_ids()
  )
  for (method in high_dimensional_methods) {
    method_provenance <- sb_high_dimensional_method_provenance(method)
    for (field in names(method_provenance)) {
      result$provenance[[field]] <- method_provenance[[field]]
    }
  }
  result$glmnet_active_set_stability <- glmnet_active_set_stability
  class(result) <- c(
    "stormbreaker_binary_cv_result",
    "stormbreaker_cv_result"
  )
  result
}

#' @export
print.stormbreaker_binary_cv_result <- function(x, ...) {
  cat("<stormbreaker_binary_cv_result:", x$status, ">\n")
  cat("Pipeline:", x$pipeline$id, "| Folds:", nrow(x$runs), "\n")
  fold_label <- if (identical(
    x$provenance$fold_source,
    "preassigned_patient_id"
  )) {
    "preassigned patient-ID folds"
  } else {
    "generated stratified folds"
  }
  cat("Outer resampling:", fold_label, "\n")
  cat(
    "Patients:", nrow(x$predictions),
    "| Successful folds:", sum(x$runs$status == "success"), "\n"
  )
  cat("Complete OOF:", x$complete_oof, "\n")
  metric <- x$metrics$roc_auc
  if (identical(metric$status, "success")) {
    cat(
      "Fold-stratified ROC AUC:",
      format(round(metric$estimate, 3), nsmall = 3),
      "| Positive-negative pairs:", metric$comparable_pairs, "\n"
    )
  } else {
    cat("ROC AUC: unavailable (incomplete)\n")
  }
  average_precision <- x$metrics$average_precision
  if (identical(average_precision$status, "success")) {
    cat(
      "Fold-stratified average precision:",
      format(round(average_precision$estimate, 3), nsmall = 3), "\n"
    )
  } else {
    cat("Fold-stratified average precision: unavailable (incomplete)\n")
  }
  brier <- x$metrics$brier_score
  if (identical(brier$status, "success")) {
    cat(
      "OOF Brier score:",
      format(round(brier$estimate, 3), nsmall = 3), "\n"
    )
  } else {
    cat("OOF Brier score: unavailable (incomplete)\n")
  }
  skill <- brier$skill_score
  if (identical(skill$status, "success")) {
    cat(
      "Brier skill vs analysis prevalence:",
      format(round(skill$estimate, 3), nsmall = 3), "\n"
    )
  } else {
    cat("Brier skill: unavailable (incomplete)\n")
  }
  log_loss <- x$metrics$log_loss
  if (identical(log_loss$status, "success")) {
    cat(
      "OOF log loss:",
      format(round(log_loss$estimate, 3), nsmall = 3), "\n"
    )
  } else {
    cat("OOF log loss: unavailable (incomplete)\n")
  }
  log_loss_skill <- log_loss$skill_score
  if (identical(log_loss_skill$status, "success")) {
    cat(
      "Log-loss skill vs analysis prevalence:",
      format(round(log_loss_skill$estimate, 3), nsmall = 3), "\n"
    )
  } else {
    cat("Log-loss skill: unavailable (incomplete)\n")
  }
  calibration <- x$metrics$calibration
  if (identical(calibration$status, "success")) {
    cat(
      "OOF calibration intercept:",
      format(round(calibration$calibration_intercept, 3), nsmall = 3),
      "| Slope:",
      format(round(calibration$calibration_slope, 3), nsmall = 3), "\n"
    )
  } else {
    cat("OOF calibration: unavailable (incomplete)\n")
  }
  moderate_calibration <- x$metrics$moderate_calibration
  if (identical(moderate_calibration$status, "success")) {
    cat(
      "OOF binary moderate calibration ICI:",
      format(round(moderate_calibration$ici, 3), nsmall = 3),
      "| E50:", format(round(moderate_calibration$e50, 3), nsmall = 3),
      "| E90:", format(round(moderate_calibration$e90, 3), nsmall = 3),
      "\n"
    )
  } else {
    cat("OOF binary moderate calibration: unavailable\n")
  }
  if ("decision_curve" %in% names(x$metrics)) {
    decision_curve <- x$metrics$decision_curve
    if (identical(decision_curve$status, "success")) {
      cat(
        "OOF binary decision curve:",
        decision_curve$n_thresholds,
        "declared thresholds\n"
      )
    } else {
      cat("OOF binary decision curve: unavailable\n")
    }
  }
  selection_stability <- x$selection_stability
  if (identical(selection_stability$status, "success")) {
    cat(
      "Selection stability (mean pairwise Jaccard):",
      format(
        round(selection_stability$mean_pairwise_jaccard, 3),
        nsmall = 3
      ), "\n"
    )
  } else if (identical(selection_stability$status, "incomplete")) {
    cat("Selection stability: unavailable (incomplete)\n")
  }
  active_set_stability <- x$glmnet_active_set_stability
  if (identical(active_set_stability$status, "success")) {
    cat(
      "Glmnet active-set stability (mean pairwise Jaccard):",
      format(
        round(active_set_stability$mean_pairwise_jaccard, 3),
        nsmall = 3
      ), "\n"
    )
  } else if (identical(active_set_stability$status, "incomplete")) {
    cat("Glmnet active-set stability: unavailable (incomplete)\n")
  } else if (identical(active_set_stability$status, "not_evaluable")) {
    cat("Glmnet active-set stability: not evaluable (empty union)\n")
  }
  invisible(x)
}
