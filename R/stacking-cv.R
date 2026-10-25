sb_stack_outer_member_matrix <- function(comparison, patient_ids, stack) {
  task_type <- stack$task_type
  value_name <- if (identical(task_type, "binary")) {
    "probability"
  } else {
    "risk"
  }
  columns <- lapply(stack$pipeline_ids, function(pipeline_id) {
    result <- comparison$results[[pipeline_id]]
    predictions <- if (is.list(result)) result$predictions else NULL
    valid_frame <- is.data.frame(predictions) &&
      all(c("patient_id", "status", value_name) %in% names(predictions))
    if (!valid_frame) {
      sb_incomplete_stack_member_oof(pipeline_id)
    }
    prediction_ids <- as.character(predictions$patient_id)
    positions <- match(patient_ids, prediction_ids)
    status <- predictions$status[positions]
    value <- predictions[[value_name]][positions]
    valid <- !anyNA(positions) && !anyDuplicated(prediction_ids) &&
      is.character(status) && !anyNA(status) && all(status == "success") &&
      is.numeric(value) && !anyNA(value) && all(is.finite(value))
    if (identical(task_type, "binary")) {
      valid <- valid && all(value >= 0 & value <= 1)
    }
    if (!valid) {
      sb_incomplete_stack_member_oof(pipeline_id)
    }
    unname(value)
  })
  names(columns) <- stack$pipeline_ids
  data.frame(
    columns,
    row.names = patient_ids,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

sb_stack_fold_prediction_rows <- function(task, stack, split, values = NULL) {
  patient_ids <- as.character(task$data[[task$id_col]][split$assessment_rows])
  success <- !is.null(values)
  if (inherits(task, "stormbreaker_binary_task")) {
    return(data.frame(
      patient_id = patient_ids,
      fold = split$fold,
      outcome = as.integer(task$data[[task$outcome_col]][split$assessment_rows]),
      positive_class = task$positive_class,
      probability = if (success) values$probability else NA_real_,
      status = if (success) "success" else "failed",
      pipeline_id = stack$id,
      learner = stack$meta_learner,
      stringsAsFactors = FALSE
    ))
  }

  rows <- data.frame(
    patient_id = patient_ids,
    fold = split$fold,
    time = task$data[[task$time_col]][split$assessment_rows],
    event = as.integer(task$data[[task$event_col]][split$assessment_rows]),
    time_unit = task$time_unit,
    risk = if (success) values$risk else NA_real_,
    status = if (success) "success" else "failed",
    pipeline_id = stack$id,
    learner = stack$meta_learner,
    stringsAsFactors = FALSE
  )
  if ("prediction_horizon" %in% names(task)) {
    rows$prediction_horizon <- task$prediction_horizon
    rows$survival_probability <- if (success) {
      values$survival_probability
    } else {
      NA_real_
    }
    rows$survival_probability_status <- if (success) "success" else "failed"
  }
  rows
}

sb_run_one_stack_fold <- function(task, stack, outer_members, split,
                                  outer_seed) {
  inner_comparison <- NULL
  member_matrix <- NULL
  meta <- NULL
  tryCatch(
    {
      analysis_task <- sb_subset_task(task, split$analysis_rows)
      inner_comparison <- sb_compare_cv(
        analysis_task,
        stack$pipelines,
        folds = stack$inner_folds,
        seed = stack$inner_seed
      )
      member_matrix <- sb_stack_member_matrix(
        inner_comparison,
        analysis_task,
        stack
      )
      meta <- sb_fit_stack_meta(analysis_task, stack, member_matrix)
      assessment_matrix <- sb_stack_outer_member_matrix(
        outer_members,
        split$assessment_ids,
        stack
      )
      values <- sb_predict_stack_meta(meta$state, assessment_matrix)
      predictions <- sb_stack_fold_prediction_rows(
        task,
        stack,
        split,
        values
      )
      run <- data.frame(
        fold = split$fold,
        seed = outer_seed,
        stack_id = stack$id,
        meta_learner = stack$meta_learner,
        n_members = length(stack$pipelines),
        n_analysis = length(split$analysis_rows),
        n_assessment = length(split$assessment_rows),
        inner_folds = stack$inner_folds,
        inner_seed = stack$inner_seed,
        meta_lambda = meta$tuning$lambda,
        status = "success",
        error_class = NA_character_,
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
      if (inherits(task, "stormbreaker_binary_task")) {
        run$n_analysis_positive <- sum(
          task$data[[task$outcome_col]][split$analysis_rows]
        )
        run <- run[c(
          "fold", "seed", "stack_id", "meta_learner", "n_members",
          "n_analysis", "n_analysis_positive", "n_assessment",
          "inner_folds", "inner_seed", "meta_lambda", "status",
          "error_class", "error_message"
        )]
      }
      list(
        predictions = predictions,
        run = run,
        audit = list(
          fold = split$fold,
          status = "success",
          comparison = inner_comparison,
          member_matrix = member_matrix,
          meta = meta,
          error_class = NA_character_,
          error_message = NA_character_
        )
      )
    },
    error = function(error) {
      predictions <- sb_stack_fold_prediction_rows(task, stack, split)
      run <- data.frame(
        fold = split$fold,
        seed = outer_seed,
        stack_id = stack$id,
        meta_learner = stack$meta_learner,
        n_members = length(stack$pipelines),
        n_analysis = length(split$analysis_rows),
        n_assessment = length(split$assessment_rows),
        inner_folds = stack$inner_folds,
        inner_seed = stack$inner_seed,
        meta_lambda = NA_real_,
        status = "failed",
        error_class = class(error)[[1L]],
        error_message = conditionMessage(error),
        stringsAsFactors = FALSE
      )
      if (inherits(task, "stormbreaker_binary_task")) {
        run$n_analysis_positive <- sum(
          task$data[[task$outcome_col]][split$analysis_rows]
        )
        run <- run[c(
          "fold", "seed", "stack_id", "meta_learner", "n_members",
          "n_analysis", "n_analysis_positive", "n_assessment",
          "inner_folds", "inner_seed", "meta_lambda", "status",
          "error_class", "error_message"
        )]
      }
      list(
        predictions = predictions,
        run = run,
        audit = list(
          fold = split$fold,
          status = "failed",
          comparison = inner_comparison,
          member_matrix = member_matrix,
          meta = meta,
          error_class = class(error)[[1L]],
          error_message = conditionMessage(error)
        )
      )
    }
  )
}

sb_stack_binary_metrics <- function(task, predictions, runs) {
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
  metrics
}

sb_stack_survival_metrics <- function(task, predictions, splits) {
  metrics <- list(harrell_c = sb_harrell_c_metric(predictions))
  if (!"prediction_horizon" %in% names(task)) {
    return(metrics)
  }
  analysis_inputs <- lapply(splits, function(split) {
    list(
      fold = split$fold,
      time = task$data[[task$time_col]][split$analysis_rows],
      event = as.integer(task$data[[task$event_col]][split$analysis_rows])
    )
  })
  horizon <- task$prediction_horizon
  metrics$ipcw_brier <- sb_survival_ipcw_brier_metric(
    predictions,
    analysis_inputs,
    horizon
  )
  metrics$time_auc <- sb_survival_time_auc_metric(
    predictions,
    analysis_inputs,
    horizon,
    metrics$ipcw_brier$censoring_curves
  )
  metrics$uno_c <- sb_survival_uno_c_metric(
    predictions,
    analysis_inputs,
    horizon,
    metrics$ipcw_brier$censoring_curves
  )
  calibration_status <- ifelse(
    predictions$status == "success" &
      predictions$survival_probability_status == "success",
    "success",
    "failed"
  )
  metrics$calibration <- tryCatch(
    sb_survival_calibration_metric(
      patient_id = predictions$patient_id,
      time = predictions$time,
      event = predictions$event,
      survival_probability = predictions$survival_probability,
      prediction_status = calibration_status,
      prediction_horizon = horizon,
      time_unit = task$time_unit,
      scope = "pooled_patient_level_out_of_fold",
      fold = predictions$fold
    ),
    error = function(error) {
      sb_survival_calibration_error_metric(
        error = error,
        prediction_horizon = horizon,
        time_unit = task$time_unit,
        scope = "pooled_patient_level_out_of_fold",
        n = nrow(predictions),
        fold = predictions$fold
      )
    }
  )
  metrics$weak_calibration <- tryCatch(
    sb_survival_weak_calibration_metric(
      patient_id = predictions$patient_id,
      time = predictions$time,
      event = predictions$event,
      survival_probability = predictions$survival_probability,
      prediction_status = calibration_status,
      prediction_horizon = horizon,
      time_unit = task$time_unit,
      scope = "pooled_patient_level_out_of_fold",
      fold = predictions$fold
    ),
    error = function(error) {
      sb_survival_weak_calibration_error_metric(
        error = error,
        prediction_horizon = horizon,
        time_unit = task$time_unit,
        scope = "pooled_patient_level_out_of_fold",
        n = nrow(predictions),
        fold = predictions$fold
      )
    }
  )
  if ("decision_thresholds" %in% names(task)) {
    metrics$decision_curve <- sb_survival_decision_curve_metric(
      predictions,
      task$decision_thresholds,
      horizon,
      task$time_unit,
      scope = "pooled_patient_level_out_of_fold"
    )
  }
  metrics
}

sb_stack_comparison_row <- function(result) {
  n_folds <- as.integer(nrow(result$runs))
  n_successful <- as.integer(sum(result$runs$status == "success"))
  prediction_variation <- sb_prediction_variation(
    result$predictions, result$task_type
  )$status
  if (identical(result$task_type, "binary")) {
    metrics <- result$metrics
    return(data.frame(
      pipeline_id = result$stack$id,
      selector = "stack_members",
      learner = result$stack$meta_learner,
      status = result$status,
      complete_oof = result$complete_oof,
      prediction_variation = prediction_variation,
      roc_auc_status = metrics$roc_auc$status,
      roc_auc = metrics$roc_auc$estimate,
      comparable_pairs = metrics$roc_auc$comparable_pairs,
      average_precision_status = metrics$average_precision$status,
      average_precision = metrics$average_precision$estimate,
      brier_status = metrics$brier_score$status,
      brier_score = metrics$brier_score$estimate,
      n_scored = metrics$brier_score$n_scored,
      reference_brier_status = metrics$brier_score$reference$status,
      reference_brier = metrics$brier_score$reference$estimate,
      brier_skill_status = metrics$brier_score$skill_score$status,
      brier_skill_score = metrics$brier_score$skill_score$estimate,
      log_loss_status = metrics$log_loss$status,
      log_loss = metrics$log_loss$estimate,
      reference_log_loss_status = metrics$log_loss$reference$status,
      reference_log_loss = metrics$log_loss$reference$estimate,
      log_loss_skill_status = metrics$log_loss$skill_score$status,
      log_loss_skill_score = metrics$log_loss$skill_score$estimate,
      calibration_status = metrics$calibration$status,
      calibration_intercept = metrics$calibration$calibration_intercept,
      calibration_slope = metrics$calibration$calibration_slope,
      moderate_calibration_status = metrics$moderate_calibration$status,
      moderate_calibration_ici = metrics$moderate_calibration$ici,
      moderate_calibration_e50 = metrics$moderate_calibration$e50,
      moderate_calibration_e90 = metrics$moderate_calibration$e90,
      selection_stability_status = "not_applicable",
      mean_pairwise_jaccard = NA_real_,
      n_successful_selection_folds = 0L,
      n_folds = n_folds,
      n_successful_folds = n_successful,
      glmnet_active_set_stability_status = "not_applicable",
      glmnet_active_set_mean_pairwise_jaccard = NA_real_,
      n_successful_glmnet_active_set_folds = 0L,
      stringsAsFactors = FALSE
    ))
  }

  metrics <- result$metrics
  row <- data.frame(
    pipeline_id = result$stack$id,
    selector = "stack_members",
    learner = result$stack$meta_learner,
    status = result$status,
    complete_oof = result$complete_oof,
    prediction_variation = prediction_variation,
    harrell_c_status = metrics$harrell_c$status,
    harrell_c = metrics$harrell_c$estimate,
    comparable_pairs = metrics$harrell_c$comparable_pairs,
    selection_stability_status = "not_applicable",
    mean_pairwise_jaccard = NA_real_,
    n_successful_selection_folds = 0L,
    n_folds = n_folds,
    n_successful_folds = n_successful,
    glmnet_active_set_stability_status = "not_applicable",
    glmnet_active_set_mean_pairwise_jaccard = NA_real_,
    n_successful_glmnet_active_set_folds = 0L,
    stringsAsFactors = FALSE
  )
  if (!"ipcw_brier" %in% names(metrics)) {
    return(row)
  }
  data.frame(
    pipeline_id = row$pipeline_id,
    selector = row$selector,
    learner = row$learner,
    status = row$status,
    complete_oof = row$complete_oof,
    prediction_variation = row$prediction_variation,
    complete_survival_probability_oof =
      result$complete_survival_probability_oof,
    prediction_horizon = metrics$ipcw_brier$prediction_horizon,
    harrell_c_status = row$harrell_c_status,
    harrell_c = row$harrell_c,
    comparable_pairs = row$comparable_pairs,
    ipcw_brier_status = metrics$ipcw_brier$status,
    ipcw_brier = metrics$ipcw_brier$estimate,
    reference_ipcw_brier_status = metrics$ipcw_brier$reference$status,
    reference_ipcw_brier = metrics$ipcw_brier$reference$estimate,
    ipcw_brier_skill_status = metrics$ipcw_brier$skill_score$status,
    ipcw_brier_skill_score = metrics$ipcw_brier$skill_score$estimate,
    time_auc_status = metrics$time_auc$status,
    time_auc = metrics$time_auc$estimate,
    time_auc_pair_denominator = metrics$time_auc$pair_denominator,
    uno_c_status = metrics$uno_c$status,
    uno_c = metrics$uno_c$estimate,
    uno_c_weighted_pair_denominator = metrics$uno_c$weighted_pair_denominator,
    calibration_status = metrics$calibration$status,
    calibration_ici = metrics$calibration$ici,
    calibration_e50 = metrics$calibration$e50,
    calibration_e90 = metrics$calibration$e90,
    mean_calibration_status = metrics$weak_calibration$mean_calibration$status,
    mean_calibration_ratio = metrics$weak_calibration$mean_calibration_ratio,
    calibration_slope_status = metrics$weak_calibration$slope_fit$status,
    calibration_slope = metrics$weak_calibration$calibration_slope,
    selection_stability_status = row$selection_stability_status,
    mean_pairwise_jaccard = row$mean_pairwise_jaccard,
    n_successful_selection_folds = row$n_successful_selection_folds,
    n_folds = row$n_folds,
    n_successful_folds = row$n_successful_folds,
    glmnet_active_set_stability_status =
      row$glmnet_active_set_stability_status,
    glmnet_active_set_mean_pairwise_jaccard =
      row$glmnet_active_set_mean_pairwise_jaccard,
    n_successful_glmnet_active_set_folds =
      row$n_successful_glmnet_active_set_folds,
    stringsAsFactors = FALSE
  )
}

#' Evaluate one predeclared stack with leakage-safe nested resampling
#'
#' `sb_run_stack_cv()` compares every stack member on one shared outer fold
#' assignment and evaluates their ridge stack on those same held-out patients.
#' Inside each outer analysis fold, members are rerun under an inner
#' event- or outcome-stratified assignment. Only those inner out-of-fold
#' predictions and analysis outcomes may train or tune the meta-learner. Outer
#' assessment outcomes are used only after stacked predictions are fixed.
#'
#' A member with incomplete inner or outer predictions causes the corresponding
#' stack fold to fail. The function never removes or replaces a member. The
#' returned comparison is descriptive and does not select a winner.
#'
#' @param task A task returned by [sb_surv_task()] or [sb_binary_task()].
#' @param stack A stack specification returned by [sb_stack()].
#' @param folds Either one outer-fold count or an exact patient-ID-named fold
#'   assignment accepted by [sb_run_cv()].
#' @param seed Positive integer used when `folds` is a scalar count.
#'
#' @return A `stormbreaker_stack_cv_result` with member results, stacked OOF
#'   predictions and metrics, per-outer-fold inner OOF/meta audits, and one
#'   descriptive member-plus-stack comparison table.
#' @export
#'
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   d <- sb_example_binary(n = 48, seed = 1)
#'   task <- sb_binary_task(
#'     d, "patient_id", "outcome", c("gene_a", "gene_b")
#'   )
#'   stack <- sb_stack(list(
#'     sb_pipeline(learner = "logistic"),
#'     sb_pipeline(learner = "gaussian_naive_bayes")
#'   ), inner_folds = 3)
#'   sb_run_stack_cv(task, stack, folds = 3, seed = 2026)
#' }
sb_run_stack_cv <- function(task, stack, folds = 5L, seed = 1L) {
  contract <- sb_stack_task_contract(task, stack)
  task <- contract$task
  stack <- contract$stack
  task_type <- contract$task_type
  split_object <- sb_make_outer_folds(task, folds, seed)
  outer_members <- sb_compare_cv(
    task,
    stack$pipelines,
    folds = folds,
    seed = seed
  )
  fold_results <- lapply(split_object$splits, function(split) {
    sb_run_one_stack_fold(
      task,
      stack,
      outer_members,
      split,
      split_object$seed
    )
  })
  predictions <- do.call(
    rbind,
    lapply(fold_results, function(result) result$predictions)
  )
  runs <- do.call(rbind, lapply(fold_results, function(result) result$run))
  inner_folds <- lapply(fold_results, function(result) result$audit)
  rownames(predictions) <- NULL
  rownames(runs) <- NULL
  patient_ids <- as.character(task$data[[task$id_col]])
  predictions <- predictions[
    match(patient_ids, predictions$patient_id),
    ,
    drop = FALSE
  ]
  rownames(predictions) <- NULL
  exactly_once <- nrow(predictions) == task$n &&
    !anyNA(predictions$patient_id) &&
    anyDuplicated(predictions$patient_id) == 0L &&
    setequal(predictions$patient_id, patient_ids)
  value_name <- if (identical(task_type, "binary")) "probability" else "risk"
  complete_oof <- exactly_once && all(predictions$status == "success") &&
    all(is.finite(predictions[[value_name]]))
  if (identical(task_type, "binary")) {
    complete_oof <- complete_oof &&
      all(predictions$probability >= 0 & predictions$probability <= 1)
    metrics <- sb_stack_binary_metrics(task, predictions, runs)
    requested_metrics_complete <- TRUE
    complete_probability <- NULL
  } else {
    complete_probability <- if ("prediction_horizon" %in% names(task)) {
      exactly_once &&
        all(predictions$survival_probability_status == "success") &&
        all(is.finite(predictions$survival_probability)) &&
        all(
          predictions$survival_probability >= 0 &
            predictions$survival_probability <= 1
        )
    } else {
      TRUE
    }
    metrics <- sb_stack_survival_metrics(
      task,
      predictions,
      split_object$splits
    )
    decision_curve_complete <-
      !("decision_thresholds" %in% names(task)) ||
        identical(metrics$decision_curve$status, "success")
    requested_metrics_complete <-
      !("prediction_horizon" %in% names(task)) ||
        (
          identical(metrics$ipcw_brier$status, "success") &&
            identical(metrics$time_auc$status, "success") &&
            identical(metrics$uno_c$status, "success") &&
            decision_curve_complete
        )
  }
  status <- if (all(runs$status == "success") && complete_oof &&
    (is.null(complete_probability) || complete_probability) &&
    requested_metrics_complete) {
    "success"
  } else {
    "failed"
  }
  fold_assignment <- data.frame(
    patient_id = patient_ids,
    fold = split_object$assignment,
    stringsAsFactors = FALSE
  )
  result <- list(
    status = status,
    task_type = task_type,
    complete_oof = complete_oof,
    predictions = predictions,
    metrics = metrics,
    runs = runs,
    folds = split_object$splits,
    inner_folds = inner_folds,
    stack = stack,
    members = outer_members,
    task = list(
      n = task$n,
      id_col = task$id_col,
      predictors = task$predictors,
      independent_unit = task$independent_unit
    ),
    provenance = list(
      seed = split_object$seed,
      folds = split_object$n_folds,
      fold_source = split_object$source,
      fold_assignment = fold_assignment,
      inner_folds = stack$inner_folds,
      inner_seed = stack$inner_seed,
      package = "stormbreaker",
      package_version = as.character(utils::packageVersion("stormbreaker")),
      r_version = R.version.string,
      glmnet_version = if (sb_glmnet_available()) {
        as.character(utils::packageVersion("glmnet"))
      } else {
        NA_character_
      }
    )
  )
  if (!is.null(complete_probability)) {
    result <- append(
      result,
      list(complete_survival_probability_oof = complete_probability),
      after = match("complete_oof", names(result))
    )
  }
  class(result) <- "stormbreaker_stack_cv_result"
  stack_row <- sb_stack_comparison_row(result)
  result$comparison <- rbind(
    outer_members$comparison,
    stack_row[names(outer_members$comparison)]
  )
  rownames(result$comparison) <- NULL
  if ("decision_thresholds" %in% names(task)) {
    result$decision_curves <- c(
      outer_members$decision_curves,
      stats::setNames(list(metrics$decision_curve), stack$id)
    )
  }
  result
}

#' @export
print.stormbreaker_stack_cv_result <- function(x, ...) {
  cat("<stormbreaker_stack_cv_result:", x$status, ">\n")
  cat(
    "Members:", length(x$stack$pipelines),
    "| Patients:", nrow(x$predictions),
    "| Outer folds:", nrow(x$runs), "\n"
  )
  cat("Complete stacked OOF:", x$complete_oof, "\n")
  if (identical(x$task_type, "binary")) {
    cat("ROC AUC:", format(round(x$metrics$roc_auc$estimate, 3)), "\n")
  } else {
    cat(
      "Harrell C:",
      format(round(x$metrics$harrell_c$estimate, 3)),
      "\n"
    )
  }
  cat("No member or stack is selected by this result.\n")
  invisible(x)
}
