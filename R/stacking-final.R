sb_stack_final_fit_failure <- function(result, stage, error) {
  result$status <- "failed"
  result$error_stage <- stage
  result$error_class <- class(error)[[1L]]
  result$error_message <- conditionMessage(error)
  class(result) <- "stormbreaker_stack_final_fit"
  result
}

sb_stack_final_prediction_semantics <- function(task, stack) {
  if (inherits(task, "stormbreaker_binary_task")) {
    semantics <- list(
      output = "probability",
      positive_class = 1L,
      direction = "higher_probability_of_class_1",
      meta_model = stack$meta_learner
    )
  } else {
    semantics <- list(
      output = "risk",
      risk_type = "stacked_ridge_cox_linear_predictor",
      risk_reference = "zero",
      risk_formula = "member_predictions_beta",
      risk_direction = "higher_is_higher_event_risk",
      meta_model = stack$meta_learner
    )
    if ("prediction_horizon" %in% names(task)) {
      semantics$output <- c("risk", "survival_probability")
      semantics$prediction_horizon <- task$prediction_horizon
      semantics$time_unit <- task$time_unit
    }
  }
  if ("decision_thresholds" %in% names(task)) {
    semantics$decision_thresholds <- task$decision_thresholds
  }
  semantics
}

sb_stack_member_fit_error <- function(pipeline_id, fit) {
  message <- paste0(
    "Stack member `", pipeline_id, "` failed during final fitting"
  )
  if (is.character(fit$error_message) && length(fit$error_message) == 1L &&
    !is.na(fit$error_message) && nzchar(fit$error_message)) {
    message <- paste0(message, ": ", fit$error_message)
  } else {
    message <- paste0(message, ".")
  }
  structure(
    list(message = message, call = NULL, parent = NULL),
    class = c(
      "stormbreaker_stack_member_fit_failure",
      "stormbreaker_error",
      "error",
      "condition"
    )
  )
}

#' Fit one leakage-safe stack on all development patients
#'
#' `sb_fit_stack_final()` first creates member out-of-fold predictions within
#' the complete development task and fits the stack's fixed ridge meta-model to
#' those predictions. It then fits every predeclared member pipeline once on all
#' development patients. The external cohort is not read and no member is
#' selected, dropped, or replaced.
#'
#' The returned fit stores portable member fits and portable meta coefficients.
#' For a survival task with a declared horizon, its Breslow baseline is learned
#' from development OOF stacked risks rather than external outcomes.
#'
#' @param task A task returned by [sb_surv_task()] or [sb_binary_task()].
#' @param stack A stack specification returned by [sb_stack()].
#'
#' @return A `stormbreaker_stack_final_fit` with stage-aware success or failure,
#'   inner resampling audit, meta fit, ordered member fits, required feature
#'   union, prediction semantics, and provenance.
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
#'   sb_fit_stack_final(task, stack)
#' }
sb_fit_stack_final <- function(task, stack) {
  contract <- sb_stack_task_contract(task, stack)
  task <- contract$task
  stack <- contract$stack
  task_type <- contract$task_type
  result <- list(
    status = "not_fitted",
    task_type = task_type,
    stack = stack,
    pipeline = list(id = stack$id, learner = stack$meta_learner),
    development = sb_final_fit_development(task),
    inner = list(
      status = "not_fitted",
      folds = stack$inner_folds,
      seed = stack$inner_seed,
      fold_source = NA_character_,
      fold_assignment = NULL,
      member_comparison = NULL
    ),
    meta = NULL,
    member_fits = list(),
    required_features = character(),
    prediction = sb_stack_final_prediction_semantics(task, stack),
    provenance = list(
      package = "stormbreaker",
      package_version = as.character(utils::packageVersion("stormbreaker")),
      r_version = R.version.string,
      survival_version = as.character(utils::packageVersion("survival")),
      glmnet_version = if (sb_glmnet_available()) {
        as.character(utils::packageVersion("glmnet"))
      } else {
        NA_character_
      }
    ),
    error_stage = NA_character_,
    error_class = NA_character_,
    error_message = NA_character_
  )

  comparison_attempt <- sb_final_fit_attempt(sb_compare_cv(
    task,
    stack$pipelines,
    folds = stack$inner_folds,
    seed = stack$inner_seed
  ))
  if (!is.null(comparison_attempt$error)) {
    return(sb_stack_final_fit_failure(
      result,
      "inner_resampling",
      comparison_attempt$error
    ))
  }
  comparison <- comparison_attempt$value
  result$inner <- list(
    status = "success",
    folds = comparison$folds,
    seed = comparison$seed,
    fold_source = comparison$fold_source,
    fold_assignment = comparison$fold_assignment,
    member_comparison = comparison$comparison
  )

  matrix_attempt <- sb_final_fit_attempt(
    sb_stack_member_matrix(comparison, task, stack)
  )
  if (!is.null(matrix_attempt$error)) {
    result$inner$status <- "failed"
    return(sb_stack_final_fit_failure(
      result,
      "inner_oof",
      matrix_attempt$error
    ))
  }
  member_matrix <- matrix_attempt$value

  meta_attempt <- sb_final_fit_attempt(
    sb_fit_stack_meta(task, stack, member_matrix)
  )
  if (!is.null(meta_attempt$error)) {
    return(sb_stack_final_fit_failure(
      result,
      "meta_fit",
      meta_attempt$error
    ))
  }
  result$meta <- meta_attempt$value

  member_error <- NULL
  sb_with_selector_reuse(for (index in seq_along(stack$pipelines)) {
    pipeline_id <- stack$pipeline_ids[[index]]
    member_fit <- sb_fit_final(task, stack$pipelines[[index]])
    result$member_fits[[pipeline_id]] <- member_fit
    if (!identical(member_fit$status, "success")) {
      member_error <- sb_stack_member_fit_error(pipeline_id, member_fit)
      break
    }
  })
  if (!is.null(member_error)) {
    return(sb_stack_final_fit_failure(result, "member_fit", member_error))
  }

  used_features <- unique(unlist(lapply(
    result$member_fits,
    function(fit) fit$required_features
  )))
  result$required_features <- task$predictors[
    task$predictors %in% used_features
  ]
  result$status <- "success"
  class(result) <- "stormbreaker_stack_final_fit"
  result
}

#' @export
print.stormbreaker_stack_final_fit <- function(x, ...) {
  cat("<stormbreaker_stack_final_fit:", x$status, ">\n")
  cat(
    "Task:", x$task_type,
    "| Members:", length(x$stack$pipelines),
    "| Development patients:", x$development$n, "\n"
  )
  cat("Meta-learner:", x$stack$meta_learner, "\n")
  if (identical(x$status, "failed")) {
    cat("Failed stage:", x$error_stage, "|", x$error_class, "\n")
  }
  invisible(x)
}

sb_validate_stack_final_fit_for_prediction <- function(fit) {
  required <- c(
    "status", "task_type", "stack", "pipeline", "development", "inner",
    "meta", "member_fits", "required_features", "prediction", "provenance"
  )
  valid <- inherits(fit, "stormbreaker_stack_final_fit") && is.list(fit) &&
    length(setdiff(required, names(fit))) == 0L &&
    identical(fit$status, "success")
  if (!valid) {
    sb_abort(
      "`fit` must be a successful result from `sb_fit_stack_final()`.",
      "stormbreaker_bad_stack_final_fit"
    )
  }
  stack <- sb_validate_stack_contract(fit$stack, fit$task_type)
  valid_members <- identical(names(fit$member_fits), stack$pipeline_ids) &&
    length(fit$member_fits) == length(stack$pipelines)
  if (!valid_members) {
    sb_abort(
      "The saved stack member fits are incomplete or reordered.",
      "stormbreaker_bad_stack_final_fit"
    )
  }
  for (member_fit in fit$member_fits) {
    sb_validate_final_fit_for_prediction(member_fit)
  }
  valid_meta <- inherits(fit$meta, "stormbreaker_stack_meta_fit") &&
    is.list(fit$meta$state) &&
    identical(fit$meta$state$features, stack$pipeline_ids)
  valid_features <- is.character(fit$required_features) &&
    !anyNA(fit$required_features) && !anyDuplicated(fit$required_features)
  valid_pipeline <- is.list(fit$pipeline) &&
    identical(fit$pipeline$id, stack$id) &&
    identical(fit$pipeline$learner, stack$meta_learner)
  if (!valid_meta || !valid_features || !valid_pipeline) {
    sb_abort(
      "The saved stack meta state or prediction metadata are invalid.",
      "stormbreaker_bad_stack_final_fit"
    )
  }
  fit
}

sb_stack_member_prediction_error <- function(pipeline_id, prediction) {
  detail <- if (is.list(prediction) &&
    is.character(prediction$error_message) &&
    length(prediction$error_message) == 1L &&
    !is.na(prediction$error_message) && nzchar(prediction$error_message)) {
    paste0(": ", prediction$error_message)
  } else {
    "."
  }
  sb_abort(
    paste0(
      "Stack member `", pipeline_id,
      "` failed during outcome-blind prediction", detail
    ),
    "stormbreaker_stack_member_prediction_failure"
  )
}

sb_predict_stack_members <- function(fit, new_data, id_col, patient_ids) {
  columns <- vector("list", length(fit$stack$pipelines))
  names(columns) <- fit$stack$pipeline_ids
  for (index in seq_along(fit$stack$pipelines)) {
    pipeline_id <- fit$stack$pipeline_ids[[index]]
    prediction <- sb_predict_final(
      fit$member_fits[[pipeline_id]],
      new_data,
      id_col
    )
    if (!is.list(prediction) || !identical(prediction$status, "success")) {
      sb_stack_member_prediction_error(pipeline_id, prediction)
    }
    contract <- sb_external_prediction_contract(prediction)
    if (!identical(contract$task_type, fit$task_type) ||
      !setequal(contract$patient_ids, patient_ids)) {
      sb_abort(
        paste0(
          "Stack member `", pipeline_id,
          "` returned misaligned patient predictions."
        ),
        "stormbreaker_stack_member_prediction_failure"
      )
    }
    value <- if (identical(fit$task_type, "binary")) {
      contract$probability
    } else {
      contract$risk
    }
    columns[[pipeline_id]] <- unname(value[
      match(patient_ids, contract$patient_ids)
    ])
  }
  data.frame(
    columns,
    row.names = patient_ids,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Predict new patients from a saved stack without reading outcomes
#'
#' `sb_predict_stack_final()` asks every saved member fit for its outcome-blind
#' prediction, preserves the declared member order, and applies the saved ridge
#' meta coefficients. It does not refit a member, the meta-model, or a survival
#' baseline. A failed member prediction makes the complete stack prediction
#' fail; no member is removed or replaced.
#'
#' The returned object uses the same final-prediction contract as
#' [sb_predict_final()] and can be passed directly to [sb_evaluate_external()].
#'
#' @param fit A successful result from [sb_fit_stack_final()].
#' @param new_data A non-empty data frame containing the declared ID and every
#'   numeric feature required by the saved member fits. Outcome columns are not
#'   required and, when present, are ignored.
#' @param id_col Name of the unique patient identifier in `new_data`.
#'
#' @return A `stormbreaker_final_prediction` with stacked patient predictions,
#'   the aligned member-prediction matrix, and stage-aware failure information.
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
#'   fit <- sb_fit_stack_final(task, stack)
#'   sb_predict_stack_final(
#'     fit,
#'     d[1:3, c("patient_id", "gene_a", "gene_b")],
#'     "patient_id"
#'   )
#' }
sb_predict_stack_final <- function(fit, new_data, id_col) {
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
    sb_validate_stack_final_fit_for_prediction(fit)
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

  member_attempt <- sb_final_fit_attempt(
    sb_predict_stack_members(fit, new_data, id_col, patient_ids)
  )
  if (!is.null(member_attempt$error)) {
    return(sb_final_prediction_failure(
      fit,
      id_col,
      patient_ids,
      "member_prediction",
      member_attempt$error
    ))
  }
  member_matrix <- member_attempt$value

  meta_attempt <- sb_final_fit_attempt({
    values <- sb_predict_stack_meta(fit$meta$state, member_matrix)
    sb_validate_final_prediction_values(fit, values, nrow(member_matrix))
  })
  if (!is.null(meta_attempt$error)) {
    return(sb_final_prediction_failure(
      fit,
      id_col,
      patient_ids,
      "meta_prediction",
      meta_attempt$error
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
      meta_attempt$value
    ),
    member_predictions = member_matrix,
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
