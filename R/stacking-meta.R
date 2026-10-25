sb_stack_task_contract <- function(task, stack) {
  stack <- sb_upgrade_legacy_names(stack)
  if (inherits(task, "stormbreaker_binary_task")) {
    task <- sb_validate_binary_task_contract(task)
    task_type <- "binary"
  } else if (inherits(task, "stormbreaker_surv_task")) {
    task <- sb_validate_task_contract(task)
    task_type <- "survival"
  } else {
    sb_abort(
      "`task` must be created by `sb_surv_task()` or `sb_binary_task()`.",
      "stormbreaker_bad_task"
    )
  }
  stack <- sb_validate_stack_contract(stack, task_type)
  list(task = task, stack = stack, task_type = task_type)
}

sb_incomplete_stack_member_oof <- function(pipeline_id) {
  sb_abort(
    paste0(
      "Stack member `", pipeline_id,
      "` does not have one complete finite out-of-fold prediction per patient."
    ),
    "stormbreaker_incomplete_stack_member_oof"
  )
}

sb_stack_member_matrix <- function(comparison, task, stack) {
  contract <- sb_stack_task_contract(task, stack)
  task <- contract$task
  stack <- contract$stack
  task_type <- contract$task_type
  valid_comparison <- inherits(
    comparison,
    "stormbreaker_cv_comparison"
  ) && is.list(comparison) && identical(comparison$task_type, task_type) &&
    is.list(comparison$results) &&
    setequal(names(comparison$results), stack$pipeline_ids)
  if (!valid_comparison) {
    sb_abort(
      "The member comparison does not match this task and stack.",
      "stormbreaker_bad_stack_member_comparison"
    )
  }

  patient_ids <- as.character(task$data[[task$id_col]])
  value_name <- if (identical(task_type, "binary")) {
    "probability"
  } else {
    "risk"
  }
  columns <- lapply(stack$pipeline_ids, function(pipeline_id) {
    result <- comparison$results[[pipeline_id]]
    predictions <- if (is.list(result)) result$predictions else NULL
    valid_frame <- isTRUE(result$complete_oof) &&
      is.data.frame(predictions) && nrow(predictions) == task$n &&
      all(c("patient_id", "status", value_name) %in% names(predictions))
    if (!valid_frame) {
      sb_incomplete_stack_member_oof(pipeline_id)
    }
    prediction_ids <- as.character(predictions$patient_id)
    valid_ids <- !anyNA(prediction_ids) &&
      !anyDuplicated(prediction_ids) && setequal(prediction_ids, patient_ids)
    status <- predictions$status
    value <- predictions[[value_name]]
    valid_values <- is.character(status) && length(status) == task$n &&
      !anyNA(status) && all(status == "success") &&
      is.numeric(value) && length(value) == task$n &&
      !anyNA(value) && all(is.finite(value))
    if (identical(task_type, "binary")) {
      valid_values <- valid_values && all(value >= 0 & value <= 1)
    }
    if (!valid_ids || !valid_values) {
      sb_incomplete_stack_member_oof(pipeline_id)
    }
    unname(value[match(patient_ids, prediction_ids)])
  })
  names(columns) <- stack$pipeline_ids
  data.frame(
    columns,
    row.names = patient_ids,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

sb_validate_stack_member_matrix <- function(member_matrix, task, stack) {
  patient_ids <- as.character(task$data[[task$id_col]])
  valid <- is.data.frame(member_matrix) && nrow(member_matrix) == task$n &&
    identical(names(member_matrix), stack$pipeline_ids) &&
    identical(rownames(member_matrix), patient_ids) &&
    all(vapply(member_matrix, is.numeric, logical(1))) &&
    all(vapply(
      member_matrix,
      function(column) !anyNA(column) && all(is.finite(column)),
      logical(1)
    ))
  if (valid && identical(stack$task_type, "binary")) {
    values <- as.matrix(member_matrix)
    valid <- all(values >= 0 & values <= 1)
  }
  if (!valid) {
    sb_abort(
      paste0(
        "`member_matrix` must contain one ordered finite column for every ",
        "stack member and one ordered row for every task patient."
      ),
      "stormbreaker_bad_stack_member_matrix"
    )
  }
  member_matrix
}

sb_fit_stack_meta <- function(task, stack, member_matrix) {
  contract <- sb_stack_task_contract(task, stack)
  task <- contract$task
  stack <- contract$stack
  task_type <- contract$task_type
  member_matrix <- sb_validate_stack_member_matrix(
    member_matrix,
    task,
    stack
  )
  patient_ids <- as.character(task$data[[task$id_col]])

  if (identical(task_type, "binary")) {
    outcome <- as.integer(task$data[[task$outcome_col]])
    tuned <- sb_select_glmnet_logistic_lambda(
      member_matrix,
      outcome,
      inner_folds = stack$inner_folds,
      seed = stack$inner_seed,
      alpha = stack$meta_alpha,
      standardize = stack$meta_standardize,
      selection_rule = stack$meta_selection_rule,
      measure = stack$meta_measure
    )
    names(tuned$fold_id) <- patient_ids
    adapter <- sb_glmnet_logistic_learner_adapter(
      alpha = tuned$alpha,
      lambda = tuned$lambda,
      standardize = tuned$standardize
    )
    model <- adapter$fit(member_matrix, outcome)
    coefficients <- adapter$extract_coefficients(model)
    state <- adapter$extract_state(model)
  } else {
    time <- task$data[[task$time_col]]
    event <- as.integer(task$data[[task$event_col]])
    tuned <- sb_select_glmnet_cox_lambda(
      member_matrix,
      time,
      event,
      inner_folds = stack$inner_folds,
      seed = stack$inner_seed,
      alpha = stack$meta_alpha,
      standardize = stack$meta_standardize,
      selection_rule = stack$meta_selection_rule,
      measure = stack$meta_measure
    )
    names(tuned$fold_id) <- patient_ids
    adapter <- sb_glmnet_cox_learner_adapter(
      alpha = tuned$alpha,
      lambda = tuned$lambda,
      standardize = tuned$standardize
    )
    model <- adapter$fit(member_matrix, time, event)
    coefficients <- adapter$extract_coefficients(model)
    context <- list(x = member_matrix, time = time, event = event)
    if ("prediction_horizon" %in% names(task)) {
      context$prediction_horizon <- task$prediction_horizon
      context$time_unit <- task$time_unit
    }
    state <- adapter$extract_state(model, context)
  }

  structure(
    list(
      task_type = task_type,
      meta_learner = stack$meta_learner,
      tuning = tuned,
      coefficients = coefficients,
      state = state
    ),
    class = "stormbreaker_stack_meta_fit"
  )
}

sb_predict_stack_meta <- function(meta_state, member_matrix) {
  member_matrix <- sb_validate_predictor_frame(
    member_matrix,
    "member_matrix"
  )
  sb_predict_final_state(meta_state, member_matrix)
}
