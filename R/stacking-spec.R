#' Define a leakage-safe stack of predeclared pipelines
#'
#' `sb_stack()` records an ordered set of existing STORMBREAKER pipelines for
#' stacked generalization. It does not fit models or inspect outcomes. Binary
#' stacks use a ridge logistic meta-learner; survival stacks use a ridge Cox
#' meta-learner. In both cases the ridge penalty is selected from member
#' out-of-fold predictions inside the current analysis data.
#'
#' Members are never searched, dropped, replaced, or reordered. At least two
#' distinct pipelines with one common task type are required.
#'
#' @param pipelines An ordered list containing at least two distinct pipeline
#'   specifications returned by [sb_pipeline()].
#' @param inner_folds Number of event- or outcome-stratified folds used to
#'   create member out-of-fold predictions and select the ridge penalty inside
#'   the current analysis data. Must be at least three.
#' @param inner_seed Positive integer seed for the inner fold assignment.
#'
#' @return A `stormbreaker_stack` specification. No model is fitted.
#' @export
#'
#' @examples
#' binary_stack <- sb_stack(list(
#'   sb_pipeline(learner = "logistic"),
#'   sb_pipeline(learner = "gaussian_naive_bayes")
#' ), inner_folds = 3, inner_seed = 2026)
#' binary_stack
sb_stack <- function(pipelines, inner_folds = 5L, inner_seed = 1L) {
  valid_collection <- is.list(pipelines) && length(pipelines) >= 2L &&
    all(vapply(
      pipelines,
      function(pipeline) inherits(pipeline, "stormbreaker_pipeline"),
      logical(1)
    ))
  if (!valid_collection) {
    sb_abort(
      paste0(
        "`pipelines` must contain at least two specifications returned by ",
        "`sb_pipeline()`."
      ),
      "stormbreaker_bad_stack_pipelines"
    )
  }

  task_types <- vapply(
    pipelines,
    function(pipeline) pipeline$task_type,
    character(1)
  )
  if (length(unique(task_types)) != 1L) {
    sb_abort(
      "Every stack member must have the same task type.",
      "stormbreaker_stack_task_mismatch"
    )
  }
  task_type <- task_types[[1L]]
  pipelines <- lapply(
    pipelines,
    sb_validate_pipeline_contract,
    task_type = task_type
  )
  pipeline_ids <- unname(vapply(
    pipelines,
    function(pipeline) pipeline$id,
    character(1)
  ))
  if (anyDuplicated(pipeline_ids) > 0L) {
    sb_abort(
      "Every stack member must have a distinct pipeline ID.",
      "stormbreaker_duplicate_stack_pipeline_id"
    )
  }

  inner_folds <- sb_check_positive_integer(
    inner_folds,
    "inner_folds",
    minimum = 3L
  )
  inner_seed <- sb_check_positive_integer(inner_seed, "inner_seed")
  meta_learner <- if (identical(task_type, "binary")) {
    "ridge_logistic"
  } else {
    "ridge_cox"
  }

  structure(
    list(
      id = paste0(
        "ridge_stack__",
        paste(pipeline_ids, collapse = "__")
      ),
      task_type = task_type,
      pipelines = pipelines,
      pipeline_ids = pipeline_ids,
      meta_learner = meta_learner,
      meta_alpha = 0,
      meta_standardize = TRUE,
      meta_selection_rule = "lambda.min",
      meta_measure = "deviance",
      inner_folds = inner_folds,
      inner_seed = inner_seed
    ),
    class = "stormbreaker_stack"
  )
}

sb_validate_stack_contract <- function(stack, task_type = NULL) {
  required <- c(
    "id", "task_type", "pipelines", "pipeline_ids", "meta_learner",
    "meta_alpha", "meta_standardize", "meta_selection_rule",
    "meta_measure", "inner_folds", "inner_seed"
  )
  if (!inherits(stack, "stormbreaker_stack") || !is.list(stack) ||
    length(setdiff(required, names(stack))) > 0L) {
    sb_abort(
      "`stack` must be created by `sb_stack()`.",
      "stormbreaker_bad_stack"
    )
  }

  if (!is.null(task_type) && !identical(stack$task_type, task_type)) {
    sb_abort(
      "The stack does not match the supplied task type.",
      "stormbreaker_stack_task_mismatch"
    )
  }
  stack
}

sb_subset_task <- function(task, rows) {
  valid_rows <- is.numeric(rows) && length(rows) >= 2L &&
    !anyNA(rows) && all(is.finite(rows)) &&
    all(rows == as.integer(rows)) && !anyDuplicated(rows)
  if (!valid_rows) {
    sb_abort(
      "Task row indices must be distinct finite whole numbers.",
      "stormbreaker_bad_task_rows"
    )
  }
  rows <- as.integer(rows)

  if (inherits(task, "stormbreaker_binary_task")) {
    task <- sb_validate_binary_task_contract(task)
    if (any(rows < 1L | rows > task$n)) {
      sb_abort(
        "Task row indices fall outside the available patients.",
        "stormbreaker_bad_task_rows"
      )
    }
    thresholds <- if ("decision_thresholds" %in% names(task)) {
      task$decision_thresholds
    } else {
      NULL
    }
    return(sb_binary_task(
      task$data[rows, , drop = FALSE],
      task$id_col,
      task$outcome_col,
      task$predictors,
      decision_thresholds = thresholds,
      group_col = task$group_col
    ))
  }

  if (inherits(task, "stormbreaker_surv_task")) {
    task <- sb_validate_task_contract(task)
    if (any(rows < 1L | rows > task$n)) {
      sb_abort(
        "Task row indices fall outside the available patients.",
        "stormbreaker_bad_task_rows"
      )
    }
    horizon <- if ("prediction_horizon" %in% names(task)) {
      task$prediction_horizon
    } else {
      NULL
    }
    thresholds <- if ("decision_thresholds" %in% names(task)) {
      task$decision_thresholds
    } else {
      NULL
    }
    return(sb_surv_task(
      task$data[rows, , drop = FALSE],
      task$id_col,
      task$time_col,
      task$event_col,
      task$predictors,
      task$time_unit,
      prediction_horizon = horizon,
      decision_thresholds = thresholds,
      group_col = task$group_col
    ))
  }

  sb_abort(
    "`task` must be created by `sb_surv_task()` or `sb_binary_task()`.",
    "stormbreaker_bad_task"
  )
}

#' @export
print.stormbreaker_stack <- function(x, ...) {
  cat("<stormbreaker_stack>\n")
  cat("Task:", x$task_type, "| Members:", length(x$pipelines), "\n")
  cat("Meta-learner:", x$meta_learner, "| Inner folds:", x$inner_folds, "\n")
  invisible(x)
}
