sb_iecv_optional_common <- function(tasks, field, description) {
  present <- vapply(tasks, function(task) field %in% names(task), logical(1L))
  if (!all(present == present[[1L]])) {
    sb_abort(
      paste0("All cohorts must use the same ", description, "."),
      "stormbreaker_internal_external_target_mismatch"
    )
  }
  if (!present[[1L]]) {
    return(invisible(NULL))
  }
  reference <- tasks[[1L]][[field]]
  same <- vapply(tasks, function(task) {
    identical(task[[field]], reference)
  }, logical(1L))
  if (!all(same)) {
    sb_abort(
      paste0("All cohorts must use the same ", description, "."),
      "stormbreaker_internal_external_target_mismatch"
    )
  }
  invisible(reference)
}

sb_validate_iecv_cohorts <- function(cohorts) {
  valid_collection <- is.list(cohorts) && length(cohorts) >= 2L &&
    !inherits(cohorts, c(
      "stormbreaker_binary_task", "stormbreaker_surv_task"
    ))
  cohort_ids <- names(cohorts)
  valid_ids <- !is.null(cohort_ids) &&
    length(cohort_ids) == length(cohorts) &&
    !anyNA(cohort_ids) && all(nzchar(trimws(cohort_ids))) &&
    anyDuplicated(cohort_ids) == 0L
  if (!valid_collection || !valid_ids) {
    sb_abort(
      paste0(
        "`cohorts` must be a named list of at least two tasks with unique, ",
        "non-blank cohort IDs."
      ),
      "stormbreaker_bad_internal_external_cohorts"
    )
  }

  task_types <- vapply(cohorts, function(task) {
    if (inherits(task, "stormbreaker_binary_task")) {
      return("binary")
    }
    if (inherits(task, "stormbreaker_surv_task")) {
      return("survival")
    }
    sb_abort(
      paste0(
        "Every cohort must be created by `sb_binary_task()` or ",
        "`sb_surv_task()`."
      ),
      "stormbreaker_bad_internal_external_cohorts"
    )
  }, character(1L))
  if (length(unique(task_types)) != 1L) {
    sb_abort(
      "All internal-external validation cohorts must use the same task type.",
      "stormbreaker_internal_external_task_mismatch"
    )
  }

  task_type <- task_types[[1L]]
  tasks <- lapply(cohorts, function(task) {
    if (identical(task_type, "binary")) {
      sb_validate_binary_task_contract(task)
    } else {
      sb_validate_task_contract(task)
    }
  })
  names(tasks) <- cohort_ids

  reference_predictors <- tasks[[1L]]$predictors
  same_predictors <- vapply(tasks, function(task) {
    length(task$predictors) == length(reference_predictors) &&
      setequal(task$predictors, reference_predictors)
  }, logical(1L))
  if (!all(same_predictors)) {
    sb_abort(
      paste0(
        "All cohorts must declare the same predictor-name set; ",
        "STORMBREAKER does not intersect or harmonize features."
      ),
      "stormbreaker_internal_external_predictor_mismatch"
    )
  }

  if (identical(task_type, "survival")) {
    time_units <- vapply(tasks, `[[`, character(1L), "time_unit")
    if (length(unique(time_units)) != 1L) {
      sb_abort(
        "All survival cohorts must use the same time unit.",
        "stormbreaker_internal_external_target_mismatch"
      )
    }
    sb_iecv_optional_common(
      tasks, "prediction_horizon", "prediction-horizon declaration"
    )
  }
  sb_iecv_optional_common(
    tasks, "decision_thresholds", "decision-threshold declaration"
  )
  grouped <- vapply(tasks, function(task) !is.null(task$group_col), logical(1))
  if (any(grouped) && !all(grouped)) {
    sb_abort(
      "All cohorts must declare grouping when any cohort declares grouping.",
      "stormbreaker_internal_external_target_mismatch"
    )
  }
  if (all(grouped)) {
    groups <- unlist(lapply(tasks, function(task) {
      unique(sb_task_group_ids(task))
    }), use.names = FALSE)
    if (anyDuplicated(groups)) {
      sb_abort(
        "Group IDs cannot occur in both development and held-out cohorts.",
        "stormbreaker_internal_external_patient_overlap"
      )
    }
  }

  patient_ids <- lapply(tasks, function(task) {
    as.character(task$data[[task$id_col]])
  })
  all_ids <- unlist(patient_ids, use.names = FALSE)
  repeated <- unique(all_ids[duplicated(all_ids)])
  if (length(repeated) > 0L) {
    examples <- paste(utils::head(repeated, 3L), collapse = ", ")
    sb_abort(
      paste0(
        "Patient IDs must be globally unique across cohorts so a patient ",
        "cannot enter both development and held-out data. Repeated IDs: ",
        examples, ". Namespace cohort-local IDs before creating the tasks."
      ),
      "stormbreaker_internal_external_patient_overlap"
    )
  }

  list(
    tasks = tasks,
    task_type = task_type,
    cohort_ids = cohort_ids,
    patient_ids = patient_ids,
    predictors = reference_predictors
  )
}

sb_validate_iecv_pipelines <- function(pipelines, task_type) {
  valid_collection <- is.list(pipelines) && length(pipelines) > 0L &&
    !inherits(pipelines, "stormbreaker_pipeline")
  if (!valid_collection) {
    sb_abort(
      paste0(
        "`pipelines` must be a non-empty list of explicit specifications ",
        "created by `sb_pipeline()`."
      ),
      "stormbreaker_bad_internal_external_pipelines"
    )
  }
  pipelines <- lapply(pipelines, function(pipeline) {
    sb_validate_pipeline_contract(pipeline, task_type)
  })
  pipeline_ids <- vapply(pipelines, `[[`, character(1L), "id")
  if (anyDuplicated(pipeline_ids)) {
    sb_abort(
      "Every internal-external validation pipeline must have a distinct ID.",
      "stormbreaker_duplicate_pipeline_id"
    )
  }
  names(pipelines) <- pipeline_ids
  pipelines
}

sb_validate_iecv_strategy <- function(pipelines, task_type) {
  if (inherits(pipelines, "stormbreaker_stack")) {
    stack <- sb_validate_stack_contract(pipelines, task_type)
    return(list(
      type = "stack",
      specification = stack,
      ids = stack$id
    ))
  }
  validated <- sb_validate_iecv_pipelines(pipelines, task_type)
  list(
    type = "pipelines",
    specification = validated,
    ids = names(validated)
  )
}

sb_fit_iecv_strategy <- function(task, strategy) {
  if (identical(strategy$type, "stack")) {
    return(stats::setNames(
      list(sb_fit_stack_final(task, strategy$specification)),
      strategy$ids
    ))
  }
  sb_with_selector_reuse(lapply(strategy$specification, function(pipeline) {
    sb_fit_final(task, pipeline)
  }))
}

sb_predict_iecv_strategy <- function(fits, strategy, new_data, id_col) {
  if (identical(strategy$type, "stack")) {
    return(stats::setNames(
      list(sb_predict_stack_final(fits[[1L]], new_data, id_col)),
      strategy$ids
    ))
  }
  lapply(fits, function(fit) {
    sb_predict_final(fit, new_data, id_col)
  })
}

sb_iecv_task_rows <- function(task, reference, task_type) {
  role_names <- if (identical(task_type, "binary")) {
    c(task$id_col, task$outcome_col)
  } else {
    c(task$id_col, task$time_col, task$event_col)
  }
  reference_roles <- if (identical(task_type, "binary")) {
    c(reference$id_col, reference$outcome_col)
  } else {
    c(reference$id_col, reference$time_col, reference$event_col)
  }
  if (!is.null(reference$group_col)) {
    role_names <- c(role_names, task$group_col)
    reference_roles <- c(reference_roles, reference$group_col)
  }
  rows <- task$data[c(role_names, reference$predictors)]
  names(rows)[seq_along(reference_roles)] <- reference_roles
  rows
}

sb_iecv_development_task <- function(tasks, reference, task_type) {
  rows <- lapply(tasks, function(task) {
    sb_iecv_task_rows(task, reference, task_type)
  })
  data <- do.call(rbind, rows)
  rownames(data) <- NULL

  decision_thresholds <- if ("decision_thresholds" %in% names(reference)) {
    reference$decision_thresholds
  } else {
    NULL
  }
  if (identical(task_type, "binary")) {
    return(sb_binary_task(
      data,
      reference$id_col,
      reference$outcome_col,
      reference$predictors,
      decision_thresholds = decision_thresholds,
      group_col = reference$group_col
    ))
  }

  prediction_horizon <- if ("prediction_horizon" %in% names(reference)) {
    reference$prediction_horizon
  } else {
    NULL
  }
  sb_surv_task(
    data,
    reference$id_col,
    reference$time_col,
    reference$event_col,
    reference$predictors,
    reference$time_unit,
    prediction_horizon = prediction_horizon,
    decision_thresholds = decision_thresholds,
    group_col = reference$group_col
  )
}

sb_iecv_endpoint <- function(task, task_type) {
  if (identical(task_type, "binary")) {
    return(sb_binary_endpoint(task$data, task$id_col, task$outcome_col))
  }
  horizon <- if ("prediction_horizon" %in% names(task)) {
    task$prediction_horizon
  } else {
    NULL
  }
  sb_surv_endpoint(
    task$data,
    task$id_col,
    task$time_col,
    task$event_col,
    task$time_unit,
    prediction_horizon = horizon
  )
}

sb_iecv_cohort_table <- function(tasks, task_type) {
  cohort_ids <- names(tasks)
  table <- data.frame(
    cohort_id = cohort_ids,
    n_patients = vapply(tasks, `[[`, integer(1L), "n"),
    stringsAsFactors = FALSE
  )
  if (identical(task_type, "binary")) {
    table$n_positive <- vapply(tasks, `[[`, integer(1L), "n_positive")
    table$n_negative <- vapply(tasks, `[[`, integer(1L), "n_negative")
    return(table)
  }
  table$n_events <- vapply(tasks, `[[`, integer(1L), "n_events")
  table$n_censored <- vapply(tasks, `[[`, integer(1L), "n_censored")
  table$time_unit <- vapply(tasks, `[[`, character(1L), "time_unit")
  if ("prediction_horizon" %in% names(tasks[[1L]])) {
    table$prediction_horizon <- vapply(
      tasks, `[[`, numeric(1L), "prediction_horizon"
    )
  }
  table
}

#' Run leave-one-cohort-out internal-external cross-validation
#'
#' `sb_run_internal_external_cv()` performs internal-external cross-validation
#' (IECV) over a named collection of compatible cohorts. In each cycle, one
#' complete cohort is excluded, either every predeclared pipeline or one
#' predeclared stack is fitted using only patients from the remaining cohorts,
#' and endpoint-free predictions are made for the excluded cohort. Only after
#' prediction is its outcome-only endpoint constructed and passed to
#' [sb_compare_external()]. Every cohort is excluded exactly once.
#'
#' This function evaluates the supplied modeling strategy; it does not invent
#' one. Development rows from the remaining cohorts are pooled without adding
#' cohort-specific intercepts, strata, random effects, or batch correction.
#' Cohorts must therefore declare the same numeric predictor-name set, time
#' target, and decision thresholds. Predictor order is normalized to the first
#' cohort, but feature intersection, measurement harmonization, eligibility,
#' outcome-definition harmonization, and outcome-blind upstream preprocessing
#' remain the caller's responsibility.
#'
#' A different fitted model is evaluated in every held-out cohort, and with
#' three or more cohorts their development sets overlap. Accordingly, the
#' cycle estimates are not independent evaluations of one locked model and
#' must not be passed to STORMBREAKER's external-performance meta-analysis
#' functions. No metric is pooled, no pipeline is ranked or selected, and no
#' final all-cohort model is fitted. If cycle results are inspected and used to
#' revise a pipeline, they are method-development evidence rather than final
#' independent external validation; that claim requires untouched data. Two
#' cohorts are executable, although the method literature has tentatively
#' recommended at least four studies for a useful heterogeneity appraisal.
#'
#' Patient IDs must be globally unique across cohorts. This detects direct
#' patient overlap between a development set and its held-out cohort. If two
#' independent cohorts reuse local labels such as `"P001"`, namespace those
#' labels before creating the tasks rather than allowing the function to hide
#' the ambiguity by prefixing them automatically.
#'
#' @param cohorts A named list of at least two [sb_binary_task()] objects or at
#'   least two [sb_surv_task()] objects. Names are unique cohort IDs. All tasks
#'   must have one task type and a common predictor-name set. Survival targets
#'   and any decision-threshold declarations must also match exactly.
#' @param pipelines Either a non-empty list of distinct, task-matched pipeline
#'   specifications returned by [sb_pipeline()] or one task-matched stack from
#'   [sb_stack()]. Specifications are fixed before the cycles; the function
#'   does not tune across held-out cohorts. A stack is evaluated as one modeling
#'   strategy; its members are not emitted as standalone IECV candidates.
#'
#' @return A `stormbreaker_internal_external_cv` containing overall execution
#'   `status`, the validation design, cohort and strategy summaries, one named
#'   `fold` per held-out cohort with its development IDs, final fits,
#'   endpoint-free predictions, and complete external comparison, plus a long
#'   cohort-by-pipeline `results` table. Failed fits and predictions remain
#'   visible and do not trigger fallback; later predeclared pipelines and
#'   held-out-cohort cycles continue. Cohort, task, and strategy-contract
#'   errors are validated before this collection behavior applies.
#'
#' @references Royston P, Parmar MKB, Sylvester R (2004). Construction and
#'   validation of a prognostic model across several studies, with an
#'   application in superficial bladder cancer. *Statistics in Medicine*,
#'   23, 907-926. \doi{10.1002/sim.1691}.
#'
#' Debray TPA, Moons KGM, Ahmed I, Koffijberg H, Riley RD (2013). A framework
#'   for developing, implementing, and evaluating clinical prediction models
#'   in an individual participant data meta-analysis. *Statistics in Medicine*,
#'   32, 3158-3180. \doi{10.1002/sim.5732}.
#'
#' Steyerberg EW, Harrell FE Jr (2016). Prediction models need appropriate
#'   internal, internal-external, and external validation. *Journal of Clinical
#'   Epidemiology*, 69, 245-247. \doi{10.1016/j.jclinepi.2015.04.005}.
#' @export
#'
#' @examples
#' a <- sb_example_binary(n = 24, seed = 1)
#' b <- sb_example_binary(n = 24, seed = 2)
#' a$patient_id <- paste0("a_", a$patient_id)
#' b$patient_id <- paste0("b_", b$patient_id)
#' predictors <- c("gene_a", "gene_b")
#' cohorts <- list(
#'   cohort_a = sb_binary_task(
#'     a, "patient_id", "outcome", predictors
#'   ),
#'   cohort_b = sb_binary_task(
#'     b, "patient_id", "outcome", predictors
#'   )
#' )
#' sb_run_internal_external_cv(
#'   cohorts,
#'   list(sb_pipeline(learner = "logistic"))
#' )
sb_run_internal_external_cv <- function(cohorts, pipelines) {
  cohort_contract <- sb_validate_iecv_cohorts(cohorts)
  tasks <- cohort_contract$tasks
  task_type <- cohort_contract$task_type
  cohort_ids <- cohort_contract$cohort_ids
  reference <- tasks[[1L]]
  strategy <- sb_validate_iecv_strategy(pipelines, task_type)
  pipeline_ids <- strategy$ids

  folds <- lapply(seq_along(tasks), function(heldout_index) {
    heldout_id <- cohort_ids[[heldout_index]]
    development_ids <- cohort_ids[-heldout_index]
    development_task <- sb_iecv_development_task(
      tasks[-heldout_index], reference, task_type
    )
    fits <- sb_fit_iecv_strategy(development_task, strategy)

    heldout <- tasks[[heldout_index]]
    predictor_data <- heldout$data[c(
      heldout$id_col, cohort_contract$predictors
    )]
    predictions <- sb_predict_iecv_strategy(
      fits,
      strategy,
      predictor_data,
      heldout$id_col
    )

    endpoint <- sb_iecv_endpoint(heldout, task_type)
    comparison <- sb_compare_external(predictions, endpoint)
    list(
      status = comparison$status,
      heldout_cohort_id = heldout_id,
      development_cohort_ids = development_ids,
      development_patient_ids = as.character(
        development_task$data[[development_task$id_col]]
      ),
      heldout_patient_ids = cohort_contract$patient_ids[[heldout_index]],
      fits = fits,
      predictions = predictions,
      comparison = comparison
    )
  })
  names(folds) <- cohort_ids

  results <- lapply(folds, function(fold) {
    table <- fold$comparison$comparison
    cbind(
      data.frame(
        heldout_cohort_id = rep(fold$heldout_cohort_id, nrow(table)),
        n_development_cohorts = rep(
          as.integer(length(fold$development_cohort_ids)), nrow(table)
        ),
        n_development_patients = rep(
          as.integer(length(fold$development_patient_ids)), nrow(table)
        ),
        n_heldout_patients = rep(
          as.integer(length(fold$heldout_patient_ids)), nrow(table)
        ),
        stringsAsFactors = FALSE
      ),
      table
    )
  })
  results <- do.call(rbind, results)
  rownames(results) <- NULL

  complete <- all(vapply(
    folds, function(fold) identical(fold$status, "success"), logical(1L)
  ))
  result <- list(
    status = if (complete) "success" else "incomplete",
    task_type = task_type,
    validation_scheme = "leave_one_cohort_out_internal_external_cv",
    design = list(
      heldout_unit = "cohort",
      development_rule = "all_other_cohorts",
      development_strategy =
        "pooled_remaining_cohort_patients_without_cohort_effects",
      heldout_outcomes = "evaluation_only",
      fold_training_sets_overlap = length(tasks) > 2L,
      final_independent_external_validation = FALSE
    ),
    pipeline_ids = pipeline_ids,
    cohorts = sb_iecv_cohort_table(tasks, task_type),
    results = results,
    folds = folds
  )
  if (identical(strategy$type, "stack")) {
    result$strategy_type <- "stack"
  }
  structure(result, class = "stormbreaker_internal_external_cv")
}

#' @export
print.stormbreaker_internal_external_cv <- function(x, ...) {
  cat("<stormbreaker_internal_external_cv:", x$status, ">\n")
  cat("leave-one-cohort-out internal-external cross-validation\n")
  count_label <- if (identical(x$strategy_type, "stack")) {
    "Stacks:"
  } else {
    "Pipelines:"
  }
  cat(
    "Cohorts:", nrow(x$cohorts),
    "|", count_label, length(x$pipeline_ids),
    "| Task:", x$task_type, "\n"
  )
  metric_columns <- if (identical(x$task_type, "binary")) {
    c("roc_auc", "brier_score")
  } else {
    c("harrell_c", "ipcw_brier", "time_auc", "uno_c")
  }
  columns <- c(
    "heldout_cohort_id", "pipeline_id", "status",
    intersect(metric_columns, names(x$results))
  )
  print(x$results[columns], row.names = FALSE, digits = 3)
  cat(
    paste0(
      "Each held-out cohort uses a different fitted model; no metric is ",
      "pooled, ranked, or used to select a winner.\n"
    )
  )
  cat(
    paste0(
      "This is method-development IECV, not final independent external ",
      "validation after pipeline choices are revised.\n"
    )
  )
  invisible(x)
}
