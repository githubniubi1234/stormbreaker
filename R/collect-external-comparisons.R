sb_external_collection_member <- function(comparison, cohort_id) {
  if (!inherits(comparison, "stormbreaker_external_comparison") ||
      !is.list(comparison)) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` is not a result from ",
        "`sb_compare_external()`."
      ),
      "stormbreaker_bad_external_collection"
    )
  }

  task_type <- comparison[["task_type"]]
  status <- comparison[["status"]]
  n_patients <- comparison[["n_patients"]]
  table <- comparison[["comparison"]]
  valid_core <- is.character(task_type) && length(task_type) == 1L &&
    !is.na(task_type) && task_type %in% c("binary", "survival") &&
    is.character(status) && length(status) == 1L && !is.na(status) &&
    status %in% c("success", "incomplete") &&
    is.numeric(n_patients) && length(n_patients) == 1L &&
    !is.na(n_patients) && is.finite(n_patients) && n_patients >= 1 &&
    n_patients == as.integer(n_patients) &&
    is.data.frame(table) && nrow(table) >= 1L
  if (!valid_core) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` has incomplete external-comparison metadata."
      ),
      "stormbreaker_bad_external_collection"
    )
  }

  pipeline_ids <- table[["pipeline_id"]]
  learners <- table[["learner"]]
  valid_pipelines <- is.character(pipeline_ids) &&
    length(pipeline_ids) == nrow(table) && !anyNA(pipeline_ids) &&
    all(nzchar(trimws(pipeline_ids))) && anyDuplicated(pipeline_ids) == 0L &&
    is.character(learners) && length(learners) == nrow(table) &&
    !anyNA(learners) && all(nzchar(trimws(learners)))
  if (!valid_pipelines) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` does not retain unique pipeline IDs and learner labels."
      ),
      "stormbreaker_bad_external_collection"
    )
  }

  time_unit <- NULL
  prediction_horizon <- NULL
  if (identical(task_type, "survival")) {
    target <- comparison[["prediction_target"]]
    if (!is.list(target)) {
      sb_abort(
        paste0(
          "Cohort `", cohort_id,
          "` does not retain a usable survival evaluation target."
        ),
        "stormbreaker_bad_external_collection"
      )
    }
    target_status <- target[["target_status"]]
    valid_target_status <- is.character(target_status) &&
      length(target_status) == 1L && !is.na(target_status) &&
      target_status %in% c("complete", "partial", "unknown")
    if (!valid_target_status) {
      sb_abort(
        paste0(
          "Cohort `", cohort_id,
          "` does not retain a usable survival evaluation target."
        ),
        "stormbreaker_bad_external_collection"
      )
    }
    if (!identical(target_status, "complete")) {
      sb_abort(
        paste0(
          "Cohort `", cohort_id,
          "` has an unknown or partially recoverable survival target."
        ),
        "stormbreaker_external_target_mismatch"
      )
    }
    has_horizon <- target[["has_horizon"]]
    time_unit <- target[["time_unit"]]
    prediction_horizon <- target[["prediction_horizon"]]
    valid_has_horizon <- is.logical(has_horizon) &&
      length(has_horizon) == 1L && !is.na(has_horizon)
    valid_time_unit <- is.character(time_unit) && length(time_unit) == 1L &&
      !is.na(time_unit) && nzchar(trimws(time_unit))
    valid_horizon <- valid_has_horizon && (
      (!has_horizon && is.null(prediction_horizon)) ||
      (has_horizon && is.numeric(prediction_horizon) &&
        length(prediction_horizon) == 1L &&
        !is.na(prediction_horizon) && is.finite(prediction_horizon) &&
        prediction_horizon > 0)
    )
    if (!valid_time_unit || !valid_horizon) {
      sb_abort(
        paste0(
          "Cohort `", cohort_id,
          "` does not retain a usable survival evaluation target."
        ),
        "stormbreaker_bad_external_collection"
      )
    }
  }

  list(
    status = status,
    task_type = task_type,
    n_patients = as.integer(n_patients),
    table = table,
    pipeline_ids = pipeline_ids,
    learners = learners,
    time_unit = time_unit,
    prediction_horizon = prediction_horizon
  )
}

#' Collect compatible external comparisons without pooling cohorts
#'
#' `sb_collect_external_comparisons()` keeps several completed
#' [sb_compare_external()] results in one ordered object. The outer list names
#' are cohort IDs. Each source comparison and every failed model row are
#' retained; the long `results` table only adds `cohort_id` and normalizes
#' pipeline order to the first cohort.
#'
#' Cohorts must contain the same task, pipeline IDs, and learner labels. For a
#' survival task they must also use the same time unit and either the same
#' effective prediction horizon or no horizon in every cohort. The effective
#' target is the endpoint-declared horizon, or the common locked prediction
#' horizon when the endpoint omitted it. These checks prevent unlike evaluation
#' targets from being placed in one table. A comparison whose target remained
#' unknown or only partially recoverable cannot be collected. These checks
#' compare only the labels stored in the supplied objects: they cannot prove
#' cohort
#' independence, common eligibility or outcome definitions, common upstream
#' preprocessing, or that the same pipeline ID represents the same locked
#' fitted state or published formula. Those study-design conditions remain the
#' caller's responsibility.
#'
#' This function does not pool or average metrics, perform association or
#' performance meta-analysis, impute a missing pipeline, rank models, select a
#' winner, or use one cohort to change another cohort's predictions. Cohort-
#' specific uncertainty also remains in the original comparison objects.
#'
#' @param comparisons A non-empty named list of
#'   `stormbreaker_external_comparison` objects returned by
#'   [sb_compare_external()]. Names must be unique, non-blank cohort IDs. A
#'   single comparison must still be wrapped, for example
#'   `list(cohort_a = comparison)`.
#'
#' @return A `stormbreaker_external_collection` containing overall `status`,
#'   common `task_type`, ordered `pipeline_ids`, a one-row-per-cohort `cohorts`
#'   table, a cohort-specific long `results` table, and the unchanged named
#'   source `comparisons`.
#' @export
#'
#' @examples
#' data <- data.frame(
#'   patient_id = paste0("p", 1:12),
#'   outcome = rep(c(0L, 1L), 6),
#'   gene_a = seq(-1, 1, length.out = 12)
#' )
#' signature <- sb_linear_signature(
#'   "published_example", "binary", c(gene_a = 0.5), intercept = 0
#' )
#' prediction <- sb_predict_signature(
#'   signature, data[c("patient_id", "gene_a")], "patient_id"
#' )
#' comparison <- sb_compare_external(
#'   list(prediction),
#'   sb_binary_endpoint(data, "patient_id", "outcome")
#' )
#' sb_collect_external_comparisons(list(cohort_a = comparison))
sb_collect_external_comparisons <- function(comparisons) {
  valid_collection <- is.list(comparisons) && length(comparisons) > 0L &&
    !inherits(comparisons, "stormbreaker_external_comparison") &&
    all(vapply(
      comparisons,
      function(comparison) {
        inherits(comparison, "stormbreaker_external_comparison")
      },
      logical(1L)
    ))
  if (!valid_collection) {
    sb_abort(
      paste0(
        "`comparisons` must be a non-empty named list returned by ",
        "`sb_compare_external()`; wrap one object as ",
        "`list(cohort_id = comparison)`."
      ),
      "stormbreaker_bad_external_collection"
    )
  }

  cohort_ids <- names(comparisons)
  valid_ids <- !is.null(cohort_ids) &&
    length(cohort_ids) == length(comparisons) &&
    !anyNA(cohort_ids) && all(nzchar(trimws(cohort_ids))) &&
    anyDuplicated(cohort_ids) == 0L
  if (!valid_ids) {
    sb_abort(
      "Outer list names must be unique, non-blank cohort IDs.",
      "stormbreaker_bad_external_collection"
    )
  }

  members <- lapply(seq_along(comparisons), function(index) {
    sb_external_collection_member(comparisons[[index]], cohort_ids[[index]])
  })
  task_types <- vapply(members, `[[`, character(1L), "task_type")
  if (length(unique(task_types)) != 1L) {
    sb_abort(
      "All collected external comparisons must use the same task type.",
      "stormbreaker_external_task_mismatch"
    )
  }

  reference_ids <- members[[1L]]$pipeline_ids
  reference_learners <- members[[1L]]$learners
  ordered_tables <- lapply(seq_along(members), function(index) {
    member <- members[[index]]
    same_ids <- length(member$pipeline_ids) == length(reference_ids) &&
      setequal(member$pipeline_ids, reference_ids)
    if (!same_ids) {
      sb_abort(
        paste0(
          "Cohort `", cohort_ids[[index]],
          "` does not contain the same pipeline ID set as the first cohort."
        ),
        "stormbreaker_external_pipeline_mismatch"
      )
    }
    positions <- match(reference_ids, member$pipeline_ids)
    if (!identical(member$learners[positions], reference_learners)) {
      sb_abort(
        paste0(
          "Cohort `", cohort_ids[[index]],
          "` does not retain the same learner labels as the first cohort."
        ),
        "stormbreaker_external_pipeline_mismatch"
      )
    }
    member$table[positions, , drop = FALSE]
  })

  task_type <- task_types[[1L]]
  if (identical(task_type, "survival")) {
    reference_unit <- members[[1L]]$time_unit
    reference_horizon <- members[[1L]]$prediction_horizon
    compatible_target <- vapply(members, function(member) {
      same_unit <- identical(member$time_unit, reference_unit)
      same_horizon_presence <- identical(
        is.null(member$prediction_horizon),
        is.null(reference_horizon)
      )
      same_horizon <- same_horizon_presence &&
        (is.null(reference_horizon) || identical(
          as.numeric(member$prediction_horizon),
          as.numeric(reference_horizon)
        ))
      same_unit && same_horizon
    }, logical(1L))
    if (!all(compatible_target)) {
      sb_abort(
        paste0(
          "Survival comparisons must use the same time unit and prediction ",
          "horizon."
        ),
        "stormbreaker_external_target_mismatch"
      )
    }
  }

  results <- lapply(seq_along(ordered_tables), function(index) {
    cbind(
      data.frame(
        cohort_id = rep(cohort_ids[[index]], nrow(ordered_tables[[index]])),
        stringsAsFactors = FALSE
      ),
      ordered_tables[[index]]
    )
  })
  results <- do.call(rbind, results)
  rownames(results) <- NULL

  cohorts <- data.frame(
    cohort_id = cohort_ids,
    status = vapply(members, `[[`, character(1L), "status"),
    n_patients = vapply(members, `[[`, integer(1L), "n_patients"),
    stringsAsFactors = FALSE
  )
  if (identical(task_type, "survival")) {
    cohorts$time_unit <- vapply(members, `[[`, character(1L), "time_unit")
    if (!is.null(members[[1L]]$prediction_horizon)) {
      cohorts$prediction_horizon <- vapply(
        members,
        `[[`,
        numeric(1L),
        "prediction_horizon"
      )
    }
  }

  complete <- all(cohorts$status == "success")
  structure(
    list(
      status = if (complete) "success" else "incomplete",
      task_type = task_type,
      pipeline_ids = reference_ids,
      cohorts = cohorts,
      results = results,
      comparisons = comparisons
    ),
    class = "stormbreaker_external_collection"
  )
}

#' @export
print.stormbreaker_external_collection <- function(x, ...) {
  cat("<stormbreaker_external_collection:", x$status, ">\n")
  cat(
    "Cohorts:", nrow(x$cohorts),
    "| Pipelines:", length(x$pipeline_ids),
    "| Task:", x$task_type, "\n"
  )
  print(x$cohorts, row.names = FALSE)
  cat(
    paste0(
      "No metric is pooled; no pipeline is ranked or selected from these ",
      "cohorts.\n"
    )
  )
  invisible(x)
}
