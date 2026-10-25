#' Define a right-censored survival task
#'
#' `sb_surv_task()` is STORMBREAKER's data contract. It records column roles
#' and rejects ambiguous or incomplete endpoint data. The function never
#' removes patients, imputes values, converts time units, or chooses a
#' prediction horizon or decision threshold from the observed data.
#'
#' @param data A data frame with one row per observation. Without `group_col`,
#'   each row is treated as an independent patient.
#' @param id_col Name of the unique observation identifier column.
#' @param time_col Name of the non-negative survival-time column.
#' @param event_col Name of the event column, coded as `0` or `1`, where `1`
#'   means the event was observed.
#' @param predictors Character vector naming numeric predictor columns.
#' @param time_unit Explicit unit label: `"day"`, `"month"`, or `"year"`.
#' @param prediction_horizon Optional single positive time point for survival
#'   probability prediction and IPCW Brier scoring. Its unit is `time_unit`.
#' @param decision_thresholds Optional unique, strictly increasing finite
#'   event-risk probabilities between zero and one. A prediction horizon is
#'   required. These should be clinically meaningful action thresholds
#'   declared before evaluation; no default grid is supplied.
#' @param group_col Optional patient/group identifier column for repeated
#'   observations. Row IDs remain unique. Generated evaluation folds and glmnet
#'   learner tuning keep each group together; supplied folds must do so too.
#'   This does not change metric weighting or bootstrap sampling, and does not
#'   configure resampling internal to other learner or selector engines.
#'
#' @return A `stormbreaker_surv_task` object. Its `data` element preserves the
#'   input rows and its remaining elements record endpoint and predictor roles,
#'   plus the prediction horizon and event-risk thresholds when declared.
#'
#' @references Vickers AJ, Cronin AM, Elkin EB, Gonen M (2008). Extensions to
#'   decision curve analysis. *BMC Medical Informatics and Decision Making*,
#'   8, 53. \doi{10.1186/1472-6947-8-53}.
#' @export
#'
#' @examples
#' d <- sb_example_survival(n = 30, seed = 1)
#' task <- sb_surv_task(
#'   d, "patient_id", "time", "event",
#'   predictors = c("gene_a", "gene_b"),
#'   time_unit = "month"
#' )
#' task
sb_surv_task <- function(data, id_col, time_col, event_col, predictors,
                         time_unit, prediction_horizon = NULL,
                         decision_thresholds = NULL, group_col = NULL) {
  if (!is.data.frame(data)) {
    sb_abort("`data` must be a data.frame.", "stormbreaker_bad_data")
  }
  if (nrow(data) < 2L) {
    sb_abort(
      "`data` must contain at least two patients.",
      "stormbreaker_too_few_patients"
    )
  }
  if (anyDuplicated(names(data))) {
    sb_abort(
      "The data frame contains duplicate column names.",
      "stormbreaker_duplicate_column_names"
    )
  }

  sb_check_name(id_col, "id_col")
  sb_check_name(time_col, "time_col")
  sb_check_name(event_col, "event_col")
  endpoint_cols <- c(id_col, time_col, event_col)
  if (anyDuplicated(endpoint_cols)) {
    sb_abort(
      "ID, time, and event must use three different columns.",
      "stormbreaker_overlapping_endpoint_roles"
    )
  }

  missing_endpoint <- setdiff(endpoint_cols, names(data))
  if (length(missing_endpoint) > 0L) {
    sb_abort(
      paste0(
        "Missing endpoint columns: ",
        paste(missing_endpoint, collapse = ", "),
        "."
      ),
      "stormbreaker_missing_endpoint_columns"
    )
  }

  allowed_units <- c("day", "month", "year")
  if (!is.character(time_unit) || length(time_unit) != 1L ||
      is.na(time_unit) || !time_unit %in% allowed_units) {
    sb_abort(
      "`time_unit` must be exactly 'day', 'month', or 'year'.",
      "stormbreaker_bad_time_unit"
    )
  }

  if (!is.null(prediction_horizon)) {
    valid_horizon <- is.numeric(prediction_horizon) &&
      length(prediction_horizon) == 1L &&
      !is.na(prediction_horizon) &&
      is.finite(prediction_horizon) &&
      prediction_horizon > 0
    if (!valid_horizon) {
      sb_abort(
        "`prediction_horizon` must be one finite positive number.",
        "stormbreaker_bad_prediction_horizon"
      )
    }
    prediction_horizon <- as.numeric(prediction_horizon)
  }
  if (!is.null(decision_thresholds) && is.null(prediction_horizon)) {
    sb_abort(
      paste0(
        "`decision_thresholds` require a declared ",
        "`prediction_horizon`."
      ),
      "stormbreaker_decision_thresholds_require_prediction_horizon"
    )
  }
  if (!is.null(decision_thresholds)) {
    decision_thresholds <- sb_validate_decision_thresholds(
      decision_thresholds
    )
  }

  id <- data[[id_col]]
  if (anyNA(id) || any(!nzchar(trimws(as.character(id))))) {
    sb_abort(
      "Patient identifiers contain missing or blank values.",
      "stormbreaker_bad_id"
    )
  }
  if (anyDuplicated(id)) {
    sb_abort(
      "Patient identifiers must be unique.",
      "stormbreaker_duplicate_id"
    )
  }

  time <- data[[time_col]]
  if (!is.numeric(time) || anyNA(time) || any(!is.finite(time)) ||
      any(time < 0)) {
    sb_abort(
      "Survival time must contain finite, non-negative numeric values.",
      "stormbreaker_bad_time"
    )
  }

  event <- data[[event_col]]
  if ((!is.numeric(event) && !is.logical(event)) || anyNA(event) ||
      (is.numeric(event) && any(!is.finite(event))) ||
      !all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only finite 0/1 values.",
      "stormbreaker_bad_event"
    )
  }
  if (sum(event == 1) < 1L) {
    sb_abort(
      "At least one observed event is required.",
      "stormbreaker_no_events"
    )
  }
  if (sum(event == 0) < 1L) {
    sb_abort(
      "At least one censored patient is required.",
      "stormbreaker_no_censoring"
    )
  }

  if (!is.character(predictors) || length(predictors) < 1L ||
      anyNA(predictors) || any(!nzchar(predictors))) {
    sb_abort(
      "`predictors` must contain at least one non-empty column name.",
      "stormbreaker_bad_predictors"
    )
  }
  predictors <- unique(predictors)
  missing_predictors <- setdiff(predictors, names(data))
  if (length(missing_predictors) > 0L) {
    sb_abort(
      paste0(
        "Missing predictor columns: ",
        paste(missing_predictors, collapse = ", "),
        "."
      ),
      "stormbreaker_missing_predictors"
    )
  }
  if (length(intersect(predictors, endpoint_cols)) > 0L) {
    sb_abort(
      "Predictors cannot reuse the ID, time, or event columns.",
      "stormbreaker_reserved_predictor"
    )
  }
  numeric_predictors <- vapply(data[predictors], is.numeric, logical(1))
  if (!all(numeric_predictors)) {
    sb_abort(
      paste0(
        "All minimum-core predictors must be numeric: ",
        paste(names(numeric_predictors)[!numeric_predictors],
              collapse = ", "),
        "."
      ),
      "stormbreaker_bad_predictor_type"
    )
  }
  finite_predictors <- vapply(
    data[predictors],
    function(x) !anyNA(x) && all(is.finite(x)),
    logical(1)
  )
  if (!all(finite_predictors)) {
    sb_abort(
      paste0(
        "Predictors contain missing or non-finite values: ",
        paste(names(finite_predictors)[!finite_predictors],
              collapse = ", "),
        "."
      ),
      "stormbreaker_bad_predictor_values"
    )
  }

  group_col <- sb_validate_group_column(
    data, group_col, c(id_col, time_col, event_col), predictors
  )
  task <- list(
    data = data,
    id_col = id_col,
    time_col = time_col,
    event_col = event_col,
    predictors = predictors,
    time_unit = time_unit,
    n = nrow(data),
    n_events = as.integer(sum(event == 1)),
    n_censored = as.integer(sum(event == 0)),
    independent_unit = "patient"
  )
  if (!is.null(prediction_horizon)) {
    task <- append(
      task,
      list(prediction_horizon = prediction_horizon),
      after = match("time_unit", names(task))
    )
  }
  if (!is.null(decision_thresholds)) {
    task <- append(
      task,
      list(decision_thresholds = decision_thresholds),
      after = match("prediction_horizon", names(task))
    )
  }
  if (!is.null(group_col)) task$group_col <- group_col
  class(task) <- "stormbreaker_surv_task"
  task
}

sb_validate_task_contract <- function(task) {
  metadata_fields <- c(
    "id_col", "time_col", "event_col", "predictors", "time_unit", "n",
    "n_events", "n_censored", "independent_unit"
  )
  if ("group_col" %in% names(task)) {
    metadata_fields <- c(metadata_fields, "group_col")
  }
  if ("prediction_horizon" %in% names(task)) {
    metadata_fields <- append(
      metadata_fields,
      "prediction_horizon",
      after = match("time_unit", metadata_fields)
    )
  }
  if ("decision_thresholds" %in% names(task)) {
    if (!"prediction_horizon" %in% names(task)) {
      sb_abort(
        paste0(
          "`decision_thresholds` require a declared ",
          "`prediction_horizon`."
        ),
        "stormbreaker_decision_thresholds_require_prediction_horizon"
      )
    }
    metadata_fields <- append(
      metadata_fields,
      "decision_thresholds",
      after = match("prediction_horizon", metadata_fields)
    )
  }
  required_fields <- c("data", metadata_fields)
  if (!is.list(task) || length(setdiff(required_fields, names(task))) > 0L) {
    sb_abort(
      "The survival task is incomplete; recreate it with `sb_surv_task()`.",
      "stormbreaker_bad_task"
    )
  }
  prediction_horizon <- if ("prediction_horizon" %in% names(task)) {
    task[["prediction_horizon"]]
  } else {
    NULL
  }
  decision_thresholds <- if ("decision_thresholds" %in% names(task)) {
    task[["decision_thresholds"]]
  } else {
    NULL
  }

  validated <- sb_surv_task(
    data = task$data,
    id_col = task$id_col,
    time_col = task$time_col,
    event_col = task$event_col,
    predictors = task$predictors,
    time_unit = task$time_unit,
    prediction_horizon = prediction_horizon,
    decision_thresholds = decision_thresholds,
    group_col = task$group_col
  )
  metadata_ok <- identical(
    unclass(task)[metadata_fields],
    unclass(validated)[metadata_fields]
  )
  if (!metadata_ok) {
    sb_abort(
      "The survival task metadata changed; recreate it with `sb_surv_task()`.",
      "stormbreaker_mutated_task"
    )
  }

  validated
}

#' @export
print.stormbreaker_surv_task <- function(x, ...) {
  cat("<stormbreaker_surv_task>\n")
  cat(if (is.null(x$group_col)) "Patients:" else "Observations:",
      x$n, "| Events:", x$n_events,
      "| Censored:", x$n_censored, "\n")
  if (!is.null(x$group_col)) {
    cat("Groups:", length(unique(sb_task_group_ids(x))), "| Column:",
        x$group_col, "\n")
  }
  cat("Time:", x$time_col, "[", x$time_unit, "]",
      "| Event:", x$event_col, "\n")
  prediction_horizon <- if ("prediction_horizon" %in% names(x)) {
    x[["prediction_horizon"]]
  } else {
    NULL
  }
  if (!is.null(prediction_horizon)) {
    cat(
      "Prediction horizon:", format(prediction_horizon),
      "[", x$time_unit, "]\n"
    )
  }
  if ("decision_thresholds" %in% names(x)) {
    cat(
      "Decision thresholds:", length(x$decision_thresholds),
      "| Event risk by horizon\n"
    )
  }
  cat("Predictors:", length(x$predictors), "\n")
  invisible(x)
}
