sb_select_endpoint_columns <- function(data, roles) {
  if (!is.data.frame(data)) {
    sb_abort("`data` must be a data.frame.", "stormbreaker_bad_data")
  }
  if (nrow(data) < 1L) {
    sb_abort(
      "An evaluation endpoint must contain at least one patient.",
      "stormbreaker_empty_endpoint"
    )
  }

  role_positions <- lapply(roles, function(role) which(names(data) == role))
  missing <- names(roles)[lengths(role_positions) == 0L]
  if (length(missing) > 0L) {
    sb_abort(
      paste0(
        "Missing endpoint columns: ",
        paste(unname(roles[missing]), collapse = ", "),
        "."
      ),
      "stormbreaker_missing_endpoint_columns"
    )
  }
  duplicated <- names(roles)[lengths(role_positions) > 1L]
  if (length(duplicated) > 0L) {
    sb_abort(
      paste0(
        "Endpoint columns are duplicated: ",
        paste(unname(roles[duplicated]), collapse = ", "),
        "."
      ),
      "stormbreaker_duplicate_endpoint_column"
    )
  }

  positions <- vapply(role_positions, `[[`, integer(1), 1L)
  selected <- data[positions]
  names(selected) <- unname(roles)
  selected
}

sb_validate_endpoint_ids <- function(id) {
  patient_ids <- as.character(id)
  if (anyNA(patient_ids) || any(!nzchar(trimws(patient_ids)))) {
    sb_abort(
      "Patient identifiers contain missing or blank values.",
      "stormbreaker_bad_id"
    )
  }
  if (anyDuplicated(patient_ids)) {
    sb_abort(
      "Patient identifiers must be unique.",
      "stormbreaker_duplicate_id"
    )
  }
  invisible(TRUE)
}

sb_validate_endpoint_time_unit <- function(time_unit) {
  valid <- is.character(time_unit) && length(time_unit) == 1L &&
    !is.na(time_unit) && time_unit %in% c("day", "month", "year")
  if (!valid) {
    sb_abort(
      "`time_unit` must be exactly 'day', 'month', or 'year'.",
      "stormbreaker_bad_time_unit"
    )
  }
  time_unit
}

sb_validate_endpoint_horizon <- function(prediction_horizon) {
  if (is.null(prediction_horizon)) {
    return(NULL)
  }
  valid <- is.numeric(prediction_horizon) &&
    length(prediction_horizon) == 1L &&
    !is.na(prediction_horizon) && is.finite(prediction_horizon) &&
    prediction_horizon > 0
  if (!valid) {
    sb_abort(
      "`prediction_horizon` must be one finite positive number.",
      "stormbreaker_bad_prediction_horizon"
    )
  }
  as.numeric(prediction_horizon)
}

#' Define a binary endpoint for external evaluation
#'
#' `sb_binary_endpoint()` records patient IDs and an observed `0/1` outcome
#' without accepting or retaining predictors. It is an outcome-only input for
#' [sb_evaluate_external()], not a modeling task. One patient or one observed
#' class is allowed; metrics that require both classes remain unavailable.
#'
#' @param data A data frame containing at least the ID and outcome columns.
#'   Unrequested columns are ignored and are not stored.
#' @param id_col Name of the unique patient identifier column.
#' @param outcome_col Name of the numeric or logical `0/1` outcome column.
#'
#' @return A `stormbreaker_binary_endpoint` containing only patient IDs,
#'   outcomes, class counts, and endpoint-role metadata.
#' @export
#'
#' @examples
#' endpoint <- sb_binary_endpoint(
#'   data.frame(patient_id = c("e1", "e2"), outcome = c(0, 1)),
#'   "patient_id",
#'   "outcome"
#' )
#' endpoint
sb_binary_endpoint <- function(data, id_col, outcome_col) {
  sb_check_name(id_col, "id_col")
  sb_check_name(outcome_col, "outcome_col")
  if (identical(id_col, outcome_col)) {
    sb_abort(
      "ID and outcome must use different columns.",
      "stormbreaker_overlapping_endpoint_roles"
    )
  }

  selected <- sb_select_endpoint_columns(
    data,
    c(id = id_col, outcome = outcome_col)
  )
  sb_validate_endpoint_ids(selected[[id_col]])
  outcome <- selected[[outcome_col]]
  valid_outcome <- (is.numeric(outcome) || is.logical(outcome)) &&
    !anyNA(outcome) &&
    (!is.numeric(outcome) || all(is.finite(outcome))) &&
    all(outcome %in% c(0, 1))
  if (!valid_outcome) {
    sb_abort(
      "`outcome` must contain only finite 0/1 values.",
      "stormbreaker_bad_binary_outcome"
    )
  }

  endpoint <- list(
    data = selected,
    id_col = id_col,
    outcome_col = outcome_col,
    positive_class = 1L,
    negative_class = 0L,
    n = as.integer(nrow(selected)),
    n_positive = as.integer(sum(outcome == 1)),
    n_negative = as.integer(sum(outcome == 0)),
    independent_unit = "patient"
  )
  class(endpoint) <- "stormbreaker_binary_endpoint"
  endpoint
}

sb_validate_binary_endpoint_contract <- function(endpoint) {
  metadata_fields <- c(
    "id_col", "outcome_col", "positive_class", "negative_class",
    "n", "n_positive", "n_negative", "independent_unit"
  )
  required_fields <- c("data", metadata_fields)
  valid_object <- inherits(endpoint, "stormbreaker_binary_endpoint") &&
    is.list(endpoint) &&
    length(setdiff(required_fields, names(endpoint))) == 0L
  if (!valid_object) {
    sb_abort(
      paste0(
        "The binary endpoint is incomplete; recreate it with ",
        "`sb_binary_endpoint()`."
      ),
      "stormbreaker_bad_endpoint"
    )
  }

  validated <- sb_binary_endpoint(
    endpoint$data,
    endpoint$id_col,
    endpoint$outcome_col
  )
  metadata_ok <- identical(
    unclass(endpoint)[metadata_fields],
    unclass(validated)[metadata_fields]
  )
  if (!metadata_ok) {
    sb_abort(
      paste0(
        "The binary endpoint metadata changed; recreate it with ",
        "`sb_binary_endpoint()`."
      ),
      "stormbreaker_mutated_endpoint"
    )
  }
  validated
}

#' @export
print.stormbreaker_binary_endpoint <- function(x, ...) {
  cat("<stormbreaker_binary_endpoint>\n")
  cat("outcome-only evaluation endpoint\n")
  cat(
    "Patients:", x$n,
    "| Positive (1):", x$n_positive,
    "| Negative (0):", x$n_negative, "\n"
  )
  cat("Outcome:", x$outcome_col, "\n")
  invisible(x)
}

#' Define a survival endpoint for external evaluation
#'
#' `sb_surv_endpoint()` records patient IDs, observed follow-up, and event
#' status without accepting or retaining predictors. It is an outcome-only
#' input for [sb_evaluate_external()], not a modeling task. One patient, all
#' events, or all censoring is allowed; metrics report when the endpoint does
#' not contain enough evaluable information.
#'
#' @param data A data frame containing at least the ID, time, and event columns.
#'   Unrequested columns are ignored and are not stored.
#' @param id_col Name of the unique patient identifier column.
#' @param time_col Name of the finite non-negative follow-up column.
#' @param event_col Name of the numeric or logical `0/1` event column, where
#'   `1` means the event was observed.
#' @param time_unit Explicit unit label: `"day"`, `"month"`, or `"year"`.
#' @param prediction_horizon Optional finite positive evaluation horizon in
#'   `time_unit`. A final prediction with a saved horizon controls scoring.
#'
#' @return A `stormbreaker_surv_endpoint` containing only patient IDs, time,
#'   event, counts, and time-semantics metadata.
#' @export
#'
#' @examples
#' endpoint <- sb_surv_endpoint(
#'   data.frame(
#'     patient_id = c("e1", "e2"),
#'     time = c(3, 5),
#'     event = c(1, 0)
#'   ),
#'   "patient_id",
#'   "time",
#'   "event",
#'   "month"
#' )
#' endpoint
sb_surv_endpoint <- function(data, id_col, time_col, event_col,
                             time_unit, prediction_horizon = NULL) {
  sb_check_name(id_col, "id_col")
  sb_check_name(time_col, "time_col")
  sb_check_name(event_col, "event_col")
  roles <- c(id_col, time_col, event_col)
  if (anyDuplicated(roles)) {
    sb_abort(
      "ID, time, and event must use three different columns.",
      "stormbreaker_overlapping_endpoint_roles"
    )
  }
  time_unit <- sb_validate_endpoint_time_unit(time_unit)
  prediction_horizon <- sb_validate_endpoint_horizon(prediction_horizon)

  selected <- sb_select_endpoint_columns(
    data,
    c(id = id_col, time = time_col, event = event_col)
  )
  sb_validate_endpoint_ids(selected[[id_col]])
  time <- selected[[time_col]]
  if (!is.numeric(time) || anyNA(time) || any(!is.finite(time)) ||
      any(time < 0)) {
    sb_abort(
      "Survival time must contain finite, non-negative numeric values.",
      "stormbreaker_bad_time"
    )
  }
  event <- selected[[event_col]]
  valid_event <- (is.numeric(event) || is.logical(event)) &&
    !anyNA(event) &&
    (!is.numeric(event) || all(is.finite(event))) &&
    all(event %in% c(0, 1))
  if (!valid_event) {
    sb_abort(
      "`event` must contain only finite 0/1 values.",
      "stormbreaker_bad_event"
    )
  }

  endpoint <- list(
    data = selected,
    id_col = id_col,
    time_col = time_col,
    event_col = event_col,
    time_unit = time_unit,
    n = as.integer(nrow(selected)),
    n_events = as.integer(sum(event == 1)),
    n_censored = as.integer(sum(event == 0)),
    independent_unit = "patient"
  )
  if (!is.null(prediction_horizon)) {
    endpoint <- append(
      endpoint,
      list(prediction_horizon = prediction_horizon),
      after = match("time_unit", names(endpoint))
    )
  }
  class(endpoint) <- "stormbreaker_surv_endpoint"
  endpoint
}

sb_validate_surv_endpoint_contract <- function(endpoint) {
  metadata_fields <- c(
    "id_col", "time_col", "event_col", "time_unit", "n",
    "n_events", "n_censored", "independent_unit"
  )
  if ("prediction_horizon" %in% names(endpoint)) {
    metadata_fields <- append(
      metadata_fields,
      "prediction_horizon",
      after = match("time_unit", metadata_fields)
    )
  }
  required_fields <- c("data", metadata_fields)
  valid_object <- inherits(endpoint, "stormbreaker_surv_endpoint") &&
    is.list(endpoint) &&
    length(setdiff(required_fields, names(endpoint))) == 0L
  if (!valid_object) {
    sb_abort(
      paste0(
        "The survival endpoint is incomplete; recreate it with ",
        "`sb_surv_endpoint()`."
      ),
      "stormbreaker_bad_endpoint"
    )
  }

  prediction_horizon <- if ("prediction_horizon" %in% names(endpoint)) {
    endpoint[["prediction_horizon"]]
  } else {
    NULL
  }
  validated <- sb_surv_endpoint(
    endpoint$data,
    endpoint$id_col,
    endpoint$time_col,
    endpoint$event_col,
    endpoint$time_unit,
    prediction_horizon
  )
  metadata_ok <- identical(
    unclass(endpoint)[metadata_fields],
    unclass(validated)[metadata_fields]
  )
  if (!metadata_ok) {
    sb_abort(
      paste0(
        "The survival endpoint metadata changed; recreate it with ",
        "`sb_surv_endpoint()`."
      ),
      "stormbreaker_mutated_endpoint"
    )
  }
  validated
}

#' @export
print.stormbreaker_surv_endpoint <- function(x, ...) {
  cat("<stormbreaker_surv_endpoint>\n")
  cat("outcome-only evaluation endpoint\n")
  cat(
    "Patients:", x$n,
    "| Events:", x$n_events,
    "| Censored:", x$n_censored, "\n"
  )
  cat(
    "Time:", x$time_col, "[", x$time_unit, "]",
    "| Event:", x$event_col, "\n"
  )
  if ("prediction_horizon" %in% names(x)) {
    cat(
      "Prediction horizon:", format(x$prediction_horizon),
      "[", x$time_unit, "]\n"
    )
  }
  invisible(x)
}
