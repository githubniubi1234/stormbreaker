sb_valid_decision_thresholds <- function(value) {
  is.numeric(value) && length(value) >= 1L &&
    !anyNA(value) && all(is.finite(value)) &&
    all(value > 0 & value < 1) &&
    !anyDuplicated(value) &&
    (length(value) == 1L || all(diff(value) > 0))
}

sb_validate_decision_thresholds <- function(value) {
  if (!sb_valid_decision_thresholds(value)) {
    sb_abort(
      paste0(
        "`decision_thresholds` must be unique, strictly increasing finite ",
        "probabilities between 0 and 1."
      ),
      "stormbreaker_bad_decision_thresholds"
    )
  }
  unname(as.numeric(value))
}

#' Define a binary classification task
#'
#' `sb_binary_task()` records an explicit patient-level `0/1` outcome. Class
#' `1` is the positive class. The function never drops patients, imputes
#' values, converts labels, or guesses the positive class. Optional decision
#' thresholds are supplied by the researcher; the package does not derive or
#' optimize them from outcomes or model performance.
#'
#' @param data A data frame with one row per observation. Without `group_col`,
#'   each row is treated as an independent patient.
#' @param id_col Name of the unique observation identifier column.
#' @param outcome_col Name of the numeric or logical outcome column, coded as
#'   `0` or `1`.
#' @param predictors Character vector naming finite numeric predictor columns.
#'   The reserved model term `"(Intercept)"` cannot be a predictor name.
#' @param decision_thresholds Optional unique, strictly increasing finite
#'   probability thresholds between zero and one. These should represent
#'   clinically meaningful action thresholds declared before model evaluation.
#'   No default grid is supplied.
#' @param group_col Optional patient/group identifier column for repeated
#'   observations. Row IDs remain unique. Generated evaluation folds, glmnet
#'   learner tuning, and e1071 SVM probability calibration keep each group
#'   together; supplied folds must do so too.
#'   This does not change metric weighting or bootstrap sampling, and does not
#'   configure resampling internal to other learner or selector engines.
#'
#' @references Vickers AJ, Elkin EB (2006). Decision curve analysis: a novel
#'   method for evaluating prediction models. *Medical Decision Making*,
#'   26(6), 565-574. \doi{10.1177/0272989X06295361}.
#'
#' @return A `stormbreaker_binary_task` object preserving the input rows and
#'   recording outcome and predictor roles.
#' @export
#'
#' @examples
#' d <- sb_example_binary(n = 30, seed = 1)
#' task <- sb_binary_task(
#'   d, "patient_id", "outcome",
#'   predictors = c("gene_a", "gene_b")
#' )
#' task
sb_binary_task <- function(data, id_col, outcome_col, predictors,
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
  sb_check_name(outcome_col, "outcome_col")
  if (identical(id_col, outcome_col)) {
    sb_abort(
      "ID and outcome must use different columns.",
      "stormbreaker_overlapping_endpoint_roles"
    )
  }
  required_columns <- c(id_col, outcome_col)
  missing_required <- setdiff(required_columns, names(data))
  if (length(missing_required) > 0L) {
    sb_abort(
      paste0(
        "Missing binary task columns: ",
        paste(missing_required, collapse = ", "), "."
      ),
      "stormbreaker_missing_endpoint_columns"
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

  outcome <- data[[outcome_col]]
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
  if (sum(outcome == 1) < 1L) {
    sb_abort(
      "At least one positive patient with outcome=1 is required.",
      "stormbreaker_no_positive_class"
    )
  }
  if (sum(outcome == 0) < 1L) {
    sb_abort(
      "At least one negative patient with outcome=0 is required.",
      "stormbreaker_no_negative_class"
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
  sb_validate_binary_predictor_names(predictors)
  missing_predictors <- setdiff(predictors, names(data))
  if (length(missing_predictors) > 0L) {
    sb_abort(
      paste0(
        "Missing predictor columns: ",
        paste(missing_predictors, collapse = ", "), "."
      ),
      "stormbreaker_missing_predictors"
    )
  }
  if (length(intersect(predictors, required_columns)) > 0L) {
    sb_abort(
      "Predictors cannot reuse the ID or outcome columns.",
      "stormbreaker_reserved_predictor"
    )
  }
  numeric_predictors <- vapply(data[predictors], is.numeric, logical(1))
  if (!all(numeric_predictors)) {
    sb_abort(
      paste0(
        "All minimum-core predictors must be numeric: ",
        paste(
          names(numeric_predictors)[!numeric_predictors],
          collapse = ", "
        ), "."
      ),
      "stormbreaker_bad_predictor_type"
    )
  }
  finite_predictors <- vapply(
    data[predictors],
    function(column) !anyNA(column) && all(is.finite(column)),
    logical(1)
  )
  if (!all(finite_predictors)) {
    sb_abort(
      paste0(
        "Predictors contain missing or non-finite values: ",
        paste(
          names(finite_predictors)[!finite_predictors],
          collapse = ", "
        ), "."
      ),
      "stormbreaker_bad_predictor_values"
    )
  }

  if (!is.null(decision_thresholds)) {
    decision_thresholds <- sb_validate_decision_thresholds(
      decision_thresholds
    )
  }
  group_col <- sb_validate_group_column(
    data, group_col, c(id_col, outcome_col), predictors
  )

  task <- list(
    data = data,
    id_col = id_col,
    outcome_col = outcome_col,
    predictors = predictors,
    positive_class = 1L,
    negative_class = 0L,
    n = nrow(data),
    n_positive = as.integer(sum(outcome == 1)),
    n_negative = as.integer(sum(outcome == 0)),
    independent_unit = "patient"
  )
  if (!is.null(decision_thresholds)) {
    task <- append(
      task,
      list(decision_thresholds = decision_thresholds),
      after = match("negative_class", names(task))
    )
  }
  if (!is.null(group_col)) task$group_col <- group_col
  class(task) <- "stormbreaker_binary_task"
  task
}

sb_validate_binary_task_contract <- function(task) {
  required_metadata_fields <- c(
    "id_col", "outcome_col", "predictors", "positive_class",
    "negative_class", "n", "n_positive", "n_negative",
    "independent_unit"
  )
  required_fields <- c("data", required_metadata_fields)
  if (!is.list(task) || length(setdiff(required_fields, names(task))) > 0L) {
    sb_abort(
      "The binary task is incomplete; recreate it with `sb_binary_task()`.",
      "stormbreaker_bad_task"
    )
  }
  has_decision_thresholds <- "decision_thresholds" %in% names(task)
  decision_thresholds <- if (has_decision_thresholds) {
    task$decision_thresholds
  } else {
    NULL
  }
  validated <- sb_binary_task(
    data = task$data,
    id_col = task$id_col,
    outcome_col = task$outcome_col,
    predictors = task$predictors,
    decision_thresholds = decision_thresholds,
    group_col = task$group_col
  )
  comparison_fields <- required_metadata_fields
  if (has_decision_thresholds) {
    comparison_fields <- c(comparison_fields, "decision_thresholds")
  }
  if ("group_col" %in% names(task)) {
    comparison_fields <- c(comparison_fields, "group_col")
  }
  metadata_ok <- identical(
    unclass(task)[comparison_fields],
    unclass(validated)[comparison_fields]
  )
  if (!metadata_ok) {
    sb_abort(
      "The binary task metadata changed; recreate it with `sb_binary_task()`.",
      "stormbreaker_mutated_task"
    )
  }
  validated
}

#' @export
print.stormbreaker_binary_task <- function(x, ...) {
  cat("<stormbreaker_binary_task>\n")
  cat(
    if (is.null(x$group_col)) "Patients:" else "Observations:", x$n,
    "| Positive (1):", x$n_positive,
    "| Negative (0):", x$n_negative, "\n"
  )
  cat("Outcome:", x$outcome_col, "| Predictors:", length(x$predictors), "\n")
  if (!is.null(x$group_col)) {
    cat("Groups:", length(unique(sb_task_group_ids(x))), "| Column:",
        x$group_col, "\n")
  }
  if ("decision_thresholds" %in% names(x)) {
    cat(
      "Decision thresholds:", length(x$decision_thresholds),
      "| Range:",
      paste(range(x$decision_thresholds), collapse = " to "), "\n"
    )
  }
  invisible(x)
}
