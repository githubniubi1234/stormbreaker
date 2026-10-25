#' Define a nominal multiclass classification task
#' @param data Data frame, one row per observation.
#' @param id_col Unique observation identifier column.
#' @param outcome_col Column containing three or more observed classes.
#' @param predictors Names of finite numeric predictor columns.
#' @param classes Explicit class order. By default use observed factor levels,
#'   or sorted character labels. No ordinal relationship is assumed.
#' @param group_col Optional patient/group column. Resampling keeps all rows
#'   from a group together, including inner tuning folds.
#' @return A `stormbreaker_multiclass_task` preserving rows and original labels.
#' @export
sb_multiclass_task <- function(data, id_col, outcome_col, predictors,
                               classes = NULL, group_col = NULL) {
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data))) sb_abort(
    "Supply a nonempty data frame with unique column names.", "stormbreaker_bad_data")
  sb_check_name(id_col, "id_col")
  sb_check_name(outcome_col, "outcome_col")
  if (id_col == outcome_col || !all(c(id_col, outcome_col) %in% names(data))) sb_abort(
    "Supply distinct existing ID and outcome columns.", "stormbreaker_bad_multiclass_task")
  id <- as.character(data[[id_col]])
  if (anyNA(id) || any(!nzchar(trimws(id))) || anyDuplicated(id)) sb_abort(
    "Observation IDs must be unique and nonblank.", "stormbreaker_bad_id")
  outcome <- data[[outcome_col]]
  if (!is.atomic(outcome) || !is.null(dim(outcome)) || anyNA(outcome) ||
      any(!nzchar(trimws(as.character(outcome)))) ||
      (is.numeric(outcome) && any(!is.finite(outcome)))) sb_abort(
    "Class labels must be nonmissing and nonblank.", "stormbreaker_bad_multiclass_outcome")
  observed <- unique(as.character(outcome))
  if (is.null(classes)) classes <- if (is.factor(outcome))
    levels(outcome)[levels(outcome) %in% observed] else sort(observed)
  if (!is.character(classes) || length(classes) < 3L || anyNA(classes) ||
      anyDuplicated(classes) || !setequal(classes, observed)) sb_abort(
    "classes must name each observed class exactly once, with at least three classes.",
    "stormbreaker_bad_multiclass_outcome")
  if (!is.character(predictors) || !length(predictors) || anyNA(predictors) ||
      anyDuplicated(predictors) || !all(predictors %in% names(data)) ||
      any(predictors %in% c(id_col, outcome_col))) sb_abort(
    "Supply distinct predictor columns excluding ID and outcome.", "stormbreaker_bad_predictors")
  sb_validate_binary_predictor_names(predictors)
  sb_validate_predictor_frame(data[predictors], "predictors")
  group_col <- sb_validate_group_column(data, group_col, c(id_col, outcome_col), predictors)
  structure(list(data = data, id_col = id_col, outcome_col = outcome_col,
    predictors = predictors, classes = classes, group_col = group_col, n = nrow(data)),
    class = "stormbreaker_multiclass_task")
}

sb_validate_multiclass_task <- function(task) {
  if (!inherits(task, "stormbreaker_multiclass_task")) sb_abort(
    "Create the task with sb_multiclass_task().", "stormbreaker_bad_multiclass_task")
  sb_multiclass_task(task$data, task$id_col, task$outcome_col, task$predictors,
    task$classes, task$group_col)
}

sb_multiclass_folds <- function(task, folds, seed) {
  folds <- sb_check_positive_integer(folds, "folds", minimum = 2L)
  seed <- sb_check_positive_integer(seed, "seed", minimum = 1L)
  y <- factor(task$data[[task$outcome_col]], levels = task$classes)
  groups <- sb_task_group_ids(task)
  if (is.null(groups) || !anyDuplicated(groups)) {
    if (any(table(y) < folds)) sb_abort(
      "Each class needs at least folds observations.", "stormbreaker_infeasible_folds")
    assignment <- integer(task$n)
    withr::with_seed(seed, {
      for (level in levels(y)) {
        rows <- which(y == level)
        assignment[rows[sample.int(length(rows))]] <- rep(seq_len(folds), length.out = length(rows))
      }
    })
  } else {
    group_index <- match(groups, unique(groups))
    counts <- unclass(table(factor(group_index), y))
    if (any(colSums(counts > 0) < folds)) sb_abort(
      "Each class needs at least folds distinct groups.", "stormbreaker_infeasible_folds")
    totals <- colSums(counts)
    assigned <- integer(nrow(counts))
    balance <- matrix(0, folds, ncol(counts))
    withr::with_seed(seed, {
      order_groups <- sample.int(nrow(counts))
      order_groups <- order_groups[order(-rowSums(counts[order_groups, , drop = FALSE] > 0),
        -rowSums(counts[order_groups, , drop = FALSE]))]
      for (g in order_groups) {
        cost <- drop(balance %*% (counts[g, ] / totals^2))
        candidates <- which(cost == min(cost))
        chosen <- candidates[sample.int(length(candidates), 1L)]
        assigned[g] <- chosen
        balance[chosen, ] <- balance[chosen, ] + counts[g, ]
      }
    })
    assignment <- assigned[group_index]
  }
  counts <- table(factor(assignment, levels = seq_len(folds)), y)
  if (any(counts == 0)) sb_abort(
    "Generated folds do not all contain every class. Reduce folds or revise grouping.",
    "stormbreaker_infeasible_folds")
  as.integer(assignment)
}
