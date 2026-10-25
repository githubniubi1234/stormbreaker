sb_prediction_variation <- function(predictions, task_type) {
  value_name <- if (identical(task_type, "binary")) "probability" else "risk"
  diagnostic <- list(status = "incomplete", n = 0L, distinct_values = NA_integer_,
    constant_folds = integer())
  if (!is.data.frame(predictions) ||
      !all(c("status", value_name) %in% names(predictions))) return(diagnostic)
  diagnostic$n <- nrow(predictions)
  values <- predictions[[value_name]]
  if (!length(values) || !is.numeric(values) || anyNA(values) ||
      any(!is.finite(values)) || anyNA(predictions$status) ||
      any(predictions$status != "success")) return(diagnostic)
  diagnostic$distinct_values <- length(unique(values))
  folds <- if ("fold" %in% names(predictions)) predictions$fold else rep(1L, length(values))
  if (anyNA(folds)) return(diagnostic)
  fold_ids <- sort(unique(folds))
  sizes <- vapply(fold_ids, function(fold) sum(folds == fold), integer(1))
  constant <- vapply(fold_ids, function(fold) {
    x <- values[folds == fold]
    length(x) >= 2L && all(x == x[1L])
  }, logical(1))
  diagnostic$constant_folds <- fold_ids[constant]
  diagnostic$status <- if (all(sizes < 2L)) "insufficient_patients" else
    if (all(constant)) "constant" else if (any(constant)) "partly_constant" else "variable"
  diagnostic
}
