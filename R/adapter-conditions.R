sb_validate_predictor_frame <- function(x, argument) {
  if (!is.data.frame(x) || nrow(x) < 1L || ncol(x) < 1L) {
    sb_abort(
      paste0("`", argument, "` must be a non-empty data frame."),
      "stormbreaker_bad_adapter_predictors"
    )
  }
  if (anyDuplicated(names(x)) || any(!nzchar(names(x)))) {
    sb_abort(
      paste0("`", argument, "` must have unique, non-empty column names."),
      "stormbreaker_bad_adapter_predictor_names"
    )
  }
  numeric_columns <- vapply(x, is.numeric, logical(1))
  if (!all(numeric_columns)) {
    sb_abort(
      paste0("`", argument, "` must contain only numeric columns."),
      "stormbreaker_bad_adapter_predictor_type"
    )
  }
  finite_columns <- vapply(
    x,
    function(column) !anyNA(column) && all(is.finite(column)),
    logical(1)
  )
  if (!all(finite_columns)) {
    sb_abort(
      paste0("`", argument, "` contains missing or non-finite values."),
      "stormbreaker_bad_adapter_predictor_values"
    )
  }
  x
}

sb_validate_endpoint_vector <- function(value, argument, n) {
  if (!is.numeric(value) || length(value) != n || anyNA(value) ||
      any(!is.finite(value))) {
    sb_abort(
      paste0(
        "`", argument,
        "` must contain one finite numeric value per analysis patient."
      ),
      "stormbreaker_bad_adapter_endpoint"
    )
  }
  value
}
