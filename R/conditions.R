sb_abort <- function(message, subclass, parent = NULL) {
  condition <- structure(
    list(message = message, call = NULL, parent = parent),
    class = c(subclass, "stormbreaker_error", "error", "condition")
  )
  stop(condition)
}

sb_check_name <- function(value, argument) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
    sb_abort(
      paste0("`", argument, "` must be one non-empty column name."),
      "stormbreaker_bad_column_argument"
    )
  }
}

sb_quote_name <- function(value) {
  if (grepl("`", value, fixed = TRUE)) {
    sb_abort(
      "Column names containing backticks are not supported.",
      "stormbreaker_unsupported_column_name"
    )
  }
  paste0("`", value, "`")
}

sb_check_positive_integer <- function(value, argument, minimum = 1L) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < minimum ||
      value > .Machine$integer.max || value != as.integer(value)) {
    sb_abort(
      paste0("`", argument, "` must be an integer >= ", minimum, "."),
      "stormbreaker_bad_integer"
    )
  }
  as.integer(value)
}

