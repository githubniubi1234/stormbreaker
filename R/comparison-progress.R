sb_check_progress <- function(progress) {
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    sb_abort("`progress` must be TRUE or FALSE.", "stormbreaker_bad_progress")
  }
}

sb_progress_message <- function(progress, ...) {
  if (isTRUE(progress)) {
    message("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  }
  invisible(NULL)
}
