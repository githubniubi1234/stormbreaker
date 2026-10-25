#' Render grouped clinical views from saved results
#' @param x A result from [sb_clinician_views()].
#' @param file New HTML filename; existing outputs are refused.
#' @param title Report title.
#' @param complete Expand all topics initially instead of only the first.
#' @return Invisibly the absolute report filename. Only unavailable or partial
#'   views and errors appear in the status table; normal status rows are hidden.
#'   Model settings appear in the report body. Original per-context diagnostics remain
#'   in x. No model is refitted and no extra metrics are requested.
#' @export
sb_render_clinician_report <- function(x, file, title = "Clinical model evaluation", complete = FALSE) {
  if (!inherits(x, "stormbreaker_clinician_views") || !length(x$plots)) sb_abort(
    "Supply clinician views with at least one available plot; inspect status for unavailable results.",
    "stormbreaker_bad_report_input")
  if (!is.logical(complete) || length(complete) != 1L || is.na(complete)) sb_abort(
    "complete must be TRUE or FALSE.", "stormbreaker_bad_report_input")
  issues <- x$status[x$status$status != "available", , drop = FALSE]
  captions <- x$captions
  if (!is.null(x$models) && nrow(x$models)) {
    settings <- vapply(seq_len(nrow(x$models)), function(i) {
      m <- x$models[i, ]
      parts <- c(m$model, if (m$selector != "None") m$selector,
        m$parameters, m$selection_parameters)
      paste(parts[nzchar(parts)], collapse = "; ")
    }, character(1))
    captions[1L] <- paste(c(captions[1L], "Model settings:", settings), collapse = "\n")
  }
  sb_render_visual_report(x$plots, file, title = title, groups = x$groups,
    expanded = if (complete) unique(x$groups) else NULL, captions = captions,
    availability = if (nrow(issues)) issues else NULL, width = x$width, height = x$height)
}
