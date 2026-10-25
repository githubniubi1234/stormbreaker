#' Open the Chinese clinician getting-started page
#' @param view Open the local page in a browser; defaults to interactive().
#' @return Invisibly the path to a local HTML page with modeling and validation workflows and
#'   copyable R examples. It does not read patient data, fit models, install
#'   packages or upload anything. Examples are read from the installed package.
#' @export
sb_start <- function(view = interactive()) {
  if (!is.logical(view) || length(view) != 1L || is.na(view)) sb_abort(
    "view must be TRUE or FALSE.", "stormbreaker_bad_report_input")
  html <- paste(readLines(system.file("start.html", package = "stormbreaker",
    mustWork = TRUE), encoding = "UTF-8", warn = FALSE), collapse = "\n")
  for (action in c("compare", "validate")) {
    recipe <- system.file("examples", paste0("doctor-", action, ".R"),
      package = "stormbreaker", mustWork = TRUE)
    code <- paste(readLines(recipe, encoding = "UTF-8", warn = FALSE), collapse = "\n")
    html <- gsub(paste0("@@", toupper(action), "@@"), sb_report_escape(code), html, fixed = TRUE)
  }
  file <- tempfile("stormbreaker-start-", fileext = ".html")
  writeLines(enc2utf8(html), file, useBytes = TRUE)
  if (view) utils::browseURL(file)
  invisible(normalizePath(file, winslash = "/"))
}
