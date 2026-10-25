sb_fit_none_selector <- function(x) {
  x <- sb_validate_predictor_frame(x, "x")
  structure(
    list(
      id = "none",
      selected_features = names(x)
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_apply_none_selector <- function(fitted, x) {
  if (!inherits(fitted, "stormbreaker_selector_fit") ||
      !identical(fitted$id, "none") ||
      !is.character(fitted$selected_features) ||
      length(fitted$selected_features) < 1L) {
    sb_abort(
      "The fitted selector is not a valid `none` selector state.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$selected_features)) {
    sb_abort(
      paste0(
        "Predictor columns must exactly match the fitted selector schema: ",
        paste(fitted$selected_features, collapse = ", "),
        "."
      ),
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x
}

sb_none_selector_adapter <- function() {
  list(
    id = "none",
    params = list(),
    needs_outcome = FALSE,
    fit = sb_fit_none_selector,
    apply = sb_apply_none_selector
  )
}
