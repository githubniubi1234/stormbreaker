#' Audit an ordered inventory of user-supplied signature records
#'
#' `sb_signature_inventory()` applies [sb_lint_signature()] sequentially to an
#' ordered, named list of literature records. Every input produces one output
#' row in the same order, including incomplete or malformed records. The list
#' names become `record_id`, a caller-supplied key that remains available even
#' when a rejected record has no usable `signature_id`.
#'
#' This is an in-memory audit inventory, not a built-in or verified literature
#' registry. It does not search publications, verify citations, infer missing
#' formula details, inspect patient data, execute predictions, merge duplicate
#' signatures, rank records, or choose a model. Repeated `signature_id` values
#' are allowed because separate source records may describe the same formula;
#' unique non-empty `record_id` names keep those source records distinguishable.
#' Unexpected internal linter errors remain errors rather than being converted
#' into rejected literature records.
#'
#' @param records A non-empty named list. Each element is one argument list or
#'   `stormbreaker_linear_signature` accepted by [sb_lint_signature()]. Names
#'   must be unique and non-empty and are returned unchanged as `record_id`.
#'   A single specification or signature object must still be wrapped, for
#'   example `list(paper_a = specification)`.
#'
#' @return A data frame with one row per supplied record, in input order. It
#'   adds `record_id` before the columns returned by [sb_lint_signature()].
#' @export
#'
#' @examples
#' records <- list(
#'   paper_a = list(
#'     signature_id = "complete_formula",
#'     task_type = "binary",
#'     coefficients = c(gene_a = 0.7),
#'     intercept = -0.2,
#'     source = "Citation A"
#'   ),
#'   paper_b = list(
#'     signature_id = "intercept_not_reported",
#'     task_type = "binary",
#'     coefficients = c(gene_a = 0.5),
#'     source = "Citation B"
#'   )
#' )
#' sb_signature_inventory(records)
sb_signature_inventory <- function(records) {
  record_names <- if (is.list(records)) names(records) else NULL
  core_signature_fields <- c("signature_id", "task_type", "coefficients")
  looks_like_one_raw_record <- is.list(records) &&
    !is.data.frame(records) &&
    !is.null(record_names) &&
    all(core_signature_fields %in% record_names) &&
    is.character(records[["signature_id"]]) &&
    length(records[["signature_id"]]) == 1L &&
    is.character(records[["task_type"]]) &&
    length(records[["task_type"]]) == 1L &&
    is.numeric(records[["coefficients"]]) &&
    is.null(dim(records[["coefficients"]]))
  looks_like_one_record <- inherits(
    records,
    "stormbreaker_linear_signature"
  ) || looks_like_one_raw_record
  if (looks_like_one_record) {
    sb_abort(
      paste0(
        "records must be an outer named list of records; wrap one ",
        "specification as list(record_id = specification)."
      ),
      "stormbreaker_bad_signature_inventory"
    )
  }
  valid <- is.list(records) && !is.data.frame(records) &&
    length(records) > 0L &&
    !is.null(record_names) && !anyNA(record_names) &&
    all(nzchar(trimws(record_names))) &&
    anyDuplicated(record_names) == 0L
  if (!valid) {
    sb_abort(
      paste0(
        "records must be a non-empty list named by unique, non-empty ",
        "record IDs."
      ),
      "stormbreaker_bad_signature_inventory"
    )
  }

  rows <- lapply(seq_along(records), function(index) {
    report <- sb_lint_signature(records[[index]])
    data.frame(
      record_id = names(records)[[index]],
      report,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
