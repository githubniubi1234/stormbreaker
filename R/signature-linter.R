sb_signature_from_lint_spec <- function(spec) {
  if (inherits(spec, "stormbreaker_linear_signature")) {
    return(sb_validate_linear_signature(spec))
  }
  if (!is.list(spec)) {
    sb_abort(
      paste0(
        "spec must be one stormbreaker_linear_signature or one named ",
        "list of arguments for sb_linear_signature()."
      ),
      "stormbreaker_bad_signature_spec"
    )
  }
  if (length(spec) > 0L && (is.null(names(spec)) ||
      anyNA(names(spec)) || any(!nzchar(names(spec))) ||
      anyDuplicated(names(spec)) > 0L)) {
    sb_abort(
      "A signature argument list must have unique, non-empty names.",
      "stormbreaker_bad_signature_spec"
    )
  }
  unknown <- setdiff(names(spec), names(formals(sb_linear_signature)))
  if (length(unknown) > 0L) {
    sb_abort(
      paste0(
        "Unknown sb_linear_signature() argument",
        if (length(unknown) > 1L) "s" else "",
        ": ",
        paste0(unknown, collapse = ", "),
        "."
      ),
      "stormbreaker_unknown_signature_argument"
    )
  }
  required <- c("signature_id", "task_type", "coefficients")
  missing <- setdiff(required, names(spec))
  if (length(missing) > 0L) {
    sb_abort(
      paste0(
        "Missing required sb_linear_signature() argument",
        if (length(missing) > 1L) "s" else "",
        ": ",
        paste0(missing, collapse = ", "),
        "."
      ),
      "stormbreaker_missing_signature_argument"
    )
  }
  do.call(sb_linear_signature, spec)
}

sb_is_expected_signature_lint_error <- function(error) {
  inherits(error, c(
    "stormbreaker_bad_signature_spec",
    "stormbreaker_unknown_signature_argument",
    "stormbreaker_missing_signature_argument",
    "stormbreaker_bad_linear_signature",
    "stormbreaker_bad_signature_id",
    "stormbreaker_bad_signature_task",
    "stormbreaker_bad_signature_coefficients",
    "stormbreaker_missing_signature_intercept",
    "stormbreaker_bad_signature_intercept",
    "stormbreaker_bad_signature_metadata",
    "stormbreaker_bad_signature_horizon",
    "stormbreaker_bad_signature_baseline_survival",
    "stormbreaker_bad_prediction_horizon",
    "stormbreaker_bad_time_unit"
  ))
}

sb_signature_lint_report <- function(status, spec, signature = NULL,
                                     error = NULL) {
  source_object <- if (is.list(signature)) signature else spec
  safe_character <- function(name) {
    if (!is.list(source_object) || !name %in% names(source_object)) {
      return(NA_character_)
    }
    value <- source_object[[name]]
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(trimws(value))) {
      return(NA_character_)
    }
    value
  }
  coefficients <- if (is.list(source_object)) {
    source_object[["coefficients"]]
  } else {
    NULL
  }
  n_features <- if (is.numeric(coefficients) && is.null(dim(coefficients)) &&
      length(coefficients) > 0L) {
    as.integer(length(coefficients))
  } else {
    NA_integer_
  }
  output <- if (is.list(signature) && is.list(signature[["prediction"]])) {
    paste(signature[["prediction"]][["output"]], collapse = "+")
  } else {
    NA_character_
  }

  data.frame(
    status = status,
    signature_id = safe_character("signature_id"),
    task_type = safe_character("task_type"),
    n_features = n_features,
    output = output,
    source = safe_character("source"),
    issue_code = if (is.null(error)) NA_character_ else class(error)[[1L]],
    message = if (is.null(error)) NA_character_ else conditionMessage(error),
    stringsAsFactors = FALSE
  )
}

#' Check whether one fixed linear signature is complete enough to execute
#'
#' `sb_lint_signature()` checks one [sb_linear_signature()] object or one named
#' list of arguments for that constructor. It does not inspect patient data,
#' read outcomes, fit a model, verify a citation, or infer omitted formula or
#' preprocessing details. Expected specification problems are retained as a
#' one-row `rejected` report rather than disappearing when construction fails.
#' Unexpected internal errors remain errors.
#'
#' `admitted` means only that STORMBREAKER can represent the stated formula.
#' It does not establish that the transcription is faithful to a publication,
#' that input identifiers or scales are compatible, that a cohort is suitable,
#' or that the signature has predictive or clinical validity.
#'
#' @param spec One `stormbreaker_linear_signature` object or one named list of
#'   arguments accepted by [sb_linear_signature()].
#'
#' @return A one-row data frame containing admission status, safely reportable
#'   signature metadata, output type, and an issue code and message for a
#'   rejected specification.
#' @export
#'
#' @examples
#' sb_lint_signature(list(
#'   signature_id = "complete_binary_formula",
#'   task_type = "binary",
#'   coefficients = c(gene_a = 0.7),
#'   intercept = -0.2
#' ))
#'
#' sb_lint_signature(list(
#'   signature_id = "intercept_not_reported",
#'   task_type = "binary",
#'   coefficients = c(gene_a = 0.7)
#' ))
sb_lint_signature <- function(spec) {
  signature <- NULL
  result <- tryCatch(
    {
      signature <- sb_signature_from_lint_spec(spec)
      signature
    },
    stormbreaker_error = function(error) {
      if (!sb_is_expected_signature_lint_error(error)) {
        stop(error)
      }
      error
    }
  )

  if (inherits(result, "stormbreaker_error")) {
    return(sb_signature_lint_report(
      status = "rejected",
      spec = spec,
      signature = signature,
      error = result
    ))
  }
  sb_signature_lint_report(
    status = "admitted",
    spec = spec,
    signature = result
  )
}
