sb_read_cohort_csv <- function(path, id_col, label) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !file.exists(path) || dir.exists(path)) {
    sb_abort(paste0("Supply an existing ", label, " CSV file."),
      "stormbreaker_bad_cohort_file")
  }
  data <- tryCatch(utils::read.csv(path, colClasses = "character",
    check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM",
    na.strings = c("", "NA"), fill = FALSE), error = function(error) {
      sb_abort(paste0("Cannot read ", label, " CSV: ", conditionMessage(error)),
        "stormbreaker_bad_cohort_file")
    })
  if (anyDuplicated(names(data))) {
    sb_abort(paste0(label, " CSV has duplicate columns: ",
      paste(unique(names(data)[duplicated(names(data))]), collapse = ", "), "."),
      "stormbreaker_duplicate_column_names")
  }
  if (!id_col %in% names(data)) {
    sb_abort(paste0(label, " CSV is missing ID column '", id_col, "'."),
      "stormbreaker_missing_endpoint_columns")
  }
  if (!nrow(data) || anyNA(data[[id_col]]) || any(!nzchar(trimws(data[[id_col]])))) {
    sb_abort(paste0(label, " CSV has no samples or contains missing/blank IDs."),
      "stormbreaker_bad_id")
  }
  duplicated_ids <- unique(data[[id_col]][duplicated(data[[id_col]])])
  if (length(duplicated_ids)) {
    sb_abort(paste0(label, " CSV contains duplicate IDs: ",
      paste(utils::head(duplicated_ids, 10L), collapse = ", "), "."),
      "stormbreaker_duplicate_id")
  }
  for (column in setdiff(names(data), id_col)) {
    data[[column]] <- utils::type.convert(data[[column]], as.is = TRUE)
  }
  data
}

#' Read expression and clinical CSV files with an exact sample-ID join
#'
#' Read a sample-by-gene expression CSV and optionally join a clinical CSV.
#' IDs are preserved as text, including leading zeroes. Clinical rows are
#' aligned to expression order by ID, never by row position. Neither missing
#' measurements nor outcome labels are repaired or recoded.
#'
#' @param expression_file Path to a UTF-8 comma-separated file with one row
#'   per sample, an ID column and named gene columns. A combined file may also
#'   contain clinical columns; assign predictor roles explicitly in the task.
#' @param clinical_file Optional path to a UTF-8 clinical CSV with exactly the
#'   same unique sample IDs. NULL reads only the first file, suitable for
#'   prediction without outcomes.
#' @param id_col ID column present in both files; defaults to `"patient_id"`.
#' @return A data frame in expression-file row order. Pass explicit outcome,
#'   time and predictor names to [sb_binary_task()] or [sb_surv_task()].
#' @details Duplicate IDs or column names, missing/extra clinical IDs, and
#'   overlapping non-ID columns are errors. Numeric columns are read as
#'   numbers, while non-numeric labels remain text. The function does not
#'   normalize expression, map probes, impute values or infer the endpoint.
#'   Use study-specific processing before these modeling-ready CSV files.
#' @export
#' @examples
#' example_file <- system.file("extdata", "demo-binary-development.csv",
#'   package = "stormbreaker")
#' data <- sb_read_cohort(example_file)
#' sb_binary_task(data, "patient_id", "outcome", c("gene_a", "gene_b"))
sb_read_cohort <- function(expression_file, clinical_file = NULL,
                           id_col = "patient_id") {
  sb_check_name(id_col, "id_col")
  expression <- sb_read_cohort_csv(expression_file, id_col, "Expression")
  if (is.null(clinical_file)) return(expression)
  clinical <- sb_read_cohort_csv(clinical_file, id_col, "Clinical")
  missing <- setdiff(expression[[id_col]], clinical[[id_col]])
  extra <- setdiff(clinical[[id_col]], expression[[id_col]])
  if (length(missing) || length(extra)) {
    show_ids <- function(x) {
      if (length(x)) paste(utils::head(x, 10L), collapse = ", ") else "none"
    }
    sb_abort(paste0("Sample IDs do not match. Missing clinical IDs (", length(missing),
      "): ", show_ids(missing), "; extra clinical IDs (", length(extra), "): ",
      show_ids(extra), ". No rows were dropped."), "stormbreaker_unmatched_cohort_ids")
  }
  overlap <- intersect(setdiff(names(expression), id_col), setdiff(names(clinical), id_col))
  if (length(overlap)) {
    sb_abort(paste0("Both CSV files contain non-ID columns: ",
      paste(overlap, collapse = ", "), ". Keep each column in one file."),
      "stormbreaker_overlapping_cohort_columns")
  }
  clinical <- clinical[match(expression[[id_col]], clinical[[id_col]]), , drop = FALSE]
  cbind(expression, clinical[setdiff(names(clinical), id_col)])
}
