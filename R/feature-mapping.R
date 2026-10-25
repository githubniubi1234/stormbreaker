sb_validate_feature_map_targets <- function(target_features) {
  valid <- is.character(target_features) && length(target_features) > 0L &&
    !anyNA(target_features) && all(nzchar(trimws(target_features))) &&
    !anyDuplicated(target_features)
  if (!valid) {
    sb_abort(
      "`target_features` must be unique, non-empty character names.",
      "stormbreaker_bad_feature_map_targets"
    )
  }
  target_features
}

sb_validate_explicit_feature_map <- function(feature_map) {
  required <- c("source_feature", "target_feature")
  if (!is.data.frame(feature_map) || nrow(feature_map) < 1L) {
    sb_abort(
      "`feature_map` must be a non-empty data frame.",
      "stormbreaker_bad_feature_map"
    )
  }
  positions <- lapply(required, function(column) {
    which(names(feature_map) == column)
  })
  if (any(lengths(positions) != 1L)) {
    sb_abort(
      paste0(
        "`feature_map` must contain exactly one `source_feature` and one ",
        "`target_feature` column."
      ),
      "stormbreaker_bad_feature_map"
    )
  }
  source <- feature_map[[positions[[1L]]]]
  target <- feature_map[[positions[[2L]]]]
  valid_values <- is.character(source) && is.character(target) &&
    !anyNA(source) && !anyNA(target) &&
    all(nzchar(trimws(source))) && all(nzchar(trimws(target)))
  if (!valid_values) {
    sb_abort(
      "Feature-map source and target identifiers must be non-empty strings.",
      "stormbreaker_bad_feature_map"
    )
  }
  if (anyDuplicated(source) || anyDuplicated(target)) {
    sb_abort(
      paste0(
        "`feature_map` must be one-to-one; resolve one-to-many or many-to-one ",
        "biology before applying it."
      ),
      "stormbreaker_non_bijective_feature_map"
    )
  }
  feature_map
}

#' Apply an explicit one-to-one feature map without reading outcomes
#'
#' `sb_apply_feature_map()` renames and reorders numeric new-patient features
#' using a caller-supplied mapping table. It performs no database lookup or
#' annotation inference. Only the requested patient ID and mapped source
#' columns are read; outcome and unrelated columns are omitted from the returned
#' data.
#'
#' `feature_map` must contain one unique `source_feature` for every unique
#' `target_feature`, and its target set must exactly equal `target_features`.
#' The order of `target_features` defines the returned predictor order. Extra
#' mapping-table columns are retained as caller-supplied provenance.
#'
#' One-to-many and many-to-one mappings are deliberately rejected. The function
#' never chooses a first annotation, aggregates probes or transcripts, strips
#' identifier versions, changes case, normalizes expression, imputes values, or
#' fills an absent feature with zero. Those decisions require an explicit,
#' source-specific biological justification before this function is called.
#'
#' A successful result can supply its `data` element directly to
#' [sb_predict_final()] or [sb_predict_signature()]. Acceptance confirms only
#' that the supplied map is structurally one-to-one and executable on these
#' columns. It does not verify the species, annotation release, platform,
#' biological equivalence, expression scale, or cited source.
#'
#' @param new_data A non-empty data frame containing one unique patient ID and
#'   every source feature named in `feature_map`. Other columns are ignored.
#' @param id_col Name of the unique patient identifier column.
#' @param feature_map A non-empty data frame with exactly one column named
#'   `source_feature` and one named `target_feature`. Both must contain unique,
#'   non-empty character identifiers. Additional provenance columns are kept.
#' @param target_features Unique, non-empty target feature names. Their set must
#'   exactly match `feature_map$target_feature`, and their order controls the
#'   returned predictor order.
#'
#' @return A `stormbreaker_mapped_features` object containing the reordered
#'   mapping table and an outcome-free `data` frame with the ID followed by the
#'   mapped numeric target features.
#' @export
#'
#' @examples
#' new_data <- data.frame(
#'   patient_id = c("p1", "p2"),
#'   ENSG_B = c(2, 3),
#'   ENSG_A = c(0.1, 0.4)
#' )
#' feature_map <- data.frame(
#'   source_feature = c("ENSG_B", "ENSG_A"),
#'   target_feature = c("gene_b", "gene_a"),
#'   annotation_release = "caller_declared"
#' )
#' sb_apply_feature_map(
#'   new_data,
#'   "patient_id",
#'   feature_map,
#'   c("gene_a", "gene_b")
#' )
sb_apply_feature_map <- function(
    new_data,
    id_col,
    feature_map,
    target_features) {
  patient_ids <- sb_extract_final_prediction_ids(new_data, id_col)
  target_features <- sb_validate_feature_map_targets(target_features)
  feature_map <- sb_validate_explicit_feature_map(feature_map)

  source <- feature_map$source_feature
  target <- feature_map$target_feature
  if (id_col %in% c(source, target)) {
    sb_abort(
      "The patient ID column cannot be a mapped source or target feature.",
      "stormbreaker_overlapping_feature_map_roles"
    )
  }
  missing_targets <- setdiff(target_features, target)
  extra_targets <- setdiff(target, target_features)
  if (length(missing_targets) > 0L || length(extra_targets) > 0L) {
    sb_abort(
      paste0(
        "The feature-map targets must exactly match `target_features`; ",
        "missing: ",
        if (length(missing_targets) == 0L) "none" else {
          paste(missing_targets, collapse = ", ")
        },
        "; extra: ",
        if (length(extra_targets) == 0L) "none" else {
          paste(extra_targets, collapse = ", ")
        },
        "."
      ),
      "stormbreaker_feature_map_target_mismatch"
    )
  }

  order <- match(target_features, target)
  ordered_map <- feature_map[order, , drop = FALSE]
  source <- ordered_map$source_feature
  source_counts <- vapply(source, function(feature) {
    sum(names(new_data) == feature)
  }, integer(1L))
  if (any(source_counts == 0L)) {
    sb_abort(
      paste0(
        "Prediction data are missing mapped source features: ",
        paste(source[source_counts == 0L], collapse = ", "),
        "."
      ),
      "stormbreaker_missing_feature_map_source"
    )
  }
  if (any(source_counts > 1L)) {
    sb_abort(
      paste0(
        "Prediction data duplicate mapped source features: ",
        paste(source[source_counts > 1L], collapse = ", "),
        "."
      ),
      "stormbreaker_duplicate_feature_map_source"
    )
  }

  source_positions <- vapply(source, function(feature) {
    which(names(new_data) == feature)
  }, integer(1L))
  x <- new_data[source_positions]
  names(x) <- target_features
  x <- sb_validate_predictor_frame(x, "mapped source features")
  id_position <- which(names(new_data) == id_col)
  mapped_data <- cbind(new_data[id_position], x)
  names(mapped_data) <- c(id_col, target_features)
  rownames(mapped_data) <- NULL
  rownames(ordered_map) <- NULL

  structure(
    list(
      status = "success",
      method = "explicit_one_to_one_feature_mapping",
      id_col = id_col,
      n = as.integer(length(patient_ids)),
      mapping = ordered_map,
      data = mapped_data
    ),
    class = "stormbreaker_mapped_features"
  )
}

#' @export
print.stormbreaker_mapped_features <- function(x, ...) {
  cat("<stormbreaker_mapped_features:", x$status, ">\n")
  cat("explicit one-to-one feature mapping\n")
  cat("Patients:", x$n, "| Mapped features:", nrow(x$mapping), "\n")
  cat("No feature was inferred, aggregated, imputed, or filled with zero.\n")
  invisible(x)
}
