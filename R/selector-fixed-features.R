sb_validate_fixed_features <- function(features) {
  valid <- is.character(features) && length(features) > 0L &&
    !anyNA(features) && all(nzchar(trimws(features))) &&
    anyDuplicated(features) == 0L
  if (!valid) {
    sb_abort(
      "`selector_features` must contain unique, non-empty feature names.",
      "stormbreaker_bad_fixed_features"
    )
  }
  features
}

sb_validate_feature_set_id <- function(feature_set_id) {
  valid <- is.character(feature_set_id) && length(feature_set_id) == 1L &&
    !is.na(feature_set_id) && nzchar(trimws(feature_set_id))
  if (!valid) {
    sb_abort(
      "`feature_set_id` must be one non-empty character value.",
      "stormbreaker_bad_feature_set_id"
    )
  }
  feature_set_id
}

sb_fit_fixed_features <- function(x, features, feature_set_id) {
  x <- sb_validate_predictor_frame(x, "x")
  features <- sb_validate_fixed_features(features)
  feature_set_id <- sb_validate_feature_set_id(feature_set_id)
  missing <- setdiff(features, names(x))
  if (length(missing) > 0L) {
    sb_abort(
      paste0(
        "The predeclared feature set is missing from the task predictors: ",
        paste(missing, collapse = ", "),
        "."
      ),
      "stormbreaker_missing_fixed_feature"
    )
  }
  structure(
    list(
      id = "fixed_features",
      feature_set_id = feature_set_id,
      source_features = names(x),
      selected_features = features
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_apply_fixed_features <- function(fitted, x) {
  valid_fit <- inherits(fitted, "stormbreaker_selector_fit") &&
    identical(fitted$id, "fixed_features") &&
    is.character(fitted$feature_set_id) &&
    length(fitted$feature_set_id) == 1L &&
    is.character(fitted$source_features) &&
    is.character(fitted$selected_features) &&
    length(fitted$selected_features) > 0L
  if (!valid_fit) {
    sb_abort(
      "The fitted selector is not a valid fixed-feature state.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$source_features)) {
    sb_abort(
      paste0(
        "Predictor columns must exactly match the fitted selector schema: ",
        paste(fitted$source_features, collapse = ", "),
        "."
      ),
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x[fitted$selected_features]
}

sb_fixed_features_selector_adapter <- function(features, feature_set_id) {
  features <- sb_validate_fixed_features(features)
  feature_set_id <- sb_validate_feature_set_id(feature_set_id)
  list(
    id = "fixed_features",
    params = list(
      features = features,
      feature_set_id = feature_set_id,
      n_features = as.integer(length(features))
    ),
    needs_outcome = FALSE,
    fit = function(x) {
      sb_fit_fixed_features(x, features, feature_set_id)
    },
    apply = sb_apply_fixed_features
  )
}

sb_selector_is_learned <- function(selector) {
  selector %in% c(
    "ordered_chain",
    sb_survival_weight_selector_ids(),
    sb_survival_boosting_selector_ids(),
    sb_pamr_selector_ids(),
    sb_binary_model_selector_ids(),
    "univ_cox_top_k",
    "univ_binary_auc_top_k",
    "stepcox",
    "rsf_minimal_depth",
    "coxboost_nonzero",
    "lasso_nonzero",
    "enet_cox_nonzero",
    "rsf_then_stepcox",
    "coxboost_then_stepcox",
    "lasso_then_stepcox",
    sb_abess_selector_ids(),
    sb_forest_selector_ids(),
    sb_high_dimensional_selector_ids()
  )
}
