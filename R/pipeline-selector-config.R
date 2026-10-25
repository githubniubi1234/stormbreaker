sb_pipeline_selector_config <- function(selector, admission, uses_inner_cv, options) {
  is_mime_selector <- selector %in% sb_mime_survival_selector_ids()
  is_mime_rsf_selector <- selector %in% sb_mime_rsf_selector_ids()
  is_mime_tuned_selector <- selector %in% sb_mime_tuned_selector_ids()
  is_mime_composite_selector <- selector %in%
    sb_mime_composite_selector_ids()
  requires_top_k <- identical(
    admission$selector_parameter[[1L]],
    "selector_n_features"
  )
  requires_fixed_features <- identical(
    admission$selector_parameter[[1L]],
    "selector_features_and_feature_set_id"
  )
  requires_stepcox <- identical(
    admission$selector_parameter[[1L]],
    "stepcox_direction"
  ) || is_mime_composite_selector
  if (!requires_top_k && !is.null(options$selector_n_features)) {
    sb_abort(
      "`selector_n_features` applies only to a Top-K selector.",
      "stormbreaker_irrelevant_selector_parameter"
    )
  }
  if (requires_top_k && is.null(options$selector_n_features)) {
    sb_abort(
      paste0(
        "`selector_n_features` is required when ",
        "a Top-K selector is requested."
      ),
      "stormbreaker_missing_selector_parameter"
    )
  }
  has_fixed_parameter <- !is.null(options$selector_features) ||
    !is.null(options$feature_set_id)
  if (!requires_fixed_features && has_fixed_parameter) {
    sb_abort(
      paste0(
        "`selector_features` and `feature_set_id` apply only when ",
        "`selector = 'fixed_features'`."
      ),
      "stormbreaker_irrelevant_selector_parameter"
    )
  }
  if (requires_fixed_features &&
    (is.null(options$selector_features) || is.null(options$feature_set_id))) {
    sb_abort(
      paste0(
        "`selector_features` and `feature_set_id` are both required when ",
        "`selector = 'fixed_features'`."
      ),
      "stormbreaker_missing_selector_parameter"
    )
  }
  if (!requires_stepcox && !is.null(options$stepcox_direction)) {
    sb_abort(
      paste0(
        "`stepcox_direction` applies only to `stepcox` or a sequential ",
        "MIME selector ending in StepCox."
      ),
      "stormbreaker_irrelevant_selector_parameter"
    )
  }
  if (requires_stepcox && is.null(options$stepcox_direction)) {
    sb_abort(
      "`stepcox_direction` is required when `selector = 'stepcox'`.",
      "stormbreaker_missing_selector_parameter"
    )
  }
  has_mime_selector_control <- !is.null(options$selector_inner_folds) ||
    !is.null(options$selector_seed) || !is.null(options$selector_rsf_nodesize)
  if (!is_mime_selector && has_mime_selector_control) {
    sb_abort(
      paste0(
        "`selector_inner_folds`, `selector_seed`, and ",
        "`selector_rsf_nodesize` apply only to MIME survival selectors."
      ),
      "stormbreaker_irrelevant_selector_parameter"
    )
  }
  if (is_mime_selector && is.null(options$selector_seed)) {
    sb_abort(
      "`selector_seed` is required for a MIME survival selector.",
      "stormbreaker_missing_selector_parameter"
    )
  }
  if (is_mime_rsf_selector && is.null(options$selector_rsf_nodesize)) {
    sb_abort(
      "`selector_rsf_nodesize` is required for an RSF selector.",
      "stormbreaker_missing_selector_parameter"
    )
  }
  if (!is_mime_rsf_selector && !is.null(options$selector_rsf_nodesize)) {
    sb_abort(
      "`selector_rsf_nodesize` applies only to an RSF selector.",
      "stormbreaker_irrelevant_selector_parameter"
    )
  }
  if (is_mime_tuned_selector && is.null(options$selector_inner_folds)) {
    sb_abort(
      paste0(
        "`selector_inner_folds` is required for a CoxBoost or Lasso ",
        "selector."
      ),
      "stormbreaker_missing_selector_parameter"
    )
  }
  if (!is_mime_tuned_selector && !is.null(options$selector_inner_folds)) {
    sb_abort(
      paste0(
        "`selector_inner_folds` applies only to a CoxBoost or Lasso ",
        "selector."
      ),
      "stormbreaker_irrelevant_selector_parameter"
    )
  }
  selector_params <- list()
  budget_options <- list()
  if (selector %in% sb_pamr_selector_ids()) selector_params <- sb_pamr_selector_params(selector)
  if (selector %in% sb_survival_weight_selector_ids()) selector_params <- sb_survival_weight_selector_params(selector)
  if (selector %in% sb_survival_boosting_selector_ids()) selector_params <- sb_survival_boosting_selector_params(selector)
  if (selector %in% sb_binary_model_selector_ids()) {
    selector_params <- sb_binary_model_selector_options(selector, options$selector_options)
  } else if (selector %in% sb_boruta_selector_ids()) {
    budget_options <- sb_selector_budget_options(selector, options$selector_options)
    selector_params <- utils::modifyList(selector_params, budget_options)
  } else if (!is.null(options$selector_options)) sb_abort(
    "selector_options applies only to binary model and Boruta selectors.", "stormbreaker_bad_selector_options")
  if (requires_top_k) {
    selector_params$n_features <- sb_check_positive_integer(
      options$selector_n_features,
      "selector_n_features",
      minimum = 1L
    )
  }
  if (requires_fixed_features) {
    selector_params <- list(
      features = sb_validate_fixed_features(options$selector_features),
      feature_set_id = sb_validate_feature_set_id(options$feature_set_id),
      n_features = as.integer(length(options$selector_features))
    )
  }
  if (requires_stepcox && !is_mime_selector) {
    selector_params <- list(
      direction = sb_validate_stepcox_direction(options$stepcox_direction)
    )
  }
  if (is_mime_rsf_selector) {
    selector_params <- list(
      nodesize = sb_check_positive_integer(
        options$selector_rsf_nodesize,
        "selector_rsf_nodesize",
        minimum = 1L
      ),
      seed = sb_check_positive_integer(
        options$selector_seed,
        "selector_seed",
        minimum = 1L
      )
    )
  } else if (is_mime_tuned_selector) {
    selector_params <- list(
      inner_folds = sb_check_positive_integer(
        options$selector_inner_folds,
        "selector_inner_folds",
        minimum = 3L
      ),
      seed = sb_check_positive_integer(
        options$selector_seed,
        "selector_seed",
        minimum = 1L
      )
    )
  }
  if (is_mime_composite_selector) {
    selector_params$direction <- sb_validate_stepcox_direction(
      options$stepcox_direction
    )
  }
  if (uses_inner_cv && "n_features" %in% names(selector_params) &&
    selector_params$n_features < 2L) {
    sb_abort(
      paste0(
        "A glmnet pipeline must retain at least 2 features; ",
        "glmnet cannot fit a path with only one variable."
      ),
      "stormbreaker_glmnet_too_few_selected_features"
    )
  }
  selector_suffix <- if (length(budget_options)) {
    paste0("__", paste(paste0(names(budget_options), unlist(budget_options)), collapse = "__"))
  } else if (selector %in% sb_binary_model_selector_ids()) {
    paste0("__", paste(vapply(names(selector_params), function(n) {
      v <- selector_params[[n]]
      paste0(n, if (is.numeric(v)) sb_format_roundtrip_double(as.numeric(v)) else v)
    }, character(1)), collapse = "__"))
  } else if (requires_top_k) {
    paste0("__k", selector_params$n_features)
  } else if (requires_fixed_features) {
    paste0("__set", selector_params$feature_set_id)
  } else if (is_mime_rsf_selector) {
    paste0(
      "__nodesize", selector_params$nodesize,
      "__seed", selector_params$seed,
      if (is_mime_composite_selector) {
        paste0("__direction", selector_params$direction)
      } else {
        ""
      }
    )
  } else if (is_mime_tuned_selector) {
    paste0(
      "__inner", selector_params$inner_folds,
      "__seed", selector_params$seed,
      if (is_mime_composite_selector) {
        paste0("__direction", selector_params$direction)
      } else {
        ""
      }
    )
  } else if (requires_stepcox) {
    paste0("__direction", selector_params$direction)
  } else {
    ""
  }
  list(params = selector_params, suffix = selector_suffix)
}
