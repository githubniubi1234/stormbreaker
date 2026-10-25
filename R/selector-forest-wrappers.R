sb_boruta_selector_ids <- function() {
  c("boruta_binary_roughfix", "boruta_survival_roughfix")
}

sb_varselrf_selector_ids <- function() {
  "varselrf_binary_oob_1se"
}

sb_forest_selector_ids <- function() {
  c(sb_boruta_selector_ids(), sb_varselrf_selector_ids())
}

sb_selector_budget_options <- function(id, options = NULL) {
  defaults <- sb_forest_selector_params(id)[c("num_trees", "maximum_runs", "seed")]
  if (is.null(options) || identical(options, list())) return(list())
  if (!is.list(options) || is.null(names(options)) || anyDuplicated(names(options)) ||
      anyNA(names(options)) || any(!names(options) %in% names(defaults))) {
    sb_abort(paste0("selector_options for ", id, " accepts only ",
      paste(names(defaults), collapse = ", "), "."), "stormbreaker_bad_selector_options")
  }
  for (name in names(options)) {
    options[[name]] <- sb_check_positive_integer(options[[name]], name,
      minimum = if (name == "maximum_runs") 11L else 1L)
  }
  # Canonical order and omission of explicit defaults preserve existing IDs.
  options <- options[intersect(names(defaults), names(options))]
  options[!vapply(names(options), function(n) identical(options[[n]], defaults[[n]]), logical(1))]
}

sb_forest_selector_params <- function(selector) {
  selector <- sub("_multiclass_", "_binary_", selector, fixed = TRUE)
  if (selector %in% sb_boruta_selector_ids()) {
    return(list(
      method = "Boruta",
      importance = "ranger_normalized_permutation",
      p_value = 0.01,
      multiple_comparison_adjustment = TRUE,
      maximum_runs = 100L,
      num_trees = 500L,
      num_threads = 1L,
      tentative_rule = "TentativeRoughFix",
      seed = 1L
    ))
  }
  if (selector %in% sb_varselrf_selector_ids()) {
    return(list(
      method = "varSelRF",
      importance = "randomForest_mean_decrease_accuracy",
      c_sd = 1,
      mtry_factor = 1,
      initial_trees = 5000L,
      iteration_trees = 2000L,
      drop_fraction = 0.2,
      whole_range = TRUE,
      recompute_importance = FALSE,
      seed = 1L
    ))
  }
  sb_abort(
    paste0("Unknown forest selector `", selector, "`."),
    "stormbreaker_unsupported_selector"
  )
}

sb_forest_selector_dependency <- function(selector) {
  selector <- sub("_multiclass_", "_binary_", selector, fixed = TRUE)
  if (selector %in% sb_boruta_selector_ids()) {
    return(list(
      package = "Boruta",
      engine = "ranger",
      available = requireNamespace("Boruta", quietly = TRUE) &&
        requireNamespace("ranger", quietly = TRUE)
    ))
  }
  if (selector %in% sb_varselrf_selector_ids()) {
    return(list(
      package = "varSelRF",
      engine = "randomForest",
      available = requireNamespace("varSelRF", quietly = TRUE) &&
        requireNamespace("randomForest", quietly = TRUE)
    ))
  }
  sb_abort(
    paste0("Unknown forest selector `", selector, "`."),
    "stormbreaker_unsupported_selector"
  )
}

sb_forest_selector_provenance <- function(selector) {
  packages <- if (selector %in% sb_boruta_selector_ids()) {
    c("Boruta", "ranger")
  } else if (selector %in% sb_varselrf_selector_ids()) {
    c("varSelRF", "randomForest")
  } else {
    sb_abort(
      paste0("Unknown forest selector `", selector, "`."),
      "stormbreaker_unsupported_selector"
    )
  }
  fields <- paste0(tolower(packages), "_version")
  versions <- lapply(packages, function(package) {
    if (requireNamespace(package, quietly = TRUE)) {
      as.character(utils::packageVersion(package))
    } else {
      NA_character_
    }
  })
  stats::setNames(versions, fields)
}

sb_require_forest_selector_dependency <- function(selector) {
  dependency <- sb_forest_selector_dependency(selector)
  if (!dependency$available) {
    sb_abort(
      paste0(
        "Packages `", dependency$package, "` and `", dependency$engine,
        "` are required for selector `", selector, "`."
      ),
      "stormbreaker_missing_forest_selector_dependency"
    )
  }
  invisible(dependency)
}

sb_forest_selector_stop_on_warning <- function(code, selector) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      sb_abort(
        paste0(
          "Selector `", selector, "` raised a warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_forest_selector_warning",
        parent = warning
      )
    }
  )
}

sb_validate_forest_selector_inputs <- function(x, endpoint, selector) {
  if (selector %in% c("boruta_multiclass_roughfix", "varselrf_multiclass_oob_1se")) {
    return(sb_multiclass_selector_inputs(x, endpoint$outcome))
  }
  if (identical(selector, "boruta_survival_roughfix")) {
    inputs <- sb_validate_classical_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      selector
    )
    return(list(
      x = inputs$x,
      endpoint = list(time = inputs$time, event = inputs$event)
    ))
  }
  inputs <- sb_validate_classical_binary_inputs(
    x,
    endpoint$outcome,
    selector
  )
  list(x = inputs$x, endpoint = list(outcome = inputs$outcome))
}

sb_new_forest_selector_fit <- function(selector, source_features,
                                        selected_features, params, scores,
                                        selection_history = NULL,
                                        audit = list()) {
  valid <- is.character(source_features) && length(source_features) > 0L &&
    !anyNA(source_features) && !anyDuplicated(source_features) &&
    is.character(selected_features) &&
    (length(selected_features) > 0L || identical(selector, "boruta_multiclass_roughfix")) &&
    !anyNA(selected_features) && !anyDuplicated(selected_features) &&
    all(selected_features %in% source_features) && is.data.frame(scores) &&
    identical(scores$feature, source_features) &&
    identical(scores$selected, source_features %in% selected_features)
  if (!valid) {
    sb_abort(
      "The native forest selector did not return a valid feature subset.",
      "stormbreaker_invalid_forest_selector"
    )
  }
  structure(
    list(
      id = selector,
      source_features = source_features,
      selected_features = selected_features,
      params = params,
      scores = scores,
      selection_history = selection_history,
      audit = audit
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_call_boruta <- function(x, response, params) {
  Boruta::Boruta(
    x,
    response,
    getImp = Boruta::getImpRfZ,
    pValue = params$p_value,
    mcAdj = params$multiple_comparison_adjustment,
    maxRuns = params$maximum_runs,
    doTrace = 0L,
    holdHistory = TRUE,
    num.trees = params$num_trees,
    num.threads = params$num_threads
  )
}

sb_fit_boruta_selector <- function(x, endpoint, selector, capture_history = FALSE,
                                  options = list()) {
  sb_require_forest_selector_dependency(selector)
  params <- utils::modifyList(sb_forest_selector_params(selector), options)
  inputs <- sb_validate_forest_selector_inputs(x, endpoint, selector)
  source_features <- names(inputs$x)
  response <- if (identical(selector, "boruta_survival_roughfix")) {
    survival::Surv(inputs$endpoint$time, inputs$endpoint$event)
  } else if (identical(selector, "boruta_multiclass_roughfix")) {
    inputs$endpoint$outcome
  } else {
    factor(inputs$endpoint$outcome, levels = c(0L, 1L))
  }

  fitted <- sb_forest_selector_stop_on_warning(
    withr::with_seed(
      params$seed,
      sb_call_boruta(inputs$x, response, params)
    ),
    selector
  )
  runs <- nrow(fitted$ImpHistory)
  initial_decision <- fitted$finalDecision
  rough_fixed <- any(fitted$finalDecision == "Tentative")
  if (rough_fixed) {
    fitted <- sb_forest_selector_stop_on_warning(
      Boruta::TentativeRoughFix(fitted),
      selector
    )
  }
  statistics <- Boruta::attStats(fitted)
  valid_statistics <- is.data.frame(statistics) &&
    identical(rownames(statistics), source_features) &&
    all(c(
      "meanImp", "medianImp", "minImp", "maxImp", "normHits", "decision"
    ) %in% names(statistics)) &&
    !anyNA(statistics$decision) &&
    all(as.character(statistics$decision) %in% c("Confirmed", "Rejected"))
  if (!valid_statistics) {
    sb_abort(
      "Boruta did not return complete rough-fixed feature decisions.",
      "stormbreaker_invalid_boruta_selector"
    )
  }
  native_selected <- Boruta::getSelectedAttributes(
    fitted,
    withTentative = FALSE
  )
  selected_features <- source_features[source_features %in% native_selected]
  if (length(selected_features) < 1L && !identical(selector, "boruta_multiclass_roughfix")) {
    sb_abort(
      "Boruta did not confirm any predictor after its rough fix.",
      "stormbreaker_boruta_no_selected_features"
    )
  }
  scores <- data.frame(
    feature = source_features,
    mean_importance = as.numeric(statistics$meanImp),
    median_importance = as.numeric(statistics$medianImp),
    minimum_importance = as.numeric(statistics$minImp),
    maximum_importance = as.numeric(statistics$maxImp),
    normalized_hits = as.numeric(statistics$normHits),
    decision = as.character(statistics$decision),
    selected = source_features %in% selected_features,
    stringsAsFactors = FALSE
  )
  if (anyNA(scores[-c(1L, 7L, 8L)]) ||
    any(!is.finite(as.matrix(scores[c(
      "mean_importance", "median_importance", "minimum_importance",
      "maximum_importance", "normalized_hits"
    )])))) {
    sb_abort(
      "Boruta returned non-finite feature-importance statistics.",
      "stormbreaker_invalid_boruta_selector"
    )
  }
  decision_counts <- table(factor(
    scores$decision,
    levels = c("Confirmed", "Rejected")
  ))
  result <- sb_new_forest_selector_fit(
    selector,
    source_features,
    selected_features,
    params,
    scores,
    audit = list(
      completed_runs = as.integer(runs),
      tentative_rough_fix_applied = rough_fixed,
      confirmed = as.integer(decision_counts[["Confirmed"]]),
      rejected = as.integer(decision_counts[["Rejected"]])
    )
  )
  if (capture_history) result$diagnostic_history <- list(
    importance = fitted$ImpHistory,
    decision_before_rough_fix = initial_decision,
    decision_after_rough_fix = fitted$finalDecision)
  result
}

sb_fit_varselrf_selector <- function(x, outcome, selector) {
  sb_require_forest_selector_dependency(selector)
  params <- sb_forest_selector_params(selector)
  inputs <- sb_validate_forest_selector_inputs(
    x,
    list(outcome = outcome),
    selector
  )
  source_features <- names(inputs$x)
  if (length(source_features) < 2L) {
    sb_abort(
      "varSelRF requires at least two predictors.",
      "stormbreaker_varselrf_too_few_predictors"
    )
  }
  fitted <- sb_forest_selector_stop_on_warning(
    withr::with_seed(params$seed, varSelRF::varSelRF(
      xdata = inputs$x,
      Class = if (identical(selector, "varselrf_multiclass_oob_1se"))
        inputs$endpoint$outcome else factor(inputs$endpoint$outcome, levels = c(0L, 1L)),
      c.sd = params$c_sd,
      mtryFactor = params$mtry_factor,
      ntree = params$initial_trees,
      ntreeIterat = params$iteration_trees,
      vars.drop.num = NULL,
      vars.drop.frac = params$drop_fraction,
      whole.range = params$whole_range,
      recompute.var.imp = params$recompute_importance,
      verbose = FALSE,
      returnFirstForest = FALSE,
      keep.forest = FALSE
    )),
    selector
  )
  selected_features <- fitted$selected.vars
  valid_selected <- is.character(selected_features) &&
    length(selected_features) >= 2L && !anyNA(selected_features) &&
    !anyDuplicated(selected_features) &&
    all(selected_features %in% source_features)
  importance <- as.matrix(fitted$initialImportances)
  valid_importance <- ncol(importance) == 1L &&
    identical(rownames(importance), source_features) &&
    !anyNA(importance) && all(is.finite(importance))
  history <- fitted$selec.history
  required_history <- c("Number.Variables", "OOB", "sd.OOB")
  valid_history <- is.data.frame(history) && nrow(history) > 0L &&
    all(required_history %in% names(history)) &&
    all(is.finite(as.matrix(history[required_history])))
  if (!valid_selected || !valid_importance || !valid_history) {
    sb_abort(
      "varSelRF did not return a valid OOB-selected feature subset.",
      "stormbreaker_invalid_varselrf_selector"
    )
  }
  scores <- data.frame(
    feature = source_features,
    mean_decrease_accuracy = as.numeric(importance[, 1L]),
    selected = source_features %in% selected_features,
    stringsAsFactors = FALSE
  )
  selection_history <- data.frame(
    number_of_variables = as.integer(history$Number.Variables),
    oob_error = as.numeric(history$OOB),
    oob_standard_error = as.numeric(history$sd.OOB),
    stringsAsFactors = FALSE
  )
  sb_new_forest_selector_fit(
    selector,
    source_features,
    selected_features,
    params,
    scores,
    selection_history = selection_history,
    audit = list(
      selection_rule = "smallest_subset_within_one_sd_of_minimum_oob_error",
      best_model_n_features = as.integer(fitted$best.model.nvars)
    )
  )
}

sb_apply_forest_selector <- function(fitted, x) {
  valid <- inherits(fitted, "stormbreaker_selector_fit") &&
    fitted$id %in% sb_forest_selector_ids() &&
    is.character(fitted$source_features) &&
    is.character(fitted$selected_features) &&
    length(fitted$selected_features) > 0L &&
    !anyNA(fitted$selected_features) &&
    !anyDuplicated(fitted$selected_features) &&
    all(fitted$selected_features %in% fitted$source_features)
  if (!valid) {
    sb_abort(
      "The fitted selector is not a valid forest-selector state.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$source_features)) {
    sb_abort(
      "Predictor columns do not match the fitted forest selector schema.",
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x[fitted$selected_features]
}

sb_forest_selector_adapter <- function(selector, capture_history = FALSE, options = list()) {
  if (is.null(options)) options <- list()
  if (!selector %in% sb_forest_selector_ids()) {
    sb_abort(
      paste0("Unknown forest selector `", selector, "`."),
      "stormbreaker_unsupported_selector"
    )
  }
  fit <- if (identical(selector, "boruta_survival_roughfix")) {
    function(x, time, event) {
      sb_fit_boruta_selector(x, list(time = time, event = event), selector, capture_history, options)
    }
  } else if (identical(selector, "boruta_binary_roughfix")) {
    function(x, outcome) {
      sb_fit_boruta_selector(x, list(outcome = outcome), selector, capture_history, options)
    }
  } else {
    function(x, outcome) {
      sb_fit_varselrf_selector(x, outcome, selector)
    }
  }
  list(
    id = selector,
    params = c(options, if (capture_history) list(capture_history = TRUE) else list()),
    needs_outcome = TRUE,
    fit = fit,
    apply = sb_apply_forest_selector
  )
}
