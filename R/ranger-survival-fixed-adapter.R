sb_ranger_binary_learner_ids <- function() {
  c(
    "ranger_binary",
    "ranger_extratrees_binary",
    "ranger_hellinger_binary"
  )
}

sb_ranger_survival_learner_ids <- function() {
  c(
    "ranger_survival",
    "ranger_extratrees_survival",
    "ranger_cindex_survival",
    "ranger_maxstat_survival"
  )
}

sb_ranger_learner_ids <- function() {
  c(sb_ranger_binary_learner_ids(), sb_ranger_survival_learner_ids())
}

sb_ranger_splitrule_spec <- function(learner) {
  switch(learner,
    ranger_binary = list(splitrule = "gini", replace = TRUE),
    ranger_extratrees_binary = list(
      splitrule = "extratrees",
      replace = FALSE,
      num_random_splits = 1L
    ),
    ranger_hellinger_binary = list(
      splitrule = "hellinger",
      replace = TRUE
    ),
    ranger_survival = list(splitrule = "logrank", replace = TRUE),
    ranger_extratrees_survival = list(
      splitrule = "extratrees",
      replace = FALSE,
      num_random_splits = 1L
    ),
    ranger_cindex_survival = list(splitrule = "C", replace = TRUE),
    ranger_maxstat_survival = list(
      splitrule = "maxstat",
      replace = TRUE,
      alpha = 0.5,
      minprop = 0.1
    ),
    sb_abort(
      paste0("Unknown fixed ranger learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_ranger_parameter_semantics <- function(learner) {
  switch(learner,
    ranger_binary = paste0(
      "Gini probability-forest trees, mtry rule, node size, bootstrap ",
      "sampling, and seed fixed before resampling"
    ),
    ranger_extratrees_binary = paste0(
      "Extremely randomized probability-forest split locations, one random ",
      "split per candidate variable, full no-replacement sampling, mtry rule, ",
      "node size, and seed fixed before resampling"
    ),
    ranger_hellinger_binary = paste0(
      "Hellinger-distance probability-forest trees, mtry rule, node size, ",
      "bootstrap sampling, and seed fixed before resampling"
    ),
    ranger_survival = paste0(
      "logrank survival-forest trees, mtry rule, node size, bootstrap ",
      "sampling, and seed fixed before resampling"
    ),
    ranger_extratrees_survival = paste0(
      "extremely randomized survival-forest split locations, one random ",
      "split per candidate variable, full no-replacement sampling, mtry rule, ",
      "node size, and seed fixed before resampling"
    ),
    ranger_cindex_survival = paste0(
      "C-index survival-forest split rule, mtry rule, node size, bootstrap ",
      "sampling, and seed fixed before resampling"
    ),
    ranger_maxstat_survival = paste0(
      "maximally selected rank-statistic survival-forest split rule, alpha ",
      "0.5, minprop 0.1, mtry rule, node size, bootstrap sampling, and seed ",
      "fixed before resampling"
    ),
    sb_abort(
      paste0("Unknown fixed ranger learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_ranger_available <- function() {
  requireNamespace("ranger", quietly = TRUE)
}

sb_validate_ranger_mtry <- function(mtry) {
  if (is.null(mtry)) {
    return(NULL)
  }
  sb_check_positive_integer(mtry, "ranger_mtry", minimum = 1L)
}

sb_ranger_survival_params <- function(
    num_trees = 500L,
    mtry = NULL,
    min_node_size = 3L,
    seed = 1L,
    learner = "ranger_survival") {
  if (!learner %in% sb_ranger_survival_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed ranger survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  mtry <- sb_validate_ranger_mtry(mtry)
  split_spec <- sb_ranger_splitrule_spec(learner)
  params <- list(
    num_trees = sb_check_positive_integer(
      num_trees,
      "ranger_num_trees",
      minimum = 1L
    ),
    mtry = mtry,
    mtry_rule = if (is.null(mtry)) "floor_sqrt_p" else "fixed",
    min_node_size = sb_check_positive_integer(
      min_node_size,
      "ranger_min_node_size",
      minimum = 1L
    ),
    seed = sb_check_positive_integer(
      seed,
      "ranger_seed",
      minimum = 1L
    ),
    splitrule = split_spec$splitrule,
    importance = "impurity",
    replace = split_spec$replace,
    sample_fraction = 1,
    max_depth = 0L,
    num_threads = 1L,
    oob_error = FALSE,
    write_forest = TRUE,
    na_action = "na.fail"
  )
  for (name in setdiff(names(split_spec), c("splitrule", "replace"))) {
    params[[name]] <- split_spec[[name]]
  }
  params
}

sb_fit_ranger_survival <- function(
  x, time, event, params, learner = "ranger_survival") {
  if (!learner %in% sb_ranger_survival_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed ranger survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time < 0)) {
    sb_abort(
      "`time` must contain non-negative values.",
      "stormbreaker_bad_adapter_time"
    )
  }
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }
  if (!any(event == 1L)) {
    sb_abort(
      "A fixed ranger survival analysis fold must contain at least one event.",
      "stormbreaker_ranger_survival_no_event"
    )
  }
  if (!sb_ranger_available()) {
    sb_abort(
      "Package `ranger` is required for the fixed survival forest learner.",
      "stormbreaker_missing_ranger"
    )
  }

  effective_mtry <- if (is.null(params$mtry)) {
    max(1L, floor(sqrt(ncol(x))))
  } else {
    params$mtry
  }
  if (effective_mtry > ncol(x)) {
    sb_abort(
      "`ranger_mtry` cannot exceed the number of fitted features.",
      "stormbreaker_ranger_mtry_exceeds_features"
    )
  }

  if (identical(learner, "ranger_cindex_survival")) {
    threads <- getOption("stormbreaker.ranger_cindex_threads", params$num_threads)
    if (!is.numeric(threads) || length(threads) != 1L || is.na(threads) ||
      !threads %in% 1:4) {
      sb_abort("`stormbreaker.ranger_cindex_threads` must be an integer from 1 to 4.",
        "stormbreaker_bad_ranger_threads")
    }
    params$num_threads <- as.integer(threads)
  }
  fit_arguments <- list(
    x = as.matrix(x),
    y = survival::Surv(time, event),
    num.trees = params$num_trees,
    mtry = effective_mtry,
    importance = params$importance,
    write.forest = params$write_forest,
    min.node.size = params$min_node_size,
    max.depth = params$max_depth,
    replace = params$replace,
    sample.fraction = params$sample_fraction,
    splitrule = params$splitrule,
    oob.error = params$oob_error,
    num.threads = params$num_threads,
    seed = params$seed,
    na.action = params$na_action,
    verbose = FALSE
  )
  if ("num_random_splits" %in% names(params)) {
    fit_arguments$num.random.splits <- params$num_random_splits
  }
  if ("alpha" %in% names(params)) {
    fit_arguments$alpha <- params$alpha
  }
  if ("minprop" %in% names(params)) {
    fit_arguments$minprop <- params$minprop
  }
  fitted <- withCallingHandlers(
    do.call(ranger::ranger, fit_arguments),
    warning = function(warning) {
      message <- conditionMessage(warning)
      if (!startsWith(message, "Unused arguments:") && !identical(message,
        "Argument 'num.random.splits' ignored if splitrule is not 'extratrees'.")) {
        return(invisible(NULL))
      }
      sb_abort(
        paste0(
          "The fixed ranger survival learner raised a warning: ",
          message
        ),
        "stormbreaker_ranger_survival_warning",
        parent = warning
      )
    }
  )
  valid_forest <- inherits(fitted$forest, "ranger.forest") &&
    identical(as.integer(fitted$num.trees), params$num_trees) &&
    identical(fitted$treetype, "Survival")
  if (!valid_forest) {
    sb_abort(
      "Fixed ranger survival did not return the declared forest.",
      "stormbreaker_incomplete_ranger_fit"
    )
  }
  importance <- fitted$variable.importance
  valid_importance <- is.numeric(importance) &&
    length(importance) == ncol(x) && !anyNA(importance) &&
    all(is.finite(importance)) &&
    setequal(names(importance), names(x))
  if (!valid_importance) {
    sb_abort(
      "Fixed ranger survival did not return one finite importance per feature.",
      "stormbreaker_invalid_ranger_importance"
    )
  }
  importance <- unname(importance[names(x)])

  structure(
    list(
      learner = learner,
      features = names(x),
      params = params,
      num_trees = params$num_trees,
      effective_mtry = as.integer(effective_mtry),
      forest = fitted$forest,
      feature_importance = importance
    ),
    class = "stormbreaker_ranger_survival"
  )
}

sb_predict_ranger_survival_curves <- function(model, x) {
  if (!inherits(model, "stormbreaker_ranger_survival")) {
    sb_abort(
      "The supplied model is not a fitted fixed ranger survival learner.",
      "stormbreaker_bad_ranger_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fixed ranger survival schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  if (!sb_ranger_available()) {
    sb_abort(
      "Package `ranger` is required for fixed survival forest prediction.",
      "stormbreaker_missing_ranger"
    )
  }

  prediction <- stats::predict(
    model$forest,
    data = x,
    num.trees = model$num_trees,
    type = "response",
    seed = 0L,
    num.threads = 1L,
    verbose = FALSE
  )
  chf <- if (is.null(dim(prediction$chf))) {
    matrix(prediction$chf, nrow = nrow(x), byrow = TRUE)
  } else {
    as.matrix(prediction$chf)
  }
  survival_probability <- if (is.null(dim(prediction$survival))) {
    matrix(prediction$survival, nrow = nrow(x), byrow = TRUE)
  } else {
    as.matrix(prediction$survival)
  }
  times <- as.numeric(prediction$unique.death.times)
  valid <- nrow(chf) == nrow(x) && nrow(survival_probability) == nrow(x) &&
    ncol(chf) == length(times) &&
    identical(dim(chf), dim(survival_probability)) &&
    length(times) > 0L && !anyNA(times) && all(is.finite(times)) &&
    !anyNA(chf) && all(is.finite(chf)) &&
    !anyNA(survival_probability) && all(is.finite(survival_probability)) &&
    all(survival_probability >= 0 & survival_probability <= 1)
  if (!valid) {
    sb_abort(
      "Fixed ranger survival produced invalid prediction curves.",
      "stormbreaker_invalid_ranger_prediction"
    )
  }
  risk <- unname(rowSums(chf))
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "Fixed ranger survival produced invalid mortality predictions.",
      "stormbreaker_invalid_ranger_prediction"
    )
  }
  list(
    risk = risk,
    time = times,
    survival = survival_probability
  )
}

sb_ranger_survival_probability_at_horizon <- function(curves,
                                                       prediction_horizon) {
  eligible <- which(curves$time <= prediction_horizon)
  probability <- if (length(eligible) == 0L) {
    rep(1, length(curves$risk))
  } else {
    unname(curves$survival[, max(eligible)])
  }
  if (length(probability) != length(curves$risk) || anyNA(probability) ||
      any(!is.finite(probability)) ||
      any(probability < 0 | probability > 1)) {
    sb_abort(
      "Fixed ranger survival produced invalid horizon probabilities.",
      "stormbreaker_invalid_ranger_prediction"
    )
  }
  probability
}

sb_predict_ranger_survival <- function(model, x) {
  sb_predict_ranger_survival_curves(model, x)$risk
}

sb_predict_ranger_survival_survival <- function(
    model, x_analysis, time, event, x_assessment, prediction_horizon) {
  x_analysis <- sb_validate_predictor_frame(x_analysis, "x_analysis")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x_analysis))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x_analysis))
  if (!identical(names(x_analysis), model$features)) {
    sb_abort(
      "Analysis columns do not match the fixed ranger survival schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  if (any(time < 0)) {
    sb_abort(
      "`time` must contain non-negative values.",
      "stormbreaker_bad_adapter_time"
    )
  }
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    time
  )
  curves <- sb_predict_ranger_survival_curves(model, x_assessment)
  list(
    prediction_horizon = prediction_horizon,
    source = "ranger_survival_forest_direct",
    ties = "not_applicable",
    probability = sb_ranger_survival_probability_at_horizon(
      curves,
      prediction_horizon
    )
  )
}

sb_extract_ranger_survival_tree_model <- function(model) {
  if (!inherits(model, "stormbreaker_ranger_survival")) {
    sb_abort(
      "The supplied model is not a fitted fixed ranger survival learner.",
      "stormbreaker_bad_ranger_fit"
    )
  }
  list(
    num_trees = model$num_trees,
    effective_mtry = model$effective_mtry,
    features = model$features,
    params = model$params,
    importance_type = model$params$importance,
    model_serialization = "R_xdr_v3_ranger.forest",
    model_raw = serialize(
      model$forest,
      connection = NULL,
      ascii = FALSE,
      xdr = TRUE,
      version = 3L
    ),
    feature_importance = data.frame(
      feature = model$features,
      importance = model$feature_importance,
      stringsAsFactors = FALSE
    )
  )
}

sb_ranger_survival_learner_adapter <- function(
    num_trees = 500L,
    mtry = NULL,
    min_node_size = 3L,
    seed = 1L,
    learner = "ranger_survival") {
  params <- sb_ranger_survival_params(
    num_trees = num_trees,
    mtry = mtry,
    min_node_size = min_node_size,
    seed = seed,
    learner = learner
  )
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_ranger_survival(x, time, event, params, learner)
    },
    predict = sb_predict_ranger_survival,
    predict_survival = sb_predict_ranger_survival_survival,
    extract_tree_model = sb_extract_ranger_survival_tree_model,
    extract_state = function(model, context = list()) {
      sb_extract_ranger_survival_learner_state(model, context)
    },
    predict_state = sb_predict_ranger_survival_learner_state
  )
}
