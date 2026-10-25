sb_ranger_binary_params <- function(
    num_trees = 500L,
    mtry = NULL,
    min_node_size = 10L,
    seed = 1L,
    learner = "ranger_binary") {
  if (!learner %in% sb_ranger_binary_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed ranger binary learner `", learner, "`."),
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
    probability = TRUE,
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

sb_fit_ranger_binary <- function(
  x, outcome, params, learner = "ranger_binary") {
  if (!learner %in% sb_ranger_binary_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed ranger binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "A fixed ranger binary analysis fold must contain both classes.",
      "stormbreaker_ranger_binary_single_class"
    )
  }
  if (!sb_ranger_available()) {
    sb_abort(
      "Package `ranger` is required for the fixed binary forest learner.",
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

  fit_arguments <- list(
    x = as.matrix(x),
    y = factor(as.character(outcome), levels = c("0", "1")),
    num.trees = params$num_trees,
    mtry = effective_mtry,
    probability = params$probability,
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
          "The fixed ranger binary learner raised a warning: ",
          message
        ),
        "stormbreaker_ranger_binary_warning",
        parent = warning
      )
    }
  )
  valid_forest <- inherits(fitted$forest, "ranger.forest") &&
    identical(as.integer(fitted$num.trees), params$num_trees) &&
    identical(fitted$treetype, "Probability estimation") &&
    identical(fitted$forest$levels, c("0", "1"))
  if (!valid_forest) {
    sb_abort(
      "Fixed ranger binary did not return the declared probability forest.",
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
      "Fixed ranger binary did not return one finite importance per feature.",
      "stormbreaker_invalid_ranger_importance"
    )
  }

  structure(
    list(
      learner = learner,
      features = names(x),
      params = params,
      num_trees = params$num_trees,
      effective_mtry = as.integer(effective_mtry),
      forest = fitted$forest,
      feature_importance = unname(importance[names(x)])
    ),
    class = "stormbreaker_ranger_binary"
  )
}

sb_predict_ranger_binary <- function(model, x) {
  if (!inherits(model, "stormbreaker_ranger_binary")) {
    sb_abort(
      "The supplied model is not a fitted fixed ranger binary learner.",
      "stormbreaker_bad_ranger_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fixed ranger binary schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  if (!sb_ranger_available()) {
    sb_abort(
      "Package `ranger` is required for fixed binary forest prediction.",
      "stormbreaker_missing_ranger"
    )
  }

  predictions <- stats::predict(
    model$forest,
    data = x,
    num.trees = model$num_trees,
    type = "response",
    seed = 0L,
    num.threads = 1L,
    verbose = FALSE
  )$predictions
  valid_matrix <- is.matrix(predictions) && nrow(predictions) == nrow(x) &&
    identical(colnames(predictions), c("0", "1")) &&
    !anyNA(predictions) && all(is.finite(predictions))
  if (!valid_matrix) {
    sb_abort(
      "Fixed ranger binary produced invalid class probability predictions.",
      "stormbreaker_invalid_ranger_prediction"
    )
  }
  probability <- unname(as.numeric(predictions[, "1"]))
  if (any(probability < 0 | probability > 1)) {
    sb_abort(
      "Fixed ranger binary produced invalid class-one probabilities.",
      "stormbreaker_invalid_ranger_prediction"
    )
  }
  probability
}

sb_extract_ranger_binary_tree_model <- function(model) {
  if (!inherits(model, "stormbreaker_ranger_binary")) {
    sb_abort(
      "The supplied model is not a fitted fixed ranger binary learner.",
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

sb_ranger_binary_learner_adapter <- function(
    num_trees = 500L,
    mtry = NULL,
    min_node_size = 10L,
    seed = 1L,
    learner = "ranger_binary") {
  params <- sb_ranger_binary_params(
    num_trees = num_trees,
    mtry = mtry,
    min_node_size = min_node_size,
    seed = seed,
    learner = learner
  )
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_ranger_binary(x, outcome, params, learner)
    },
    predict = sb_predict_ranger_binary,
    extract_tree_model = sb_extract_ranger_binary_tree_model,
    extract_state = function(model, context = list()) {
      sb_extract_ranger_binary_learner_state(model)
    },
    predict_state = sb_predict_ranger_binary_learner_state
  )
}
