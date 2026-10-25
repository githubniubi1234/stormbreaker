sb_lightgbm_binary_learner_ids <- function() {
  c(
    "lightgbm_binary",
    "lightgbm_binary_dart",
    "lightgbm_binary_goss",
    "lightgbm_binary_rf"
  )
}

sb_lightgbm_binary_variant <- function(learner) {
  variant <- switch(learner,
    lightgbm_binary = "gbdt",
    lightgbm_binary_dart = "dart",
    lightgbm_binary_goss = "goss",
    lightgbm_binary_rf = "rf",
    NULL
  )
  if (is.null(variant)) {
    sb_abort(
      paste0("Unknown fixed LightGBM learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  variant
}

sb_lightgbm_available <- function() {
  requireNamespace("lightgbm", quietly = TRUE)
}

sb_validate_lightgbm_learning_rate <- function(learning_rate) {
  valid <- is.numeric(learning_rate) && length(learning_rate) == 1L &&
    !is.na(learning_rate) && is.finite(learning_rate) && learning_rate > 0
  if (!valid) {
    sb_abort(
      "`lightgbm_learning_rate` must be one finite number greater than 0.",
      "stormbreaker_bad_lightgbm_learning_rate"
    )
  }
  as.numeric(learning_rate)
}

sb_lightgbm_binary_params <- function(
    nrounds = 100L,
    learning_rate = 0.1,
    num_leaves = 31L,
    min_data_in_leaf = 20L,
    seed = 1L,
    learner = "lightgbm_binary") {
  variant <- sb_lightgbm_binary_variant(learner)
  params <- list(
    nrounds = sb_check_positive_integer(
      nrounds,
      "lightgbm_nrounds",
      minimum = 1L
    ),
    learning_rate = sb_validate_lightgbm_learning_rate(learning_rate),
    num_leaves = sb_check_positive_integer(
      num_leaves,
      "lightgbm_num_leaves",
      minimum = 2L
    ),
    min_data_in_leaf = sb_check_positive_integer(
      min_data_in_leaf,
      "lightgbm_min_data_in_leaf",
      minimum = 0L
    ),
    seed = sb_check_positive_integer(
      seed,
      "lightgbm_seed",
      minimum = 1L
    ),
    objective = "binary",
    metric = "binary_logloss",
    boosting = if (identical(variant, "dart")) {
      "dart"
    } else if (identical(variant, "rf")) {
      "rf"
    } else {
      "gbdt"
    },
    device_type = "cpu",
    num_threads = 1L,
    deterministic = TRUE,
    force_col_wise = TRUE,
    feature_fraction = if (identical(variant, "rf")) 0.8 else 1,
    bagging_fraction = if (identical(variant, "rf")) 0.8 else 1,
    bagging_freq = if (identical(variant, "rf")) 1L else 0L,
    max_depth = -1L,
    lambda_l1 = 0,
    lambda_l2 = 0,
    boost_from_average = TRUE,
    verbosity = -1L
  )
  if (identical(variant, "dart")) {
    return(c(params, list(
      drop_rate = 0.1,
      max_drop = 50L,
      skip_drop = 0,
      uniform_drop = FALSE,
      xgboost_dart_mode = FALSE
    )))
  }
  if (identical(variant, "goss")) {
    return(c(params, list(
      data_sample_strategy = "goss",
      top_rate = 0.2,
      other_rate = 0.1
    )))
  }
  params
}

sb_fit_lightgbm_binary <- function(x, outcome, params, learner) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "A fixed LightGBM analysis fold must contain both classes.",
      "stormbreaker_lightgbm_single_class"
    )
  }
  if (!sb_lightgbm_available()) {
    sb_abort(
      "Package `lightgbm` is required for the fixed LightGBM learner.",
      "stormbreaker_missing_lightgbm"
    )
  }

  matrix_x <- as.matrix(x)
  training_params <- params[setdiff(names(params), "nrounds")]
  # Native alias avoids BiocGenerics format.list truncating binary_logloss.
  training_params$metric <- "binary"
  dataset <- lightgbm::lgb.Dataset(
    data = matrix_x,
    label = outcome,
    colnames = names(x),
    free_raw_data = TRUE
  )
  booster <- lightgbm::lgb.train(
    params = training_params,
    data = dataset,
    nrounds = params$nrounds,
    verbose = -1L,
    record = FALSE,
    early_stopping_rounds = NULL,
    serializable = TRUE
  )

  structure(
    list(
      learner = learner,
      features = names(x),
      params = params,
      nrounds = params$nrounds,
      booster = booster
    ),
    class = "stormbreaker_lightgbm_binary"
  )
}

sb_predict_lightgbm_binary <- function(model, x) {
  if (!inherits(model, "stormbreaker_lightgbm_binary")) {
    sb_abort(
      "The supplied model is not a fitted fixed LightGBM learner.",
      "stormbreaker_bad_lightgbm_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fixed LightGBM schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  if (!sb_lightgbm_available()) {
    sb_abort(
      "Package `lightgbm` is required for fixed LightGBM prediction.",
      "stormbreaker_missing_lightgbm"
    )
  }

  probability <- as.numeric(stats::predict(
    model$booster,
    as.matrix(x),
    type = "response",
    num_iteration = model$nrounds,
    params = list(num_threads = 1L)
  ))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      "Fixed LightGBM produced invalid class-one probabilities.",
      "stormbreaker_invalid_lightgbm_prediction"
    )
  }
  probability
}

sb_extract_lightgbm_tree_model <- function(model) {
  if (!inherits(model, "stormbreaker_lightgbm_binary")) {
    sb_abort(
      "The supplied model is not a fitted fixed LightGBM learner.",
      "stormbreaker_bad_lightgbm_fit"
    )
  }
  if (!sb_lightgbm_available()) {
    sb_abort(
      "Package `lightgbm` is required to audit a LightGBM model.",
      "stormbreaker_missing_lightgbm"
    )
  }

  model_string <- model$booster$save_model_to_string(
    num_iteration = model$nrounds
  )
  feature_line <- grep(
    "^feature_names=",
    strsplit(model_string, "\n", fixed = TRUE)[[1L]],
    value = TRUE
  )
  native_features <- strsplit(
    sub("^feature_names=", "", feature_line[[1L]]),
    " ",
    fixed = TRUE
  )[[1L]]

  observed <- if (!is.null(model$reused_native_importance) &&
    identical(model_string, model$reused_native_importance$model_string)) {
    model$reused_native_importance$observed
  } else as.data.frame(
    lightgbm::lgb.importance(model$booster, percentage = FALSE),
    stringsAsFactors = FALSE
  )
  importance <- data.frame(
    feature = model$features,
    gain = numeric(length(model$features)),
    cover = numeric(length(model$features)),
    frequency = integer(length(model$features)),
    stringsAsFactors = FALSE
  )
  if (nrow(observed) > 0L) {
    index <- match(observed$Feature, native_features)
    if (anyNA(index)) {
      sb_abort(
        "LightGBM importance names could not be mapped to the fitted features.",
        "stormbreaker_lightgbm_importance_mismatch"
      )
    }
    importance$gain[index] <- as.numeric(observed$Gain)
    importance$cover[index] <- as.numeric(observed$Cover)
    importance$frequency[index] <- as.integer(observed$Frequency)
  }

  list(
    nrounds = model$nrounds,
    features = model$features,
    params = model$params,
    model_string = model_string,
    feature_importance = importance
  )
}

sb_lightgbm_binary_learner_adapter <- function(
    nrounds = 100L,
    learning_rate = 0.1,
    num_leaves = 31L,
    min_data_in_leaf = 20L,
    seed = 1L,
    learner = "lightgbm_binary") {
  params <- sb_lightgbm_binary_params(
    nrounds = nrounds,
    learning_rate = learning_rate,
    num_leaves = num_leaves,
    min_data_in_leaf = min_data_in_leaf,
    seed = seed,
    learner = learner
  )
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_lightgbm_binary(x, outcome, params, learner)
    },
    predict = sb_predict_lightgbm_binary,
    extract_tree_model = sb_extract_lightgbm_tree_model,
    extract_state = function(model, context = list()) {
      sb_extract_lightgbm_learner_state(model)
    },
    predict_state = sb_predict_lightgbm_learner_state
  )
}
