sb_catboost_available <- function() {
  requireNamespace("catboost", quietly = TRUE)
}

sb_validate_catboost_learning_rate <- function(learning_rate) {
  valid <- is.numeric(learning_rate) && length(learning_rate) == 1L &&
    !is.na(learning_rate) && is.finite(learning_rate) && learning_rate > 0
  if (!valid) {
    sb_abort(
      "`catboost_learning_rate` must be one finite number greater than 0.",
      "stormbreaker_bad_catboost_learning_rate"
    )
  }
  as.numeric(learning_rate)
}

sb_validate_catboost_depth <- function(depth) {
  depth <- sb_check_positive_integer(depth, "catboost_depth", minimum = 1L)
  if (depth > 16L) {
    sb_abort(
      "`catboost_depth` must be an integer between 1 and 16 for CPU training.",
      "stormbreaker_bad_catboost_depth"
    )
  }
  depth
}

sb_validate_catboost_l2_leaf_reg <- function(l2_leaf_reg) {
  valid <- is.numeric(l2_leaf_reg) && length(l2_leaf_reg) == 1L &&
    !is.na(l2_leaf_reg) && is.finite(l2_leaf_reg) && l2_leaf_reg > 0
  if (!valid) {
    sb_abort(
      "`catboost_l2_leaf_reg` must be one finite number greater than 0.",
      "stormbreaker_bad_catboost_l2_leaf_reg"
    )
  }
  as.numeric(l2_leaf_reg)
}

sb_catboost_binary_params <- function(
    iterations = 100L,
    learning_rate = 0.1,
    depth = 6L,
    l2_leaf_reg = 3,
    seed = 1L) {
  list(
    iterations = sb_check_positive_integer(
      iterations,
      "catboost_iterations",
      minimum = 1L
    ),
    learning_rate = sb_validate_catboost_learning_rate(learning_rate),
    depth = sb_validate_catboost_depth(depth),
    l2_leaf_reg = sb_validate_catboost_l2_leaf_reg(l2_leaf_reg),
    random_seed = sb_check_positive_integer(
      seed,
      "catboost_seed",
      minimum = 1L
    ),
    loss_function = "Logloss",
    eval_metric = "Logloss",
    task_type = "CPU",
    thread_count = 1L,
    bootstrap_type = "Bayesian",
    bagging_temperature = 1,
    random_strength = 1,
    rsm = 1,
    grow_policy = "SymmetricTree",
    boosting_type = "Plain",
    use_best_model = FALSE,
    boost_from_average = TRUE,
    allow_writing_files = FALSE,
    save_snapshot = FALSE,
    logging_level = "Silent"
  )
}

sb_fit_catboost_binary <- function(x, outcome, params) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "A fixed CatBoost analysis fold must contain both classes.",
      "stormbreaker_catboost_single_class"
    )
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for the fixed CatBoost learner.",
      "stormbreaker_missing_catboost"
    )
  }

  training_pool <- catboost::catboost.load_pool(
    data = as.matrix(x),
    label = outcome
  )
  booster <- catboost::catboost.train(
    learn_pool = training_pool,
    test_pool = NULL,
    params = params
  )
  tree_count <- as.integer(booster$tree_count)
  if (!identical(tree_count, params$iterations)) {
    sb_abort(
      "Fixed CatBoost did not fit the declared number of iterations.",
      "stormbreaker_incomplete_catboost_fit"
    )
  }

  structure(
    list(
      features = names(x),
      params = params,
      iterations = params$iterations,
      booster = booster
    ),
    class = "stormbreaker_catboost_binary"
  )
}

sb_predict_catboost_binary <- function(model, x) {
  if (!inherits(model, "stormbreaker_catboost_binary")) {
    sb_abort(
      "The supplied model is not a fitted fixed CatBoost learner.",
      "stormbreaker_bad_catboost_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fixed CatBoost schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for fixed CatBoost prediction.",
      "stormbreaker_missing_catboost"
    )
  }

  prediction_pool <- catboost::catboost.load_pool(as.matrix(x))
  probability <- as.numeric(catboost::catboost.predict(
    model$booster,
    prediction_pool,
    prediction_type = "Probability",
    ntree_start = 0L,
    ntree_end = model$iterations,
    thread_count = 1L
  ))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      "Fixed CatBoost produced invalid class-one probabilities.",
      "stormbreaker_invalid_catboost_prediction"
    )
  }
  probability
}

sb_extract_catboost_tree_model <- function(model) {
  if (!inherits(model, c(
    "stormbreaker_catboost_binary",
    "stormbreaker_catboost_cox",
    "stormbreaker_catboost_aft_survival"
  ))) {
    sb_abort(
      "The supplied model is not a fitted fixed CatBoost learner.",
      "stormbreaker_bad_catboost_fit"
    )
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required to audit a CatBoost model.",
      "stormbreaker_missing_catboost"
    )
  }

  importance <- as.numeric(catboost::catboost.get_feature_importance(
    model$booster,
    type = "FeatureImportance",
    thread_count = 1L
  ))
  valid_importance <- length(importance) == length(model$features) &&
    !anyNA(importance) && all(is.finite(importance))
  if (!valid_importance) {
    sb_abort(
      "Fixed CatBoost did not return one finite importance per feature.",
      "stormbreaker_invalid_catboost_importance"
    )
  }

  list(
    iterations = model$iterations,
    tree_count = as.integer(model$booster$tree_count),
    features = model$features,
    params = model$params,
    model_serialization = "R_xdr_v3_catboost.Model",
    model_raw = serialize(
      model$booster,
      connection = NULL,
      ascii = FALSE,
      xdr = TRUE,
      version = 3L
    ),
    feature_importance = data.frame(
      feature = model$features,
      importance = importance,
      stringsAsFactors = FALSE
    )
  )
}

sb_catboost_binary_learner_adapter <- function(
    iterations = 100L,
    learning_rate = 0.1,
    depth = 6L,
    l2_leaf_reg = 3,
    seed = 1L) {
  params <- sb_catboost_binary_params(
    iterations = iterations,
    learning_rate = learning_rate,
    depth = depth,
    l2_leaf_reg = l2_leaf_reg,
    seed = seed
  )
  list(
    id = "catboost_binary",
    params = params,
    fit = function(x, outcome) {
      sb_fit_catboost_binary(x, outcome, params)
    },
    predict = sb_predict_catboost_binary,
    extract_tree_model = sb_extract_catboost_tree_model,
    extract_state = function(model, context = list()) {
      sb_extract_catboost_learner_state(model)
    },
    predict_state = sb_predict_catboost_learner_state
  )
}
