sb_classical_binary_boosting_learner_ids <- function() {
  c("adaboost_binary", "logitboost_binary")
}

sb_fastadaboost_available <- function() {
  requireNamespace("fastAdaboost", quietly = TRUE)
}

sb_catools_available <- function() {
  requireNamespace("caTools", quietly = TRUE)
}

sb_classical_binary_boosting_dependency <- function(learner) {
  switch(learner,
    adaboost_binary = list(
      field = "fastadaboost_version",
      package = "fastAdaboost",
      available = sb_fastadaboost_available()
    ),
    logitboost_binary = list(
      field = "catools_version",
      package = "caTools",
      available = sb_catools_available()
    ),
    NULL
  )
}

sb_classical_binary_boosting_params <- function(learner, iterations = NULL) {
  if (!learner %in% sb_classical_binary_boosting_learner_ids()) {
    sb_abort(
      "The requested classical binary boosting learner is not supported.",
      "stormbreaker_bad_classical_binary_boosting_learner"
    )
  }
  if (is.null(iterations)) {
    iterations <- if (identical(learner, "adaboost_binary")) 50L else 51L
  }
  iterations <- sb_check_positive_integer(
    iterations,
    if (identical(learner, "adaboost_binary")) {
      "adaboost_iterations"
    } else {
      "logitboost_iterations"
    },
    minimum = 1L
  )
  if (identical(learner, "adaboost_binary")) {
    return(list(
      iterations = iterations,
      method = "AdaBoost.M1",
      base_learner = "rpart_classification_tree",
      rpart_cp = 0,
      iteration_semantics = "maximum_until_native_convergence",
      tuning = "none"
    ))
  }
  list(
    iterations = iterations,
    method = "LogitBoost",
    base_learner = "decision_stump",
    iteration_semantics = "exact",
    tuning = "none"
  )
}



sb_fastadaboost_fit <- function(formula, data, iterations) {
  fastAdaboost::adaboost(formula, data, nIter = iterations)
}

sb_catools_logitboost_fit <- function(x, outcome, iterations) {
  caTools::LogitBoost(x, outcome, nIter = iterations)
}

sb_validate_classical_binary_boosting_engine <- function(engine, learner,
                                                         iterations) {
  valid <- if (identical(learner, "adaboost_binary")) {
    inherits(engine, "adaboost") &&
      length(engine$trees) >= 1L && length(engine$trees) <= iterations &&
      length(engine$weights) == length(engine$trees) &&
      all(vapply(engine$trees, inherits, logical(1), "rpart")) &&
      all(is.finite(engine$weights)) &&
      identical(unname(engine$classnames), c("0", "1"))
  } else {
    inherits(engine, "LogitBoost") && is.matrix(engine$Stump) &&
      nrow(engine$Stump) == iterations && ncol(engine$Stump) == 3L &&
      !anyNA(engine$Stump) && all(is.finite(engine$Stump)) &&
      identical(as.character(engine$lablist), c("0", "1"))
  }
  if (!valid) {
    sb_abort(
      "The classical binary boosting fit is incomplete or invalid.",
      "stormbreaker_invalid_classical_binary_boosting_fit"
    )
  }
  engine
}

sb_fit_classical_binary_boosting <- function(x, outcome, learner,
                                                   params) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "A classical binary boosting analysis fold must contain both classes.",
      "stormbreaker_classical_binary_boosting_single_class"
    )
  }
  params <- sb_classical_binary_boosting_params(learner, params$iterations)

  if (identical(learner, "adaboost_binary")) {
    if (!sb_fastadaboost_available()) {
      sb_abort(
        "Package `fastAdaboost` is required for fixed AdaBoost.M1.",
        "stormbreaker_missing_fastadaboost"
      )
    }
    response <- ".stormbreaker_outcome"
    while (response %in% names(x)) {
      response <- paste0(".", response)
    }
    training <- data.frame(
      factor(outcome, levels = c(0L, 1L)),
      x,
      check.names = FALSE
    )
    names(training)[[1L]] <- response
    formula <- stats::as.formula(
      paste0("`", response, "` ~ ."),
      env = baseenv()
    )
    engine <- sb_fastadaboost_fit(formula, training, params$iterations)
  } else {
    if (!sb_catools_available()) {
      sb_abort(
        "Package `caTools` is required for fixed LogitBoost.",
        "stormbreaker_missing_catools"
      )
    }
    engine <- sb_catools_logitboost_fit(
      as.matrix(x),
      outcome,
      params$iterations
    )
  }
  engine <- sb_validate_classical_binary_boosting_engine(
    engine,
    learner,
    params$iterations
  )

  structure(
    list(
      learner = learner,
      features = names(x),
      params = params,
      iterations = params$iterations,
      fitted_iterations = if (identical(
        learner,
        "adaboost_binary"
      )) {
        as.integer(length(engine$trees))
      } else {
        as.integer(nrow(engine$Stump))
      },
      engine = engine
    ),
    class = "stormbreaker_classical_binary_boosting_fixed"
  )
}

sb_predict_classical_binary_boosting <- function(model, x) {
  if (!inherits(model, "stormbreaker_classical_binary_boosting_fixed") ||
    !model$learner %in% sb_classical_binary_boosting_learner_ids()) {
    sb_abort(
      "The supplied model is not a fitted classical binary boosting learner.",
      "stormbreaker_bad_classical_binary_boosting_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the binary boosting schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }

  if (identical(model$learner, "adaboost_binary")) {
    if (!sb_fastadaboost_available()) {
      sb_abort(
        "Package `fastAdaboost` is required for fixed AdaBoost.M1 prediction.",
        "stormbreaker_missing_fastadaboost"
      )
    }
    prediction <- stats::predict(model$engine, newdata = x)
    positive_column <- match("1", unname(model$engine$classnames))
    probability <- prediction$prob[, positive_column]
  } else {
    if (!sb_catools_available()) {
      sb_abort(
        "Package `caTools` is required for fixed LogitBoost prediction.",
        "stormbreaker_missing_catools"
      )
    }
    prediction <- stats::predict(
      model$engine,
      as.matrix(x),
      type = "raw",
      nIter = model$iterations
    )
    positive_column <- match("1", as.character(model$engine$lablist))
    probability <- prediction[, positive_column]
  }
  probability <- unname(as.numeric(probability))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      "Classical binary boosting produced invalid class-one probabilities.",
      "stormbreaker_invalid_classical_binary_boosting_prediction"
    )
  }
  probability
}

sb_classical_binary_boosting_feature_usage <- function(model) {
  split_features <- if (identical(model$learner, "adaboost_binary")) {
    unlist(lapply(model$engine$trees, function(tree) {
      variables <- as.character(tree$frame$var)
      variables[variables != "<leaf>"]
    }), use.names = FALSE)
  } else {
    indices <- as.integer(model$engine$Stump[, "feature"])
    if (any(indices < 1L | indices > length(model$features))) {
      sb_abort(
        "LogitBoost returned a stump outside the fitted feature schema.",
        "stormbreaker_invalid_classical_binary_boosting_fit"
      )
    }
    model$features[indices]
  }
  if (length(setdiff(split_features, model$features)) > 0L) {
    sb_abort(
      "Binary boosting returned a split outside the fitted feature schema.",
      "stormbreaker_invalid_classical_binary_boosting_fit"
    )
  }
  as.numeric(tabulate(
    match(split_features, model$features),
    nbins = length(model$features)
  ))
}

sb_extract_classical_binary_boosting_model <- function(model) {
  if (!inherits(model, "stormbreaker_classical_binary_boosting_fixed")) {
    sb_abort(
      "The supplied model is not a fitted classical binary boosting learner.",
      "stormbreaker_bad_classical_binary_boosting_fit"
    )
  }
  list(
    iterations = model$iterations,
    fitted_iterations = model$fitted_iterations,
    features = model$features,
    params = model$params,
    importance_type = "split_count",
    model_serialization = paste0(
      "R_xdr_v3_",
      if (identical(model$learner, "adaboost_binary")) {
        "adaboost"
      } else {
        "LogitBoost"
      }
    ),
    model_raw = serialize(
      model$engine,
      connection = NULL,
      ascii = FALSE,
      xdr = TRUE,
      version = 3L
    ),
    feature_importance = data.frame(
      feature = model$features,
      importance = sb_classical_binary_boosting_feature_usage(model),
      stringsAsFactors = FALSE
    )
  )
}

sb_classical_binary_boosting_learner_adapter <- function(learner,
                                                         iterations = NULL) {
  params <- sb_classical_binary_boosting_params(learner, iterations)
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_classical_binary_boosting(x, outcome, learner, params)
    },
    predict = sb_predict_classical_binary_boosting,
    extract_tree_model = sb_extract_classical_binary_boosting_model,
    extract_state = function(model, context = list()) {
      sb_extract_classical_binary_boosting_learner_state(model)
    },
    predict_state = sb_predict_classical_binary_boosting_learner_state
  )
}
