sb_bad_learner_state <- function(message) {
  sb_abort(message, "stormbreaker_bad_learner_state")
}

sb_validate_learner_state_x <- function(state, x, learner) {
  valid_state <- inherits(state, "stormbreaker_learner_state") &&
    is.list(state) && identical(state[["learner"]], learner) &&
    is.character(state[["features"]]) &&
    length(state[["features"]]) > 0L &&
    !anyNA(state[["features"]]) &&
    !anyDuplicated(state[["features"]])
  if (!valid_state) {
    sb_bad_learner_state("The saved learner state is incomplete or invalid.")
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), state[["features"]])) {
    sb_abort(
      "Prediction columns do not match the saved learner state.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  x
}

sb_extract_cox_learner_state <- function(model, learner, coefficients,
                                         predict_survival, context,
                                         params = list(), audit = list()) {
  valid_coefficients <- is.data.frame(coefficients) &&
    all(c("feature", "coefficient") %in% names(coefficients)) &&
    nrow(coefficients) > 0L &&
    is.character(coefficients[["feature"]]) &&
    !anyNA(coefficients[["feature"]]) &&
    !anyDuplicated(coefficients[["feature"]]) &&
    is.numeric(coefficients[["coefficient"]]) &&
    !anyNA(coefficients[["coefficient"]]) &&
    all(is.finite(coefficients[["coefficient"]]))
  if (!valid_coefficients) {
    sb_bad_learner_state(
      "The fitted Cox learner did not yield portable coefficients."
    )
  }

  baseline_hazard <- NULL
  prediction_horizon <- if ("prediction_horizon" %in% names(context)) {
    context[["prediction_horizon"]]
  } else {
    NULL
  }
  if (!is.null(prediction_horizon)) {
    x <- context[["x"]]
    probability_result <- predict_survival(
      model,
      x,
      context[["time"]],
      context[["event"]],
      x[1L, , drop = FALSE],
      prediction_horizon
    )
    baseline_hazard <- list(
      prediction_horizon = probability_result$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context[["time_unit"]]
      } else {
        NA_character_
      },
      cumulative_baseline_hazard =
        probability_result$cumulative_baseline_hazard,
      baseline_survival = probability_result$baseline_survival,
      source = probability_result$source,
      ties = probability_result$ties
    )
  }

  structure(
    list(
      learner = learner,
      features = coefficients[["feature"]],
      coefficient = unname(coefficients[["coefficient"]]),
      risk_reference = "zero",
      risk_formula = "x_beta",
      params = params,
      audit = audit,
      baseline_hazard = baseline_hazard
    ),
    class = c("stormbreaker_cox_learner_state", "stormbreaker_learner_state")
  )
}

sb_predict_cox_learner_state <- function(state, x) {
  valid_learner <- is.list(state) && state[["learner"]] %in%
    sb_linear_survival_learner_ids()
  if (!isTRUE(valid_learner)) {
    sb_bad_learner_state("The saved Cox learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, state[["learner"]])
  coefficient <- state[["coefficient"]]
  valid_coefficient <- is.numeric(coefficient) &&
    length(coefficient) == ncol(x) && !anyNA(coefficient) &&
    all(is.finite(coefficient))
  if (!valid_coefficient) {
    sb_bad_learner_state("The saved Cox coefficients are invalid.")
  }
  risk <- unname(as.vector(as.matrix(x) %*% coefficient))
  if (anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "The saved Cox state produced invalid risk predictions.",
      "stormbreaker_invalid_survival_prediction"
    )
  }
  result <- list(risk = risk)
  baseline <- state[["baseline_hazard"]]
  if (!is.null(baseline)) {
    valid_baseline <- is.list(baseline) &&
      is.numeric(baseline[["prediction_horizon"]]) &&
      length(baseline[["prediction_horizon"]]) == 1L &&
      !is.na(baseline[["prediction_horizon"]]) &&
      is.finite(baseline[["prediction_horizon"]]) &&
      baseline[["prediction_horizon"]] > 0
    if (!valid_baseline) {
      sb_bad_learner_state("The saved Cox baseline-hazard state is invalid.")
    }
    result$prediction_horizon <- baseline[["prediction_horizon"]]
    result$time_unit <- baseline[["time_unit"]]
    result$survival_probability <- sb_survival_probability_from_hazard(
      risk,
      baseline[["cumulative_baseline_hazard"]]
    )
  }
  result
}

sb_extract_logistic_learner_state <- function(learner, coefficients,
                                              params = list()) {
  valid_coefficients <- is.data.frame(coefficients) &&
    all(c("term", "coefficient") %in% names(coefficients)) &&
    nrow(coefficients) >= 2L &&
    identical(coefficients[["term"]][1L], "(Intercept)") &&
    is.character(coefficients[["term"]]) &&
    !anyNA(coefficients[["term"]]) &&
    !anyDuplicated(coefficients[["term"]]) &&
    is.numeric(coefficients[["coefficient"]]) &&
    !anyNA(coefficients[["coefficient"]]) &&
    all(is.finite(coefficients[["coefficient"]]))
  if (!valid_coefficients) {
    sb_bad_learner_state(
      "The fitted logistic learner did not yield portable coefficients."
    )
  }
  structure(
    list(
      learner = learner,
      features = coefficients[["term"]][-1L],
      intercept = unname(coefficients[["coefficient"]][1L]),
      coefficient = unname(coefficients[["coefficient"]][-1L]),
      positive_class = 1L,
      params = params
    ),
    class = c(
      "stormbreaker_logistic_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_logistic_learner_state <- function(state, x) {
  valid_learner <- is.list(state) && state[["learner"]] %in%
    c(
      "logistic", "glmnet_logistic", sb_abess_binary_learner_ids(),
      sb_ncvreg_binary_learner_ids()
    )
  if (!isTRUE(valid_learner)) {
    sb_bad_learner_state("The saved logistic learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, state[["learner"]])
  intercept <- state[["intercept"]]
  coefficient <- state[["coefficient"]]
  valid_coefficients <- is.numeric(intercept) && length(intercept) == 1L &&
    !is.na(intercept) && is.finite(intercept) &&
    is.numeric(coefficient) && length(coefficient) == ncol(x) &&
    !anyNA(coefficient) && all(is.finite(coefficient))
  if (!valid_coefficients) {
    sb_bad_learner_state("The saved logistic coefficients are invalid.")
  }
  probability <- unname(stats::plogis(
    intercept + as.vector(as.matrix(x) %*% coefficient)
  ))
  valid_probability <- length(probability) == nrow(x) &&
    !anyNA(probability) && all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid_probability) {
    sb_abort(
      "The saved logistic state produced invalid class-one probabilities.",
      "stormbreaker_invalid_logistic_prediction"
    )
  }
  list(probability = probability)
}

sb_extract_naive_bayes_learner_state <- function(model) {
  if (!inherits(model, "stormbreaker_gaussian_naive_bayes")) {
    sb_bad_learner_state("The fitted Gaussian Naive Bayes state is invalid.")
  }
  structure(
    list(
      learner = "gaussian_naive_bayes",
      features = model$features,
      payload = unclass(model)
    ),
    class = c(
      "stormbreaker_naive_bayes_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_naive_bayes_learner_state <- function(state, x) {
  x <- sb_validate_learner_state_x(state, x, "gaussian_naive_bayes")
  payload <- state[["payload"]]
  required <- c(
    "features", "classes", "var_smoothing", "epsilon", "counts",
    "priors", "means", "variances"
  )
  if (!is.list(payload) || length(setdiff(required, names(payload))) > 0L) {
    sb_bad_learner_state("The saved Gaussian Naive Bayes state is invalid.")
  }
  model <- structure(payload, class = "stormbreaker_gaussian_naive_bayes")
  list(probability = sb_predict_gaussian_naive_bayes(model, x))
}

sb_extract_lightgbm_learner_state <- function(model) {
  extracted <- sb_extract_lightgbm_tree_model(model)
  structure(
    c(list(learner = model$learner), extracted),
    class = c(
      "stormbreaker_lightgbm_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_lightgbm_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_lightgbm_binary_learner_ids()) {
    sb_bad_learner_state("The saved LightGBM learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  valid_state <- is.character(state[["model_string"]]) &&
    length(state[["model_string"]]) == 1L &&
    !is.na(state[["model_string"]]) &&
    nzchar(state[["model_string"]]) &&
    is.numeric(state[["nrounds"]]) && length(state[["nrounds"]]) == 1L &&
    !is.na(state[["nrounds"]]) && state[["nrounds"]] >= 1
  if (!valid_state) {
    sb_bad_learner_state("The saved LightGBM learner state is invalid.")
  }
  if (!sb_lightgbm_available()) {
    sb_abort(
      "Package `lightgbm` is required for fixed LightGBM prediction.",
      "stormbreaker_missing_lightgbm"
    )
  }
  booster <- tryCatch(
    lightgbm::lgb.load(model_str = state[["model_string"]]),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved LightGBM state could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_lightgbm_state",
        parent = error
      )
    }
  )
  model <- structure(
    list(
      learner = learner,
      features = state[["features"]],
      params = state[["params"]],
      nrounds = as.integer(state[["nrounds"]]),
      booster = booster
    ),
    class = "stormbreaker_lightgbm_binary"
  )
  list(probability = sb_predict_lightgbm_binary(model, x))
}

sb_extract_catboost_learner_state <- function(model) {
  extracted <- sb_extract_catboost_tree_model(model)
  structure(
    c(list(learner = "catboost_binary"), extracted),
    class = c(
      "stormbreaker_catboost_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_catboost_learner_state <- function(state, x) {
  x <- sb_validate_learner_state_x(state, x, "catboost_binary")
  valid_state <- is.raw(state[["model_raw"]]) &&
    length(state[["model_raw"]]) > 0L &&
    is.numeric(state[["iterations"]]) &&
    length(state[["iterations"]]) == 1L &&
    !is.na(state[["iterations"]]) && state[["iterations"]] >= 1
  if (!valid_state) {
    sb_bad_learner_state("The saved CatBoost learner state is invalid.")
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for fixed CatBoost prediction.",
      "stormbreaker_missing_catboost"
    )
  }
  booster <- tryCatch(
    unserialize(state[["model_raw"]]),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved CatBoost state could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_catboost_state",
        parent = error
      )
    }
  )
  model <- structure(
    list(
      features = state[["features"]],
      params = state[["params"]],
      iterations = as.integer(state[["iterations"]]),
      booster = booster
    ),
    class = "stormbreaker_catboost_binary"
  )
  list(probability = sb_predict_catboost_binary(model, x))
}

sb_extract_ranger_binary_learner_state <- function(model) {
  learner <- model[["learner"]]
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_ranger_binary_learner_ids()) {
    sb_bad_learner_state("The fitted ranger binary learner state is invalid.")
  }
  extracted <- sb_extract_ranger_binary_tree_model(model)
  structure(
    c(list(learner = learner), extracted),
    class = c(
      "stormbreaker_ranger_binary_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_ranger_binary_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_ranger_binary_learner_ids()) {
    sb_bad_learner_state("The saved ranger binary learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  valid_state <- is.raw(state[["model_raw"]]) &&
    length(state[["model_raw"]]) > 0L &&
    is.numeric(state[["num_trees"]]) &&
    length(state[["num_trees"]]) == 1L &&
    !is.na(state[["num_trees"]]) && state[["num_trees"]] >= 1 &&
    is.numeric(state[["effective_mtry"]]) &&
    length(state[["effective_mtry"]]) == 1L &&
    !is.na(state[["effective_mtry"]]) && state[["effective_mtry"]] >= 1 &&
    is.list(state[["params"]])
  if (!valid_state) {
    sb_bad_learner_state("The saved ranger binary learner state is invalid.")
  }
  if (!sb_ranger_available()) {
    sb_abort(
      "Package `ranger` is required for fixed binary forest prediction.",
      "stormbreaker_missing_ranger"
    )
  }
  forest <- tryCatch(
    unserialize(state[["model_raw"]]),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved ranger binary state could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_ranger_state",
        parent = error
      )
    }
  )
  if (!inherits(forest, "ranger.forest") ||
    !identical(forest$levels, c("0", "1"))) {
    sb_bad_learner_state("The saved ranger binary forest is invalid.")
  }
  model <- structure(
    list(
      learner = learner,
      features = state[["features"]],
      params = state[["params"]],
      num_trees = as.integer(state[["num_trees"]]),
      effective_mtry = as.integer(state[["effective_mtry"]]),
      forest = forest
    ),
    class = "stormbreaker_ranger_binary"
  )
  list(probability = sb_predict_ranger_binary(model, x))
}

sb_extract_classical_binary_boosting_learner_state <- function(model) {
  extracted <- sb_extract_classical_binary_boosting_model(model)
  structure(
    c(list(learner = model$learner), extracted),
    class = c(
      "stormbreaker_classical_binary_boosting_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_classical_binary_boosting_learner_state <- function(state, x) {
  learner <- state[["learner"]]
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_classical_binary_boosting_learner_ids()) {
    sb_bad_learner_state("The saved binary boosting learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  valid_state <- is.raw(state[["model_raw"]]) &&
    length(state[["model_raw"]]) > 0L &&
    is.numeric(state[["iterations"]]) &&
    length(state[["iterations"]]) == 1L &&
    !is.na(state[["iterations"]]) && state[["iterations"]] >= 1L &&
    is.numeric(state[["fitted_iterations"]]) &&
    length(state[["fitted_iterations"]]) == 1L &&
    !is.na(state[["fitted_iterations"]]) &&
    state[["fitted_iterations"]] >= 1L &&
    state[["fitted_iterations"]] <= state[["iterations"]] &&
    is.list(state[["params"]])
  if (!valid_state) {
    sb_bad_learner_state("The saved binary boosting learner state is invalid.")
  }
  dependency <- sb_classical_binary_boosting_dependency(learner)
  if (!dependency$available) {
    sb_abort(
      paste0(
        "Package `", dependency$package,
        "` is required for fixed binary boosting prediction."
      ),
      if (identical(learner, "adaboost_binary")) {
        "stormbreaker_missing_fastadaboost"
      } else {
        "stormbreaker_missing_catools"
      }
    )
  }
  engine <- tryCatch(
    unserialize(state[["model_raw"]]),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved binary boosting state could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_classical_binary_boosting_state",
        parent = error
      )
    }
  )
  engine <- sb_validate_classical_binary_boosting_engine(
    engine,
    learner,
    as.integer(state[["iterations"]])
  )
  fitted_iterations <- if (identical(learner, "adaboost_binary")) {
    length(engine$trees)
  } else {
    nrow(engine$Stump)
  }
  if (!identical(
    as.integer(state[["fitted_iterations"]]),
    as.integer(fitted_iterations)
  )) {
    sb_bad_learner_state("The saved binary boosting learner state is invalid.")
  }
  model <- structure(
    list(
      learner = learner,
      features = state[["features"]],
      params = state[["params"]],
      iterations = as.integer(state[["iterations"]]),
      fitted_iterations = as.integer(state[["fitted_iterations"]]),
      engine = engine
    ),
    class = "stormbreaker_classical_binary_boosting_fixed"
  )
  list(probability = sb_predict_classical_binary_boosting(model, x))
}

sb_extract_catboost_cox_learner_state <- function(model, context = list()) {
  extracted <- sb_extract_catboost_tree_model(model)
  baseline_hazard <- NULL
  prediction_horizon <- if ("prediction_horizon" %in% names(context)) {
    context[["prediction_horizon"]]
  } else {
    NULL
  }
  if (!is.null(prediction_horizon)) {
    probability_result <- sb_predict_catboost_cox_survival(
      model,
      context[["x"]],
      context[["time"]],
      context[["event"]],
      context[["x"]][1L, , drop = FALSE],
      prediction_horizon
    )
    baseline_hazard <- list(
      prediction_horizon = probability_result$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context[["time_unit"]]
      } else {
        NA_character_
      },
      cumulative_baseline_hazard =
        probability_result$cumulative_baseline_hazard,
      baseline_survival = probability_result$baseline_survival,
      source = probability_result$source,
      ties = probability_result$ties
    )
  }
  structure(
    c(
      list(
        learner = "catboost_cox",
        risk_reference = "catboost_native",
        risk_formula = "catboost_raw_formula_value"
      ),
      extracted,
      list(baseline_hazard = baseline_hazard)
    ),
    class = c(
      "stormbreaker_catboost_cox_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_catboost_cox_learner_state <- function(state, x) {
  x <- sb_validate_learner_state_x(state, x, "catboost_cox")
  valid_state <- is.raw(state[["model_raw"]]) &&
    length(state[["model_raw"]]) > 0L &&
    is.numeric(state[["iterations"]]) &&
    length(state[["iterations"]]) == 1L &&
    !is.na(state[["iterations"]]) && state[["iterations"]] >= 1
  if (!valid_state) {
    sb_bad_learner_state("The saved CatBoost Cox learner state is invalid.")
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for fixed CatBoost Cox prediction.",
      "stormbreaker_missing_catboost"
    )
  }
  booster <- tryCatch(
    unserialize(state[["model_raw"]]),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved CatBoost Cox state could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_catboost_state",
        parent = error
      )
    }
  )
  model <- structure(
    list(
      features = state[["features"]],
      params = state[["params"]],
      iterations = as.integer(state[["iterations"]]),
      booster = booster
    ),
    class = "stormbreaker_catboost_cox"
  )
  result <- list(risk = sb_predict_catboost_cox(model, x))
  baseline <- state[["baseline_hazard"]]
  if (!is.null(baseline)) {
    valid_baseline <- is.list(baseline) &&
      is.numeric(baseline[["prediction_horizon"]]) &&
      length(baseline[["prediction_horizon"]]) == 1L &&
      !is.na(baseline[["prediction_horizon"]]) &&
      is.finite(baseline[["prediction_horizon"]]) &&
      baseline[["prediction_horizon"]] > 0 &&
      is.character(baseline[["time_unit"]]) &&
      length(baseline[["time_unit"]]) == 1L &&
      !is.na(baseline[["time_unit"]]) && nzchar(baseline[["time_unit"]])
    if (!valid_baseline) {
      sb_bad_learner_state(
        "The saved CatBoost Cox baseline-hazard state is invalid."
      )
    }
    result$prediction_horizon <- baseline[["prediction_horizon"]]
    result$time_unit <- baseline[["time_unit"]]
    result$survival_probability <- sb_survival_probability_from_hazard(
      result$risk,
      baseline[["cumulative_baseline_hazard"]]
    )
  }
  result
}

sb_extract_gbm_cox_learner_state <- function(model, context = list()) {
  extracted <- sb_extract_gbm_cox_tree_model(model)
  learner <- if ("learner" %in% names(model)) {
    model[["learner"]]
  } else {
    "gbm_cox"
  }
  baseline_hazard <- NULL
  prediction_horizon <- if ("prediction_horizon" %in% names(context)) {
    context[["prediction_horizon"]]
  } else {
    NULL
  }
  if (!is.null(prediction_horizon)) {
    probability_result <- sb_predict_gbm_cox_survival(
      model,
      context[["x"]],
      context[["time"]],
      context[["event"]],
      context[["x"]][1L, , drop = FALSE],
      prediction_horizon
    )
    baseline_hazard <- list(
      prediction_horizon = probability_result$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context[["time_unit"]]
      } else {
        NA_character_
      },
      cumulative_baseline_hazard =
        probability_result$cumulative_baseline_hazard,
      baseline_survival = probability_result$baseline_survival,
      source = probability_result$source,
      ties = probability_result$ties
    )
  }
  structure(
    c(
      list(
        learner = learner,
        risk_reference = "gbm_native",
        risk_formula = "gbm_link_prediction"
      ),
      extracted,
      list(baseline_hazard = baseline_hazard)
    ),
    class = c(
      "stormbreaker_gbm_cox_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_gbm_cox_learner_state <- function(state, x) {
  learner <- state[["learner"]]
  if (!is.character(learner) || length(learner) != 1L ||
    !learner %in% c("gbm_cox", "gbm_cox_cv")) {
    sb_bad_learner_state("The saved GBM-Cox learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  valid_state <- is.raw(state[["model_raw"]]) &&
    length(state[["model_raw"]]) > 0L &&
    is.numeric(state[["n_trees"]]) &&
    length(state[["n_trees"]]) == 1L &&
    !is.na(state[["n_trees"]]) && state[["n_trees"]] >= 1 &&
    is.list(state[["params"]])
  if (!valid_state) {
    sb_bad_learner_state("The saved GBM-Cox learner state is invalid.")
  }
  if (!sb_gbm_available()) {
    sb_abort(
      "Package `gbm` is required for fixed GBM-Cox prediction.",
      "stormbreaker_missing_gbm"
    )
  }
  booster <- tryCatch(
    unserialize(state[["model_raw"]]),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved GBM-Cox state could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_gbm_state",
        parent = error
      )
    }
  )
  if (!inherits(booster, "gbm")) {
    sb_bad_learner_state("The saved GBM-Cox model is invalid.")
  }
  model <- structure(
    list(
      learner = learner,
      features = state[["features"]],
      params = state[["params"]],
      n_trees = as.integer(state[["n_trees"]]),
      booster = booster
    ),
    class = "stormbreaker_gbm_cox"
  )
  result <- list(risk = sb_predict_gbm_cox(model, x))
  baseline <- state[["baseline_hazard"]]
  if (!is.null(baseline)) {
    valid_baseline <- is.list(baseline) &&
      is.numeric(baseline[["prediction_horizon"]]) &&
      length(baseline[["prediction_horizon"]]) == 1L &&
      !is.na(baseline[["prediction_horizon"]]) &&
      is.finite(baseline[["prediction_horizon"]]) &&
      baseline[["prediction_horizon"]] > 0 &&
      is.character(baseline[["time_unit"]]) &&
      length(baseline[["time_unit"]]) == 1L &&
      !is.na(baseline[["time_unit"]]) && nzchar(baseline[["time_unit"]])
    if (!valid_baseline) {
      sb_bad_learner_state(
        "The saved GBM-Cox baseline-hazard state is invalid."
      )
    }
    result$prediction_horizon <- baseline[["prediction_horizon"]]
    result$time_unit <- baseline[["time_unit"]]
    result$survival_probability <- sb_survival_probability_from_hazard(
      result$risk,
      baseline[["cumulative_baseline_hazard"]]
    )
  }
  result
}

sb_extract_ranger_survival_learner_state <- function(model,
                                                     context = list()) {
  learner <- model[["learner"]]
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_ranger_survival_learner_ids()) {
    sb_bad_learner_state("The fitted ranger survival learner state is invalid.")
  }
  extracted <- sb_extract_ranger_survival_tree_model(model)
  survival_probability <- NULL
  prediction_horizon <- if ("prediction_horizon" %in% names(context)) {
    context[["prediction_horizon"]]
  } else {
    NULL
  }
  if (!is.null(prediction_horizon)) {
    probability_result <- sb_predict_ranger_survival_survival(
      model,
      context[["x"]],
      context[["time"]],
      context[["event"]],
      context[["x"]][1L, , drop = FALSE],
      prediction_horizon
    )
    survival_probability <- list(
      prediction_horizon = probability_result$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context[["time_unit"]]
      } else {
        NA_character_
      },
      source = probability_result$source,
      ties = probability_result$ties
    )
  }
  structure(
    c(
      list(
        learner = learner,
        risk_reference = "ranger_native",
        risk_formula = "row_sums_cumulative_hazard"
      ),
      extracted,
      list(
        baseline_hazard = NULL,
        survival_probability = survival_probability
      )
    ),
    class = c(
      "stormbreaker_ranger_survival_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_ranger_survival_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_ranger_survival_learner_ids()) {
    sb_bad_learner_state("The saved ranger survival learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  valid_state <- is.raw(state[["model_raw"]]) &&
    length(state[["model_raw"]]) > 0L &&
    is.numeric(state[["num_trees"]]) &&
    length(state[["num_trees"]]) == 1L &&
    !is.na(state[["num_trees"]]) && state[["num_trees"]] >= 1 &&
    is.numeric(state[["effective_mtry"]]) &&
    length(state[["effective_mtry"]]) == 1L &&
    !is.na(state[["effective_mtry"]]) && state[["effective_mtry"]] >= 1 &&
    is.list(state[["params"]])
  if (!valid_state) {
    sb_bad_learner_state("The saved ranger survival learner state is invalid.")
  }
  if (!sb_ranger_available()) {
    sb_abort(
      "Package `ranger` is required for fixed survival forest prediction.",
      "stormbreaker_missing_ranger"
    )
  }
  forest <- tryCatch(
    unserialize(state[["model_raw"]]),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved ranger survival state could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_ranger_state",
        parent = error
      )
    }
  )
  if (!inherits(forest, "ranger.forest")) {
    sb_bad_learner_state("The saved ranger survival forest is invalid.")
  }
  model <- structure(
    list(
      learner = learner,
      features = state[["features"]],
      params = state[["params"]],
      num_trees = as.integer(state[["num_trees"]]),
      effective_mtry = as.integer(state[["effective_mtry"]]),
      forest = forest
    ),
    class = "stormbreaker_ranger_survival"
  )
  curves <- sb_predict_ranger_survival_curves(model, x)
  result <- list(risk = curves$risk)
  probability_state <- state[["survival_probability"]]
  if (!is.null(probability_state)) {
    valid_probability_state <- is.list(probability_state) &&
      is.numeric(probability_state[["prediction_horizon"]]) &&
      length(probability_state[["prediction_horizon"]]) == 1L &&
      !is.na(probability_state[["prediction_horizon"]]) &&
      is.finite(probability_state[["prediction_horizon"]]) &&
      probability_state[["prediction_horizon"]] > 0 &&
      is.character(probability_state[["time_unit"]]) &&
      length(probability_state[["time_unit"]]) == 1L &&
      !is.na(probability_state[["time_unit"]]) &&
      nzchar(probability_state[["time_unit"]]) &&
      identical(
        probability_state[["source"]],
        "ranger_survival_forest_direct"
      ) &&
      identical(probability_state[["ties"]], "not_applicable")
    if (!valid_probability_state) {
      sb_bad_learner_state(
        "The saved ranger survival-probability state is invalid."
      )
    }
    result$prediction_horizon <-
      probability_state[["prediction_horizon"]]
    result$time_unit <- probability_state[["time_unit"]]
    result$survival_probability <-
      sb_ranger_survival_probability_at_horizon(
        curves,
        probability_state[["prediction_horizon"]]
      )
  }
  result
}
