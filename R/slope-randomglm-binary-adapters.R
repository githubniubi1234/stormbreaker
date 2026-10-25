
sb_slope_binary_learner_ids <- function() character()

sb_random_glm_binary_learner_ids <- function() {
  "random_glm_binary"
}

sb_slope_randomglm_binary_learner_ids <- function() {
  c(sb_slope_binary_learner_ids(), sb_random_glm_binary_learner_ids())
}

sb_slope_randomglm_binary_params <- function(learner) {
  switch(learner,
    slope_binary = list(
      q = c(0.05, 0.1, 0.2),
      gamma = 0,
      inner_folds = 5L,
      repeats = 1L,
      measure = "deviance",
      center = "mean",
      scale = "l2",
      path_length = 30L,
      threads = 1L,
      seed = 1L
    ),
    random_glm_binary = list(
      n_bags = 100L,
      n_candidate_covariates = 3L,
      max_interaction_order = 1L,
      replace = TRUE,
      threads = 1L,
      seed = 1L
    ),
    sb_abort(
      paste0("Unknown SLOPE/randomGLM binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_slope_randomglm_binary_dependency <- function(learner) {
  package <- switch(learner,
    slope_binary = "SLOPE",
    random_glm_binary = "randomGLM",
    sb_abort(
      paste0("Unknown SLOPE/randomGLM binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
  list(
    package = package,
    available = requireNamespace(package, quietly = TRUE)
  )
}

sb_require_slope_randomglm_binary_dependency <- function(learner) {
  dependency <- sb_slope_randomglm_binary_dependency(learner)
  if (!dependency$available) {
    sb_abort(
      paste0(
        "Package `", dependency$package, "` is required for learner `",
        learner, "`."
      ),
      if (identical(learner, "slope_binary")) {
        "stormbreaker_missing_slope"
      } else {
        "stormbreaker_missing_randomglm"
      }
    )
  }
  invisible(TRUE)
}

sb_slope_randomglm_binary_parameter_semantics <- function(learner) {
  params <- sb_slope_randomglm_binary_params(learner)
  switch(learner,
    slope_binary = paste0(
      "SLOPE binomial path with q in {",
      paste(params$q, collapse = ","), "}, gamma=", params$gamma,
      ", mean centering, L2 scaling, and ", params$path_length,
      " alpha values; q and alpha selected by native ",
      params$inner_folds, "-fold deviance CV inside each analysis fold with ",
      "seed ", params$seed
    ),
    random_glm_binary = paste0(
      "randomGLM ensemble with ", params$n_bags,
      " bootstrap bags, the native random-subspace rule, at most ",
      params$n_candidate_covariates,
      " forward-selection candidates per bag, no interactions, serial ",
      "execution, and seed ", params$seed, " fixed before resampling"
    )
  )
}

sb_new_slope_binary_model <- function(features, params, intercept, coefficient,
                                      tuning) {
  structure(
    list(
      learner = "slope_binary",
      features = features,
      params = params,
      intercept = intercept,
      coefficient = coefficient,
      tuning = tuning
    ),
    class = "stormbreaker_slope_binary_model"
  )
}

sb_fit_slope_binary <- function(inputs, params) {
  sb_reject_retired_learner("slope_binary")
}


sb_compact_random_glm_binary <- function(fit, features, params) {
  intercepts <- unname(as.numeric(fit$interceptOfForwardRegression))
  bag_coefficients <- unname(fit$coefOfForwardRegression)
  bag_features <- lapply(fit$featuresInForwardRegression, function(indices) {
    as.integer(indices)
  })
  valid <- length(intercepts) == params$n_bags &&
    length(bag_coefficients) == params$n_bags &&
    length(bag_features) == params$n_bags &&
    !anyNA(intercepts) && all(is.finite(intercepts))
  if (valid) {
    valid <- all(vapply(seq_len(params$n_bags), function(index) {
      coefficients <- unname(as.numeric(bag_coefficients[[index]]))
      indices <- bag_features[[index]]
      length(coefficients) == length(indices) &&
        !anyNA(coefficients) && all(is.finite(coefficients)) &&
        !anyNA(indices) && all(indices >= 1L & indices <= length(features))
    }, logical(1)))
  }
  if (!isTRUE(valid)) {
    sb_abort(
      "randomGLM returned an invalid bagged logistic ensemble.",
      "stormbreaker_invalid_randomglm_fit"
    )
  }
  structure(
    list(
      learner = "random_glm_binary",
      features = features,
      params = params,
      intercepts = intercepts,
      bag_features = bag_features,
      bag_coefficients = lapply(bag_coefficients, function(x) {
        unname(as.numeric(x))
      })
    ),
    class = "stormbreaker_random_glm_binary_model"
  )
}

sb_fit_random_glm_binary <- function(inputs, params) {
  fit <- withr::with_seed(
    params$seed,
    randomGLM::randomGLM(
      x = as.matrix(inputs$x),
      y = factor(inputs$outcome, levels = c(0L, 1L)),
      type = "binary",
      maxInteractionOrder = params$max_interaction_order,
      nBags = params$n_bags,
      replace = params$replace,
      nCandidateCovariates = min(
        params$n_candidate_covariates,
        ncol(inputs$x)
      ),
      keepModels = FALSE,
      randomSeed = params$seed,
      nThreads = params$threads,
      verbose = 0
    )
  )
  sb_compact_random_glm_binary(fit, names(inputs$x), params)
}

sb_fit_slope_randomglm_binary <- function(x, outcome, learner, params) {
  sb_reject_retired_learner(learner)
  sb_require_slope_randomglm_binary_dependency(learner)
  inputs <- sb_validate_classical_binary_inputs(x, outcome, learner)
  if (identical(learner, "slope_binary")) {
    return(sb_fit_slope_binary(inputs, params))
  }
  sb_fit_random_glm_binary(inputs, params)
}

sb_validate_slope_randomglm_binary_x <- function(model, x) {
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      paste0("Prediction columns do not match learner `", model$learner, "`."),
      "stormbreaker_learner_schema_mismatch"
    )
  }
  x
}

sb_predict_slope_randomglm_binary <- function(model, x) {
  x <- sb_validate_slope_randomglm_binary_x(model, x)
  matrix_x <- as.matrix(x)
  probability <- if (inherits(model, "stormbreaker_slope_binary_model")) {
    stats::plogis(model$intercept + matrix_x %*% model$coefficient)
  } else if (inherits(model, "stormbreaker_random_glm_binary_model")) {
    bag_probabilities <- lapply(seq_along(model$intercepts), function(index) {
      selected <- model$bag_features[[index]]
      linear_predictor <- rep(model$intercepts[[index]], nrow(matrix_x))
      if (length(selected) > 0L) {
        linear_predictor <- linear_predictor + as.numeric(
          matrix_x[, selected, drop = FALSE] %*%
            model$bag_coefficients[[index]]
        )
      }
      stats::plogis(linear_predictor)
    })
    Reduce(`+`, bag_probabilities) / length(bag_probabilities)
  } else {
    sb_abort(
      "The supplied SLOPE/randomGLM binary model is invalid.",
      "stormbreaker_bad_slope_randomglm_fit"
    )
  }
  sb_class_one_probability(probability, nrow(x), model$learner)
}

sb_extract_slope_binary_coefficients <- function(model) {
  data.frame(
    term = c("(Intercept)", model$features),
    coefficient = c(model$intercept, unname(model$coefficient)),
    nonzero = c(model$intercept, unname(model$coefficient)) != 0,
    stringsAsFactors = FALSE
  )
}

sb_extract_slope_randomglm_binary_state <- function(model, context = list()) {
  state <- unclass(model)
  state$diagnostic_path <- NULL
  class(state) <- c(
    "stormbreaker_slope_randomglm_binary_state",
    "stormbreaker_learner_state"
  )
  state
}

sb_restore_slope_randomglm_binary_model <- function(state) {
  if (identical(state$learner, "slope_binary")) {
    return(sb_new_slope_binary_model(
      state$features,
      state$params,
      state$intercept,
      state$coefficient,
      state$tuning
    ))
  }
  if (identical(state$learner, "random_glm_binary")) {
    model <- state
    class(model) <- "stormbreaker_random_glm_binary_model"
    return(model)
  }
  sb_bad_learner_state("The saved SLOPE/randomGLM binary state is invalid.")
}

sb_predict_slope_randomglm_binary_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_slope_randomglm_binary_learner_ids()) {
    sb_bad_learner_state("The saved SLOPE/randomGLM binary state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_restore_slope_randomglm_binary_model(state)
  list(probability = sb_predict_slope_randomglm_binary(model, x))
}

sb_slope_randomglm_binary_learner_adapter <- function(
  learner,
  params = sb_slope_randomglm_binary_params(learner)
) {
  adapter <- list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_slope_randomglm_binary(x, outcome, learner, params)
    },
    predict = sb_predict_slope_randomglm_binary,
    extract_state = sb_extract_slope_randomglm_binary_state,
    predict_state = sb_predict_slope_randomglm_binary_state
  )
  if (identical(learner, "slope_binary")) {
    adapter$extract_coefficients <- sb_extract_slope_binary_coefficients
  }
  adapter
}
