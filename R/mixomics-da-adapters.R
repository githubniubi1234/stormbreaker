sb_mixomics_da_learner_ids <- function() {
  c("plsda_binary", "splsda_binary")
}

sb_mixomics_da_params <- function(learner) {
  if (!learner %in% sb_mixomics_da_learner_ids()) {
    sb_abort(
      paste0("Unknown mixOmics discriminant learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    method = if (identical(learner, "plsda_binary")) {
      "PLS-DA"
    } else {
      "sPLS-DA"
    },
    requested_components = 2L,
    component_rule = "min_requested_features_analysis_rows_minus_one",
    scale = TRUE,
    tolerance = 1e-6,
    maximum_iterations = 100L,
    near_zero_variance_filter = FALSE,
    logratio = "none",
    keep_x_rule = if (identical(learner, "plsda_binary")) {
      "all_features"
    } else {
      "min_10_and_p_per_component"
    },
    probability_mapping = "binary_softmax_of_native_dummy_scores"
  )
}

sb_require_mixomics <- function() {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    sb_abort(
      "Package `mixOmics` is required for PLS-DA learners.",
      "stormbreaker_missing_mixomics"
    )
  }
  invisible(TRUE)
}

sb_mixomics_da_active_spec <- function(x, params) {
  components <- min(
    params$requested_components,
    ncol(x),
    nrow(x) - 1L
  )
  list(
    components = as.integer(components),
    keep_x = if (identical(params$method, "sPLS-DA")) {
      rep.int(as.integer(min(10L, ncol(x))), components)
    } else {
      NULL
    }
  )
}

sb_mixomics_da_native_fit <- function(x, grouping, params, active, learner) {
  fit <- function() {
    arguments <- list(
      X = data.matrix(x),
      Y = grouping,
      ncomp = active$components,
      scale = params$scale,
      tol = params$tolerance,
      max.iter = params$maximum_iterations,
      near.zero.var = params$near_zero_variance_filter,
      logratio = params$logratio,
      all.outputs = TRUE
    )
    if (identical(params$method, "sPLS-DA")) {
      arguments$keepX <- active$keep_x
      return(do.call(mixOmics::splsda, arguments))
    }
    do.call(mixOmics::plsda, arguments)
  }

  tryCatch(
    withCallingHandlers(
      fit(),
      warning = function(warning) {
        sb_abort(
          paste0(
            "Learner `", learner, "` raised a fitting warning: ",
            conditionMessage(warning)
          ),
          "stormbreaker_mixomics_da_fit_warning",
          parent = warning
        )
      }
    ),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) {
        stop(error)
      }
      sb_abort(
        paste0(
          "Learner `", learner, "` could not fit its mixOmics model: ",
          conditionMessage(error)
        ),
        "stormbreaker_mixomics_da_fit_error",
        parent = error
      )
    }
  )
}

sb_mixomics_da_projection <- function(native, features, components, learner) {
  x <- native$X
  y <- native$ind.mat
  variates <- native$variates$X
  loadings <- native$loadings$X
  index <- seq_len(components)
  p_matrix <- crossprod(x, variates)
  c_matrix <- crossprod(y, variates)
  coefficients <- loadings[, index, drop = FALSE] %*%
    solve(
      t(p_matrix[, index, drop = FALSE]) %*%
        loadings[, index, drop = FALSE]
    ) %*%
    t(c_matrix)[index, , drop = FALSE]
  rownames(coefficients) <- features

  engine <- list(
    coefficients = coefficients,
    x_center = unname(attr(x, "scaled:center")),
    x_scale = unname(attr(x, "scaled:scale")),
    y_center = unname(attr(y, "scaled:center")),
    y_scale = unname(attr(y, "scaled:scale"))
  )
  valid <- identical(dim(engine$coefficients), c(length(features), 2L)) &&
    all(vapply(engine, is.numeric, logical(1))) &&
    all(vapply(engine, function(value) {
      !anyNA(value) && all(is.finite(value))
    }, logical(1))) &&
    length(engine$x_center) == length(features) &&
    length(engine$x_scale) == length(features) &&
    length(engine$y_center) == 2L && length(engine$y_scale) == 2L &&
    all(engine$x_scale > 0) && all(engine$y_scale > 0)
  if (!valid) {
    sb_abort(
      paste0("Learner `", learner, "` returned an incomplete projection."),
      "stormbreaker_incomplete_mixomics_da_fit"
    )
  }
  engine
}

sb_fit_mixomics_da <- function(x, outcome, learner, params) {
  sb_require_mixomics()
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      paste0("A `", learner, "` analysis fold must contain both classes."),
      "stormbreaker_mixomics_da_single_class"
    )
  }

  active <- sb_mixomics_da_active_spec(x, params)
  native <- sb_mixomics_da_native_fit(
    x,
    factor(outcome, levels = c(0L, 1L)),
    params,
    active,
    learner
  )
  engine <- tryCatch(
    sb_mixomics_da_projection(native, names(x), active$components, learner),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) {
        stop(error)
      }
      sb_abort(
        paste0(
          "Learner `", learner, "` could not form its prediction state: ",
          conditionMessage(error)
        ),
        "stormbreaker_mixomics_da_projection_error",
        parent = error
      )
    }
  )
  component_features <- if (identical(params$method, "sPLS-DA")) {
    lapply(seq_len(active$components), function(component) {
      names(x)[native$loadings$X[, component] != 0]
    })
  } else {
    list()
  }
  active_params <- params
  active_params$effective_components <- active$components
  active_params$active_keep_x <- active$keep_x

  structure(
    list(
      learner = learner,
      features = names(x),
      params = active_params,
      engine = engine,
      diagnostic_components = list(loadings = native$loadings$X,
        scores = native$variates$X, loading_definition = "Native mixOmics X component weights"),
      audit = list(
        effective_components = active$components,
        active_keep_x = active$keep_x,
        native_iterations = as.integer(native$iter),
        component_features = component_features,
        probability_mapping = params$probability_mapping
      )
    ),
    class = "stormbreaker_mixomics_da"
  )
}

sb_predict_mixomics_da <- function(model, x) {
  if (!inherits(model, "stormbreaker_mixomics_da")) {
    sb_abort(
      "The fitted mixOmics discriminant model is invalid.",
      "stormbreaker_bad_mixomics_da_model"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted mixOmics model.",
      "stormbreaker_mixomics_da_schema_mismatch"
    )
  }

  matrix_x <- data.matrix(x)
  matrix_x <- sweep(matrix_x, 2L, model$engine$x_center, "-")
  matrix_x <- sweep(matrix_x, 2L, model$engine$x_scale, "/")
  scores <- matrix_x %*% model$engine$coefficients
  scores <- sweep(scores, 2L, model$engine$y_scale, "*")
  scores <- sweep(scores, 2L, model$engine$y_center, "+")
  probability <- stats::plogis(scores[, "1"] - scores[, "0"])
  probability <- unname(as.numeric(probability))
  if (length(probability) != nrow(x) || anyNA(probability) ||
      !all(is.finite(probability)) ||
      !all(probability >= 0 & probability <= 1)) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid probabilities."),
      "stormbreaker_invalid_mixomics_da_prediction"
    )
  }
  probability
}

sb_extract_mixomics_da_learner_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_mixomics_da")) {
    sb_bad_learner_state("The fitted mixOmics discriminant state is invalid.")
  }
  portable <- unclass(model)
  portable$diagnostic_components <- NULL
  structure(
    portable,
    class = c(
      "stormbreaker_mixomics_da_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_mixomics_da_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
      !learner %in% sb_mixomics_da_learner_ids()) {
    sb_bad_learner_state("The saved mixOmics discriminant state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- structure(unclass(state), class = "stormbreaker_mixomics_da")
  list(probability = sb_predict_mixomics_da(model, x))
}

sb_mixomics_da_learner_adapter <- function(
    learner,
    params = sb_mixomics_da_params(learner)) {
  if (!learner %in% sb_mixomics_da_learner_ids()) {
    sb_abort(
      paste0("Unknown mixOmics discriminant learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_mixomics_da(x, outcome, learner, params)
    },
    predict = sb_predict_mixomics_da,
    extract_state = sb_extract_mixomics_da_learner_state,
    predict_state = sb_predict_mixomics_da_learner_state
  )
}
