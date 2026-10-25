sb_bagged_knn_survival_learner_ids <- function() {
  "bagged_knn_survival"
}

sb_bagged_knn_survival_available <- function() {
  requireNamespace("bnnSurvival", quietly = TRUE)
}

sb_require_bagged_knn_survival_dependency <- function() {
  if (!sb_bagged_knn_survival_available()) {
    sb_abort(
      paste0(
        "Package `bnnSurvival` is required for bagged k-nearest-neighbor ",
        "survival prediction."
      ),
      "stormbreaker_missing_bnn_survival"
    )
  }
  invisible(TRUE)
}

sb_bagged_knn_survival_params <- function(
  learner = sb_bagged_knn_survival_learner_ids()
) {
  if (!learner %in% sb_bagged_knn_survival_learner_ids()) {
    sb_abort(
      paste0("Unknown bagged kNN survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    num_base_learners = 50L,
    neighbor_rule = "max(1,floor(n_analysis/10))",
    feature_rule = "all_selected_features",
    metric = "mahalanobis",
    distance_weighting = "constant_one",
    replace = TRUE,
    sample_fraction = 1,
    cores = 1L,
    seed = 1L
  )
}

sb_bagged_knn_survival_parameter_semantics <- function(learner) {
  params <- sb_bagged_knn_survival_params(learner)
  paste0(
    "bagged k-nearest-neighbor survival probability with ",
    params$num_base_learners, " bootstrap learners, k=max(1,floor(n/10)), ",
    "all selected features, Mahalanobis distance, constant neighbor weights, ",
    "one-core execution, and seed ", params$seed,
    " fixed before resampling"
  )
}



sb_with_bagged_knn_single_core <- function(code) {
  previous <- getOption("mc.cores")
  on.exit(options(mc.cores = previous), add = TRUE)
  options(mc.cores = 1L)
  force(code)
}

sb_bagged_knn_constant_weight <- function(distance) {
  distance * 0 + 1
}

sb_bagged_knn_internal_names <- function(p) {
  paste0(".stormbreaker_bknn_x", seq_len(p))
}

sb_bagged_knn_internal_data <- function(x) {
  internal <- x
  internal_features <- sb_bagged_knn_internal_names(ncol(x))
  names(internal) <- internal_features
  uses_dummy <- ncol(x) == 1L
  if (uses_dummy) {
    internal$.stormbreaker_bknn_zero <- 0
  }
  list(
    data = internal,
    internal_features = names(internal),
    uses_zero_variance_bridge = uses_dummy
  )
}

sb_bagged_knn_formula_data <- function(inputs) {
  prepared <- sb_bagged_knn_internal_data(inputs$x)
  analysis <- prepared$data
  rownames(analysis) <- NULL
  analysis$.stormbreaker_bknn_time <- inputs$time
  analysis$.stormbreaker_bknn_event <- inputs$event
  response <- paste0(
    "Surv(.stormbreaker_bknn_time, .stormbreaker_bknn_event)"
  )
  formula <- stats::as.formula(
    paste(response, "~", paste(prepared$internal_features, collapse = " + ")),
    env = asNamespace("survival")
  )
  c(prepared, list(analysis = analysis, formula = formula))
}

sb_validate_bagged_knn_engine <- function(engine, internal_features,
                                          effective_k, num_base_learners) {
  valid <- methods::is(engine, "bnnSurvivalEnsemble") &&
    identical(
      as.integer(methods::slot(engine, "num_base_learners")),
      as.integer(num_base_learners)
    ) &&
    identical(as.integer(methods::slot(engine, "k")), as.integer(effective_k)) &&
    identical(
      as.integer(methods::slot(engine, "num_features_per_base_learner")),
      as.integer(length(internal_features))
    ) &&
    identical(
      colnames(methods::slot(engine, "train_data"))[-c(1L, 2L)],
      internal_features
    )
  timepoints <- if (isTRUE(valid)) {
    as.numeric(methods::slot(engine, "timepoints"))
  } else {
    numeric()
  }
  valid <- isTRUE(valid) && length(timepoints) > 0L && !anyNA(timepoints) &&
    all(is.finite(timepoints)) && all(timepoints >= 0) &&
    all(diff(timepoints) > 0)
  if (!isTRUE(valid)) {
    sb_abort(
      "The bagged kNN survival learner returned an invalid native fit.",
      "stormbreaker_invalid_bagged_knn_survival_fit"
    )
  }
  engine
}

sb_validate_bagged_knn_model <- function(model) {
  valid <- inherits(model, "stormbreaker_bagged_knn_survival_model") &&
    is.list(model) &&
    identical(model$learner, sb_bagged_knn_survival_learner_ids()) &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.character(model$internal_features) &&
    length(model$internal_features) %in% c(
      length(model$features),
      length(model$features) + 1L
    ) &&
    is.logical(model$uses_zero_variance_bridge) &&
    length(model$uses_zero_variance_bridge) == 1L &&
    is.numeric(model$effective_k) && length(model$effective_k) == 1L &&
    !is.na(model$effective_k) && model$effective_k >= 1L &&
    is.list(model$params)
  if (!isTRUE(valid)) {
    sb_abort(
      "The supplied model is not a valid fitted bagged kNN survival learner.",
      "stormbreaker_bad_bagged_knn_survival_fit"
    )
  }
  sb_validate_bagged_knn_engine(
    model$engine,
    model$internal_features,
    model$effective_k,
    model$params$num_base_learners
  )
  model
}

sb_fit_bagged_knn_survival <- function(x, time, event, learner, params) {
  sb_require_bagged_knn_survival_dependency()
  inputs <- sb_validate_classical_survival_inputs(x, time, event, learner)
  prepared <- sb_bagged_knn_formula_data(inputs)
  effective_k <- max(1L, floor(nrow(inputs$x) / 10L))
  engine <- sb_with_bagged_knn_single_core(
    withr::with_seed(
      params$seed,
      bnnSurvival::bnnSurvival(
        formula = prepared$formula,
        data = prepared$analysis,
        k = effective_k,
        num_base_learners = params$num_base_learners,
        num_features_per_base_learner = length(prepared$internal_features),
        metric = params$metric,
        weighting_function = sb_bagged_knn_constant_weight,
        replace = params$replace,
        sample_fraction = params$sample_fraction
      )
    )
  )
  model <- structure(
    list(
      learner = learner,
      features = names(inputs$x),
      internal_features = prepared$internal_features,
      uses_zero_variance_bridge = prepared$uses_zero_variance_bridge,
      effective_k = as.integer(effective_k),
      params = params,
      engine = engine
    ),
    class = "stormbreaker_bagged_knn_survival_model"
  )
  sb_validate_bagged_knn_model(model)
}

sb_bagged_knn_prediction_data <- function(model, x) {
  model <- sb_validate_bagged_knn_model(model)
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      paste0(
        "Prediction columns do not match the fitted bagged kNN survival ",
        "learner schema."
      ),
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prepared <- sb_bagged_knn_internal_data(x)
  if (!identical(prepared$internal_features, model$internal_features) ||
    !identical(
      prepared$uses_zero_variance_bridge,
      model$uses_zero_variance_bridge
    )) {
    sb_abort(
      "The internal bagged kNN survival predictor schema is invalid.",
      "stormbreaker_bad_bagged_knn_survival_fit"
    )
  }
  list(model = model, x = x, internal = prepared$data)
}

sb_bagged_knn_distance_function <- function(k) {
  inverse <- product <- transpose_product <- diagonal <- NULL
  bound_scale <- train_magnitude <- NULL
  function(x, center, cov, ...) {
    x <- as.matrix(x)
    if (is.null(inverse)) {
      inverse <<- solve(cov, ...)
      product <<- x %*% inverse
      transpose_product <<- x %*% t(inverse)
      diagonal <<- rowSums(product * x)
      train_magnitude <<- max(abs(x))
      bound_scale <<- 128 * .Machine$double.eps * ncol(x)^3 *
        max(abs(inverse))
    }
    # Expanded quadratic products screen neighbors only. Retain both cross
    # terms: a numerically inverted covariance need not be exactly symmetric.
    values <- diagonal - as.vector(product %*% center) -
      as.vector(transpose_product %*% center) +
      sum((center %*% inverse) * center)
    # Cancellation can change near-ties. This conservative absolute rounding
    # envelope covers both the expanded and centered quadratic evaluations;
    # the normal-minimum term also covers underflow. Recompute every possible
    # top-k neighbor with the original centered native expression below.
    bound <- bound_scale * (train_magnitude + max(abs(center)))^2 +
      128 * .Machine$double.xmin * ncol(x)^3
    exact <- function(rows = NULL) {
      xx <- if (is.null(rows)) x else x[rows, , drop = FALSE]
      stats::mahalanobis(xx, center, inverse, inverted = TRUE)
    }
    if (any(!is.finite(values)) || !is.finite(bound) || k >= nrow(x)) {
      return(exact())
    }
    index <- order(values)
    selected <- which(values <= values[index[k]] + 2 * bound)
    values[selected] <- exact(selected)
    stats::setNames(values, rownames(x))
  }
}

sb_predict_bagged_knn_engine <- function(object, test_data) {
  ns <- asNamespace("bnnSurvival")
  # The native 0.1.5 S4 methods store their implementation as a .local closure.
  # Use local copies of those closures, not namespace/S4-table replacement or
  # a second implementation of native sorting, KM and ensemble aggregation.
  base_predict <- body(methods::getMethod("predict",
    "bnnSurvivalBaseLearner", where = ns))[[2L]][[3L]]
  predict_base <- function(object, train_data, test_data, timepoints,
                           metric, weighting_function, k) {
    # A fresh closure belongs only to this bootstrap learner and this call.
    distance <- sb_bagged_knn_distance_function(k)
    local_base <- base_predict
    environment(local_base) <- list2env(list(mahalanobis = distance),
      parent = environment(base_predict))
    local_base(object, train_data, test_data, timepoints,
      metric, weighting_function, k)
  }
  native_ensemble <- body(methods::getMethod("predict",
    "bnnSurvivalEnsemble", where = ns))[[2L]][[3L]]
  local_ensemble <- native_ensemble
  environment(local_ensemble) <- list2env(list(predict = predict_base),
    parent = environment(native_ensemble))
  local_ensemble(object, test_data)
}

sb_predict_bagged_knn_survival_curves <- function(model, x) {
  prepared <- sb_bagged_knn_prediction_data(model, x)
  prediction_object <- sb_with_bagged_knn_single_core(
    sb_predict_bagged_knn_engine(prepared$model$engine, prepared$internal)
  )
  survival <- as.matrix(bnnSurvival::predictions(prediction_object))
  timepoints <- as.numeric(bnnSurvival::timepoints(prediction_object))
  valid <- nrow(survival) == nrow(prepared$x) &&
    ncol(survival) == length(timepoints) && length(timepoints) > 0L &&
    !anyNA(timepoints) && all(is.finite(timepoints)) &&
    all(timepoints >= 0) && all(diff(timepoints) > 0) &&
    !anyNA(survival) && all(is.finite(survival)) &&
    all(survival >= 0 & survival <= 1) &&
    all(apply(survival, 1L, function(values) {
      all(diff(values) <= 1e-10)
    }))
  if (!isTRUE(valid)) {
    sb_abort(
      "The bagged kNN learner produced invalid survival curves.",
      "stormbreaker_invalid_bagged_knn_survival_prediction"
    )
  }
  interval <- diff(c(0, timepoints))
  previous <- if (ncol(survival) == 1L) {
    matrix(1, nrow = nrow(survival), ncol = 1L)
  } else {
    cbind(1, survival[, -ncol(survival), drop = FALSE])
  }
  risk <- -as.vector(previous %*% interval)
  if (length(risk) != nrow(prepared$x) || anyNA(risk) ||
    any(!is.finite(risk))) {
    sb_abort(
      "The bagged kNN learner produced invalid restricted-mean risk scores.",
      "stormbreaker_invalid_bagged_knn_survival_prediction"
    )
  }
  list(
    risk = unname(risk),
    time = timepoints,
    survival = survival
  )
}

sb_predict_bagged_knn_survival_risk <- function(model, x) {
  sb_predict_bagged_knn_survival_curves(model, x)$risk
}

sb_bagged_knn_direct_survival_source <- function() {
  "bnn_survival_bagged_knn_direct"
}

sb_predict_bagged_knn_survival_probability <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon,
  prediction_details = NULL
) {
  inputs <- sb_validate_classical_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      paste0(
        "Analysis columns do not match the fitted bagged kNN survival ",
        "learner schema."
      ),
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  curves <- if (is.null(prediction_details)) {
    sb_predict_bagged_knn_survival_curves(model, x_assessment)
  } else {
    prediction_details
  }
  eligible <- which(curves$time <= prediction_horizon)
  probability <- if (length(eligible) == 0L) {
    rep(1, length(curves$risk))
  } else {
    unname(curves$survival[, max(eligible)])
  }
  list(
    prediction_horizon = prediction_horizon,
    source = sb_bagged_knn_direct_survival_source(),
    ties = "not_applicable",
    probability = probability
  )
}

sb_extract_bagged_knn_survival_state <- function(model, context = list()) {
  model <- sb_validate_bagged_knn_model(model)
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    horizon <- sb_validate_prediction_horizon(
      context$prediction_horizon,
      context$time
    )
    survival_probability <- list(
      prediction_horizon = horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      source = sb_bagged_knn_direct_survival_source(),
      ties = "not_applicable"
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      internal_features = model$internal_features,
      uses_zero_variance_bridge = model$uses_zero_variance_bridge,
      effective_k = model$effective_k,
      params = model$params,
      engine_serialization = "R_xdr_v3_bnnSurvival_prediction_engine",
      engine_raw = serialize(
        model$engine,
        connection = NULL,
        ascii = FALSE,
        xdr = TRUE,
        version = 3L
      ),
      training_data_retention =
        "required_by_bagged_knn_prediction_state",
      risk_reference = "bnnSurvival_native_training_time_grid",
      risk_formula = "negative_step_RMST_over_native_training_time_grid",
      baseline_hazard = NULL,
      survival_probability = survival_probability
    ),
    class = c(
      "stormbreaker_bagged_knn_survival_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_bagged_knn_survival_model <- function(state) {
  if (!is.raw(state$engine_raw) || length(state$engine_raw) == 0L) {
    sb_bad_learner_state(
      "The saved bagged kNN survival prediction engine is invalid."
    )
  }
  engine <- tryCatch(
    unserialize(state$engine_raw),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved bagged kNN survival engine could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_bagged_knn_survival_state",
        parent = error
      )
    }
  )
  model <- structure(
    list(
      learner = state$learner,
      features = state$features,
      internal_features = state$internal_features,
      uses_zero_variance_bridge = state$uses_zero_variance_bridge,
      effective_k = state$effective_k,
      params = state$params,
      engine = engine
    ),
    class = "stormbreaker_bagged_knn_survival_model"
  )
  sb_validate_bagged_knn_model(model)
}

sb_validate_bagged_knn_probability_state <- function(state) {
  probability <- state$survival_probability
  valid <- is.list(probability) &&
    is.numeric(probability$prediction_horizon) &&
    length(probability$prediction_horizon) == 1L &&
    !is.na(probability$prediction_horizon) &&
    is.finite(probability$prediction_horizon) &&
    probability$prediction_horizon > 0 &&
    is.character(probability$time_unit) &&
    length(probability$time_unit) == 1L &&
    !is.na(probability$time_unit) && nzchar(probability$time_unit) &&
    identical(
      probability$source,
      sb_bagged_knn_direct_survival_source()
    ) && identical(probability$ties, "not_applicable")
  if (!isTRUE(valid)) {
    sb_bad_learner_state(
      "The saved bagged kNN direct survival-probability state is invalid."
    )
  }
  probability
}

sb_predict_bagged_knn_survival_state <- function(state, x) {
  valid <- inherits(state, "stormbreaker_bagged_knn_survival_state") &&
    is.list(state) &&
    identical(state$learner, sb_bagged_knn_survival_learner_ids()) &&
    identical(
      state$training_data_retention,
      "required_by_bagged_knn_prediction_state"
    )
  if (!isTRUE(valid)) {
    sb_bad_learner_state("The saved bagged kNN survival state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, state$learner)
  sb_require_bagged_knn_survival_dependency()
  model <- sb_restore_bagged_knn_survival_model(state)
  curves <- sb_predict_bagged_knn_survival_curves(model, x)
  result <- list(risk = curves$risk)
  if (!is.null(state$survival_probability)) {
    probability <- sb_validate_bagged_knn_probability_state(state)
    eligible <- which(curves$time <= probability$prediction_horizon)
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- if (length(eligible) == 0L) {
      rep(1, length(curves$risk))
    } else {
      unname(curves$survival[, max(eligible)])
    }
  }
  result
}

sb_bagged_knn_survival_learner_adapter <- function(
  learner,
  params = sb_bagged_knn_survival_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_bagged_knn_survival(x, time, event, learner, params)
    },
    predict = sb_predict_bagged_knn_survival_risk,
    predict_details = sb_predict_bagged_knn_survival_curves,
    predict_survival = sb_predict_bagged_knn_survival_probability,
    extract_state = sb_extract_bagged_knn_survival_state,
    predict_state = sb_predict_bagged_knn_survival_state
  )
}
