sb_flexsurvspline_learner_ids <- function() {
  c(
    "flexsurvspline_hazard_k1",
    "flexsurvspline_odds_k1",
    "flexsurvspline_normal_k1"
  )
}

sb_flexsurvspline_params <- function(learner) {
  scale <- switch(learner,
    flexsurvspline_hazard_k1 = "hazard",
    flexsurvspline_odds_k1 = "odds",
    flexsurvspline_normal_k1 = "normal",
    sb_abort(
      paste0("Unknown fixed flexsurvspline learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
  list(
    scale = scale,
    k = 1L,
    timescale = "log",
    spline = "rp",
    knot_rule = "analysis_event_log_time_quantiles"
  )
}

sb_flexsurvspline_parameter_semantics <- function(learner) {
  params <- sb_flexsurvspline_params(learner)
  scale <- switch(params$scale,
    hazard = "log cumulative-hazard",
    odds = "log cumulative-odds",
    normal = "normal survival-transform"
  )
  paste0(
    "Royston-Parmar ", scale,
    " scale with one internal analysis-event log-time knot; ",
    "spline coefficients and knots estimated in each analysis fold"
  )
}

sb_flexsurvspline_dependency <- function(learner) {
  if (!learner %in% sb_flexsurvspline_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed flexsurvspline learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    package = "flexsurv",
    packages = "flexsurv",
    available = requireNamespace("flexsurv", quietly = TRUE)
  )
}

sb_require_flexsurvspline_dependency <- function(learner) {
  dependency <- sb_flexsurvspline_dependency(learner)
  if (!dependency$available) {
    sb_abort(
      paste0("Package `flexsurv` is required for learner `", learner, "`."),
      "stormbreaker_missing_flexsurvspline_dependency"
    )
  }
  invisible(dependency)
}

sb_flexsurvspline_gamma_names <- function(params) {
  paste0("gamma", 0:(params$k + 1L))
}

sb_new_flexsurvspline_model <- function(
  learner, features, engine_features, params, prediction, engine = NULL
) {
  structure(
    list(
      learner = learner,
      features = features,
      engine_features = engine_features,
      params = params,
      prediction = prediction,
      engine = engine
    ),
    class = "stormbreaker_flexsurvspline_fixed"
  )
}

sb_fit_flexsurvspline <- function(x, time, event, learner, params) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_flexsurvspline_dependency(learner)

  numerical <- sb_flexsurv_numerical_settings(inputs$x)
  engine_x <- sb_flexsurv_numerical_predictors(inputs$x, numerical)
  engine_features <- sprintf(".stormbreaker_x%04d", seq_len(ncol(engine_x)))
  names(engine_x) <- engine_features
  formula_data <- sb_additional_survival_formula_data(
    engine_x,
    inputs$time,
    inputs$event
  )
  critical_warnings <- character()
  engine <- tryCatch(
    withCallingHandlers(
      sb_flexsurvspline_numerical_fit(
        formula = formula_data$formula,
        data = formula_data$data,
        k = params$k,
        scale = params$scale,
        timescale = params$timescale,
        spline = params$spline,
        hessian = TRUE,
        maxit = numerical$maxit
      ),
      warning = function(warning) {
        message <- conditionMessage(warning)
        critical <- grepl(
          paste0(
            "Hessian not positive definite|Hessian is not finite|",
            "not converged to the maximum likelihood"
          ),
          message,
          ignore.case = TRUE
        )
        if (critical) {
          critical_warnings <<- c(critical_warnings, message)
          invokeRestart("muffleWarning")
        }
      }
    ),
    error = function(error) {
      sb_abort(
        paste0(
          "Fixed flexsurvspline ", params$scale, "-scale fitting failed: ",
          conditionMessage(error)
        ),
        "stormbreaker_flexsurvspline_fit_failure",
        parent = error
      )
    }
  )
  if (length(critical_warnings) > 0L) {
    sb_abort(
      paste0(
        "Fixed flexsurvspline ", params$scale,
        "-scale fit reported an unreliable likelihood fit: ",
        paste(unique(critical_warnings), collapse = " | ")
      ),
      "stormbreaker_flexsurvspline_nonconvergence"
    )
  }
  convergence <- if (is.list(engine$opt)) engine$opt$convergence else NULL
  covariance_valid <- is.matrix(engine$cov) &&
    !anyNA(engine$cov) && all(is.finite(engine$cov))
  if (length(convergence) != 1L || is.na(convergence) || convergence != 0) {
    sb_abort(
      paste0(
        "Fixed flexsurvspline ", params$scale,
        "-scale fit did not converge (optimizer code ",
        if (length(convergence) == 1L && !is.na(convergence)) {
          convergence
        } else {
          "unavailable"
        },
        ")."
      ),
      "stormbreaker_flexsurvspline_nonconvergence"
    )
  }
  if (!covariance_valid) {
    sb_abort(
      paste0(
        "Fixed flexsurvspline ", params$scale,
        "-scale fit did not return a finite Hessian covariance matrix."
      ),
      "stormbreaker_flexsurvspline_nonconvergence"
    )
  }

  coefficients <- stats::coef(engine)
  gamma_names <- sb_flexsurvspline_gamma_names(params)
  gamma <- unname(coefficients[gamma_names])
  names(gamma) <- gamma_names
  coefficient <- unname(coefficients[engine_features])
  knots <- unname(engine$knots)
  valid <- is.numeric(gamma) && length(gamma) == length(gamma_names) &&
    !anyNA(gamma) && all(is.finite(gamma)) &&
    is.numeric(coefficient) && length(coefficient) == length(engine_features) &&
    !anyNA(coefficient) && all(is.finite(coefficient)) &&
    is.numeric(knots) && length(knots) == params$k + 2L &&
    !anyNA(knots) && all(is.finite(knots)) &&
    length(engine$k) == 1L && isTRUE(engine$k == params$k) &&
    identical(engine$scale, params$scale) &&
    identical(engine$timescale, params$timescale) &&
    identical(engine$spline, params$spline)
  if (!valid) {
    sb_abort(
      paste0(
        "Fixed flexsurvspline ", params$scale,
        "-scale fit did not return finite prediction parameters."
      ),
      "stormbreaker_incomplete_flexsurvspline_fit"
    )
  }
  gamma[[1L]] <- gamma[[1L]] - sum(coefficient / numerical$scale * numerical$center)
  coefficient <- unname(coefficient / numerical$scale)
  numerical$initialization <- "training_intercept_only_spline"
  prediction <- list(
    gamma = gamma,
    coefficient = coefficient,
    knots = knots
  )
  model <- sb_new_flexsurvspline_model(
    learner,
    names(inputs$x),
    engine_features,
    params,
    prediction,
    engine
  )
  model$numerical_fitting <- numerical
  sb_predict_flexsurvspline(model, inputs$x)
  model
}

sb_validate_flexsurvspline_model_x <- function(model, x) {
  if (!inherits(model, "stormbreaker_flexsurvspline_fixed") ||
    !model$learner %in% sb_flexsurvspline_learner_ids()) {
    sb_abort(
      "The supplied model is not a fitted flexsurvspline learner.",
      "stormbreaker_bad_flexsurvspline_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the flexsurvspline learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  sb_require_flexsurvspline_dependency(model$learner)
  x
}

sb_flexsurvspline_location <- function(model, x) {
  model$prediction$gamma[[1L]] +
    as.vector(as.matrix(x) %*% model$prediction$coefficient)
}

sb_predict_flexsurvspline <- function(model, x) {
  x <- sb_validate_flexsurvspline_model_x(model, x)
  risk <- unname(as.numeric(sb_flexsurvspline_location(model, x)))
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid risk scores."),
      "stormbreaker_invalid_flexsurvspline_prediction"
    )
  }
  risk
}

sb_flexsurvspline_direct_probability <- function(model, x, horizon) {
  x <- sb_validate_flexsurvspline_model_x(model, x)
  gamma <- model$prediction$gamma
  gamma_matrix <- matrix(
    rep(gamma, each = nrow(x)),
    nrow = nrow(x),
    ncol = length(gamma)
  )
  gamma_matrix[, 1L] <- sb_flexsurvspline_location(model, x)
  probability <- flexsurv::psurvspline(
    q = rep(horizon, nrow(x)),
    gamma = gamma_matrix,
    knots = model$prediction$knots,
    scale = model$params$scale,
    timescale = model$params$timescale,
    spline = model$params$spline,
    lower.tail = FALSE
  )
  probability <- unname(as.numeric(probability))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      paste0(
        "Learner `", model$learner,
        "` produced invalid horizon survival probabilities."
      ),
      "stormbreaker_invalid_flexsurvspline_prediction"
    )
  }
  probability
}

sb_flexsurvspline_direct_survival_source <- function(learner) {
  params <- sb_flexsurvspline_params(learner)
  paste0(
    "flexsurvspline_", params$scale, "_k", params$k, "_",
    params$spline, "_distribution_direct"
  )
}

sb_predict_flexsurvspline_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  inputs <- sb_validate_additional_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the flexsurvspline learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  list(
    prediction_horizon = prediction_horizon,
    source = sb_flexsurvspline_direct_survival_source(model$learner),
    ties = "not_applicable",
    probability = sb_flexsurvspline_direct_probability(
      model,
      x_assessment,
      prediction_horizon
    )
  )
}

sb_flexsurvspline_risk_reference <- function(learner) {
  params <- sb_flexsurvspline_params(learner)
  paste0(
    "flexsurvspline_", params$scale, "_k", params$k, "_",
    params$spline, "_native_location"
  )
}

sb_flexsurvspline_risk_formula <- function(learner) {
  sb_flexsurvspline_params(learner)
  "gamma0_plus_x_beta_native_location"
}

sb_extract_flexsurvspline_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_flexsurvspline_fixed")) {
    sb_bad_learner_state("The fitted flexsurvspline learner state is invalid.")
  }
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_flexsurvspline_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    survival_probability <- list(
      prediction_horizon = probability$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      source = probability$source,
      ties = probability$ties
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      gamma = model$prediction$gamma,
      coefficient = model$prediction$coefficient,
      knots = model$prediction$knots,
      numerical_fitting = model$numerical_fitting,
      risk_reference = sb_flexsurvspline_risk_reference(model$learner),
      risk_formula = sb_flexsurvspline_risk_formula(model$learner),
      survival_probability = survival_probability
    ),
    class = c(
      "stormbreaker_flexsurvspline_fixed_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_flexsurvspline_model <- function(state, learner) {
  params <- sb_flexsurvspline_params(learner)
  gamma_names <- sb_flexsurvspline_gamma_names(params)
  valid <- sb_pipeline_parameters_equal(state$params, params) &&
    is.numeric(state$gamma) &&
    identical(names(state$gamma), gamma_names) &&
    !anyNA(state$gamma) && all(is.finite(state$gamma)) &&
    is.numeric(state$coefficient) &&
    length(state$coefficient) == length(state$features) &&
    !anyNA(state$coefficient) && all(is.finite(state$coefficient)) &&
    is.numeric(state$knots) && length(state$knots) == params$k + 2L &&
    !anyNA(state$knots) && all(is.finite(state$knots))
  if (!valid) {
    sb_bad_learner_state(
      "The saved flexsurvspline prediction state is invalid."
    )
  }
  sb_new_flexsurvspline_model(
    learner,
    state$features,
    sprintf(".stormbreaker_x%04d", seq_along(state$features)),
    params,
    list(
      gamma = state$gamma,
      coefficient = state$coefficient,
      knots = state$knots
    )
  )
}

sb_predict_flexsurvspline_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_flexsurvspline_learner_ids()) {
    sb_bad_learner_state("The saved flexsurvspline learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  sb_require_flexsurvspline_dependency(learner)
  model <- sb_restore_flexsurvspline_model(state, learner)
  result <- list(risk = sb_predict_flexsurvspline(model, x))
  if (!is.null(state$survival_probability)) {
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
        sb_flexsurvspline_direct_survival_source(learner)
      ) &&
      identical(probability$ties, "not_applicable")
    if (!valid) {
      sb_bad_learner_state(
        "The saved flexsurvspline survival-probability state is invalid."
      )
    }
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- sb_flexsurvspline_direct_probability(
      model,
      x,
      probability$prediction_horizon
    )
  }
  result
}

sb_flexsurvspline_learner_adapter <- function(
  learner,
  params = sb_flexsurvspline_params(learner)
) {
  if (!learner %in% sb_flexsurvspline_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed flexsurvspline learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_flexsurvspline(x, time, event, learner, params)
    },
    predict = sb_predict_flexsurvspline,
    predict_survival = sb_predict_flexsurvspline_survival,
    extract_state = sb_extract_flexsurvspline_state,
    predict_state = sb_predict_flexsurvspline_state
  )
}
