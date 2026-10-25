sb_flexsurv_parametric_learner_ids <- function() {
  c(
    "flexsurv_gamma",
    "flexsurv_gompertz",
    "flexsurv_gengamma"
  )
}

sb_flexsurv_parametric_params <- function(learner) {
  switch(learner,
    flexsurv_gamma = list(
      distribution = "gamma",
      parameterization = "shape_rate",
      location_parameter = "rate",
      risk_transform = "log_location"
    ),
    flexsurv_gompertz = list(
      distribution = "gompertz",
      parameterization = "shape_rate",
      location_parameter = "rate",
      risk_transform = "log_location"
    ),
    flexsurv_gengamma = list(
      distribution = "gengamma",
      parameterization = "stable_mu_sigma_Q",
      location_parameter = "mu",
      risk_transform = "negative_location"
    ),
    sb_abort(
      paste0("Unknown fixed flexsurv learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_flexsurv_parametric_parameter_semantics <- function(learner) {
  params <- sb_flexsurv_parametric_params(learner)
  switch(params$distribution,
    gamma = paste0(
      "Gamma shape-rate distribution and native rate regression fixed; ",
      "parameters estimated in each analysis fold"
    ),
    gompertz = paste0(
      "Gompertz shape-rate proportional-hazards distribution fixed; ",
      "parameters estimated in each analysis fold"
    ),
    gengamma = paste0(
      "stable generalized-gamma mu-sigma-Q parameterization fixed; ",
      "parameters estimated in each analysis fold"
    )
  )
}

sb_flexsurv_parametric_dependency <- function(learner) {
  if (!learner %in% sb_flexsurv_parametric_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed flexsurv learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    package = "flexsurv",
    packages = "flexsurv",
    available = requireNamespace("flexsurv", quietly = TRUE)
  )
}

sb_require_flexsurv_parametric_dependency <- function(learner) {
  dependency <- sb_flexsurv_parametric_dependency(learner)
  if (!dependency$available) {
    sb_abort(
      paste0("Package `flexsurv` is required for learner `", learner, "`."),
      "stormbreaker_missing_flexsurv_dependency"
    )
  }
  invisible(dependency)
}

sb_flexsurv_parametric_ancillary_names <- function(distribution) {
  switch(distribution,
    gamma = "shape",
    gompertz = "shape",
    gengamma = c("sigma", "Q")
  )
}

sb_new_flexsurv_parametric_model <- function(
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
    class = "stormbreaker_flexsurv_parametric_fixed"
  )
}

sb_fit_flexsurv_parametric <- function(
  x, time, event, learner, params
) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_flexsurv_parametric_dependency(learner)

  numerical <- sb_flexsurv_numerical_settings(inputs$x)
  engine_x <- sb_flexsurv_numerical_predictors(inputs$x, numerical)
  engine_features <- sprintf(".stormbreaker_x%04d", seq_len(ncol(engine_x)))
  names(engine_x) <- engine_features
  formula_data <- sb_additional_survival_formula_data(
    engine_x,
    inputs$time,
    inputs$event
  )
  native_arguments <- list(
    formula = formula_data$formula,
    data = formula_data$data,
    dist = params$distribution,
    hessian = TRUE,
    control = list(maxit = numerical$maxit)
  )
  if (identical(params$distribution, "gompertz")) {
    # Use native normalized-time starting values and optimizer units, while
    # keeping the actual fit and all returned parameters in the input units.
    # Otherwise a long censored time can round its initial CDF to exactly one.
    time_scale <- mean(inputs$time)
    native_arguments$inits <- c(0.001 / time_scale, 1 / time_scale,
                               rep(0, ncol(engine_x)))
    native_arguments$control <- list(
      parscale = c(1 / time_scale, rep(1, ncol(engine_x) + 1L)),
      reltol = 1e-10,
      maxit = numerical$maxit
    )
  }
  critical_warnings <- character()
  engine <- tryCatch(
    withCallingHandlers(
      do.call(flexsurv::flexsurvreg, native_arguments),
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
          "Fixed flexsurv ", params$distribution, " fitting failed: ",
          conditionMessage(error)
        ),
        "stormbreaker_flexsurv_fit_failure",
        parent = error
      )
    }
  )
  if (length(critical_warnings) > 0L) {
    sb_abort(
      paste0(
        "Fixed flexsurv ", params$distribution,
        " reported an unreliable likelihood fit: ",
        paste(unique(critical_warnings), collapse = " | ")
      ),
      "stormbreaker_flexsurv_nonconvergence"
    )
  }
  convergence <- if (is.list(engine$opt)) engine$opt$convergence else NULL
  covariance_valid <- is.matrix(engine$cov) &&
    !anyNA(engine$cov) && all(is.finite(engine$cov))
  if (length(convergence) != 1L || is.na(convergence) || convergence != 0) {
    sb_abort(
      paste0(
        "Fixed flexsurv ", params$distribution,
        " did not converge (optimizer code ",
        if (length(convergence) == 1L && !is.na(convergence)) {
          convergence
        } else {
          "unavailable"
        },
        ")."
      ),
      "stormbreaker_flexsurv_nonconvergence"
    )
  }
  if (!covariance_valid) {
    sb_abort(
      paste0(
        "Fixed flexsurv ", params$distribution,
        " did not return a finite Hessian covariance matrix."
      ),
      "stormbreaker_flexsurv_nonconvergence"
    )
  }

  coefficients <- stats::coef(engine)
  ancillary_names <- sb_flexsurv_parametric_ancillary_names(
    params$distribution
  )
  ancillary <- unname(engine$res[ancillary_names, "est"])
  names(ancillary) <- ancillary_names
  valid <- is.numeric(coefficients) &&
    params$location_parameter %in% names(coefficients) &&
    all(engine_features %in% names(coefficients)) &&
    !anyNA(coefficients) && all(is.finite(coefficients)) &&
    is.numeric(ancillary) && length(ancillary) == length(ancillary_names) &&
    !anyNA(ancillary) && all(is.finite(ancillary))
  if (valid && params$distribution %in% c("gamma", "gengamma")) {
    valid <- ancillary[[1L]] > 0
  }
  if (!valid) {
    sb_abort(
      paste0(
        "Fixed flexsurv ", params$distribution,
        " did not return finite prediction parameters."
      ),
      "stormbreaker_incomplete_flexsurv_fit"
    )
  }
  prediction <- list(
    intercept = unname(coefficients[[params$location_parameter]]) -
      sum(unname(coefficients[engine_features]) / numerical$scale * numerical$center),
    coefficient = unname(coefficients[engine_features] / numerical$scale),
    ancillary = ancillary
  )
  model <- sb_new_flexsurv_parametric_model(
    learner,
    names(inputs$x),
    engine_features,
    params,
    prediction,
    engine
  )
  model$numerical_fitting <- numerical
  sb_predict_flexsurv_parametric(model, inputs$x)
  model
}

sb_validate_flexsurv_parametric_model_x <- function(model, x) {
  if (!inherits(model, "stormbreaker_flexsurv_parametric_fixed") ||
    !model$learner %in% sb_flexsurv_parametric_learner_ids()) {
    sb_abort(
      "The supplied model is not a fitted flexsurv parametric learner.",
      "stormbreaker_bad_flexsurv_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the flexsurv learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  sb_require_flexsurv_parametric_dependency(model$learner)
  x
}

sb_flexsurv_parametric_location <- function(model, x) {
  model$prediction$intercept +
    as.vector(as.matrix(x) %*% model$prediction$coefficient)
}

sb_predict_flexsurv_parametric <- function(model, x) {
  x <- sb_validate_flexsurv_parametric_model_x(model, x)
  location <- sb_flexsurv_parametric_location(model, x)
  risk <- if (identical(model$params$risk_transform, "log_location")) {
    location
  } else {
    -location
  }
  risk <- unname(as.numeric(risk))
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid risk scores."),
      "stormbreaker_invalid_flexsurv_prediction"
    )
  }
  risk
}

sb_flexsurv_parametric_direct_probability <- function(model, x, horizon) {
  x <- sb_validate_flexsurv_parametric_model_x(model, x)
  location <- sb_flexsurv_parametric_location(model, x)
  ancillary <- model$prediction$ancillary
  probability <- switch(model$params$distribution,
    gamma = stats::pgamma(
      horizon,
      shape = ancillary[["shape"]],
      rate = exp(location),
      lower.tail = FALSE
    ),
    gompertz = flexsurv::pgompertz(
      horizon,
      shape = ancillary[["shape"]],
      rate = exp(location),
      lower.tail = FALSE
    ),
    gengamma = flexsurv::pgengamma(
      horizon,
      mu = location,
      sigma = ancillary[["sigma"]],
      Q = ancillary[["Q"]],
      lower.tail = FALSE
    )
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
      "stormbreaker_invalid_flexsurv_prediction"
    )
  }
  probability
}

sb_flexsurv_parametric_direct_survival_source <- function(learner) {
  params <- sb_flexsurv_parametric_params(learner)
  paste0("flexsurv_", params$distribution, "_distribution_direct")
}

sb_predict_flexsurv_parametric_survival <- function(
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
      "Analysis columns do not match the flexsurv learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  list(
    prediction_horizon = prediction_horizon,
    source = sb_flexsurv_parametric_direct_survival_source(model$learner),
    ties = "not_applicable",
    probability = sb_flexsurv_parametric_direct_probability(
      model,
      x_assessment,
      prediction_horizon
    )
  )
}

sb_flexsurv_parametric_risk_reference <- function(learner) {
  params <- sb_flexsurv_parametric_params(learner)
  paste0(
    "flexsurv_", params$distribution, "_", params$parameterization,
    "_native"
  )
}

sb_flexsurv_parametric_risk_formula <- function(learner) {
  params <- sb_flexsurv_parametric_params(learner)
  if (identical(params$risk_transform, "log_location")) {
    "log_rate_location"
  } else {
    "negative_mu_log_time_location"
  }
}

sb_extract_flexsurv_parametric_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_flexsurv_parametric_fixed")) {
    sb_bad_learner_state("The fitted flexsurv learner state is invalid.")
  }
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_flexsurv_parametric_survival(
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
      intercept = model$prediction$intercept,
      coefficient = model$prediction$coefficient,
      ancillary = model$prediction$ancillary,
      numerical_fitting = model$numerical_fitting,
      risk_reference = sb_flexsurv_parametric_risk_reference(model$learner),
      risk_formula = sb_flexsurv_parametric_risk_formula(model$learner),
      survival_probability = survival_probability
    ),
    class = c(
      "stormbreaker_flexsurv_parametric_fixed_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_flexsurv_parametric_model <- function(state, learner) {
  params <- sb_flexsurv_parametric_params(learner)
  ancillary_names <- sb_flexsurv_parametric_ancillary_names(
    params$distribution
  )
  valid <- sb_pipeline_parameters_equal(state$params, params) &&
    is.numeric(state$intercept) && length(state$intercept) == 1L &&
    !is.na(state$intercept) && is.finite(state$intercept) &&
    is.numeric(state$coefficient) &&
    length(state$coefficient) == length(state$features) &&
    !anyNA(state$coefficient) && all(is.finite(state$coefficient)) &&
    is.numeric(state$ancillary) &&
    identical(names(state$ancillary), ancillary_names) &&
    !anyNA(state$ancillary) && all(is.finite(state$ancillary))
  if (valid && params$distribution %in% c("gamma", "gengamma")) {
    valid <- state$ancillary[[1L]] > 0
  }
  if (!valid) {
    sb_bad_learner_state("The saved flexsurv prediction state is invalid.")
  }
  sb_new_flexsurv_parametric_model(
    learner,
    state$features,
    sprintf(".stormbreaker_x%04d", seq_along(state$features)),
    params,
    list(
      intercept = state$intercept,
      coefficient = state$coefficient,
      ancillary = state$ancillary
    )
  )
}

sb_predict_flexsurv_parametric_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_flexsurv_parametric_learner_ids()) {
    sb_bad_learner_state("The saved flexsurv learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  sb_require_flexsurv_parametric_dependency(learner)
  model <- sb_restore_flexsurv_parametric_model(state, learner)
  result <- list(risk = sb_predict_flexsurv_parametric(model, x))
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
        sb_flexsurv_parametric_direct_survival_source(learner)
      ) &&
      identical(probability$ties, "not_applicable")
    if (!valid) {
      sb_bad_learner_state(
        "The saved flexsurv survival-probability state is invalid."
      )
    }
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <-
      sb_flexsurv_parametric_direct_probability(
        model,
        x,
        probability$prediction_horizon
      )
  }
  result
}

sb_flexsurv_parametric_learner_adapter <- function(
  learner,
  params = sb_flexsurv_parametric_params(learner)
) {
  if (!learner %in% sb_flexsurv_parametric_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed flexsurv learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_flexsurv_parametric(x, time, event, learner, params)
    },
    predict = sb_predict_flexsurv_parametric,
    predict_survival = sb_predict_flexsurv_parametric_survival,
    extract_state = sb_extract_flexsurv_parametric_state,
    predict_state = sb_predict_flexsurv_parametric_state
  )
}
