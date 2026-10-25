sb_ncvreg_binary_learner_ids <- function() {
  c("ncvreg_mcp_binary_bic", "ncvreg_scad_binary_bic")
}

sb_ncvreg_survival_learner_ids <- function() {
  c("ncvreg_mcp_cox_bic", "ncvreg_scad_cox_bic")
}

sb_ncvreg_learner_ids <- function() {
  c(sb_ncvreg_binary_learner_ids(), sb_ncvreg_survival_learner_ids())
}

sb_ncvreg_available <- function() {
  requireNamespace("ncvreg", quietly = TRUE)
}

sb_ncvreg_dependency <- function() {
  list(package = "ncvreg", available = sb_ncvreg_available())
}

sb_ncvreg_learner_params <- function(learner) {
  if (!learner %in% sb_ncvreg_learner_ids()) {
    sb_abort(
      paste0("Unknown ncvreg learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  penalty <- if (grepl("_mcp_", learner, fixed = TRUE)) {
    "MCP"
  } else {
    "SCAD"
  }
  list(
    family = if (learner %in% sb_ncvreg_binary_learner_ids()) {
      "binomial"
    } else {
      "cox"
    },
    penalty = penalty,
    gamma = if (identical(penalty, "MCP")) 3 else 3.7,
    alpha = 1,
    nlambda = 100L,
    lambda_min_ratio_n_gt_p = 0.001,
    lambda_min_ratio_n_le_p = 0.05,
    eps = 1e-4,
    max_iter = 10000L,
    convex_diagnostics = TRUE,
    standardization = "native_internal",
    selection = "analysis_fold_bic",
    bic_tie_rule = "first_minimum_on_descending_lambda_path",
    warning_policy = "fail_fit",
    return_x = FALSE
  )
}

sb_ncvreg_parameter_semantics <- function(learner) {
  params <- sb_ncvreg_learner_params(learner)
  paste0(
    params$penalty,
    " gamma ",
    sb_format_roundtrip_double(params$gamma),
    ", alpha 1, native standardization, 100-lambda path, and documented ",
    "lambda-min ratio rule fixed; lambda selected by BIC inside each ",
    "analysis fold"
  )
}

sb_require_ncvreg_dependency <- function() {
  if (!sb_ncvreg_available()) {
    sb_abort(
      "Package `ncvreg` is required for MCP and SCAD learners.",
      "stormbreaker_missing_ncvreg"
    )
  }
  invisible(TRUE)
}

sb_ncvreg_stop_on_warning <- function(code, learner, stage) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      sb_abort(
        paste0(
          "Learner `", learner, "` raised a ", stage, " warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_ncvreg_fit_warning",
        parent = warning
      )
    }
  )
}

sb_ncvreg_lambda_min_ratio <- function(params, n, p) {
  if (n > p) {
    params$lambda_min_ratio_n_gt_p
  } else {
    params$lambda_min_ratio_n_le_p
  }
}

sb_ncvreg_selected_coefficients <- function(engine, features, learner,
                                             selected_lambda) {
  coefficient <- drop(stats::coef(engine, lambda = selected_lambda))
  expected <- if (learner %in% sb_ncvreg_binary_learner_ids()) {
    c("(Intercept)", features)
  } else {
    features
  }
  if (is.numeric(coefficient) && length(coefficient) == length(expected) &&
    is.null(names(coefficient))) {
    names(coefficient) <- expected
  }
  valid <- is.numeric(coefficient) &&
    identical(names(coefficient), expected) &&
    !anyNA(coefficient) && all(is.finite(coefficient))
  if (!valid) {
    sb_abort(
      paste0(
        "Learner `", learner,
        "` did not return finite coefficients in the fitted feature order."
      ),
      "stormbreaker_invalid_ncvreg_coefficients"
    )
  }
  coefficient
}

sb_fit_ncvreg_engine <- function(inputs, learner, params) {
  sb_require_ncvreg_dependency()
  matrix_x <- as.matrix(inputs$x)
  storage.mode(matrix_x) <- "double"
  lambda_min <- sb_ncvreg_lambda_min_ratio(
    params,
    nrow(matrix_x),
    ncol(matrix_x)
  )
  dfmax <- ncol(matrix_x) + if (identical(params$family, "binomial")) 1L else 0L

  engine <- sb_ncvreg_stop_on_warning(
    if (identical(params$family, "binomial")) {
      ncvreg::ncvreg(
        X = matrix_x,
        y = as.numeric(inputs$outcome),
        family = "binomial",
        penalty = params$penalty,
        gamma = params$gamma,
        alpha = params$alpha,
        lambda.min = lambda_min,
        nlambda = params$nlambda,
        eps = params$eps,
        max.iter = params$max_iter,
        convex = params$convex_diagnostics,
        dfmax = dfmax,
        penalty.factor = rep(1, ncol(matrix_x)),
        warn = TRUE,
        returnX = params$return_x
      )
    } else {
      response <- cbind(
        time = as.numeric(inputs$time),
        event = as.numeric(inputs$event)
      )
      ncvreg::ncvsurv(
        X = matrix_x,
        y = response,
        penalty = params$penalty,
        gamma = params$gamma,
        alpha = params$alpha,
        lambda.min = lambda_min,
        nlambda = params$nlambda,
        eps = params$eps,
        max.iter = params$max_iter,
        convex = params$convex_diagnostics,
        dfmax = dfmax,
        penalty.factor = rep(1, ncol(matrix_x)),
        warn = TRUE,
        returnX = params$return_x
      )
    },
    learner,
    "fitting"
  )
  expected_class <- if (identical(params$family, "binomial")) {
    "ncvreg"
  } else {
    "ncvsurv"
  }
  if (!inherits(engine, expected_class)) {
    sb_abort(
      paste0("Learner `", learner, "` returned an invalid native fit."),
      "stormbreaker_invalid_ncvreg_fit"
    )
  }

  bic <- sb_ncvreg_stop_on_warning(
    stats::BIC(engine),
    learner,
    "BIC calculation"
  )
  valid_bic <- is.numeric(bic) && length(bic) == length(engine$lambda) &&
    length(bic) > 0L && !anyNA(bic) && all(is.finite(bic)) &&
    is.numeric(engine$lambda) && all(is.finite(engine$lambda)) &&
    all(engine$lambda > 0)
  if (!valid_bic) {
    sb_abort(
      paste0("Learner `", learner, "` returned an invalid BIC path."),
      "stormbreaker_invalid_ncvreg_bic"
    )
  }
  selected_index <- which.min(bic)
  selected_lambda <- unname(engine$lambda[[selected_index]])
  coefficient <- sb_ncvreg_selected_coefficients(
    engine,
    names(inputs$x),
    learner,
    selected_lambda
  )
  convex_min <- engine$convex.min
  if (!is.numeric(convex_min) || length(convex_min) != 1L ||
    is.na(convex_min) || !is.finite(convex_min)) {
    convex_min <- NA_integer_
  } else {
    convex_min <- as.integer(convex_min)
  }

  list(
    coefficient = coefficient,
    diagnostic_path = list(lambda = engine$lambda, beta = engine$beta, bic = bic,
      engine_class = class(engine)),
    effective = list(
      selected_index = as.integer(selected_index),
      selected_lambda = selected_lambda,
      selected_bic = unname(bic[[selected_index]]),
      path_length = as.integer(length(engine$lambda)),
      lambda_min_ratio = lambda_min,
      convex_min = convex_min
    )
  )
}

sb_new_ncvreg_model <- function(learner, features, params, fitted) {
  binary <- learner %in% sb_ncvreg_binary_learner_ids()
  coefficient <- fitted$coefficient
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      intercept = if (binary) unname(coefficient[[1L]]) else 0,
      coefficient = stats::setNames(
        unname(if (binary) coefficient[-1L] else coefficient),
        features
      ),
      effective = fitted$effective,
      diagnostic_path = fitted$diagnostic_path
    ),
    class = "stormbreaker_ncvreg_model"
  )
}

sb_fit_ncvreg_learner <- function(x, endpoint, learner, params) {
  inputs <- if (learner %in% sb_ncvreg_binary_learner_ids()) {
    sb_validate_classical_binary_inputs(
      x,
      endpoint$outcome,
      learner
    )
  } else {
    sb_validate_classical_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      learner
    )
  }
  fitted <- sb_fit_ncvreg_engine(inputs, learner, params)
  sb_new_ncvreg_model(learner, names(inputs$x), params, fitted)
}

sb_validate_ncvreg_model_x <- function(model, x) {
  valid <- inherits(model, "stormbreaker_ncvreg_model") &&
    model$learner %in% sb_ncvreg_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.numeric(model$coefficient) &&
    identical(names(model$coefficient), model$features) &&
    !anyNA(model$coefficient) && all(is.finite(model$coefficient))
  if (!valid) {
    sb_abort(
      "The supplied model is not a valid fitted ncvreg learner.",
      "stormbreaker_bad_ncvreg_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted ncvreg learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  x
}

sb_predict_ncvreg_learner <- function(model, x) {
  x <- sb_validate_ncvreg_model_x(model, x)
  score <- model$intercept + drop(as.matrix(x) %*% model$coefficient)
  prediction <- if (model$learner %in% sb_ncvreg_binary_learner_ids()) {
    stats::plogis(score)
  } else {
    score
  }
  valid <- length(prediction) == nrow(x) && !anyNA(prediction) &&
    all(is.finite(prediction))
  if (model$learner %in% sb_ncvreg_binary_learner_ids()) {
    valid <- valid && all(prediction >= 0 & prediction <= 1)
  }
  if (!valid) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid predictions."),
      "stormbreaker_invalid_ncvreg_prediction"
    )
  }
  unname(prediction)
}

sb_extract_ncvreg_coefficients <- function(model) {
  sb_validate_ncvreg_model_x(
    model,
    as.data.frame(
      matrix(
        0,
        nrow = 1L,
        ncol = length(model$features),
        dimnames = list(NULL, model$features)
      ),
      check.names = FALSE
    )
  )
  if (model$learner %in% sb_ncvreg_binary_learner_ids()) {
    values <- c(model$intercept, unname(model$coefficient))
    return(data.frame(
      term = c("(Intercept)", model$features),
      coefficient = values,
      nonzero = values != 0,
      stringsAsFactors = FALSE
    ))
  }
  values <- unname(model$coefficient)
  data.frame(
    feature = model$features,
    coefficient = values,
    nonzero = values != 0,
    direction = ifelse(
      values > 0,
      "higher_hazard",
      ifelse(values < 0, "lower_hazard", "zero")
    ),
    stringsAsFactors = FALSE
  )
}

sb_predict_ncvreg_cox_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  if (!model$learner %in% sb_ncvreg_survival_learner_ids()) {
    sb_abort(
      "A binary ncvreg learner cannot predict survival probability.",
      "stormbreaker_bad_ncvreg_fit"
    )
  }
  inputs <- sb_validate_classical_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  analysis_risk <- sb_predict_ncvreg_learner(model, inputs$x)
  assessment_risk <- sb_predict_ncvreg_learner(model, x_assessment)
  cumulative_hazard <- sb_breslow_baseline_hazard(
    inputs$time,
    inputs$event,
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = paste0(model$learner, "_breslow_native_link"),
    ties = "breslow",
    risk = assessment_risk
  )
}

sb_extract_ncvreg_state <- function(model, context = list()) {
  coefficients <- sb_extract_ncvreg_coefficients(model)
  state_params <- c(model$params, model$effective)
  if (model$learner %in% sb_ncvreg_binary_learner_ids()) {
    return(sb_extract_logistic_learner_state(
      model$learner,
      coefficients,
      params = state_params
    ))
  }
  sb_extract_cox_learner_state(
    model,
    model$learner,
    coefficients,
    sb_predict_ncvreg_cox_survival,
    context,
    params = state_params,
    audit = list(lambda_selection = "analysis_fold_bic")
  )
}

sb_predict_ncvreg_state <- function(state, x) {
  if (state$learner %in% sb_ncvreg_binary_learner_ids()) {
    return(sb_predict_logistic_learner_state(state, x))
  }
  sb_predict_cox_learner_state(state, x)
}

sb_ncvreg_learner_adapter <- function(
  learner,
  params = sb_ncvreg_learner_params(learner)
) {
  if (!learner %in% sb_ncvreg_learner_ids()) {
    sb_abort(
      paste0("Unknown ncvreg learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  adapter <- list(
    id = learner,
    params = params,
    fit = if (learner %in% sb_ncvreg_binary_learner_ids()) {
      function(x, outcome) {
        sb_fit_ncvreg_learner(
          x,
          list(outcome = outcome),
          learner,
          params
        )
      }
    } else {
      function(x, time, event) {
        sb_fit_ncvreg_learner(
          x,
          list(time = time, event = event),
          learner,
          params
        )
      }
    },
    predict = sb_predict_ncvreg_learner,
    extract_coefficients = sb_extract_ncvreg_coefficients,
    extract_state = sb_extract_ncvreg_state,
    predict_state = sb_predict_ncvreg_state
  )
  if (learner %in% sb_ncvreg_survival_learner_ids()) {
    adapter$predict_survival <- sb_predict_ncvreg_cox_survival
  }
  adapter
}
