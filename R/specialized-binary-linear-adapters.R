sb_specialized_binary_linear_learner_ids <- function() {
  c(
    "stepglm_binary",
    "glmboost_binary",
    "plsrglm_binary"
  )
}

sb_validate_stepglm_direction <- function(direction) {
  if (!is.character(direction) || length(direction) != 1L ||
    is.na(direction) || !direction %in% c("forward", "backward", "both")) {
    sb_abort(
      "`stepglm_direction` must be one of `forward`, `backward`, or `both`.",
      "stormbreaker_bad_stepglm_direction"
    )
  }
  unname(as.character(direction))
}

sb_specialized_binary_linear_params <- function(
  learner,
  stepglm_direction = NULL
) {
  if (!learner %in% sb_specialized_binary_linear_learner_ids()) {
    sb_abort(
      paste0("Unknown specialized binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  if (identical(learner, "stepglm_binary")) {
    if (is.null(stepglm_direction)) {
      stepglm_direction <- "forward"
    }
    return(list(
      direction = sb_validate_stepglm_direction(stepglm_direction),
      criterion = "AIC",
      k = 2,
      max_steps = 1000L
    ))
  }
  if (!is.null(stepglm_direction)) {
    sb_abort(
      paste0(
        "`stepglm_direction` applies only when ",
        "`learner = 'stepglm_binary'`."
      ),
      "stormbreaker_irrelevant_stepglm_parameter"
    )
  }
  if (identical(learner, "glmboost_binary")) {
    return(list(
      iterations = 100L,
      shrinkage = 0.1,
      center = FALSE,
      family = "binomial_logit"
    ))
  }
  list(
    components = 2L,
    scale_x = TRUE,
    model = "pls-glm-logistic"
  )
}

sb_specialized_binary_linear_dependency <- function(learner) {
  package <- switch(learner,
    stepglm_binary = NA_character_,
    glmboost_binary = "mboost",
    plsrglm_binary = "plsRglm",
    sb_abort(
      paste0("Unknown specialized binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
  list(
    package = package,
    available = is.na(package) || requireNamespace(package, quietly = TRUE)
  )
}

sb_require_specialized_binary_linear_dependency <- function(learner) {
  dependency <- sb_specialized_binary_linear_dependency(learner)
  if (!dependency$available) {
    sb_abort(
      paste0(
        "Package `", dependency$package, "` is required for learner `",
        learner, "`."
      ),
      "stormbreaker_missing_specialized_binary_dependency"
    )
  }
  invisible(dependency)
}

sb_specialized_binary_stop_on_warning <- function(code, learner) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      sb_abort(
        paste0(
          "Learner `", learner, "` raised a fitting warning: ",
          conditionMessage(warning)
        ),
        paste0(
          "stormbreaker_",
          sub("_binary$", "", learner),
          "_warning"
        ),
        parent = warning
      )
    }
  )
}

sb_new_specialized_binary_linear_model <- function(
  learner, features, params, intercept, coefficient, center, effective
) {
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      intercept = intercept,
      coefficient = coefficient,
      center = center,
      effective = effective
    ),
    class = "stormbreaker_specialized_binary_linear"
  )
}

sb_fit_stepglm_binary <- function(inputs, params, capture_history = FALSE) {
  analysis <- inputs$x
  outcome_name <- ".stormbreaker_outcome"
  while (outcome_name %in% names(analysis)) {
    outcome_name <- paste0(outcome_name, "_")
  }
  analysis[[outcome_name]] <- inputs$outcome
  rhs <- paste(vapply(names(inputs$x), sb_quote_name, ""), collapse = " + ")
  full_formula <- stats::as.formula(
    paste(sb_quote_name(outcome_name), "~", rhs),
    env = environment()
  )
  null_formula <- stats::as.formula(
    paste(sb_quote_name(outcome_name), "~ 1"),
    env = environment()
  )
  start_formula <- if (identical(params$direction, "forward")) {
    null_formula
  } else {
    full_formula
  }
  model <- sb_specialized_binary_stop_on_warning(
    {
      start <- stats::glm(
        start_formula,
        data = analysis,
        family = stats::binomial(link = "logit"),
        model = FALSE,
        x = FALSE,
        y = FALSE
      )
      stats::step(
        start,
        scope = list(lower = null_formula, upper = full_formula),
        direction = params$direction,
        trace = 0,
        keep = NULL,
        steps = params$max_steps,
        k = params$k
      )
    },
    "stepglm_binary"
  )
  fitted <- stats::coef(model)
  selected <- setdiff(names(fitted), "(Intercept)")
  valid <- isTRUE(model$converged) &&
    identical(names(fitted)[[1L]], "(Intercept)") &&
    all(selected %in% names(inputs$x)) &&
    identical(as.integer(model$rank), as.integer(length(fitted))) &&
    !anyNA(fitted) && all(is.finite(fitted))
  if (!valid) {
    sb_abort(
      paste0(
        "Fixed Stepwise Logistic did not return a converged, full-rank fit ",
        "with finite coefficients."
      ),
      "stormbreaker_incomplete_stepglm_fit"
    )
  }
  coefficient <- stats::setNames(rep(0, ncol(inputs$x)), names(inputs$x))
  coefficient[selected] <- fitted[selected]
  list(
    intercept = unname(fitted[["(Intercept)"]]),
    coefficient = coefficient,
    center = stats::setNames(rep(0, ncol(inputs$x)), names(inputs$x)),
    effective = list(
      selected_features = selected,
      direction = params$direction,
      aic = unname(stats::AIC(model)),
      selection_path = if (capture_history) model$anova else NULL
    )
  )
}

sb_fit_glmboost_binary <- function(inputs, params, capture_history = FALSE) {
  model <- sb_specialized_binary_stop_on_warning(
    mboost::glmboost(
      x = as.matrix(inputs$x),
      y = factor(inputs$outcome, levels = c(0L, 1L)),
      family = mboost::Binomial(type = "glm", link = "logit"),
      control = mboost::boost_control(
        mstop = params$iterations,
        nu = params$shrinkage,
        risk = "inbag",
        stopintern = FALSE,
        trace = FALSE
      ),
      center = params$center
    ),
    "glmboost_binary"
  )
  fitted <- stats::coef(model)
  offset <- attr(fitted, "offset", exact = TRUE)
  selected <- names(fitted)
  center <- model$center
  valid <- is.numeric(fitted) && all(selected %in% names(inputs$x)) &&
    !anyNA(fitted) && all(is.finite(fitted)) &&
    is.numeric(offset) && length(offset) == 1L &&
    !is.na(offset) && is.finite(offset) &&
    is.numeric(center) && length(center) == ncol(inputs$x) &&
    !anyNA(center) && all(is.finite(center)) &&
    identical(as.integer(mboost::mstop(model)), params$iterations)
  if (!valid) {
    sb_abort(
      "Fixed glmBoost did not return finite linear prediction parameters.",
      "stormbreaker_incomplete_glmboost_fit"
    )
  }
  coefficient <- stats::setNames(rep(0, ncol(inputs$x)), names(inputs$x))
  coefficient[selected] <- fitted
  center <- stats::setNames(as.numeric(center), names(inputs$x))
  list(
    intercept = unname(offset),
    coefficient = coefficient,
    center = center,
    effective = list(
      selected_features = selected,
      iterations = params$iterations,
      training_loss = if (capture_history) data.frame(step = seq_along(mboost::risk(model)) - 1L,
        loss = as.numeric(mboost::risk(model))) else NULL
    )
  )
}

sb_fit_plsrglm_binary <- function(inputs, params) {
  if ("Intercept" %in% names(inputs$x)) {
    sb_abort(
      "Predictor name `Intercept` is reserved by the plsRglm engine.",
      "stormbreaker_reserved_plsrglm_intercept_name"
    )
  }
  model <- sb_specialized_binary_stop_on_warning(
    plsRglm::PLS_glm(
      dataY = inputs$outcome,
      dataX = inputs$x,
      nt = params$components,
      modele = params$model,
      typeVC = "none",
      scaleX = params$scale_x,
      pvals.expli = FALSE,
      MClassed = FALSE,
      sparse = FALSE,
      naive = FALSE,
      verbose = FALSE
    ),
    "plsrglm_binary"
  )
  fitted <- if (is.list(model)) drop(model$Coeffs) else NULL
  expected_terms <- c("Intercept", names(inputs$x))
  effective_components <- if (is.list(model)) model$computed_nt else NULL
  valid <- is.numeric(fitted) && identical(names(fitted), expected_terms) &&
    !anyNA(fitted) && all(is.finite(fitted)) &&
    is.numeric(effective_components) && length(effective_components) == 1L &&
    !is.na(effective_components) && is.finite(effective_components) &&
    effective_components >= 1L &&
    effective_components <= min(params$components, ncol(inputs$x)) &&
    is.list(model$FinalModel) && isTRUE(model$FinalModel$converged) &&
    !isTRUE(model$FinalModel$boundary)
  if (!valid) {
    sb_abort(
      "Fixed plsRglm did not return finite logistic prediction parameters.",
      "stormbreaker_incomplete_plsrglm_fit"
    )
  }
  list(
    intercept = unname(fitted[["Intercept"]]),
    coefficient = fitted[names(inputs$x)],
    center = stats::setNames(rep(0, ncol(inputs$x)), names(inputs$x)),
    effective = list(components = as.integer(effective_components)),
    diagnostic_components = list(loadings = model$pp, scores = model$tt,
      loading_definition = "Native plsRglm X reconstruction loadings")
  )
}

sb_fit_specialized_binary_linear <- function(
  x, outcome, learner, params
) {
  inputs <- sb_validate_classical_binary_inputs(x, outcome, learner)
  sb_validate_binary_predictor_names(names(inputs$x))
  sb_require_specialized_binary_linear_dependency(learner)
  fitted <- switch(learner,
    stepglm_binary = sb_fit_stepglm_binary(inputs, params),
    glmboost_binary = sb_fit_glmboost_binary(inputs, params),
    plsrglm_binary = sb_fit_plsrglm_binary(inputs, params)
  )
  model <- sb_new_specialized_binary_linear_model(
    learner = learner,
    features = names(inputs$x),
    params = params,
    intercept = fitted$intercept,
    coefficient = fitted$coefficient,
    center = fitted$center,
    effective = fitted$effective
  )
  if (!is.null(fitted$diagnostic_components)) model$diagnostic_components <- fitted$diagnostic_components
  model
}

sb_predict_specialized_binary_linear <- function(model, x) {
  if (!inherits(model, "stormbreaker_specialized_binary_linear")) {
    sb_abort(
      "The supplied model is not a fitted specialized binary learner.",
      "stormbreaker_bad_specialized_binary_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the specialized binary learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  matrix_x <- as.matrix(x)
  centered <- sweep(matrix_x, 2L, model$center, "-")
  probability <- stats::plogis(
    model$intercept + as.vector(centered %*% model$coefficient)
  )
  sb_class_one_probability(probability, nrow(x), model$learner)
}

sb_extract_specialized_binary_linear_coefficients <- function(model) {
  original_intercept <- model$intercept -
    sum(model$center * model$coefficient)
  data.frame(
    term = c("(Intercept)", model$features),
    coefficient = unname(c(original_intercept, model$coefficient)),
    stringsAsFactors = FALSE
  )
}

sb_extract_specialized_binary_linear_state <- function(model) {
  if (!inherits(model, "stormbreaker_specialized_binary_linear")) {
    sb_bad_learner_state(
      "The fitted specialized binary learner state is invalid."
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      intercept = model$intercept,
      coefficient = model$coefficient,
      center = model$center,
      effective = model$effective
    ),
    class = c(
      "stormbreaker_specialized_binary_linear_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_specialized_binary_linear_model <- function(state, learner) {
  valid <- is.list(state) && identical(state$learner, learner) &&
    is.character(state$features) && length(state$features) > 0L &&
    !anyNA(state$features) && !anyDuplicated(state$features) &&
    is.list(state$params) &&
    is.numeric(state$intercept) && length(state$intercept) == 1L &&
    !is.na(state$intercept) && is.finite(state$intercept) &&
    is.numeric(state$coefficient) &&
    identical(names(state$coefficient), state$features) &&
    !anyNA(state$coefficient) && all(is.finite(state$coefficient)) &&
    is.numeric(state$center) &&
    identical(names(state$center), state$features) &&
    !anyNA(state$center) && all(is.finite(state$center)) &&
    is.list(state$effective)
  if (!valid) {
    sb_bad_learner_state("The saved specialized binary learner state is invalid.")
  }
  sb_new_specialized_binary_linear_model(
    learner = learner,
    features = state$features,
    params = state$params,
    intercept = state$intercept,
    coefficient = state$coefficient,
    center = state$center,
    effective = state$effective
  )
}

sb_predict_specialized_binary_linear_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_specialized_binary_linear_learner_ids()) {
    sb_bad_learner_state("The saved specialized binary learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_restore_specialized_binary_linear_model(state, learner)
  list(probability = sb_predict_specialized_binary_linear(model, x))
}

sb_specialized_binary_linear_learner_adapter <- function(
  learner,
  params = sb_specialized_binary_linear_params(learner)
) {
  if (!learner %in% sb_specialized_binary_linear_learner_ids()) {
    sb_abort(
      paste0("Unknown specialized binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_specialized_binary_linear(x, outcome, learner, params)
    },
    predict = sb_predict_specialized_binary_linear,
    extract_coefficients = sb_extract_specialized_binary_linear_coefficients,
    extract_state = function(model, context = list()) {
      sb_extract_specialized_binary_linear_state(model)
    },
    predict_state = sb_predict_specialized_binary_linear_state
  )
}
