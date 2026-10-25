sb_coxphf_learner_ids <- function() {
  "coxphf_firth"
}

sb_coxphf_available <- function() {
  requireNamespace("coxphf", quietly = TRUE)
}

sb_require_coxphf_dependency <- function() {
  if (!sb_coxphf_available()) {
    sb_abort(
      "Package `coxphf` is required for Firth-corrected Cox regression.",
      "stormbreaker_missing_coxphf"
    )
  }
  invisible(TRUE)
}

sb_coxphf_learner_params <- function(
  learner = sb_coxphf_learner_ids()
) {
  if (!learner %in% sb_coxphf_learner_ids()) {
    sb_abort(
      paste0("Unknown Firth Cox learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    firth = TRUE,
    penalty = 0.5,
    profile_likelihood = FALSE,
    max_iterations = 100L,
    max_step_halvings = 5L,
    convergence_epsilon = 1e-6,
    gradient_convergence = 1e-4,
    maximum_standardized_step = 0.5,
    ties = "breslow"
  )
}

sb_coxphf_parameter_semantics <- function(learner) {
  params <- sb_coxphf_learner_params(learner)
  paste0(
    "Firth bias-reduced Cox regression with penalty ", params$penalty,
    ", Breslow ties, Wald summaries without profile-likelihood fitting, and ",
    params$max_iterations, " maximum iterations fixed before resampling"
  )
}

sb_coxphf_stop_on_warning <- function(code, stage) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      sb_abort(
        paste0(
          "The Firth Cox learner raised a ", stage, " warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_coxphf_warning",
        parent = warning
      )
    }
  )
}

sb_coxphf_internal_names <- function(p) {
  paste0(".stormbreaker_coxphf_x", seq_len(p))
}

sb_coxphf_formula_data <- function(inputs) {
  analysis <- inputs$x
  internal_features <- sb_coxphf_internal_names(ncol(analysis))
  names(analysis) <- internal_features
  analysis$.stormbreaker_coxphf_time <- inputs$time
  analysis$.stormbreaker_coxphf_event <- inputs$event
  rownames(analysis) <- as.character(seq_len(nrow(analysis)))
  response <- paste0(
    "Surv(.stormbreaker_coxphf_time, .stormbreaker_coxphf_event)"
  )
  formula <- stats::as.formula(
    paste(response, "~", paste(internal_features, collapse = " + ")),
    env = asNamespace("coxphf")
  )
  list(
    data = analysis,
    formula = formula,
    internal_features = internal_features
  )
}

sb_validate_coxphf_model <- function(model) {
  valid <- inherits(model, "stormbreaker_coxphf_model") && is.list(model) &&
    identical(model$learner, sb_coxphf_learner_ids()) &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.character(model$internal_features) &&
    length(model$internal_features) == length(model$features) &&
    is.numeric(model$coefficient) &&
    length(model$coefficient) == length(model$features) &&
    !anyNA(model$coefficient) && all(is.finite(model$coefficient)) &&
    is.list(model$params)
  if (!isTRUE(valid)) {
    sb_abort(
      "The supplied model is not a valid fitted Firth Cox learner.",
      "stormbreaker_bad_coxphf_fit"
    )
  }
  model
}

sb_fit_coxphf_learner <- function(x, time, event, learner, params) {
  sb_require_coxphf_dependency()
  inputs <- sb_validate_classical_survival_inputs(
    x,
    time,
    event,
    learner
  )
  prepared <- sb_coxphf_formula_data(inputs)
  engine <- sb_coxphf_stop_on_warning(
    coxphf::coxphf(
      formula = prepared$formula,
      data = prepared$data,
      pl = params$profile_likelihood,
      maxit = params$max_iterations,
      maxhs = params$max_step_halvings,
      epsilon = params$convergence_epsilon,
      gconv = params$gradient_convergence,
      maxstep = params$maximum_standardized_step,
      firth = params$firth,
      penalty = params$penalty
    ),
    "fitting"
  )
  coefficient <- unname(engine$coefficients)
  valid <- inherits(engine, "coxphf") &&
    identical(names(engine$coefficients), prepared$internal_features) &&
    length(coefficient) == ncol(inputs$x) && !anyNA(coefficient) &&
    all(is.finite(coefficient)) && identical(engine$method.ties, params$ties)
  if (!isTRUE(valid)) {
    sb_abort(
      "The Firth Cox learner returned invalid coefficients.",
      "stormbreaker_invalid_coxphf_fit"
    )
  }
  sb_validate_coxphf_model(structure(
    list(
      learner = learner,
      features = names(inputs$x),
      internal_features = prepared$internal_features,
      coefficient = coefficient,
      params = params,
      iterations = as.integer(engine$iter)
    ),
    class = "stormbreaker_coxphf_model"
  ))
}

sb_validate_coxphf_model_x <- function(model, x) {
  model <- sb_validate_coxphf_model(model)
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted Firth Cox learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  list(model = model, x = x)
}

sb_predict_coxphf_learner <- function(model, x) {
  validated <- sb_validate_coxphf_model_x(model, x)
  risk <- unname(as.vector(
    as.matrix(validated$x) %*% validated$model$coefficient
  ))
  if (length(risk) != nrow(validated$x) || anyNA(risk) ||
    any(!is.finite(risk))) {
    sb_abort(
      "The Firth Cox learner produced invalid risk predictions.",
      "stormbreaker_invalid_coxphf_prediction"
    )
  }
  risk
}

sb_extract_coxphf_coefficients <- function(model) {
  model <- sb_validate_coxphf_model(model)
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

sb_predict_coxphf_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
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
  analysis_risk <- sb_predict_coxphf_learner(model, inputs$x)
  assessment_risk <- sb_predict_coxphf_learner(model, x_assessment)
  cumulative_hazard <- sb_breslow_baseline_hazard(
    inputs$time,
    inputs$event,
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = "coxphf_firth_breslow",
    ties = model$params$ties,
    risk = assessment_risk
  )
}

sb_extract_coxphf_state <- function(model, context = list()) {
  model <- sb_validate_coxphf_model(model)
  sb_extract_cox_learner_state(
    model,
    learner = model$learner,
    coefficients = sb_extract_coxphf_coefficients(model),
    predict_survival = sb_predict_coxphf_survival,
    context = context,
    params = model$params,
    audit = list(
      native_linear_predictors_used = FALSE,
      fitting_iterations = model$iterations,
      estimation = "Firth_penalized_partial_likelihood"
    )
  )
}

sb_coxphf_learner_adapter <- function(
  learner,
  params = sb_coxphf_learner_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_coxphf_learner(x, time, event, learner, params)
    },
    predict = sb_predict_coxphf_learner,
    predict_survival = sb_predict_coxphf_survival,
    extract_coefficients = sb_extract_coxphf_coefficients,
    extract_state = sb_extract_coxphf_state,
    predict_state = sb_predict_cox_learner_state
  )
}
