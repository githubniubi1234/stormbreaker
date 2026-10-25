sb_abess_binary_learner_ids <- function() {
  "abess_binary_gic"
}

sb_abess_survival_learner_ids <- function() {
  "abess_cox_gic"
}

sb_abess_learner_ids <- function() {
  c(sb_abess_binary_learner_ids(), sb_abess_survival_learner_ids())
}

sb_abess_selector_ids <- function() {
  c("abess_binary_gic_nonzero", "abess_cox_gic_nonzero")
}

sb_abess_available <- function() {
  requireNamespace("abess", quietly = TRUE)
}

sb_abess_dependency <- function() {
  list(
    package = "abess",
    available = sb_abess_available()
  )
}

sb_abess_learner_params <- function(learner) {
  if (!learner %in% sb_abess_learner_ids()) {
    sb_abort(
      paste0("Unknown ABESS learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  family <- if (learner %in% sb_abess_binary_learner_ids()) {
    "binomial"
  } else {
    "cox"
  }
  list(
    family = family,
    tune_path = "sequence",
    tune_type = "gic",
    minimum_support = 0L,
    maximum_support = 10L,
    normalize = if (identical(family, "binomial")) 2L else 3L,
    lambda = 0,
    c_max = 2L,
    splicing_type = 2L,
    max_splicing_iter = 20L,
    newton = "exact",
    max_newton_iter = 10L,
    early_stop = FALSE,
    warm_start = TRUE,
    important_search = "all_predictors",
    num_threads = 1L,
    seed = 1L
  )
}

sb_abess_selector_params <- function(selector) {
  learner <- switch(selector,
    abess_binary_gic_nonzero = "abess_binary_gic",
    abess_cox_gic_nonzero = "abess_cox_gic",
    sb_abort(
      paste0("Unknown ABESS selector `", selector, "`."),
      "stormbreaker_unsupported_selector"
    )
  )
  params <- sb_abess_learner_params(learner)
  params$minimum_support <- 2L
  params
}

sb_abess_stop_on_warning <- function(code, label, subclass) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      message <- conditionMessage(warning)
      changes_model <- grepl(
        "coerc|change to family|numerical result might be unstable",
        message,
        ignore.case = TRUE
      )
      if (!changes_model) return(invisible(NULL))
      sb_abort(
        paste0(label, " raised a warning: ", message),
        subclass,
        parent = warning
      )
    }
  )
}

sb_abess_support_path <- function(params, n_predictors, label) {
  maximum <- min(params$maximum_support, n_predictors)
  if (maximum < params$minimum_support) {
    sb_abort(
      paste0(
        label,
        " requires at least ",
        params$minimum_support,
        " predictors."
      ),
      "stormbreaker_abess_too_few_predictors"
    )
  }
  seq.int(params$minimum_support, maximum)
}

sb_fit_abess_engine <- function(x, endpoint, params, label) {
  if (!sb_abess_available()) {
    sb_abort(
      "Package `abess` is required for the requested ABESS component.",
      "stormbreaker_missing_abess"
    )
  }
  matrix_x <- as.matrix(x)
  storage.mode(matrix_x) <- "double"
  if (identical(params$family, "cox") && ncol(matrix_x) < 2L) {
    sb_abort(
      paste0(
        label,
        " requires at least 2 predictors because the native ABESS Cox ",
        "path does not preserve a one-column matrix."
      ),
      "stormbreaker_abess_cox_too_few_predictors"
    )
  }
  support_size <- sb_abess_support_path(params, ncol(matrix_x), label)
  important_search <- ncol(matrix_x)
  fit_intercept <- identical(params$family, "binomial")

  engine <- sb_abess_stop_on_warning(
    withr::with_seed(params$seed, {
      if (identical(params$family, "binomial")) {
        abess::abess(
          x = matrix_x,
          y = endpoint$outcome,
          family = "binomial",
          tune.path = params$tune_path,
          tune.type = params$tune_type,
          normalize = params$normalize,
          fit.intercept = fit_intercept,
          c.max = params$c_max,
          support.size = support_size,
          lambda = params$lambda,
          splicing.type = params$splicing_type,
          max.splicing.iter = params$max_splicing_iter,
          important.search = important_search,
          warm.start = params$warm_start,
          newton = params$newton,
          max.newton.iter = params$max_newton_iter,
          early.stop = params$early_stop,
          num.threads = params$num_threads,
          seed = params$seed
        )
      } else {
        abess::abess(
          x = matrix_x,
          y = survival::Surv(endpoint$time, endpoint$event),
          family = "cox",
          tune.path = params$tune_path,
          tune.type = params$tune_type,
          normalize = params$normalize,
          fit.intercept = fit_intercept,
          c.max = params$c_max,
          support.size = support_size,
          lambda = params$lambda,
          splicing.type = params$splicing_type,
          max.splicing.iter = params$max_splicing_iter,
          important.search = important_search,
          warm.start = params$warm_start,
          newton = params$newton,
          max.newton.iter = params$max_newton_iter,
          early.stop = params$early_stop,
          num.threads = params$num_threads,
          seed = params$seed
        )
      }
    }),
    paste0(label, " fit"),
    "stormbreaker_abess_fit_warning"
  )
  valid <- inherits(engine, "abess") &&
    is.numeric(engine$best.size) && length(engine$best.size) == 1L &&
    !is.na(engine$best.size) && engine$best.size %in% support_size
  if (!valid) {
    sb_abort(
      paste0(label, " did not return one valid GIC-selected support size."),
      "stormbreaker_invalid_abess_fit"
    )
  }
  engine$call <- NULL
  engine
}

sb_abess_coefficient_matrix <- function(engine) {
  coefficients <- as.matrix(stats::coef(
    engine,
    support.size = engine$best.size
  ))
  if (ncol(coefficients) != 1L || is.null(rownames(coefficients)) ||
    anyNA(coefficients) || any(!is.finite(coefficients))) {
    sb_abort(
      "The fitted ABESS model did not return portable coefficients.",
      "stormbreaker_invalid_abess_coefficients"
    )
  }
  coefficients
}

sb_new_abess_model <- function(learner, features, params, engine) {
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      selected_support_size = as.integer(engine$best.size),
      support_path = as.integer(engine$support.size),
      criterion_path = data.frame(
        support_size = as.integer(engine$support.size),
        gic = as.numeric(engine$tune.value),
        deviance = as.numeric(engine$dev),
        stringsAsFactors = FALSE
      ),
      engine = engine
    ),
    class = "stormbreaker_abess_model"
  )
}

sb_validate_abess_model_x <- function(model, x) {
  if (!inherits(model, "stormbreaker_abess_model") ||
    !model$learner %in% sb_abess_learner_ids()) {
    sb_abort(
      "The supplied model is not a fitted ABESS learner.",
      "stormbreaker_bad_abess_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted ABESS learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  x
}

sb_fit_abess_learner <- function(x, endpoint, learner, params) {
  if (learner %in% sb_abess_binary_learner_ids()) {
    inputs <- sb_validate_classical_binary_inputs(
      x,
      endpoint$outcome,
      learner
    )
    endpoint <- list(outcome = inputs$outcome)
  } else {
    inputs <- sb_validate_classical_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      learner
    )
    endpoint <- list(time = inputs$time, event = inputs$event)
  }
  engine <- sb_fit_abess_engine(inputs$x, endpoint, params, learner)
  sb_new_abess_model(learner, names(inputs$x), params, engine)
}

sb_predict_abess_learner <- function(model, x) {
  x <- sb_validate_abess_model_x(model, x)
  type <- if (model$learner %in% sb_abess_binary_learner_ids()) {
    "response"
  } else {
    "link"
  }
  prediction <- sb_abess_stop_on_warning(
    stats::predict(
      model$engine,
      as.matrix(x),
      type = type,
      support.size = model$selected_support_size
    ),
    paste0(model$learner, " prediction"),
    "stormbreaker_abess_prediction_warning"
  )
  prediction <- unname(as.numeric(prediction))
  valid <- length(prediction) == nrow(x) && !anyNA(prediction) &&
    all(is.finite(prediction))
  if (model$learner %in% sb_abess_binary_learner_ids()) {
    valid <- valid && all(prediction >= 0 & prediction <= 1)
  }
  if (!valid) {
    sb_abort(
      "The fitted ABESS learner produced invalid predictions.",
      "stormbreaker_invalid_abess_prediction"
    )
  }
  prediction
}

sb_extract_abess_coefficients <- function(model) {
  if (!inherits(model, "stormbreaker_abess_model")) {
    sb_abort(
      "The supplied model is not a fitted ABESS learner.",
      "stormbreaker_bad_abess_fit"
    )
  }
  coefficients <- sb_abess_coefficient_matrix(model$engine)
  expected_rows <- c("(intercept)", model$features)
  if (!identical(rownames(coefficients), expected_rows)) {
    sb_abort(
      "The fitted ABESS coefficient names do not match the learner schema.",
      "stormbreaker_invalid_abess_coefficients"
    )
  }
  values <- unname(as.numeric(coefficients[, 1L]))
  if (model$learner %in% sb_abess_binary_learner_ids()) {
    return(data.frame(
      term = c("(Intercept)", model$features),
      coefficient = values,
      nonzero = values != 0,
      stringsAsFactors = FALSE
    ))
  }
  predictor_values <- values[-1L]
  data.frame(
    feature = model$features,
    coefficient = predictor_values,
    nonzero = predictor_values != 0,
    direction = ifelse(
      predictor_values > 0,
      "higher_hazard",
      ifelse(predictor_values < 0, "lower_hazard", "zero")
    ),
    stringsAsFactors = FALSE
  )
}

sb_predict_abess_cox_survival <- function(
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
  analysis_risk <- sb_predict_abess_learner(model, inputs$x)
  assessment_risk <- sb_predict_abess_learner(model, x_assessment)
  cumulative_hazard <- sb_breslow_baseline_hazard(
    inputs$time,
    inputs$event,
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = "abess_cox_gic_breslow_native_link",
    ties = "breslow",
    risk = assessment_risk
  )
}

sb_extract_abess_learner_state <- function(model, context = list()) {
  coefficients <- sb_extract_abess_coefficients(model)
  state_params <- c(
    model$params,
    list(selected_support_size = model$selected_support_size)
  )
  if (model$learner %in% sb_abess_binary_learner_ids()) {
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
    sb_predict_abess_cox_survival,
    context,
    params = state_params,
    audit = list(
      support_path = model$support_path,
      criterion_path = model$criterion_path,
      support_selection = "analysis_fold_gic"
    )
  )
}

sb_predict_abess_learner_state <- function(state, x) {
  if (state$learner %in% sb_abess_binary_learner_ids()) {
    return(sb_predict_logistic_learner_state(state, x))
  }
  sb_predict_cox_learner_state(state, x)
}

sb_abess_learner_adapter <- function(
  learner,
  params = sb_abess_learner_params(learner)
) {
  if (!learner %in% sb_abess_learner_ids()) {
    sb_abort(
      paste0("Unknown ABESS learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  adapter <- list(
    id = learner,
    params = params,
    fit = if (learner %in% sb_abess_binary_learner_ids()) {
      function(x, outcome) {
        sb_fit_abess_learner(
          x,
          list(outcome = outcome),
          learner,
          params
        )
      }
    } else {
      function(x, time, event) {
        sb_fit_abess_learner(
          x,
          list(time = time, event = event),
          learner,
          params
        )
      }
    },
    predict = sb_predict_abess_learner,
    extract_coefficients = sb_extract_abess_coefficients,
    extract_state = sb_extract_abess_learner_state,
    predict_state = sb_predict_abess_learner_state
  )
  if (learner %in% sb_abess_survival_learner_ids()) {
    adapter$predict_survival <- sb_predict_abess_cox_survival
  }
  adapter
}

sb_fit_abess_selector <- function(x, endpoint, selector, capture_history = FALSE) {
  params <- sb_abess_selector_params(selector)
  learner <- if (identical(selector, "abess_binary_gic_nonzero")) {
    "abess_binary_gic"
  } else {
    "abess_cox_gic"
  }
  if (learner %in% sb_abess_binary_learner_ids()) {
    inputs <- sb_validate_classical_binary_inputs(
      x,
      endpoint$outcome,
      selector
    )
    endpoint <- list(outcome = inputs$outcome)
  } else {
    inputs <- sb_validate_classical_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      selector
    )
    endpoint <- list(time = inputs$time, event = inputs$event)
  }
  engine <- sb_fit_abess_engine(inputs$x, endpoint, params, selector)
  coefficients <- sb_abess_coefficient_matrix(engine)
  if (!identical(rownames(coefficients), c("(intercept)", names(inputs$x)))) {
    sb_abort(
      "The ABESS selector coefficient names do not match its input schema.",
      "stormbreaker_invalid_abess_selector"
    )
  }
  values <- unname(as.numeric(coefficients[-1L, 1L]))
  scores <- data.frame(
    feature = names(inputs$x),
    coefficient = values,
    absolute_coefficient = abs(values),
    selected = values != 0,
    stringsAsFactors = FALSE
  )
  selected_features <- scores$feature[scores$selected]
  if (length(selected_features) < 2L) {
    sb_abort(
      "The ABESS GIC selector retained fewer than two predictors.",
      "stormbreaker_abess_selector_too_few_features"
    )
  }
  structure(
    list(
      id = selector,
      source_features = names(inputs$x),
      selected_features = selected_features,
      selected_support_size = as.integer(engine$best.size),
      params = params,
      scores = scores,
      coefficient_path = if (capture_history) list(beta = engine$beta,
        support_size = engine$support.size) else NULL,
      criterion_path = data.frame(
        support_size = as.integer(engine$support.size),
        gic = as.numeric(engine$tune.value),
        deviance = as.numeric(engine$dev),
        stringsAsFactors = FALSE
      )
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_apply_abess_selector <- function(fitted, x) {
  valid <- inherits(fitted, "stormbreaker_selector_fit") &&
    fitted$id %in% sb_abess_selector_ids() &&
    is.character(fitted$source_features) &&
    is.character(fitted$selected_features) &&
    length(fitted$selected_features) >= 2L &&
    !anyNA(fitted$selected_features) &&
    !anyDuplicated(fitted$selected_features) &&
    all(fitted$selected_features %in% fitted$source_features)
  if (!valid) {
    sb_abort(
      "The fitted selector is not a valid ABESS GIC selector state.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$source_features)) {
    sb_abort(
      "Predictor columns do not match the fitted ABESS selector schema.",
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x[fitted$selected_features]
}

sb_abess_selector_adapter <- function(selector, capture_history = FALSE) {
  if (!selector %in% sb_abess_selector_ids()) {
    sb_abort(
      paste0("Unknown ABESS selector `", selector, "`."),
      "stormbreaker_unsupported_selector"
    )
  }
  list(
    id = selector,
    params = list(),
    capture_history = capture_history,
    needs_outcome = TRUE,
    fit = if (identical(selector, "abess_binary_gic_nonzero")) {
      function(x, outcome) {
        sb_fit_abess_selector(x, list(outcome = outcome), selector, capture_history)
      }
    } else {
      function(x, time, event) {
        sb_fit_abess_selector(
          x,
          list(time = time, event = event),
          selector, capture_history
        )
      }
    },
    apply = sb_apply_abess_selector
  )
}
