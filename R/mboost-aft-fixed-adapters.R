sb_mboost_aft_learner_ids <- function() {
  c(
    "gamboost_weibull_aft",
    "gamboost_loglogistic_aft",
    "gamboost_lognormal_aft",
    "glmboost_weibull_aft",
    "glmboost_loglogistic_aft",
    "glmboost_lognormal_aft",
    "blackboost_weibull_aft",
    "blackboost_loglogistic_aft",
    "blackboost_lognormal_aft"
  )
}

sb_mboost_aft_params <- function(learner) {
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_mboost_aft_learner_ids()) {
    sb_abort(
      "The requested fixed mboost AFT learner is not supported.",
      "stormbreaker_unsupported_learner"
    )
  }

  distribution <- if (grepl("_weibull_", learner, fixed = TRUE)) {
    list(family = "Weibull", distribution = "weibull")
  } else if (grepl("_loglogistic_", learner, fixed = TRUE)) {
    list(family = "Loglog", distribution = "loglogistic")
  } else {
    list(family = "Lognormal", distribution = "lognormal")
  }
  common <- c(
    distribution,
    list(
      mstop = 100L,
      nu = 0.1,
      risk = "inbag",
      nuisance_range = c(0.1, 100)
    )
  )
  if (startsWith(learner, "gamboost_")) {
    return(c(common, list(
      algorithm = "gamboost",
      base_learner = "bbs",
      dfbase = 4
    )))
  }
  if (startsWith(learner, "glmboost_")) {
    return(c(common, list(
      algorithm = "glmboost",
      base_learner = "bolscw",
      center = TRUE
    )))
  }
  c(common, list(
    algorithm = "blackboost",
    base_learner = "btree",
    tree_controls = list(
      teststat = "quad",
      testtype = "Teststatistic",
      mincriterion = 0,
      minsplit = 10L,
      minbucket = 4L,
      maxdepth = 2L,
      saveinfo = FALSE
    ),
    seed = 1L
  ))
}

sb_mboost_aft_parameter_semantics <- function(learner) {
  params <- sb_mboost_aft_params(learner)
  paste0(
    params$algorithm, " ", params$distribution,
    " AFT likelihood, 100 boosting iterations, shrinkage 0.1, ",
    "positive scale-search interval [0.1,100], and ",
    params$base_learner, " base learner fixed before resampling"
  )
}

sb_mboost_aft_dependency <- function(learner) {
  params <- sb_mboost_aft_params(learner)
  packages <- if (identical(params$algorithm, "blackboost")) {
    c("mboost", "partykit")
  } else {
    "mboost"
  }
  list(
    packages = packages,
    available = all(vapply(
      packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    ))
  )
}

sb_require_mboost_aft_dependency <- function(learner) {
  dependency <- sb_mboost_aft_dependency(learner)
  if (!dependency$available) {
    missing <- dependency$packages[!vapply(
      dependency$packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )]
    sb_abort(
      paste0(
        "Package", if (length(missing) > 1L) "s" else "", " `",
        paste(missing, collapse = "`, `"), "` ",
        if (length(missing) > 1L) "are" else "is",
        " required for learner `", learner, "`."
      ),
      "stormbreaker_missing_mboost_aft_dependency"
    )
  }
  invisible(dependency)
}

sb_mboost_aft_internal_feature_names <- function(n_features) {
  paste0(".stormbreaker_mboost_aft_x", seq_len(n_features))
}

sb_mboost_aft_formula_data <- function(x, time, event) {
  engine_features <- sb_mboost_aft_internal_feature_names(ncol(x))
  engine_x <- x
  names(engine_x) <- engine_features
  formula_data <- sb_additional_survival_formula_data(
    engine_x,
    time,
    event
  )
  c(formula_data, list(engine_features = engine_features))
}

sb_mboost_aft_family <- function(params) {
  switch(params$family,
    Weibull = mboost::Weibull(params$nuisance_range),
    Loglog = mboost::Loglog(params$nuisance_range),
    Lognormal = mboost::Lognormal(params$nuisance_range)
  )
}

sb_fit_mboost_aft_engine <- function(formula_data, params) {
  if (identical(params$algorithm, "gamboost")) {
    withr::local_options(mboost_useMatrix = FALSE)
  }
  control <- mboost::boost_control(
    mstop = params$mstop,
    nu = params$nu,
    risk = params$risk,
    center = TRUE,
    trace = FALSE
  )
  fit_engine <- function() {
    family <- sb_mboost_aft_family(params)
    switch(params$algorithm,
      gamboost = mboost::gamboost(
        formula_data$formula,
        data = formula_data$data,
        family = family,
        control = control,
        baselearner = mboost::bbs,
        dfbase = params$dfbase
      ),
      glmboost = mboost::glmboost(
        formula_data$formula,
        data = formula_data$data,
        family = family,
        control = control,
        center = params$center
      ),
      blackboost = mboost::blackboost(
        formula_data$formula,
        data = formula_data$data,
        family = family,
        control = control,
        tree_controls = do.call(
          partykit::ctree_control,
          params$tree_controls
        )
      )
    )
  }
  fit <- function() {
    sb_classical_survival_stop_on_warning(
      fit_engine(),
      paste0("fixed ", params$algorithm, " ", params$distribution, " AFT learner"),
      "stormbreaker_mboost_aft_fit_warning"
    )
  }
  if ("seed" %in% names(params)) {
    return(withr::with_seed(params$seed, fit()))
  }
  fit()
}

sb_new_mboost_aft_model <- function(
  learner, features, engine_features, params, nuisance_scale, engine
) {
  structure(
    list(
      learner = learner,
      features = features,
      engine_features = engine_features,
      params = params,
      nuisance_scale = nuisance_scale,
      engine = engine
    ),
    class = c(
      paste0("stormbreaker_", learner),
      "stormbreaker_mboost_aft_fixed"
    )
  )
}

sb_fit_mboost_aft <- function(x, time, event, learner, params) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_mboost_aft_dependency(learner)
  formula_data <- sb_mboost_aft_formula_data(
    inputs$x,
    inputs$time,
    inputs$event
  )
  engine <- sb_fit_mboost_aft_engine(formula_data, params)
  nuisance_scale <- as.numeric(mboost::nuisance(engine))
  valid <- inherits(engine, params$algorithm) &&
    identical(as.integer(mboost::mstop(engine)), params$mstop) &&
    length(nuisance_scale) == 1L && !is.na(nuisance_scale) &&
    is.finite(nuisance_scale) &&
    nuisance_scale >= params$nuisance_range[[1L]] &&
    nuisance_scale <= params$nuisance_range[[2L]]
  if (!valid) {
    sb_abort(
      paste0("Learner `", learner, "` did not return its requested AFT fit."),
      "stormbreaker_incomplete_mboost_aft_fit"
    )
  }

  model <- sb_new_mboost_aft_model(
    learner,
    names(inputs$x),
    formula_data$engine_features,
    params,
    nuisance_scale,
    engine
  )
  sb_predict_mboost_aft(model, inputs$x)
  model
}

sb_validate_mboost_aft_model_x <- function(model, x) {
  valid_model <- inherits(model, "stormbreaker_mboost_aft_fixed") &&
    is.character(model$learner) && length(model$learner) == 1L &&
    !is.na(model$learner) &&
    model$learner %in% sb_mboost_aft_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.character(model$engine_features) &&
    length(model$engine_features) == length(model$features) &&
    !anyNA(model$engine_features) && !anyDuplicated(model$engine_features) &&
    is.numeric(model$nuisance_scale) &&
    length(model$nuisance_scale) == 1L &&
    !is.na(model$nuisance_scale) && is.finite(model$nuisance_scale) &&
    model$nuisance_scale > 0
  if (!valid_model) {
    sb_abort(
      "The supplied model is not a fitted mboost AFT learner.",
      "stormbreaker_bad_mboost_aft_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted survival learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  sb_require_mboost_aft_dependency(model$learner)
  x
}

sb_mboost_aft_engine_x <- function(model, x) {
  engine_x <- x
  names(engine_x) <- model$engine_features
  engine_x
}

sb_mboost_aft_link <- function(model, x) {
  x <- sb_validate_mboost_aft_model_x(model, x)
  link <- withCallingHandlers(
    as.numeric(stats::predict(
      model$engine,
      newdata = sb_mboost_aft_engine_x(model, x),
      type = "link"
    )),
    warning = function(warning) {
      spline_extrapolation <- identical(
        model$params$algorithm,
        "gamboost"
      ) && identical(
        conditionMessage(warning),
        paste0(
          "Some 'x' values are beyond 'boundary.knots'; ",
          "Linear extrapolation used."
        )
      )
      if (spline_extrapolation) {
        invokeRestart("muffleWarning")
      }
      sb_abort(
        paste0(
          "The fixed ", model$learner, " prediction raised a warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_mboost_aft_prediction_warning",
        parent = warning
      )
    }
  )
  if (length(link) != nrow(x) || anyNA(link) || any(!is.finite(link))) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid AFT scores."),
      "stormbreaker_invalid_mboost_aft_prediction"
    )
  }
  unname(link)
}

sb_predict_mboost_aft <- function(model, x) {
  -sb_mboost_aft_link(model, x)
}

sb_mboost_aft_survival_probability <- function(model, x, horizon) {
  link <- sb_mboost_aft_link(model, x)
  transformed <- (log(horizon) - link) / model$nuisance_scale
  probability <- switch(model$params$distribution,
    weibull = exp(-exp(transformed)),
    loglogistic = stats::plogis(-transformed),
    lognormal = stats::pnorm(transformed, lower.tail = FALSE)
  )
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid horizon survival."),
      "stormbreaker_invalid_mboost_aft_prediction"
    )
  }
  unname(probability)
}

sb_mboost_aft_direct_survival_source <- function(learner) {
  params <- sb_mboost_aft_params(learner)
  paste0(
    "mboost_", params$algorithm, "_",
    params$distribution, "_aft_direct"
  )
}

sb_mboost_aft_survival_source <- function(model) {
  sb_mboost_aft_direct_survival_source(model$learner)
}

sb_predict_mboost_aft_survival <- function(
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
      "Analysis columns do not match the fitted survival learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  list(
    prediction_horizon = prediction_horizon,
    source = sb_mboost_aft_survival_source(model),
    ties = "not_applicable",
    probability = sb_mboost_aft_survival_probability(
      model,
      x_assessment,
      prediction_horizon
    )
  )
}

sb_extract_mboost_aft_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_mboost_aft_fixed")) {
    sb_bad_learner_state("The fitted mboost AFT learner state is invalid.")
  }
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_mboost_aft_survival(
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
      engine_features = model$engine_features,
      params = model$params,
      nuisance_scale = model$nuisance_scale,
      risk_reference = paste0(
        "mboost_", model$params$algorithm, "_",
        model$params$distribution, "_aft_native"
      ),
      risk_formula = "negative_mboost_aft_log_time_link",
      model_serialization = "R_xdr_v3_mboost",
      engine_raw = serialize(
        model$engine,
        connection = NULL,
        ascii = FALSE,
        xdr = TRUE,
        version = 3L
      ),
      training_data_retention = "required_by_mboost_prediction_state",
      survival_probability = survival_probability
    ),
    class = c(
      "stormbreaker_mboost_aft_fixed_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_mboost_aft_model <- function(state, learner) {
  valid <- is.raw(state$engine_raw) && length(state$engine_raw) > 0L &&
    is.list(state$params) &&
    is.character(state$engine_features) &&
    length(state$engine_features) == length(state$features) &&
    !anyNA(state$engine_features) && !anyDuplicated(state$engine_features) &&
    is.numeric(state$nuisance_scale) &&
    length(state$nuisance_scale) == 1L &&
    !is.na(state$nuisance_scale) && is.finite(state$nuisance_scale) &&
    state$nuisance_scale > 0
  if (!valid) {
    sb_bad_learner_state("The saved mboost AFT learner state is invalid.")
  }
  sb_require_mboost_aft_dependency(learner)
  engine <- tryCatch(
    unserialize(state$engine_raw),
    error = function(error) {
      sb_abort(
        paste0("The saved `", learner, "` model could not be loaded."),
        "stormbreaker_invalid_mboost_aft_state",
        parent = error
      )
    }
  )
  if (!inherits(engine, state$params$algorithm) ||
    !sb_pipeline_parameters_equal(as.integer(mboost::mstop(engine)), state$params$mstop)) {
    sb_bad_learner_state("The saved mboost AFT learner engine is invalid.")
  }
  state$params$mstop <- as.integer(state$params$mstop)
  sb_new_mboost_aft_model(
    learner,
    state$features,
    state$engine_features,
    state$params,
    state$nuisance_scale,
    engine
  )
}

sb_predict_mboost_aft_state <- function(state, x) {
  learner <- if (is.list(state)) state$learner else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_mboost_aft_learner_ids()) {
    sb_bad_learner_state("The saved mboost AFT learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_restore_mboost_aft_model(state, learner)
  result <- list(risk = sb_predict_mboost_aft(model, x))
  probability <- state$survival_probability
  if (!is.null(probability)) {
    valid <- is.list(probability) &&
      is.numeric(probability$prediction_horizon) &&
      length(probability$prediction_horizon) == 1L &&
      !is.na(probability$prediction_horizon) &&
      is.finite(probability$prediction_horizon) &&
      probability$prediction_horizon > 0
    if (!valid) {
      sb_bad_learner_state("The saved mboost AFT horizon state is invalid.")
    }
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- sb_mboost_aft_survival_probability(
      model,
      x,
      probability$prediction_horizon
    )
  }
  result
}

sb_mboost_aft_learner_adapter <- function(
  learner, params = sb_mboost_aft_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_mboost_aft(x, time, event, learner, params)
    },
    predict = sb_predict_mboost_aft,
    predict_survival = sb_predict_mboost_aft_survival,
    extract_state = sb_extract_mboost_aft_state,
    predict_state = sb_predict_mboost_aft_state
  )
}
