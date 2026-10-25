sb_mboost_ranking_learner_ids <- function() {
  c(
    "gamboost_gehan_aft",
    "gamboost_unoc",
    "glmboost_gehan_aft",
    "glmboost_unoc",
    "blackboost_gehan_aft",
    "blackboost_unoc"
  )
}

sb_mboost_ranking_params <- function(learner) {
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_mboost_ranking_learner_ids()) {
    sb_abort(
      "The requested fixed mboost ranking learner is not supported.",
      "stormbreaker_unsupported_learner"
    )
  }

  objective <- if (grepl("_gehan_", learner, fixed = TRUE)) {
    list(
      family = "Gehan",
      objective = "gehan_rank_aft",
      probability_method = "kaplan_meier_log_time_residual"
    )
  } else {
    list(
      family = "Cindex",
      objective = "uno_ipcw_concordance",
      sigma = 0.1,
      ipcw = 1,
      probability_method = "analysis_cox_breslow_calibration",
      calibration_ties = "breslow"
    )
  }
  common <- c(
    objective,
    list(
      mstop = 100L,
      nu = 0.1,
      risk = "inbag"
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

sb_mboost_ranking_parameter_semantics <- function(learner) {
  params <- sb_mboost_ranking_params(learner)
  objective <- if (identical(params$family, "Gehan")) {
    paste0(
      "Gehan rank-based AFT loss with analysis-fold Kaplan-Meier ",
      "log-time residual survival"
    )
  } else {
    paste0(
      "Uno IPCW concordance loss with sigma 0.1 and fold-local IPCW, ",
      "plus analysis-fold Cox/Breslow horizon calibration"
    )
  }
  paste0(
    params$algorithm, " ", objective,
    ", 100 boosting iterations, shrinkage 0.1, and ",
    params$base_learner, " base learner fixed before resampling"
  )
}

sb_mboost_ranking_dependency <- function(learner) {
  params <- sb_mboost_ranking_params(learner)
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

sb_require_mboost_ranking_dependency <- function(learner) {
  dependency <- sb_mboost_ranking_dependency(learner)
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
      "stormbreaker_missing_mboost_ranking_dependency"
    )
  }
  invisible(dependency)
}

sb_mboost_ranking_internal_feature_names <- function(n_features) {
  paste0(".stormbreaker_mboost_ranking_x", seq_len(n_features))
}

sb_mboost_ranking_formula_data <- function(x, time, event) {
  engine_features <- sb_mboost_ranking_internal_feature_names(ncol(x))
  engine_x <- x
  names(engine_x) <- engine_features
  formula_data <- sb_additional_survival_formula_data(
    engine_x,
    time,
    event
  )
  c(formula_data, list(engine_features = engine_features))
}

sb_mboost_ranking_family <- function(params) {
  if (identical(params$family, "Gehan")) {
    return(mboost::Gehan())
  }
  mboost::Cindex(sigma = params$sigma, ipcw = params$ipcw)
}

sb_fit_mboost_ranking_engine <- function(formula_data, params) {
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
    family <- sb_mboost_ranking_family(params)
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
      paste0("fixed ", params$algorithm, " ", params$objective, " learner"),
      "stormbreaker_mboost_ranking_fit_warning"
    )
  }
  if ("seed" %in% names(params)) {
    return(withr::with_seed(params$seed, fit()))
  }
  fit()
}

sb_new_mboost_ranking_model <- function(
  learner, features, engine_features, params, engine
) {
  structure(
    list(
      learner = learner,
      features = features,
      engine_features = engine_features,
      params = params,
      engine = engine
    ),
    class = c(
      paste0("stormbreaker_", learner),
      "stormbreaker_mboost_ranking_fixed"
    )
  )
}

sb_fit_mboost_ranking <- function(x, time, event, learner, params) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_mboost_ranking_dependency(learner)
  formula_data <- sb_mboost_ranking_formula_data(
    inputs$x,
    inputs$time,
    inputs$event
  )
  engine <- sb_fit_mboost_ranking_engine(formula_data, params)
  valid <- inherits(engine, params$algorithm) &&
    identical(as.integer(mboost::mstop(engine)), params$mstop)
  if (!valid) {
    sb_abort(
      paste0("Learner `", learner, "` did not return its requested fit."),
      "stormbreaker_incomplete_mboost_ranking_fit"
    )
  }

  model <- sb_new_mboost_ranking_model(
    learner,
    names(inputs$x),
    formula_data$engine_features,
    params,
    engine
  )
  sb_predict_mboost_ranking(model, inputs$x)
  model
}

sb_validate_mboost_ranking_model_x <- function(model, x) {
  valid_model <- inherits(model, "stormbreaker_mboost_ranking_fixed") &&
    is.character(model$learner) && length(model$learner) == 1L &&
    !is.na(model$learner) &&
    model$learner %in% sb_mboost_ranking_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.character(model$engine_features) &&
    length(model$engine_features) == length(model$features) &&
    !anyNA(model$engine_features) && !anyDuplicated(model$engine_features)
  if (!valid_model) {
    sb_abort(
      "The supplied model is not a fitted mboost ranking learner.",
      "stormbreaker_bad_mboost_ranking_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted survival learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  sb_require_mboost_ranking_dependency(model$learner)
  x
}

sb_mboost_ranking_engine_x <- function(model, x) {
  engine_x <- x
  names(engine_x) <- model$engine_features
  engine_x
}

sb_mboost_ranking_link <- function(model, x) {
  x <- sb_validate_mboost_ranking_model_x(model, x)
  link <- withCallingHandlers(
    as.numeric(stats::predict(
      model$engine,
      newdata = sb_mboost_ranking_engine_x(model, x),
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
        "stormbreaker_mboost_ranking_prediction_warning",
        parent = warning
      )
    }
  )
  if (length(link) != nrow(x) || anyNA(link) || any(!is.finite(link))) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid rank scores."),
      "stormbreaker_invalid_mboost_ranking_prediction"
    )
  }
  unname(link)
}

sb_predict_mboost_ranking <- function(model, x) {
  -sb_mboost_ranking_link(model, x)
}

sb_fit_mboost_gehan_residual_survival <- function(time, event, link) {
  fit <- survival::survfit(
    survival::Surv(log(time) - link, event) ~ 1
  )
  curve <- list(
    time = as.numeric(fit$time),
    survival = as.numeric(fit$surv)
  )
  valid <- length(curve$time) > 0L &&
    length(curve$time) == length(curve$survival) &&
    !anyNA(curve$time) && all(is.finite(curve$time)) &&
    !is.unsorted(curve$time) &&
    !anyNA(curve$survival) && all(is.finite(curve$survival)) &&
    all(curve$survival >= 0 & curve$survival <= 1) &&
    all(diff(curve$survival) <= 0)
  if (!valid) {
    sb_abort(
      "The Gehan fit produced an invalid residual survival curve.",
      "stormbreaker_invalid_mboost_gehan_survival"
    )
  }
  curve
}

sb_predict_mboost_gehan_residual_survival <- function(
  curve, link, prediction_horizon
) {
  valid_curve <- is.list(curve) &&
    is.numeric(curve$time) && is.numeric(curve$survival) &&
    length(curve$time) > 0L &&
    length(curve$time) == length(curve$survival) &&
    !anyNA(curve$time) && all(is.finite(curve$time)) &&
    !is.unsorted(curve$time) &&
    !anyNA(curve$survival) && all(is.finite(curve$survival)) &&
    all(curve$survival >= 0 & curve$survival <= 1) &&
    all(diff(curve$survival) <= 0)
  if (!valid_curve) {
    sb_abort(
      "The saved Gehan residual survival curve is invalid.",
      "stormbreaker_invalid_mboost_gehan_survival"
    )
  }
  index <- findInterval(log(prediction_horizon) - link, curve$time)
  probability <- rep(1, length(index))
  use <- index > 0L
  probability[use] <- curve$survival[index[use]]
  if (anyNA(probability) || any(!is.finite(probability)) ||
    any(probability < 0 | probability > 1)) {
    sb_abort(
      "The Gehan fit produced invalid horizon survival probabilities.",
      "stormbreaker_invalid_mboost_gehan_survival"
    )
  }
  unname(probability)
}

sb_fit_mboost_unoc_calibration <- function(time, event, raw_risk, ties) {
  calibration <- sb_classical_survival_stop_on_warning(
    survival::coxph(
      survival::Surv(time, event) ~ raw_risk,
      data = data.frame(time = time, event = event, raw_risk = raw_risk),
      ties = ties,
      singular.ok = FALSE
    ),
    "mboost Uno-C Cox horizon calibration",
    "stormbreaker_mboost_unoc_calibration_warning"
  )
  coefficient <- unname(as.numeric(stats::coef(calibration)))
  if (length(coefficient) != 1L || is.na(coefficient) ||
    !is.finite(coefficient) || coefficient < 0) {
    sb_abort(
      paste0(
        "The mboost Uno-C horizon calibration did not preserve the ",
        "native higher-risk direction."
      ),
      "stormbreaker_invalid_mboost_unoc_calibration"
    )
  }
  coefficient
}

sb_mboost_ranking_direct_survival_source <- function(learner) {
  params <- sb_mboost_ranking_params(learner)
  if (!identical(params$family, "Gehan")) {
    return(NULL)
  }
  paste0("mboost_", params$algorithm, "_gehan_aft_residual_kaplan_meier")
}

sb_predict_mboost_ranking_survival <- function(
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
  analysis_link <- sb_mboost_ranking_link(model, inputs$x)
  assessment_link <- sb_mboost_ranking_link(model, x_assessment)

  if (identical(model$params$family, "Gehan")) {
    curve <- sb_fit_mboost_gehan_residual_survival(
      inputs$time,
      inputs$event,
      analysis_link
    )
    return(list(
      prediction_horizon = prediction_horizon,
      source = sb_mboost_ranking_direct_survival_source(model$learner),
      ties = "not_applicable",
      probability = sb_predict_mboost_gehan_residual_survival(
        curve,
        assessment_link,
        prediction_horizon
      ),
      residual_survival = curve
    ))
  }

  raw_analysis_risk <- -analysis_link
  raw_assessment_risk <- -assessment_link
  coefficient <- sb_fit_mboost_unoc_calibration(
    inputs$time,
    inputs$event,
    raw_analysis_risk,
    model$params$calibration_ties
  )
  calibrated_analysis_risk <- coefficient * raw_analysis_risk
  cumulative_hazard <- sb_breslow_baseline_hazard(
    inputs$time,
    inputs$event,
    calibrated_analysis_risk,
    prediction_horizon
  )
  result <- sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = paste0(model$learner, "_analysis_cox_breslow_calibration"),
    ties = model$params$calibration_ties,
    risk = coefficient * raw_assessment_risk
  )
  result$calibration_coefficient <- coefficient
  result
}

sb_extract_mboost_ranking_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_mboost_ranking_fixed")) {
    sb_bad_learner_state("The fitted mboost ranking learner state is invalid.")
  }

  residual_survival <- NULL
  calibration_coefficient <- NULL
  baseline_hazard <- NULL
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_mboost_ranking_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    if (identical(model$params$family, "Gehan")) {
      residual_survival <- probability$residual_survival
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
    } else {
      calibration_coefficient <- probability$calibration_coefficient
      baseline_hazard <- list(
        prediction_horizon = probability$prediction_horizon,
        time_unit = if ("time_unit" %in% names(context)) {
          context$time_unit
        } else {
          NA_character_
        },
        cumulative_baseline_hazard = probability$cumulative_baseline_hazard,
        baseline_survival = probability$baseline_survival,
        source = probability$source,
        ties = probability$ties
      )
    }
  }

  structure(
    list(
      learner = model$learner,
      features = model$features,
      engine_features = model$engine_features,
      params = model$params,
      risk_reference = paste0(
        "mboost_", model$params$algorithm, "_",
        tolower(model$params$family), "_native"
      ),
      risk_formula = "negative_mboost_ranking_link",
      probability_method = model$params$probability_method,
      model_serialization = "R_xdr_v3_mboost",
      engine_raw = serialize(
        model$engine,
        connection = NULL,
        ascii = FALSE,
        xdr = TRUE,
        version = 3L
      ),
      training_data_retention = "required_by_mboost_prediction_state",
      residual_survival = residual_survival,
      calibration_coefficient = calibration_coefficient,
      baseline_hazard = baseline_hazard,
      survival_probability = survival_probability
    ),
    class = c(
      "stormbreaker_mboost_ranking_fixed_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_mboost_ranking_model <- function(state, learner) {
  valid <- is.raw(state$engine_raw) && length(state$engine_raw) > 0L &&
    is.list(state$params) &&
    is.character(state$engine_features) &&
    length(state$engine_features) == length(state$features) &&
    !anyNA(state$engine_features) && !anyDuplicated(state$engine_features)
  if (!valid) {
    sb_bad_learner_state("The saved mboost ranking learner state is invalid.")
  }
  sb_require_mboost_ranking_dependency(learner)
  engine <- tryCatch(
    unserialize(state$engine_raw),
    error = function(error) {
      sb_abort(
        paste0("The saved `", learner, "` model could not be loaded."),
        "stormbreaker_invalid_mboost_ranking_state",
        parent = error
      )
    }
  )
  if (!inherits(engine, state$params$algorithm) ||
    !sb_pipeline_parameters_equal(as.integer(mboost::mstop(engine)), state$params$mstop)) {
    sb_bad_learner_state("The saved mboost ranking learner engine is invalid.")
  }
  state$params$mstop <- as.integer(state$params$mstop)
  sb_new_mboost_ranking_model(
    learner,
    state$features,
    state$engine_features,
    state$params,
    engine
  )
}

sb_predict_mboost_ranking_state <- function(state, x) {
  learner <- if (is.list(state)) state$learner else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_mboost_ranking_learner_ids()) {
    sb_bad_learner_state("The saved mboost ranking learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_restore_mboost_ranking_model(state, learner)
  link <- sb_mboost_ranking_link(model, x)
  risk <- -link
  result <- list(risk = risk)

  if (identical(model$params$family, "Gehan")) {
    probability <- state$survival_probability
    if (!is.null(probability)) {
      horizon <- probability$prediction_horizon
      valid <- is.list(probability) && is.numeric(horizon) &&
        length(horizon) == 1L && !is.na(horizon) && is.finite(horizon) &&
        horizon > 0
      if (!valid) {
        sb_bad_learner_state("The saved mboost Gehan horizon state is invalid.")
      }
      result$prediction_horizon <- horizon
      result$time_unit <- probability$time_unit
      result$survival_probability <-
        sb_predict_mboost_gehan_residual_survival(
          state$residual_survival,
          link,
          horizon
        )
    }
    return(result)
  }

  baseline <- state$baseline_hazard
  if (!is.null(baseline)) {
    coefficient <- state$calibration_coefficient
    valid <- is.list(baseline) &&
      is.numeric(baseline$prediction_horizon) &&
      length(baseline$prediction_horizon) == 1L &&
      !is.na(baseline$prediction_horizon) &&
      is.finite(baseline$prediction_horizon) &&
      baseline$prediction_horizon > 0 &&
      is.numeric(coefficient) && length(coefficient) == 1L &&
      !is.na(coefficient) && is.finite(coefficient) && coefficient >= 0
    if (!valid) {
      sb_bad_learner_state("The saved mboost Uno-C horizon state is invalid.")
    }
    result$prediction_horizon <- baseline$prediction_horizon
    result$time_unit <- baseline$time_unit
    result$survival_probability <- sb_survival_probability_from_hazard(
      coefficient * risk,
      baseline$cumulative_baseline_hazard
    )
  }
  result
}

sb_mboost_ranking_learner_adapter <- function(
  learner, params = sb_mboost_ranking_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_mboost_ranking(x, time, event, learner, params)
    },
    predict = sb_predict_mboost_ranking,
    predict_survival = sb_predict_mboost_ranking_survival,
    extract_state = sb_extract_mboost_ranking_state,
    predict_state = sb_predict_mboost_ranking_state
  )
}
