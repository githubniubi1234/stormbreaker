sb_mboost_cox_learner_ids <- function() {
  c(
    "gamboost_cox",
    "glmboost_cox",
    "blackboost_cox"
  )
}

sb_mboost_cox_params <- function(learner) {
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_mboost_cox_learner_ids()) {
    sb_abort(
      "The requested fixed mboost Cox learner is not supported.",
      "stormbreaker_unsupported_learner"
    )
  }

  common <- list(
    mstop = 100L,
    nu = 0.1,
    risk = "inbag",
    family = "CoxPH"
  )
  specific <- switch(learner,
    gamboost_cox = list(
      base_learner = "bbs",
      dfbase = 4
    ),
    glmboost_cox = list(
      base_learner = "bolscw",
      center = TRUE
    ),
    blackboost_cox = list(
      base_learner = "btree",
      tree_controls = list(
        teststat = "quad",
        testtype = "Teststatistic",
        mincriterion = 0,
        minsplit = 10L,
        minbucket = 4L,
        maxdepth = 2L,
        saveinfo = FALSE
      )
    )
  )
  c(common, specific)
}

sb_mboost_cox_dependency <- function(learner) {
  packages <- if (identical(learner, "blackboost_cox")) {
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

sb_require_mboost_cox_dependency <- function(learner) {
  dependency <- sb_mboost_cox_dependency(learner)
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
      "stormbreaker_missing_mboost_cox_dependency"
    )
  }
  invisible(dependency)
}

sb_new_mboost_cox_model <- function(learner, features, params, engine) {
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      engine = engine
    ),
    class = c(
      paste0("stormbreaker_", learner),
      "stormbreaker_mboost_cox_fixed"
    )
  )
}

sb_fit_mboost_cox <- function(x, time, event, learner, params) {
  inputs <- sb_validate_classical_survival_inputs(
    x,
    time,
    event,
    learner
  )
  sb_require_mboost_cox_dependency(learner)
  if (identical(learner, "gamboost_cox")) {
    # Native dense spline algebra avoids sparse-QR augmentation warnings;
    # the spline basis, degrees of freedom and warning checks are unchanged.
    withr::local_options(mboost_useMatrix = FALSE)
  }
  formula_data <- sb_additional_survival_formula_data(
    inputs$x,
    inputs$time,
    inputs$event
  )
  control <- mboost::boost_control(
    mstop = params$mstop,
    nu = params$nu,
    risk = params$risk,
    center = TRUE,
    trace = FALSE
  )

  engine <- sb_classical_survival_stop_on_warning(
    switch(learner,
      gamboost_cox = mboost::gamboost(
        formula_data$formula,
        data = formula_data$data,
        family = mboost::CoxPH(),
        control = control,
        baselearner = mboost::bbs,
        dfbase = params$dfbase
      ),
      glmboost_cox = mboost::glmboost(
        formula_data$formula,
        data = formula_data$data,
        family = mboost::CoxPH(),
        control = control,
        center = params$center
      ),
      blackboost_cox = mboost::blackboost(
        formula_data$formula,
        data = formula_data$data,
        family = mboost::CoxPH(),
        control = control,
        tree_controls = do.call(
          partykit::ctree_control,
          params$tree_controls
        )
      )
    ),
    paste0("fixed ", learner, " learner"),
    "stormbreaker_mboost_cox_fit_warning"
  )

  expected_class <- switch(learner,
    gamboost_cox = "gamboost",
    glmboost_cox = "glmboost",
    blackboost_cox = "blackboost"
  )
  if (!inherits(engine, expected_class) ||
    !identical(as.integer(mboost::mstop(engine)), params$mstop)) {
    sb_abort(
      paste0("Learner `", learner, "` did not return its requested fit."),
      "stormbreaker_incomplete_mboost_cox_fit"
    )
  }

  model <- sb_new_mboost_cox_model(
    learner,
    names(inputs$x),
    params,
    engine
  )
  sb_predict_mboost_cox(model, inputs$x)
  model
}

sb_validate_mboost_cox_model_x <- function(model, x) {
  valid_model <- inherits(model, "stormbreaker_mboost_cox_fixed") &&
    is.character(model$learner) && length(model$learner) == 1L &&
    !is.na(model$learner) &&
    model$learner %in% sb_mboost_cox_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features)
  if (!valid_model) {
    sb_abort(
      "The supplied model is not a fitted mboost Cox learner.",
      "stormbreaker_bad_mboost_cox_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted survival learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  sb_require_mboost_cox_dependency(model$learner)
  x
}

sb_predict_mboost_cox <- function(model, x) {
  x <- sb_validate_mboost_cox_model_x(model, x)
  risk <- withCallingHandlers(
    as.numeric(stats::predict(
      model$engine,
      newdata = x,
      type = "link"
    )),
    warning = function(warning) {
      spline_extrapolation <- identical(
        model$learner,
        "gamboost_cox"
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
        "stormbreaker_mboost_cox_prediction_warning",
        parent = warning
      )
    }
  )
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid risk scores."),
      "stormbreaker_invalid_mboost_cox_prediction"
    )
  }
  unname(risk)
}

sb_predict_mboost_cox_survival <- function(
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
  analysis_risk <- sb_predict_mboost_cox(model, inputs$x)
  assessment_risk <- sb_predict_mboost_cox(model, x_assessment)
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

sb_extract_mboost_cox_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_mboost_cox_fixed")) {
    sb_bad_learner_state("The fitted mboost Cox learner state is invalid.")
  }

  baseline_hazard <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_mboost_cox_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
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

  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      risk_reference = "mboost_native",
      risk_formula = "mboost_coxph_link_prediction",
      model_serialization = "R_xdr_v3_mboost",
      engine_raw = serialize(
        model$engine,
        connection = NULL,
        ascii = FALSE,
        xdr = TRUE,
        version = 3L
      ),
      training_data_retention = "required_by_mboost_prediction_state",
      baseline_hazard = baseline_hazard
    ),
    class = c(
      "stormbreaker_mboost_cox_fixed_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_mboost_cox_model <- function(state, learner) {
  if (!is.raw(state$engine_raw) || length(state$engine_raw) == 0L ||
    !is.list(state$params)) {
    sb_bad_learner_state("The saved mboost Cox learner state is invalid.")
  }
  sb_require_mboost_cox_dependency(learner)
  engine <- tryCatch(
    unserialize(state$engine_raw),
    error = function(error) {
      sb_abort(
        paste0("The saved `", learner, "` model could not be loaded."),
        "stormbreaker_invalid_mboost_cox_state",
        parent = error
      )
    }
  )
  expected_class <- switch(learner,
    gamboost_cox = "gamboost",
    glmboost_cox = "glmboost",
    blackboost_cox = "blackboost"
  )
  if (!inherits(engine, expected_class)) {
    sb_bad_learner_state("The saved mboost Cox learner engine is invalid.")
  }
  sb_new_mboost_cox_model(
    learner,
    state$features,
    state$params,
    engine
  )
}

sb_predict_mboost_cox_state <- function(state, x) {
  learner <- if (is.list(state)) state$learner else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_mboost_cox_learner_ids()) {
    sb_bad_learner_state("The saved mboost Cox learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_restore_mboost_cox_model(state, learner)
  risk <- sb_predict_mboost_cox(model, x)
  result <- list(risk = risk)

  baseline <- state$baseline_hazard
  if (!is.null(baseline)) {
    valid_baseline <- is.list(baseline) &&
      is.numeric(baseline$prediction_horizon) &&
      length(baseline$prediction_horizon) == 1L &&
      !is.na(baseline$prediction_horizon) &&
      is.finite(baseline$prediction_horizon) &&
      baseline$prediction_horizon > 0
    if (!valid_baseline) {
      sb_bad_learner_state("The saved mboost Cox baseline state is invalid.")
    }
    result$prediction_horizon <- baseline$prediction_horizon
    result$time_unit <- baseline$time_unit
    result$survival_probability <- sb_survival_probability_from_hazard(
      risk,
      baseline$cumulative_baseline_hazard
    )
  }
  result
}

sb_mboost_cox_learner_adapter <- function(
  learner, params = sb_mboost_cox_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_mboost_cox(x, time, event, learner, params)
    },
    predict = sb_predict_mboost_cox,
    predict_survival = sb_predict_mboost_cox_survival,
    extract_state = sb_extract_mboost_cox_state,
    predict_state = sb_predict_mboost_cox_state
  )
}
