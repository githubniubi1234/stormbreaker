sb_timereg_survival_learner_ids <- function() {
  c("aalen_additive_survival", "prop_odds_survival")
}

sb_timereg_survival_params <- function(learner) {
  switch(learner,
    aalen_additive_survival = list(
      coefficient_form = "time_varying",
      robust = 0L,
      n_sim = 0L,
      weighted_test = 0L,
      covariance = 0L,
      resample_iid = 0L,
      residuals = 0L,
      seed = 1L,
      risk_time_rule = "last_analysis_event"
    ),
    prop_odds_survival = list(
      iterations = 20L,
      n_sim = 0L,
      weighted_test = 0L,
      profile = 1L,
      sym = 0L,
      baseline_variance = 1L,
      seed = 1L
    ),
    sb_abort(
      paste0("Unknown timereg survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_timereg_available <- function() {
  requireNamespace("timereg", quietly = TRUE)
}

sb_require_timereg <- function(learner) {
  if (!sb_timereg_available()) {
    sb_abort(
      paste0("Package `timereg` is required for learner `", learner, "`."),
      "stormbreaker_missing_timereg_dependency"
    )
  }
  invisible(TRUE)
}

sb_timereg_formula_data <- function(x, time, event, learner) {
  analysis <- x
  time_name <- ".stormbreaker_time"
  while (time_name %in% names(analysis)) {
    time_name <- paste0(time_name, "_")
  }
  event_name <- ".stormbreaker_event"
  while (event_name %in% c(names(analysis), time_name)) {
    event_name <- paste0(event_name, "_")
  }
  analysis[[time_name]] <- time
  analysis[[event_name]] <- event

  response <- if (identical(learner, "aalen_additive_survival")) {
    paste0(
      "survival::Surv(", sb_quote_name(time_name), ", ",
      sb_quote_name(event_name), ")"
    )
  } else {
    paste0(
      "timereg::Event(", sb_quote_name(time_name), ", ",
      sb_quote_name(event_name), ")"
    )
  }
  rhs <- paste(vapply(names(x), sb_quote_name, ""), collapse = " + ")
  list(
    data = analysis,
    formula = stats::as.formula(
      paste(response, "~", rhs),
      env = baseenv()
    )
  )
}

sb_new_timereg_survival_model <- function(learner, features, params, time,
                                           cumulative_intercept = NULL,
                                           cumulative_coefficient = NULL,
                                           risk_reference_time = NULL,
                                           baseline_cumulative_odds = NULL,
                                           coefficient = NULL) {
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      time = time,
      cumulative_intercept = cumulative_intercept,
      cumulative_coefficient = cumulative_coefficient,
      risk_reference_time = risk_reference_time,
      baseline_cumulative_odds = baseline_cumulative_odds,
      coefficient = coefficient
    ),
    class = "stormbreaker_timereg_survival_fixed"
  )
}

sb_fit_timereg_survival <- function(x, time, event, learner, params) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_timereg(learner)
  formula_data <- sb_timereg_formula_data(
    inputs$x,
    inputs$time,
    inputs$event,
    learner
  )

  engine <- if (identical(learner, "aalen_additive_survival")) {
    withr::with_seed(
      params$seed,
      timereg::aalen(
        formula = formula_data$formula,
        data = formula_data$data,
        robust = params$robust,
        residuals = params$residuals,
        n.sim = params$n_sim,
        weighted.test = params$weighted_test,
        covariance = params$covariance,
        resample.iid = params$resample_iid,
        silent = 1L
      )
    )
  } else {
    withr::with_seed(
      params$seed,
      timereg::prop.odds(
        formula = formula_data$formula,
        data = formula_data$data,
        Nit = params$iterations,
        detail = 0L,
        n.sim = params$n_sim,
        weighted.test = params$weighted_test,
        profile = params$profile,
        sym = params$sym,
        baselinevar = params$baseline_variance
      )
    )
  }

  cumulative <- as.matrix(engine$cum)
  valid_cumulative <- is.numeric(cumulative) && nrow(cumulative) > 0L &&
    !anyNA(cumulative) && all(is.finite(cumulative)) &&
    all(diff(cumulative[, 1L]) >= 0)
  if (!valid_cumulative) {
    sb_abort(
      paste0("Learner `", learner, "` did not return a finite time process."),
      "stormbreaker_incomplete_timereg_survival_fit"
    )
  }

  if (identical(learner, "aalen_additive_survival")) {
    valid <- ncol(cumulative) == ncol(inputs$x) + 2L
    if (!valid) {
      sb_abort(
        "Aalen fitting did not retain one cumulative coefficient per feature.",
        "stormbreaker_incomplete_timereg_survival_fit"
      )
    }
    coefficient <- cumulative[, -(1:2), drop = FALSE]
    colnames(coefficient) <- names(inputs$x)
    model <- sb_new_timereg_survival_model(
      learner = learner,
      features = names(inputs$x),
      params = params,
      time = unname(cumulative[, 1L]),
      cumulative_intercept = unname(cumulative[, 2L]),
      cumulative_coefficient = unname(coefficient),
      risk_reference_time = unname(cumulative[nrow(cumulative), 1L])
    )
  } else {
    gamma <- as.matrix(engine$gamma)
    valid <- ncol(cumulative) == 2L && nrow(gamma) == ncol(inputs$x) &&
      ncol(gamma) >= 1L && is.numeric(gamma) && !anyNA(gamma[, 1L]) &&
      all(is.finite(gamma[, 1L]))
    if (!valid) {
      sb_abort(
        paste0(
          "Proportional-odds fitting did not retain one finite coefficient ",
          "per feature."
        ),
        "stormbreaker_incomplete_timereg_survival_fit"
      )
    }
    model <- sb_new_timereg_survival_model(
      learner = learner,
      features = names(inputs$x),
      params = params,
      time = unname(cumulative[, 1L]),
      baseline_cumulative_odds = unname(cumulative[, 2L]),
      coefficient = unname(gamma[, 1L])
    )
  }

  sb_predict_timereg_survival(model, inputs$x)
  model
}

sb_validate_timereg_survival_model_x <- function(model, x) {
  valid <- inherits(model, "stormbreaker_timereg_survival_fixed") &&
    is.list(model) && model$learner %in% sb_timereg_survival_learner_ids()
  if (!isTRUE(valid)) {
    sb_abort(
      "The supplied model is not a fitted timereg survival learner.",
      "stormbreaker_bad_timereg_survival_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted timereg learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  x
}

sb_timereg_process_index <- function(time, horizon) {
  eligible <- which(time <= horizon)
  if (length(eligible) == 0L) 0L else max(eligible)
}

sb_predict_timereg_survival <- function(model, x) {
  x <- sb_validate_timereg_survival_model_x(model, x)
  risk <- if (identical(model$learner, "aalen_additive_survival")) {
    row <- length(model$time)
    model$cumulative_intercept[[row]] +
      as.vector(as.matrix(x) %*% model$cumulative_coefficient[row, ])
  } else {
    as.vector(as.matrix(x) %*% model$coefficient)
  }
  risk <- unname(as.numeric(risk))
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid risk scores."),
      "stormbreaker_invalid_timereg_survival_prediction"
    )
  }
  risk
}

sb_timereg_direct_survival_source <- function(learner) {
  switch(learner,
    aalen_additive_survival =
      "timereg_aalen_native_additive_survival_direct",
    prop_odds_survival =
      "timereg_proportional_odds_survival_direct",
    NULL
  )
}

sb_timereg_direct_survival_probability <- function(model, x, horizon) {
  x <- sb_validate_timereg_survival_model_x(model, x)
  valid_horizon <- is.numeric(horizon) && length(horizon) == 1L &&
    !is.na(horizon) && is.finite(horizon) && horizon > 0
  if (!valid_horizon) {
    sb_abort(
      "`horizon` must be one finite positive number.",
      "stormbreaker_bad_prediction_horizon"
    )
  }
  index <- sb_timereg_process_index(model$time, horizon)

  probability <- if (identical(
    model$learner,
    "aalen_additive_survival"
  )) {
    cumulative_hazard <- if (index == 0L) {
      rep(0, nrow(x))
    } else {
      model$cumulative_intercept[[index]] +
        as.vector(
          as.matrix(x) %*% model$cumulative_coefficient[index, ]
        )
    }
    exp(-cumulative_hazard)
  } else {
    baseline_odds <- if (index == 0L) {
      0
    } else {
      model$baseline_cumulative_odds[[index]]
    }
    if (!is.finite(baseline_odds) || baseline_odds < 0) {
      rep(NA_real_, nrow(x))
    } else if (baseline_odds == 0) {
      rep(1, nrow(x))
    } else {
      stats::plogis(-(log(baseline_odds) +
        sb_predict_timereg_survival(model, x)))
    }
  }

  probability <- unname(as.numeric(probability))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      paste0(
        "Learner `", model$learner,
        "` produced an invalid native horizon survival probability; ",
        "the estimate was not clipped or recalibrated."
      ),
      "stormbreaker_invalid_timereg_survival_prediction"
    )
  }
  probability
}

sb_predict_timereg_survival_survival <- function(
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
      "Analysis columns do not match the fitted timereg learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  list(
    prediction_horizon = prediction_horizon,
    source = sb_timereg_direct_survival_source(model$learner),
    ties = "not_applicable",
    probability = sb_timereg_direct_survival_probability(
      model,
      x_assessment,
      prediction_horizon
    )
  )
}

sb_extract_timereg_survival_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_timereg_survival_fixed")) {
    sb_bad_learner_state("The fitted timereg survival learner is invalid.")
  }
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_timereg_survival_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x,
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

  state <- unclass(model)
  state$survival_probability <- survival_probability
  state$risk_reference <- if (identical(
    model$learner,
    "aalen_additive_survival"
  )) {
    "timereg_native_last_analysis_event_time"
  } else {
    "zero"
  }
  state$risk_formula <- if (identical(
    model$learner,
    "aalen_additive_survival"
  )) {
    "cumulative_hazard_at_last_analysis_event"
  } else {
    "x_beta"
  }
  structure(
    state,
    class = c(
      "stormbreaker_timereg_survival_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_timereg_survival_model <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_timereg_survival_learner_ids()) {
    sb_bad_learner_state("The saved timereg survival state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  valid_time <- is.numeric(state$time) && length(state$time) > 0L &&
    !anyNA(state$time) && all(is.finite(state$time)) &&
    all(diff(state$time) >= 0)
  if (!valid_time) {
    sb_bad_learner_state("The saved timereg time process is invalid.")
  }

  if (identical(learner, "aalen_additive_survival")) {
    valid <- is.numeric(state$cumulative_intercept) &&
      length(state$cumulative_intercept) == length(state$time) &&
      !anyNA(state$cumulative_intercept) &&
      all(is.finite(state$cumulative_intercept)) &&
      is.matrix(state$cumulative_coefficient) &&
      nrow(state$cumulative_coefficient) == length(state$time) &&
      ncol(state$cumulative_coefficient) == ncol(x) &&
      !anyNA(state$cumulative_coefficient) &&
      all(is.finite(state$cumulative_coefficient)) &&
      is.numeric(state$risk_reference_time) &&
      length(state$risk_reference_time) == 1L &&
      identical(
        state$risk_reference_time,
        state$time[[length(state$time)]]
      )
  } else {
    valid <- is.numeric(state$baseline_cumulative_odds) &&
      length(state$baseline_cumulative_odds) == length(state$time) &&
      !anyNA(state$baseline_cumulative_odds) &&
      all(is.finite(state$baseline_cumulative_odds)) &&
      is.numeric(state$coefficient) && length(state$coefficient) == ncol(x) &&
      !anyNA(state$coefficient) && all(is.finite(state$coefficient))
  }
  if (!isTRUE(valid)) {
    sb_bad_learner_state("The saved timereg prediction parameters are invalid.")
  }

  model <- sb_new_timereg_survival_model(
    learner = learner,
    features = state$features,
    params = state$params,
    time = state$time,
    cumulative_intercept = state$cumulative_intercept,
    cumulative_coefficient = state$cumulative_coefficient,
    risk_reference_time = state$risk_reference_time,
    baseline_cumulative_odds = state$baseline_cumulative_odds,
    coefficient = state$coefficient
  )
  list(model = model, x = x)
}

sb_predict_timereg_survival_state <- function(state, x) {
  restored <- sb_restore_timereg_survival_model(state, x)
  result <- list(
    risk = sb_predict_timereg_survival(restored$model, restored$x)
  )
  probability <- state$survival_probability
  if (!is.null(probability)) {
    valid <- is.list(probability) &&
      is.numeric(probability$prediction_horizon) &&
      length(probability$prediction_horizon) == 1L &&
      !is.na(probability$prediction_horizon) &&
      is.finite(probability$prediction_horizon) &&
      probability$prediction_horizon > 0 &&
      is.character(probability$time_unit) &&
      length(probability$time_unit) == 1L && !is.na(probability$time_unit) &&
      nzchar(probability$time_unit) &&
      identical(
        probability$source,
        sb_timereg_direct_survival_source(state$learner)
      ) && identical(probability$ties, "not_applicable")
    if (!valid) {
      sb_bad_learner_state(
        "The saved timereg survival-probability state is invalid."
      )
    }
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- sb_timereg_direct_survival_probability(
      restored$model,
      restored$x,
      probability$prediction_horizon
    )
  }
  result
}

sb_timereg_survival_learner_adapter <- function(
  learner,
  params = sb_timereg_survival_params(learner)
) {
  if (!learner %in% sb_timereg_survival_learner_ids()) {
    sb_abort(
      paste0("Unknown timereg survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_timereg_survival(x, time, event, learner, params)
    },
    predict = sb_predict_timereg_survival,
    predict_survival = sb_predict_timereg_survival_survival,
    extract_state = sb_extract_timereg_survival_state,
    predict_state = sb_predict_timereg_survival_state
  )
}
