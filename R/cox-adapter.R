sb_fit_coxph <- function(x, time, event) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time < 0)) {
    sb_abort(
      "`time` must contain non-negative values.",
      "stormbreaker_bad_adapter_time"
    )
  }
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }
  if (ncol(x) >= nrow(x)) {
    sb_abort(
      paste0(
        "The unpenalized Cox learner cannot safely fit ", ncol(x),
        " predictors with ", nrow(x),
        " patients; use a fold-local selector or glmnet Cox."
      ),
      "stormbreaker_high_dimensional_cox"
    )
  }

  time_name <- ".stormbreaker_time"
  while (time_name %in% names(x)) {
    time_name <- paste0(time_name, "_")
  }
  event_name <- ".stormbreaker_event"
  while (event_name %in% c(names(x), time_name)) {
    event_name <- paste0(event_name, "_")
  }
  analysis <- x
  analysis[[time_name]] <- time
  analysis[[event_name]] <- event
  rhs <- paste(vapply(names(x), sb_quote_name, ""), collapse = " + ")
  response <- paste0(
    "survival::Surv(",
    sb_quote_name(time_name), ", ",
    sb_quote_name(event_name), ")"
  )
  formula <- stats::as.formula(
    paste(response, "~", rhs),
    env = parent.frame()
  )

  model <- withCallingHandlers(
    survival::coxph(
      formula,
      data = analysis,
      ties = "efron",
      singular.ok = FALSE,
      x = TRUE,
      y = TRUE,
      model = TRUE
    ),
    warning = function(warning) {
      sb_abort(
        paste0("The Cox learner raised a warning: ",
               conditionMessage(warning)),
        "stormbreaker_cox_warning",
        parent = warning
      )
    }
  )

  if (NROW(model$y) != nrow(analysis)) {
    sb_abort(
      "The Cox learner dropped analysis patients.",
      "stormbreaker_rows_dropped"
    )
  }
  coefficients <- stats::coef(model)
  if (length(coefficients) < 1L || anyNA(coefficients) ||
      any(!is.finite(coefficients))) {
    sb_abort(
      "The Cox learner produced missing or non-finite coefficients.",
      "stormbreaker_invalid_cox_fit"
    )
  }

  model
}

sb_predict_coxph <- function(model, x) {
  x <- sb_validate_predictor_frame(x, "x")
  expected_features <- attr(model, "stormbreaker_features", exact = TRUE)
  if (!identical(names(x), expected_features)) {
    sb_abort(
      "Prediction columns do not match the fitted Cox learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  risk <- unname(stats::predict(
    model,
    newdata = x,
    type = "lp",
    reference = "zero"
  ))
  if (length(risk) != nrow(x) || anyNA(risk) ||
      any(!is.finite(risk))) {
    sb_abort(
      "The Cox learner produced incomplete or non-finite risk predictions.",
      "stormbreaker_invalid_cox_prediction"
    )
  }
  risk
}

sb_extract_coxph_coefficients <- function(model) {
  features <- attr(model, "stormbreaker_features", exact = TRUE)
  coefficient <- unname(stats::coef(model))
  data.frame(
    feature = features,
    coefficient = coefficient,
    nonzero = coefficient != 0,
    direction = ifelse(
      coefficient > 0,
      "higher_hazard",
      ifelse(coefficient < 0, "lower_hazard", "zero")
    ),
    stringsAsFactors = FALSE
  )
}

sb_predict_coxph_survival <- function(model, x_analysis, time, event,
                                      x_assessment, prediction_horizon) {
  x_analysis <- sb_validate_predictor_frame(x_analysis, "x_analysis")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x_analysis))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x_analysis))
  expected_features <- attr(model, "stormbreaker_features", exact = TRUE)
  if (!identical(names(x_analysis), expected_features)) {
    sb_abort(
      "Analysis columns do not match the fitted Cox learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    time
  )
  risk <- sb_predict_coxph(model, x_assessment)
  baseline <- withCallingHandlers(
    survival::basehaz(model, centered = FALSE),
    warning = function(warning) {
      sb_abort(
        paste0(
          "The Cox baseline-hazard calculation raised a warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_baseline_hazard_warning",
        parent = warning
      )
    }
  )
  eligible <- baseline$time <= prediction_horizon
  cumulative_hazard <- if (any(eligible)) {
    unname(utils::tail(baseline$hazard[eligible], 1L))
  } else {
    0
  }
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = "coxph_basehaz_zero_reference",
    ties = model$method,
    risk = risk
  )
}

sb_coxph_learner_adapter <- function() {
  list(
    id = "coxph",
    fit = function(x, time, event) {
      model <- sb_fit_coxph(x, time, event)
      attr(model, "stormbreaker_features") <- names(x)
      model
    },
    predict = sb_predict_coxph,
    predict_survival = sb_predict_coxph_survival,
    extract_coefficients = sb_extract_coxph_coefficients,
    extract_state = function(model, context = list()) {
      sb_extract_cox_learner_state(
        model,
        learner = "coxph",
        coefficients = sb_extract_coxph_coefficients(model),
        predict_survival = sb_predict_coxph_survival,
        context = context
      )
    },
    predict_state = sb_predict_cox_learner_state
  )
}
