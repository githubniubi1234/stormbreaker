sb_glmnet_available <- function() {
  requireNamespace("glmnet", quietly = TRUE)
}

sb_glmnet_lambda_in_path <- function(lambda, path) {
  limits <- range(path)
  # Screening can change column order and round an otherwise identical path
  # endpoint. Accept roundoff only; retain the requested lambda/native path.
  tolerance <- 16 * .Machine$double.eps
  (limits[[1L]] - lambda) / limits[[1L]] <= tolerance &&
    (lambda - limits[[2L]]) / limits[[2L]] <= tolerance
}

sb_predict_glmnet_cox <- function(model, x) {
  x <- sb_validate_predictor_frame(x, "x")
  expected_features <- attr(
    model,
    "stormbreaker_features",
    exact = TRUE
  )
  if (!identical(names(x), expected_features)) {
    sb_abort(
      paste0(
        "Prediction columns do not match the fitted glmnet Cox ",
        "learner schema."
      ),
      "stormbreaker_learner_schema_mismatch"
    )
  }
  lambda <- attr(model, "stormbreaker_lambda", exact = TRUE)
  risk <- drop(stats::predict(
    model,
    newx = as.matrix(x),
    s = lambda,
    type = "link"
  ))
  if (length(risk) != nrow(x) || anyNA(risk) ||
      any(!is.finite(risk))) {
    sb_abort(
      paste0(
        "The glmnet Cox learner produced incomplete or non-finite ",
        "risk predictions."
      ),
      "stormbreaker_invalid_glmnet_prediction"
    )
  }
  unname(risk)
}

sb_predict_glmnet_cox_survival <- function(model, x_analysis, time, event,
                                           x_assessment,
                                           prediction_horizon) {
  x_analysis <- sb_validate_predictor_frame(x_analysis, "x_analysis")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x_analysis))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x_analysis))
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    time
  )
  analysis_risk <- sb_predict_glmnet_cox(model, x_analysis)
  assessment_risk <- sb_predict_glmnet_cox(model, x_assessment)
  cumulative_hazard <- sb_breslow_baseline_hazard(
    time,
    as.integer(event),
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = "glmnet_fixed_beta_breslow",
    ties = "breslow",
    risk = assessment_risk
  )
}

sb_glmnet_cox_learner_adapter <- function(alpha, lambda,
                                          standardize = TRUE) {
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      !is.finite(alpha) || alpha < 0 || alpha > 1) {
    sb_abort(
      "`alpha` must be one finite number between 0 and 1.",
      "stormbreaker_bad_glmnet_alpha"
    )
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) ||
      !is.finite(lambda) || lambda <= 0) {
    sb_abort(
      "`lambda` must be one finite positive number.",
      "stormbreaker_bad_glmnet_lambda"
    )
  }
  if (!is.logical(standardize) || length(standardize) != 1L ||
      is.na(standardize)) {
    sb_abort(
      "`standardize` must be `TRUE` or `FALSE`.",
      "stormbreaker_bad_glmnet_standardize"
    )
  }

  params <- list(
    alpha = as.numeric(alpha),
    lambda = as.numeric(lambda),
    standardize = standardize
  )
  list(
    id = "glmnet_cox",
    params = params,
    fit = function(x, time, event) {
      if (!sb_glmnet_available()) {
        sb_abort(
          paste0(
            "Package `glmnet` is required for the glmnet Cox learner."
          ),
          "stormbreaker_missing_glmnet"
        )
      }
      x <- sb_validate_predictor_frame(x, "x")
      time <- sb_validate_endpoint_vector(time, "time", nrow(x))
      event <- sb_validate_endpoint_vector(event, "event", nrow(x))
      if (any(time <= 0)) {
        sb_abort(
          paste0(
            "The glmnet Cox learner requires strictly positive ",
            "survival times."
          ),
          "stormbreaker_glmnet_nonpositive_time"
        )
      }
      if (!all(event %in% c(0, 1))) {
        sb_abort(
          "`event` must contain only 0/1 values.",
          "stormbreaker_bad_adapter_event"
        )
      }

      path_warnings <- character()
      model <- withCallingHandlers(
        glmnet::glmnet(
          x = as.matrix(x),
          y = survival::Surv(time, event),
          family = "cox",
          cox.ties = "breslow",
          alpha = params$alpha,
          standardize = params$standardize
        ),
        warning = function(warning) {
          if (sb_is_glmnet_tail_convergence_warning(warning)) {
            path_warnings <<- c(path_warnings, conditionMessage(warning))
            invokeRestart("muffleWarning")
          }
          sb_abort(
            paste0(
              "The glmnet Cox learner raised a warning: ",
              conditionMessage(warning)
            ),
            "stormbreaker_glmnet_fit_warning",
            parent = warning
          )
        }
      )
      lambda_min <- min(model$lambda)
      lambda_max <- max(model$lambda)
      if (!sb_glmnet_lambda_in_path(params$lambda, model$lambda)) {
        sb_abort(
          paste0(
            "The requested `lambda` is outside the fitted glmnet Cox ",
            "path [", format(lambda_min), ", ", format(lambda_max), "]."
          ),
          "stormbreaker_glmnet_lambda_outside_path"
        )
      }
      attr(model, "stormbreaker_features") <- names(x)
      attr(model, "stormbreaker_lambda") <- params$lambda
      attr(model, "stormbreaker_fit_warnings") <- path_warnings
      model
    },
    predict = sb_predict_glmnet_cox,
    predict_survival = sb_predict_glmnet_cox_survival,
    extract_coefficients = function(model) {
      features <- attr(model, "stormbreaker_features", exact = TRUE)
      lambda <- attr(model, "stormbreaker_lambda", exact = TRUE)
      coefficient_matrix <- as.matrix(stats::coef(model, s = lambda))
      coefficient <- unname(coefficient_matrix[features, 1L])
      nonzero <- coefficient != 0
      data.frame(
        feature = features,
        coefficient = coefficient,
        nonzero = nonzero,
        direction = ifelse(
          coefficient > 0,
          "higher_hazard",
          ifelse(coefficient < 0, "lower_hazard", "zero")
        ),
        stringsAsFactors = FALSE
      )
    },
    extract_state = function(model, context = list()) {
      features <- attr(model, "stormbreaker_features", exact = TRUE)
      lambda <- attr(model, "stormbreaker_lambda", exact = TRUE)
      coefficient_matrix <- as.matrix(stats::coef(model, s = lambda))
      coefficient <- unname(coefficient_matrix[features, 1L])
      coefficients <- data.frame(
        feature = features,
        coefficient = coefficient,
        stringsAsFactors = FALSE
      )
      sb_extract_cox_learner_state(
        model,
        learner = "glmnet_cox",
        coefficients = coefficients,
        predict_survival = sb_predict_glmnet_cox_survival,
        context = context,
        params = params
      )
    },
    predict_state = sb_predict_cox_learner_state
  )
}
