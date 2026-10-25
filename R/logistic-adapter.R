sb_validate_binary_predictor_names <- function(predictors) {
  if ("(Intercept)" %in% predictors) {
    sb_abort(
      paste0(
        "Predictor name `(Intercept)` is reserved for the model intercept ",
        "in binary pipelines."
      ),
      "stormbreaker_reserved_intercept_name"
    )
  }
  invisible(predictors)
}

sb_validate_binary_outcome_vector <- function(outcome, n) {
  valid <- (is.numeric(outcome) || is.logical(outcome)) &&
    length(outcome) == n && !anyNA(outcome) &&
    (!is.numeric(outcome) || all(is.finite(outcome))) &&
    all(outcome %in% c(0, 1))
  if (!valid) {
    sb_abort(
      "`outcome` must contain one finite 0/1 value per analysis patient.",
      "stormbreaker_bad_adapter_binary_outcome"
    )
  }
  as.integer(outcome)
}

sb_fit_logistic <- function(x, outcome) {
  x <- sb_validate_predictor_frame(x, "x")
  sb_validate_binary_predictor_names(names(x))
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "A logistic analysis fold must contain both outcome classes.",
      "stormbreaker_logistic_single_class"
    )
  }

  outcome_name <- ".stormbreaker_outcome"
  while (outcome_name %in% names(x)) {
    outcome_name <- paste0(outcome_name, "_")
  }
  analysis <- x
  analysis[[outcome_name]] <- outcome
  rhs <- paste(vapply(names(x), sb_quote_name, ""), collapse = " + ")
  formula <- stats::as.formula(
    paste(sb_quote_name(outcome_name), "~", rhs),
    env = parent.frame()
  )
  fit_warnings <- character()
  model <- withCallingHandlers(
    stats::glm(
      formula,
      data = analysis,
      family = stats::binomial(link = "logit"),
      x = TRUE,
      y = TRUE,
      model = TRUE
    ),
    warning = function(warning) {
      if (identical(
        conditionMessage(warning),
        "glm.fit: fitted probabilities numerically 0 or 1 occurred"
      )) {
        fit_warnings <<- c(fit_warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
      sb_abort(
        paste0(
          "The logistic learner raised a warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_logistic_warning",
        parent = warning
      )
    }
  )
  coefficients <- stats::coef(model)
  fitted_matrix <- unname(model$x[, -1L, drop = FALSE])
  expected_matrix <- unname(as.matrix(x))
  storage.mode(expected_matrix) <- "double"
  valid_fit <- isTRUE(model$converged) && NROW(model$y) == nrow(analysis) &&
    length(coefficients) == ncol(x) + 1L &&
    model$rank == ncol(x) + 1L &&
    identical(fitted_matrix, expected_matrix) &&
    !anyNA(coefficients) && all(is.finite(coefficients))
  if (!valid_fit) {
    sb_abort(
      paste0(
        "The logistic learner did not produce a converged, full-rank fit ",
        "with one finite coefficient per declared term."
      ),
      "stormbreaker_invalid_logistic_fit"
    )
  }
  if (length(fit_warnings) > 0L) {
    # A saturated probability can coexist with a finite MLE. Reuse the
    # existing post-fit separation calculation rather than ignoring a warning
    # that can also accompany an unbounded fit.
    separation <- sb_binary_moderate_calibration_separation(model$x, outcome)
    if (identical(separation$status, "failed")) {
      sb_abort(
        paste0(
          "The logistic separation calculation failed: ",
          separation$error_message
        ),
        "stormbreaker_logistic_separation_check_failed"
      )
    }
    if (isTRUE(separation$detected)) {
      sb_abort(
        paste0(
          "The logistic design has complete or quasi-complete separation; ",
          "no finite maximum-likelihood fit is reported. Use fewer predictors ",
          "or explicitly choose a regularized learner."
        ),
        "stormbreaker_logistic_separation"
      )
    }
  }
  attr(model, "stormbreaker_features") <- names(x)
  attr(model, "stormbreaker_fit_warnings") <- fit_warnings
  model
}

sb_predict_logistic <- function(model, x) {
  x <- sb_validate_predictor_frame(x, "x")
  expected_features <- attr(model, "stormbreaker_features", exact = TRUE)
  if (!identical(names(x), expected_features)) {
    sb_abort(
      "Prediction columns do not match the fitted logistic learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  probability <- unname(stats::predict(model, newdata = x, type = "response"))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      "The logistic learner produced invalid class-one probabilities.",
      "stormbreaker_invalid_logistic_prediction"
    )
  }
  probability
}

sb_logistic_learner_adapter <- function() {
  list(
    id = "logistic",
    fit = sb_fit_logistic,
    predict = sb_predict_logistic,
    extract_coefficients = function(model) {
      coefficients <- stats::coef(model)
      features <- attr(model, "stormbreaker_features", exact = TRUE)
      data.frame(
        term = c("(Intercept)", features),
        coefficient = unname(coefficients),
        stringsAsFactors = FALSE
      )
    },
    extract_state = function(model, context = list()) {
      coefficients <- stats::coef(model)
      features <- attr(model, "stormbreaker_features", exact = TRUE)
      sb_extract_logistic_learner_state(
        learner = "logistic",
        coefficients = data.frame(
          term = c("(Intercept)", features),
          coefficient = unname(coefficients),
          stringsAsFactors = FALSE
        )
      )
    },
    predict_state = sb_predict_logistic_learner_state
  )
}
