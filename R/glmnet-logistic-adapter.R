sb_glmnet_logistic_learner_adapter <- function(alpha, lambda,
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
    id = "glmnet_logistic",
    params = params,
    fit = function(x, outcome) {
      if (!sb_glmnet_available()) {
        sb_abort(
          "Package `glmnet` is required for the glmnet logistic learner.",
          "stormbreaker_missing_glmnet"
        )
      }
      x <- sb_validate_predictor_frame(x, "x")
      sb_validate_binary_predictor_names(names(x))
      outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
      if (!all(c(0L, 1L) %in% outcome)) {
        sb_abort(
          "A glmnet logistic analysis fold must contain both classes.",
          "stormbreaker_logistic_single_class"
        )
      }

      path_warnings <- character()
      model <- withCallingHandlers(
        glmnet::glmnet(
          x = as.matrix(x),
          y = outcome,
          family = "binomial",
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
              "The glmnet logistic learner raised a warning: ",
              conditionMessage(warning)
            ),
            "stormbreaker_glmnet_logistic_fit_warning",
            parent = warning
          )
        }
      )
      lambda_min <- min(model$lambda)
      lambda_max <- max(model$lambda)
      if (!sb_glmnet_lambda_in_path(params$lambda, model$lambda)) {
        sb_abort(
          paste0(
            "The requested `lambda` is outside the fitted glmnet logistic ",
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
    predict = function(model, x) {
      x <- sb_validate_predictor_frame(x, "x")
      expected_features <- attr(
        model,
        "stormbreaker_features",
        exact = TRUE
      )
      if (!identical(names(x), expected_features)) {
        sb_abort(
          paste0(
            "Prediction columns do not match the fitted glmnet logistic ",
            "learner schema."
          ),
          "stormbreaker_learner_schema_mismatch"
        )
      }
      lambda <- attr(model, "stormbreaker_lambda", exact = TRUE)
      probability <- drop(stats::predict(
        model,
        newx = as.matrix(x),
        s = lambda,
        type = "response"
      ))
      valid <- length(probability) == nrow(x) && !anyNA(probability) &&
        all(is.finite(probability)) &&
        all(probability >= 0 & probability <= 1)
      if (!valid) {
        sb_abort(
          paste0(
            "The glmnet logistic learner produced invalid class-one ",
            "probabilities."
          ),
          "stormbreaker_invalid_glmnet_logistic_prediction"
        )
      }
      unname(probability)
    },
    extract_coefficients = function(model) {
      features <- attr(model, "stormbreaker_features", exact = TRUE)
      lambda <- attr(model, "stormbreaker_lambda", exact = TRUE)
      coefficient_matrix <- as.matrix(stats::coef(model, s = lambda))
      terms <- c("(Intercept)", features)
      coefficient <- unname(coefficient_matrix[terms, 1L])
      data.frame(
        term = terms,
        coefficient = coefficient,
        nonzero = coefficient != 0,
        stringsAsFactors = FALSE
      )
    },
    extract_state = function(model, context = list()) {
      features <- attr(model, "stormbreaker_features", exact = TRUE)
      lambda <- attr(model, "stormbreaker_lambda", exact = TRUE)
      coefficient_matrix <- as.matrix(stats::coef(model, s = lambda))
      terms <- c("(Intercept)", features)
      sb_extract_logistic_learner_state(
        learner = "glmnet_logistic",
        coefficients = data.frame(
          term = terms,
          coefficient = unname(coefficient_matrix[terms, 1L]),
          stringsAsFactors = FALSE
        ),
        params = params
      )
    },
    predict_state = sb_predict_logistic_learner_state
  )
}
