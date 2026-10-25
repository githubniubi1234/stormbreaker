sb_select_glmnet_logistic_lambda <- function(x, outcome,
                                             inner_folds = 5L,
                                             seed = 1L,
                                             alpha = 1,
                                             standardize = TRUE,
                                             selector = NULL,
                                             fitted_selector = NULL,
                                             patient_ids = NULL,
                                             selection_rule = "lambda.min",
                                             measure = "deviance",
                                             group_ids = NULL,
                                             capture_path = FALSE) {
  if (!sb_glmnet_available()) {
    sb_abort(
      "Package `glmnet` is required for inner logistic tuning.",
      "stormbreaker_missing_glmnet"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  inner_folds <- sb_check_positive_integer(
    inner_folds,
    "inner_folds",
    minimum = 3L
  )
  seed <- sb_check_positive_integer(seed, "seed", minimum = 1L)
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      !is.finite(alpha) || alpha < 0 || alpha > 1) {
    sb_abort(
      "`alpha` must be one finite number between 0 and 1.",
      "stormbreaker_bad_glmnet_alpha"
    )
  }
  if (!is.logical(standardize) || length(standardize) != 1L ||
      is.na(standardize)) {
    sb_abort(
      "`standardize` must be `TRUE` or `FALSE`.",
      "stormbreaker_bad_glmnet_standardize"
    )
  }
  selection_rule <- sb_validate_glmnet_lambda_rule(selection_rule)
  measure <- sb_validate_glmnet_tuning_measure(measure, "binary")
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "Inner glmnet logistic tuning requires both outcome classes.",
      "stormbreaker_logistic_single_class"
    )
  }

  fold_id <- sb_make_stratified_fold_id(
    outcome,
    inner_folds,
    seed,
    stratum_label = "outcome",
    group_ids = group_ids
  )
  uses_fold_local_selection <- !is.null(selector) ||
    !is.null(fitted_selector) || !is.null(patient_ids)
  if (uses_fold_local_selection &&
      (is.null(selector) || is.null(fitted_selector) ||
       is.null(patient_ids))) {
    sb_abort(
      paste0(
        "Fold-local glmnet selection requires `selector`, ",
        "`fitted_selector`, and `patient_ids` together."
      ),
      "stormbreaker_bad_glmnet_selection_context"
    )
  }
  x_matrix <- as.matrix(x)
  fold_filter <- NULL
  if (uses_fold_local_selection) {
    rownames(x_matrix) <- as.character(patient_ids)
    fold_filter <- sb_make_glmnet_fold_selector_filter(
      selector,
      fitted_selector,
      patient_ids,
      fold_id,
      task_type = "binary", group_ids = group_ids
    )
  }
  cv_fit <- tryCatch(
    withCallingHandlers(
      if (uses_fold_local_selection) {
        glmnet::cv.glmnet(
          x = x_matrix,
          y = outcome,
          family = "binomial",
          alpha = as.numeric(alpha),
          foldid = fold_id,
          type.measure = measure,
          standardize = standardize,
          exclude = fold_filter$exclude,
          alignment = "lambda",
          parallel = FALSE
        )
      } else {
        glmnet::cv.glmnet(
          x = x_matrix,
          y = outcome,
          family = "binomial",
          alpha = as.numeric(alpha),
          foldid = fold_id,
          type.measure = measure,
          standardize = standardize
        )
      },
      warning = function(warning) {
        sb_abort(
          paste0(
            "Inner glmnet logistic tuning raised a warning: ",
            conditionMessage(warning)
          ),
          "stormbreaker_glmnet_logistic_tuning_warning",
          parent = warning
        )
      }
    ),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) {
        stop(error)
      }
      sb_abort(
        paste0(
          "Inner glmnet logistic tuning failed: ",
          conditionMessage(error)
        ),
        "stormbreaker_glmnet_logistic_tuning_failure",
        parent = error
      )
    }
  )
  inner_selections <- if (uses_fold_local_selection) {
    fold_filter$inner_selections()
  } else {
    list()
  }

  curve <- data.frame(
    lambda = cv_fit$lambda,
    mean_value = cv_fit$cvm,
    standard_error = cv_fit$cvsd,
    nonzero = cv_fit$nzero,
    stringsAsFactors = FALSE
  )
  names(curve)[[2L]] <- if (identical(measure, "deviance")) {
    "mean_deviance"
  } else {
    "mean_measure"
  }
  valid_curve <- nrow(curve) > 0L &&
    all(is.finite(curve$lambda)) && all(curve$lambda > 0) &&
    all(is.finite(cv_fit$cvm)) &&
    all(is.finite(curve$standard_error)) &&
    all(is.finite(curve$nonzero))
  if (!valid_curve) {
    sb_abort(
      "Inner glmnet logistic tuning produced an incomplete error curve.",
      "stormbreaker_invalid_glmnet_logistic_tuning"
    )
  }
  selected_lambda <- unname(cv_fit[[selection_rule]])
  if (length(selected_lambda) != 1L || !is.finite(selected_lambda) ||
      selected_lambda <= 0 || !selected_lambda %in% curve$lambda) {
    sb_abort(
      "Inner glmnet logistic tuning produced an invalid selected lambda.",
      "stormbreaker_invalid_glmnet_logistic_tuning"
    )
  }

  list(
    lambda = selected_lambda,
    selection_rule = selection_rule,
    measure = measure,
    alpha = as.numeric(alpha),
    standardize = standardize,
    inner_folds = inner_folds,
    seed = seed,
    fold_id = fold_id,
    alignment = "lambda",
    inner_selection_method = if (uses_fold_local_selection) {
      "glmnet_exclude_function"
    } else {
      "none"
    },
    inner_selections = inner_selections,
    curve = curve,
    coefficient_path = if (capture_path) list(lambda = cv_fit$glmnet.fit$lambda,
      beta = cv_fit$glmnet.fit$beta) else NULL,
    warnings = character()
  )
}
