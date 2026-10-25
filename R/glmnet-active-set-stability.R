sb_glmnet_active_set_stability <- function(pipeline, coefficients,
                                            task_predictors) {
  is_glmnet <- pipeline$learner %in% c("glmnet_cox", "glmnet_logistic")
  alpha <- if (is_glmnet) {
    as.numeric(pipeline$learner_params$alpha)
  } else {
    NA_real_
  }
  applicable <- is_glmnet && alpha > 0
  fixed_features <- identical(pipeline$selector, "fixed_features")
  candidate_features <- if (fixed_features) {
    pipeline$selector_params$features
  } else {
    task_predictors
  }
  active_sets <- lapply(coefficients, function(audit) {
    success <- identical(audit$status, "success")
    list(
      fold = as.integer(audit$fold),
      status = if (success) "success" else audit$status,
      selected_features = if (success) {
        as.character(audit$nonzero_features)
      } else {
        character()
      },
      error_class = NA_character_,
      error_message = NA_character_
    )
  })
  summary <- sb_selection_stability(
    if (applicable) "glmnet_active_set" else "none",
    active_sets,
    candidate_features
  )

  all_folds_fitted <- summary$n_folds >= 2L &&
    summary$n_successful_folds == summary$n_folds
  empty_union <- any(summary$pairwise$status == "not_evaluable")
  if (applicable && identical(summary$status, "incomplete") &&
      all_folds_fitted && empty_union) {
    summary$status <- "not_evaluable"
    summary$error_message <- paste0(
      "Glmnet active-set stability is not evaluable because at least one ",
      "outer-fold pair has an empty feature union."
    )
  }

  list(
    learner = pipeline$learner,
    alpha = alpha,
    status = summary$status,
    active_set_definition = if (applicable) {
      "exact_nonzero_predictor_coefficients_at_selected_lambda"
    } else {
      NA_character_
    },
    scope = "outer_analysis_fold_glmnet_active_sets",
    aggregation = summary$aggregation,
    frequency_denominator = "successful_glmnet_coefficient_folds",
    n_candidate_features = summary$n_candidate_features,
    n_folds = summary$n_folds,
    n_successful_folds = summary$n_successful_folds,
    n_pairs = summary$n_pairs,
    n_evaluable_pairs = summary$n_evaluable_pairs,
    mean_pairwise_jaccard = summary$mean_pairwise_jaccard,
    feature_frequency = summary$feature_frequency,
    per_fold = summary$per_fold,
    pairwise = summary$pairwise,
    error_message = summary$error_message
  )
}
