sb_additive_shap_linear_terms <- function(state, x, learner) {
  if (learner %in% c("lda_binary", sb_extended_logistic_learner_ids())) {
    p <- if (learner == "lda_binary") sb_discriminant_parameters(state) else sb_extended_logistic_parameters(state)
    return(list(constant = p$intercept, terms = sweep(as.matrix(x), 2L, p$coefficient, "*"),
      score_scale = "class_1_log_odds"))
  }
  if (learner %in% c(
    "logistic", "glmnet_logistic", "slope_binary", sb_abess_binary_learner_ids(),
    sb_ncvreg_binary_learner_ids()
  )) {
    return(list(
      constant = state$intercept,
      terms = sweep(as.matrix(x), 2L, state$coefficient, "*"),
      score_scale = "class_1_log_odds"
    ))
  }
  list(
    constant = 0,
    terms = sweep(as.matrix(x), 2L, state$coefficient, "*"),
    score_scale = "zero_reference_log_hazard_ratio"
  )
}

sb_additive_shap_naive_bayes_terms <- function(state, x) {
  model <- structure(
    state$payload,
    class = "stormbreaker_gaussian_naive_bayes"
  )
  matrix_x <- as.matrix(x)
  class_terms <- function(class) {
    centered <- sweep(matrix_x, 2L, model$means[class, ], "-")
    scaled_distance <- sweep(
      centered^2,
      2L,
      model$variances[class, ],
      "/"
    )
    log_variance <- matrix(
      log(2 * pi * model$variances[class, ]),
      nrow = nrow(matrix_x),
      ncol = ncol(matrix_x),
      byrow = TRUE
    )
    -0.5 * (log_variance + scaled_distance)
  }
  terms <- class_terms("1") - class_terms("0")
  dimnames(terms) <- list(NULL, names(x))
  list(
    constant = log(model$priors[["1"]]) - log(model$priors[["0"]]),
    terms = terms,
    score_scale = "class_1_log_posterior_odds"
  )
}

sb_abort_unsupported_additive_shap <- function(learner) {
  sb_abort(
    paste0(
      "Exact additive SHAP is not implemented for learner `", learner,
      "`; STORMBREAKER will not substitute a model-agnostic approximation."
    ),
    "stormbreaker_unsupported_additive_shap"
  )
}

sb_additive_shap_terms <- function(state, x) {
  learner <- state$learner
  if (learner == "earth_binary") return(sb_additive_shap_earth_terms(state, x))
  if (learner %in% c(
    "logistic", "glmnet_logistic", "lda_binary", "slope_binary", sb_abess_binary_learner_ids(),
    sb_ncvreg_binary_learner_ids(),
    sb_linear_survival_learner_ids(), sb_extended_logistic_learner_ids()
  )) {
    return(sb_additive_shap_linear_terms(state, x, learner))
  }
  if (identical(learner, "gaussian_naive_bayes")) {
    return(sb_additive_shap_naive_bayes_terms(state, x))
  }
  sb_abort_unsupported_additive_shap(learner)
}

#' Explain a final additive model with exact interventional SHAP values
#'
#' Also supports saved additive log-odds predictors from stepwise GLM, GLMBoost,
#' plsRglm, bias-reduced logistic and Cauchy-prior logistic regression, retaining
#' fitted centering and original parameter estimates.
#' LDA uses its exact original posterior log-odds equation, not a refitted GLM.
#' SLOPE retains its original selected penalized logistic coefficients.
#' Degree-one MARS uses its exact nonlinear additive hinge-basis contributions.
#'
#' `sb_additive_shap()` explains a successful final fit whose saved prediction
#' score is additive over features. It currently supports ordinary, glmnet,
#' MCP/SCAD, and ABESS best-subset logistic regression, every saved linear
#' survival learner, and Gaussian Naive Bayes. It reads only the requested ID
#' and fitted numeric
#' features from
#' `new_data`; outcomes and unrelated columns are ignored. `background_data`
#' likewise supplies only a caller-chosen empirical feature reference and need
#' not contain an ID or outcome.
#'
#' For an additive score \eqn{f(x) = c + \sum_j g_j(x_j)}, the returned value for
#' feature \eqn{j} is
#' \eqn{g_j(x_j) - mean_b(g_j(x_{bj}))}, where the mean is over the supplied
#' background rows. The reference value is
#' \eqn{c + \sum_j mean_b(g_j(x_{bj}))}. Their sum reconstructs the model score.
#' These are exact empirical interventional SHAP values for this additive score;
#' they are not a coefficient-times-input shortcut mislabeled as SHAP.
#'
#' The score is class-1 log odds for logistic regression, class-1 log posterior
#' odds for Gaussian Naive Bayes, and the saved zero-reference `x_beta`
#' log-hazard ratio for a supported linear survival learner. Contributions are
#' not changes in
#' probability, survival probability, or absolute risk. The returned final
#' prediction keeps those task-specific outputs separate.
#'
#' The background choice is part of the scientific question and changes both
#' the reference value and contributions. It should be declared and reported;
#' the function never defaults to the explained patients or reads their
#' outcomes. This interventional explanation does not preserve correlations or
#' estimate conditional SHAP values. It is descriptive, not causal, and does
#' not rank features or measure explanation stability.
#'
#' Tree learners are deliberately unsupported here because their prediction
#' functions are not additive feature scores. The function does not silently
#' replace exact model-specific SHAP with a permutation approximation.
#'
#' @param fit A successful additive-model `stormbreaker_final_fit` returned by
#'   [sb_fit_final()].
#' @param new_data A non-empty data frame containing one row per explained
#'   patient, the requested ID column, and every fitted numeric feature. Other
#'   columns are ignored.
#' @param id_col Name of the unique patient identifier in `new_data`.
#' @param background_data A non-empty data frame containing every fitted numeric
#'   feature. Rows define the empirical interventional reference distribution.
#'   An ID and outcomes are neither required nor read; unrelated columns are
#'   ignored.
#'
#' @return A `stormbreaker_additive_shap` object containing the model and score
#'   semantics, background feature and score-term means, scalar reference
#'   value, a patient-by-feature contribution matrix, reconstructed additive
#'   scores, and the unchanged state-only final prediction.
#'
#' @references Lundberg SM, Lee S-I (2017). A unified approach to interpreting
#'   model predictions. *Advances in Neural Information Processing Systems 30*,
#'   4765-4774. \url{https://arxiv.org/abs/1705.07874}.
#' @export
#'
#' @examples
#' development <- sb_example_binary(n = 40, seed = 1)
#' task <- sb_binary_task(
#'   development,
#'   "patient_id",
#'   "outcome",
#'   c("gene_a", "gene_b")
#' )
#' fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
#' sb_additive_shap(
#'   fit,
#'   development[1:3, ],
#'   "patient_id",
#'   background_data = development[4:20, ]
#' )
sb_additive_shap <- function(fit, new_data, id_col, background_data) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  learner <- fit$learner_state$learner
  supported <- learner %in% c(
    "earth_binary",
    "logistic", "glmnet_logistic", "lda_binary", "slope_binary", sb_abess_binary_learner_ids(),
    sb_ncvreg_binary_learner_ids(),
    sb_linear_survival_learner_ids(), sb_extended_logistic_learner_ids(),
    "gaussian_naive_bayes"
  )
  if (!supported) {
    sb_abort_unsupported_additive_shap(learner)
  }

  patient_ids <- sb_extract_final_prediction_ids(new_data, id_col)
  x <- sb_select_final_prediction_features(
    new_data,
    id_col,
    fit$required_features
  )
  background_x <- sb_select_final_feature_frame(
    background_data,
    fit$required_features,
    "Background data",
    "background_data fitted features"
  )

  prediction <- sb_predict_final(fit, new_data, id_col)
  if (!identical(prediction$status, "success")) {
    sb_abort(
      paste0(
        "The final prediction needed for additive SHAP failed: ",
        prediction$error_message
      ),
      "stormbreaker_additive_shap_prediction_failure"
    )
  }

  new_components <- sb_additive_shap_terms(fit$learner_state, x)
  background_components <- sb_additive_shap_terms(
    fit$learner_state,
    background_x
  )
  reference_terms <- colMeans(background_components$terms)
  reference_value <- background_components$constant + sum(reference_terms)
  contributions <- sweep(
    new_components$terms,
    2L,
    reference_terms,
    "-"
  )
  finite <- is.finite(reference_value) &&
    all(is.finite(reference_terms)) && all(is.finite(contributions))
  if (!finite) {
    sb_abort(
      "The additive SHAP calculation produced non-finite score values.",
      "stormbreaker_nonfinite_additive_shap"
    )
  }
  rownames(contributions) <- patient_ids
  colnames(contributions) <- fit$required_features
  scores <- data.frame(
    patient_id = patient_ids,
    score = unname(reference_value + rowSums(contributions)),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      status = "success",
      task_type = fit$task_type,
      pipeline_id = fit$pipeline$id,
      learner = learner,
      method = "exact_interventional_shap_for_additive_score",
      score_scale = new_components$score_scale,
      background = list(
        n = as.integer(nrow(background_x)),
        feature_summary = data.frame(
          feature = fit$required_features,
          mean_feature_value = unname(colMeans(background_x)),
          mean_score_term = unname(reference_terms),
          stringsAsFactors = FALSE
        )
      ),
      reference_value = unname(reference_value),
      contributions = contributions,
      scores = scores,
      prediction = prediction
    ),
    class = "stormbreaker_additive_shap"
  )
}

#' @export
print.stormbreaker_additive_shap <- function(x, ...) {
  cat("<stormbreaker_additive_shap:", x$status, ">\n")
  cat("exact interventional SHAP for an additive model score\n")
  cat(
    "Patients:", nrow(x$contributions),
    "| Features:", ncol(x$contributions),
    "| Background rows:", x$background$n, "\n"
  )
  cat("Pipeline:", x$pipeline_id, "| Learner:", x$learner, "\n")
  cat("Score scale:", x$score_scale, "\n")
  cat("Contributions are not probability-scale effects.\n")
  cat("No feature is ranked and no model output is changed.\n")
  invisible(x)
}
