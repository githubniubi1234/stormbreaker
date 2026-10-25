sb_tree_shap_ids <- function() c("random_forest_binary", "gbm_binary",
  "gbm_cox", "gbm_cox_cv", "rsf_survival", "decision_tree_binary",
  "rpart_survival", "cforest_binary", "cforest_survival",
  sb_ranger_learner_ids(), sb_xgboost_learner_ids(),
  sb_lightgbm_binary_learner_ids(), "catboost_binary", "catboost_cox")

#' Approximate interventional SHAP for saved tree-model predictions
#' @param fit Successful saved binary or survival tree-model final fit.
#' @param new_data Observations to explain, including ID and predictor columns.
#' @param id_col ID column in new_data.
#' @param background_data Development/reference predictors, not outcome columns.
#' @param permutations Number of random feature orders, default 32.
#' @param background_n Maximum background rows, default 20 (seeded sampling).
#' @param seed Seed for orders and background sampling.
#' @return A stormbreaker_tree_shap result accepted by existing SHAP plots.
#' @details This is Monte Carlo permutation interventional SHAP, NOT exact
#'   TreeSHAP. Every path telescopes from the same empirical background mean
#'   to the explained prediction, so additivity alone does not prove estimator
#'   accuracy. Binary explanations use class-1 probability; survival uses the
#'   saved risk score, not survival probability or causal effects. Runtime grows
#'   with features, patients, permutations and background size. No refitting.
#' @examples
#' if (requireNamespace("rpart", quietly = TRUE)) {
#'   d <- sb_example_binary(60, seed = 42)
#'   task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
#'   fit <- sb_fit_final(task, sb_pipeline(learner = "decision_tree_binary"))
#'   sb_tree_shap(fit, d[1:2, ], "patient_id", d,
#'     permutations = 2, background_n = 3)
#' }
#' @export
sb_tree_shap <- function(fit, new_data, id_col, background_data,
                         permutations = 32L, background_n = 20L, seed = 1L) {
  fit <- sb_validate_final_fit_for_prediction(sb_upgrade_legacy_names(fit))
  if (!fit$pipeline$learner %in% sb_tree_shap_ids()) stop("This learner is not in the tree-SHAP menu.")
  result <- sb_model_shap_validated(fit, new_data, id_col, background_data,
    permutations, background_n, seed)
  class(result) <- "stormbreaker_tree_shap"
  result
}

#' Approximate SHAP for a saved binary or survival model
#' @inheritParams sb_tree_shap
#' @param fit Successful saved binary or survival final fit, including non-tree models.
#' @return A stormbreaker_model_shap result accepted by the SHAP plotting functions.
#' @details Uses the same permutation interventional estimator as sb_tree_shap,
#'   with the original model's prediction interface. Supports non-tree models,
#'   including LogitBoost and RDA. Binary values explain class-1 probability;
#'   survival values explain the saved risk score. No model is refitted.
#'   Background data should come from the development cohort. This is an
#'   approximation, not exact TreeSHAP; more permutations reduce Monte Carlo noise.
#' @examples
#' d <- sb_example_binary(60, seed = 42)
#' task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
#' fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
#' s <- sb_model_shap(fit, d[1:2, ], "patient_id", d,
#'   permutations = 2, background_n = 3)
#' s$reference_value + rowSums(s$contributions)
#' @export
sb_model_shap <- function(fit, new_data, id_col, background_data,
                          permutations = 32L, background_n = 20L, seed = 1L) {
  fit <- sb_validate_final_fit_for_prediction(sb_upgrade_legacy_names(fit))
  sb_model_shap_validated(fit, new_data, id_col, background_data,
    permutations, background_n, seed)
}

sb_model_shap_validated <- function(fit, new_data, id_col, background_data,
                                    permutations, background_n, seed) {
  if (!fit$task_type %in% c("binary", "survival")) stop("Supply a binary or survival fit.")
  permutations <- sb_check_positive_integer(permutations, "permutations")
  background_n <- sb_check_positive_integer(background_n, "background_n")
  seed <- sb_check_positive_integer(seed, "seed")
  ids <- sb_extract_final_prediction_ids(new_data, id_col)
  features <- fit$required_features
  x <- sb_select_final_prediction_features(new_data, id_col, features)
  b <- sb_select_final_feature_frame(background_data, features,
    "Background data", "background_data")
  withr::local_seed(seed)
  if (nrow(b) > background_n) b <- b[sample.int(nrow(b), background_n), , drop = FALSE]
  score <- function(d) {
    values <- sb_predict_final_state(fit$learner_state, d)
    values <- sb_validate_final_prediction_values(fit, values, nrow(d))
    if (fit$task_type == "binary") values$probability else values$risk
  }
  reference <- mean(score(b))
  target <- score(x)
  if (length(target) != nrow(x) || any(!is.finite(target)) || !is.finite(reference)) stop("Saved prediction score is unavailable.")
  phi <- matrix(0, nrow(x), ncol(x), dimnames = list(ids, features))
  for (i in seq_len(nrow(x))) for (iteration in seq_len(permutations)) {
    order <- sample.int(ncol(x))
    hybrid <- b
    previous <- reference
    for (j in order) {
      hybrid[[j]] <- x[[j]][i]
      current <- mean(score(hybrid))
      phi[i, j] <- phi[i, j] + (current - previous) / permutations
      previous <- current
    }
  }
  structure(list(status = "success", task_type = fit$task_type,
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    method = "approximate_permutation_interventional_shap",
    score_scale = if (fit$task_type == "binary") "class_1_probability" else "saved_risk_score",
    reference_value = reference, contributions = phi,
    scores = data.frame(patient_id = ids, score = target),
    background = list(n = nrow(b)), permutations = permutations, seed = seed),
    class = "stormbreaker_model_shap")
}
