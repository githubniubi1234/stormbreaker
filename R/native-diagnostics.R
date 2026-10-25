sb_native_diagnostic_plan <- function(id) {
  route <- function(extractor, types) data.frame(extractor = extractor, type = types)
  if (id == "cancerclass_binary") return(route("sb_cancerclass_diagnostics", c("centroids", "score_mapping")))
  if (id == "loclda_binary") return(route("sb_loclda_neighborhood", ""))
  if (id %in% c("hda_binary", "hdda_binary_cattell")) return(route("sb_high_dimensional_diagnostics", if (id == "hda_binary") c("loadings", "densities") else c("loadings", "spectrum", "centers")))
  if (id %in% sb_mclust_da_learner_ids()) return(route("sb_mclust_da_diagnostics", c("centers", "correlations", "mixture_weights")))
  if (id %in% c("sda_binary", "slda_binary")) return(route("sb_shrinkage_discriminant_diagnostics", if (id == "sda_binary") c("coefficients", "regularization") else c("loadings", "class_means")))
  if (id %in% c("fda_binary", "pda_binary", "mda_binary")) return(route("sb_fda_diagnostics", c("centers", if (id == "mda_binary") "mixture_weights", "query_variates")))
  if (id == "gaussian_process_binary") return(route("sb_gaussian_process_diagnostics", c("basis_weights", "query_contributions")))
  if (id == "random_ferns_binary") return(route("sb_ferns_diagnostics", c("split_usage", "thresholds")))
  if(id=="rotation_forest_binary") return(route("sb_rotation_forest_diagnostics",c("loadings","tree_size","tree_predictions")))
  if (id %in% c("c50_binary","c50_rules_binary")) return(route("sb_c50_diagnostics",c("usage","splits","size")))
  if (id == "tree_binary") return(route("sb_classification_tree_diagnostics", c("structure","leaf_probabilities")))
  if (id == "evtree_binary") return(route("sb_evtree_structure", ""))
  if (id %in% c("bart_binary", "bart_survival", "bart_aft_survival")) return(route("sb_bart_posterior", ""))
  if (id == "bagged_knn_survival") return(route("sb_bagged_knn_diagnostics", c("bootstrap", "bag_curves")))
  if (id %in% c("knn_binary", "kernel_weighted_knn_binary")) return(route("sb_knn_neighborhood", ""))
  if (id == "bagging_binary") return(route("sb_bagging_diagnostics", c("split_usage", "tree_size", "importance")))
  if (id %in% c("cforest_binary", "cforest_survival")) return(route("sb_cforest_diagnostics", c("split_usage", "tree_size")))
  if (id %in% c("aorsf_binary", "aorsf_survival")) return(route("sb_aorsf_diagnostics", c("tree_size", "oob_score")))
  if (id %in% sb_rlt_learner_ids()) return(route("sb_rlt_diagnostics", c("split_usage", "tree_size")))
  if (id == "wsrf_binary") return(route("sb_wsrf_diagnostics", c("importance", "tree_nodes", "tree_oob")))
  if (id %in% c("pamr_binary", "pamr_survival")) return(route("sb_pamr_diagnostics", ""))
  if (id == "penaft_gehan_elastic_net") return(route("sb_penaft_diagnostics", c("coefficients", "cv_curve", "coefficient_path")))
  if (id == sb_bujar_linear_learner_ids()) return(route("sb_bujar_diagnostics", c("coefficients", "residual_survival")))
  if (id %in% sb_ahaz_survival_learner_ids()) return(route("sb_ahaz_diagnostics",
    if (id == "ahaz_lin_ying") "coefficients" else c("coefficients", "bic_path", "coefficient_path")))
  if (id %in% sb_timereg_survival_learner_ids()) return(route("sb_timereg_diagnostics", c("effects", "baseline")))
  if (id == "mtlr_survival") return(route("sb_mtlr_diagnostics", c("feature_weights", "biases")))
  if (id == "ktsp_binary") return(route("sb_ktsp_diagnostics", ""))
  if (id %in% sb_flexsurv_parametric_learner_ids()) return(route("sb_flexsurv_diagnostics", "location_coefficients"))
  if (id %in% sb_flexsurvspline_learner_ids()) return(route("sb_flexsurv_diagnostics", c("location_coefficients", "spline_baseline")))
  if (id %in% sb_survreg_learner_ids()) return(route("sb_aft_diagnostics", c("time_ratios", "log_time_coefficients")))
  if (id %in% c("ctree_binary", "ctree_survival")) return(route("sb_ctree_structure", ""))
  if (id == "earth_binary") return(route("sb_earth_diagnostics", c("basis_coefficients", "pruning_gcv")))
  if (id %in% c("nnet_binary", "neuralnet_multilayer_binary")) return(route("sb_neural_diagnostics", c("topology", "weights")))
  if (id == "random_glm_binary") return(route("sb_randomglm_diagnostics", c("inclusion_frequency", "bag_coefficients", "bag_size")))
  if (id == "slope_binary") return(rbind(route("sb_slope_diagnostics", ""),
    route("sb_linear_coefficients_diagnostic", "")))
  if (id %in% sb_lvq_binary_learner_ids()) return(route("sb_lvq_diagnostics", c("prototypes", "prototype_distances")))
  if (id %in% c("lda_binary", "qda_binary", "rda_binary")) return(rbind(
    if (id == "lda_binary") route("sb_linear_coefficients_diagnostic", ""),
    route("sb_discriminant_diagnostics", c("class_means", if (id != "lda_binary") "class_correlation"))))
  if (id %in% c("gaussian_naive_bayes", "kernel_naive_bayes_binary")) return(
    route("sb_naive_bayes_diagnostics", c("conditional_density", if (id == "gaussian_naive_bayes") "log_likelihood_ratio")))
  if (id %in% c("glmnet_logistic", "glmnet_cox")) return(
    route("sb_algorithm_diagnostics", c("cv_curve", "nonzero_path", "coefficient_path")))
  if (id %in% sb_ncvreg_learner_ids()) return(route("sb_algorithm_diagnostics", c("coefficient_path", "bic_path")))
  if (id %in% sb_abess_learner_ids()) return(route("sb_algorithm_diagnostics", c("support_criterion", "coefficient_path")))
  if (id %in% sb_mgcv_gam_learner_ids()) return(route("sb_gam_effects", ""))
  if (id %in% c("decision_tree_binary", "rpart_survival")) return(route("sb_tree_structure", ""))
  if (id %in% sb_pre_rule_ensemble_learner_ids()) return(route("sb_rule_ensemble", c("rule_coefficients", "rule_cv_curve")))
  if (id %in% c(sb_ranger_learner_ids(), "random_forest_binary", "rrf_binary", "rsf_survival")) return(
    route("sb_forest_diagnostics", c("importance", "tree_size", "oob_error")))
  if (id %in% c("plsrglm_binary", "plsda_binary", "splsda_binary", "plsrcox", "plsrcox_cv")) return(
    rbind(route("sb_pls_diagnostics", c("loadings", "score_distribution", if (id == "plsrcox_cv") "cv_curve")),
      if (id == "plsrglm_binary") route("sb_linear_coefficients_diagnostic", "")))
  if (id %in% c("superpc", "superpc_cv")) return(route("sb_superpc_diagnostics",
    c("feature_scores", "loadings", if (id == "superpc_cv") "cv_curve")))
  if (id %in% sb_survivalsvm_learner_ids()) return(route("sb_survival_svm_diagnostics",
    c("score_distribution", "basis_weights", "coefficients")))
  if (id %in% sb_linear_margin_diagnostic_ids()) return(route("sb_svm_diagnostics", c("coefficients", "decision_distribution")))
  if (id == "lssvm_radial_binary") return(route("sb_svm_diagnostics", "decision_distribution"))
  if (id %in% c("svm_linear_binary", "svm_radial_binary", "svm_radial_weighted_binary",
      "svm_polynomial_binary", "svm_sigmoid_binary")) return(
    route("sb_svm_diagnostics", c("support_counts", "decision_distribution", if (id == "svm_linear_binary") "coefficients")))
  if (id == "bst_hinge_binary") return(route("sb_boosting_diagnostics", c("training_loss", "component_usage", "coefficient_path")))
  if (id %in% sb_classical_binary_boosting_learner_ids()) return(route("sb_boosting_diagnostics",
    c("split_usage", if (id == "adaboost_binary") c("tree_weights", "tree_size") else "stump_sequence")))
  if (id %in% c("coxboost", "coxboost_cv")) return(route("sb_boosting_diagnostics",
    c("coefficient_path", "nonzero_path", if (id == "coxboost_cv") "cv_curve")))
  if (id %in% sb_lightgbm_binary_learner_ids()) return(route("sb_boosting_diagnostics", c("importance", "tree_size")))
  if (id %in% sb_catboost_diagnostic_ids()) return(route("sb_boosting_diagnostics", c("importance", "interaction")))
  if (id %in% sb_xgboost_learner_ids()) return(route("sb_boosting_diagnostics", "importance"))
  if (id %in% sb_gbm_diagnostic_ids()) return(route("sb_boosting_diagnostics",
    c("importance", "training_loss", if (id == "gbm_cox_cv") "cv_curve")))
  if (id %in% sb_mboost_diagnostic_ids()) return(route("sb_boosting_diagnostics", c("training_loss", "component_usage")))
  if (id %in% c("logistic", sb_extended_logistic_learner_ids(), sb_linear_survival_learner_ids())) return(
    route("sb_linear_coefficients_diagnostic", ""))
  data.frame(extractor = character(), type = character())
}

#' Detect and extract implemented diagnostics for a saved model
#' @param fit A successful final fit. Missing captured artifacts are reported;
#'   this function never enables capture, refits, tunes or creates a surrogate.
#' @param complete FALSE stops after the first available registered diagnostic;
#'   TRUE attempts every registered native view, retaining all success/failure
#'   statuses. The registry is implemented coverage, not a claim that every
#'   catalog option currently has diagnostics or every candidate is available.
#' @param newdata Optional explicit query rows for nearest-neighbor, BART,
#'   Gaussian-process, local-LDA, FDA/PDA/MDA and rotation-forest diagnostics.
#'   Other families do not use this argument. No implicit training
#'   queries are generated when it is NULL.
#' @return A `stormbreaker_native_diagnostics` list with extracted `diagnostics`
#'   and a `status` table, including reasons and error classes. Registration
#'   alone is not availability: each requested extraction actually executes.
#'   Unsupported learners get `not_implemented`; missing original artifacts get
#'   `missing_capture`; expected unavailability stays distinct from errors.
#'   Captured sample-level score diagnostics can remain sensitive.
#' @seealso [sb_clinician_views()]
#' @export
sb_native_diagnostics <- function(fit, complete = FALSE, newdata = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!is.logical(complete) || length(complete) != 1L || is.na(complete)) sb_abort(
    "complete must be TRUE or FALSE.", "stormbreaker_bad_report_input")
  plan <- sb_native_diagnostic_plan(fit$pipeline$learner)
  diagnostics <- list()
  rows <- list()
  if (!nrow(plan)) rows[[1L]] <- data.frame(view = fit$pipeline$learner, extractor = "", type = "",
    status = "not_implemented", reason = "No native diagnostic route is implemented for this learner.", error_class = "")
  for (i in seq_len(nrow(plan))) {
    extractor <- get(plan$extractor[i], envir = asNamespace("stormbreaker"))
    args <- list(fit = fit)
    if ("newdata" %in% names(formals(extractor))) args$newdata <- newdata
    if (nzchar(plan$type[i])) args$type <- plan$type[i]
    value <- tryCatch(do.call(extractor, args), error = identity)
    label <- paste(fit$pipeline$learner, if (nzchar(plan$type[i])) plan$type[i] else sub("sb_", "", plan$extractor[i]), sep = ": ")
    status <- "available"
    reason <- error_class <- ""
    if (inherits(value, "error")) {
      status <- if (inherits(value, "stormbreaker_missing_diagnostic_artifacts")) "missing_capture" else
        if (inherits(value, c("stormbreaker_unavailable_diagnostic", "stormbreaker_unsupported_diagnostic"))) "unavailable" else "error"
      reason <- conditionMessage(value)
      error_class <- class(value)[1L]
    } else diagnostics[[label]] <- value
    rows[[i]] <- data.frame(view = label, extractor = plan$extractor[i], type = plan$type[i],
      status = status, reason = reason, error_class = error_class)
    if (!complete && status == "available") break
  }
  structure(list(diagnostics = diagnostics, status = do.call(rbind, rows), learner = fit$pipeline$learner,
    pipeline_id = fit$pipeline$id), class = "stormbreaker_native_diagnostics")
}
