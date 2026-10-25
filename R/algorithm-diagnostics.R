#' Extract native regularization paths from a captured final model
#' @param fit Final model. Native paths require capture_diagnostics = TRUE;
#'   saved glmnet tuning curves do not require extra capture.
#' @param type Diagnostic type: coefficient_path for native glmnet, ABESS or ncvreg
#'   engines, support_criterion for the ABESS GIC selection path, or cv_curve
#'   and nonzero_path for the saved glmnet inner-validation sequence, or
#'   bic_path for MCP/SCAD's original BIC selection sequence.
#' @return A `stormbreaker_algorithm_diagnostic` list with the native path
#'   table and origin metadata. No fitting or parameter search occurs.
#' @export
sb_algorithm_diagnostics <- function(fit, type = "coefficient_path") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("coefficient_path", "support_criterion", "cv_curve", "nonzero_path", "bic_path")) sb_abort(
    "This diagnostic type has not been implemented.", "stormbreaker_unsupported_diagnostic")
  if (type %in% c("cv_curve", "nonzero_path")) return(sb_glmnet_tuning_diagnostic(fit, type))
  if (is.null(fit$diagnostic_artifacts)) sb_abort(
    "Native artifacts were not captured. Explicitly refit on the original development task with capture_diagnostics = TRUE if this diagnostic is needed.",
    "stormbreaker_missing_diagnostic_artifacts")
  model <- fit$diagnostic_artifacts$learner
  if (inherits(model, "stormbreaker_ncvreg_model") &&
      type %in% c("coefficient_path", "bic_path")) return(sb_ncvreg_path_diagnostic(fit, type))
  if (identical(type, "bic_path")) sb_abort(
    "BIC-path diagnostics currently require an MCP/SCAD learner.",
    "stormbreaker_unsupported_diagnostic")
  engine <- if (is.list(model) && !is.null(model$engine)) model$engine else model
  if (type == "coefficient_path" && inherits(model, "stormbreaker_abess_model")) {
    b <- as.matrix(engine$beta)
    d <- data.frame(support_size = rep(engine$support.size, each = nrow(b)),
      feature = rep(model$features, ncol(b)), coefficient = as.vector(b))
    return(structure(list(type = "abess_coefficient_path", table = d,
      selected_support_size = model$selected_support_size, learner = fit$pipeline$learner,
      pipeline_id = fit$pipeline$id, source = "native ABESS support-size path"),
      class = "stormbreaker_algorithm_diagnostic"))
  }
  if (identical(type, "support_criterion")) {
    if (!inherits(model, "stormbreaker_abess_model")) sb_abort(
      "Support criterion diagnostics require the original ABESS model.",
      "stormbreaker_unsupported_diagnostic")
    path <- model$criterion_path
    selected <- model$selected_support_size
    if (!is.data.frame(path) || !all(c("support_size", "gic", "deviance") %in% names(path)) ||
        !nrow(path) || !all(vapply(path[c("support_size", "gic", "deviance")],
          function(z) is.numeric(z) && all(is.finite(z)), logical(1))) ||
        anyDuplicated(path$support_size) || any(path$support_size < 0) ||
        any(path$support_size != floor(path$support_size)) ||
        length(selected) != 1L || is.na(selected) || !selected %in% path$support_size) sb_abort(
      "Native ABESS selection path is incomplete or invalid.",
      "stormbreaker_invalid_diagnostic")
    path$selected <- path$support_size == selected
    return(structure(list(type = type, table = path,
      learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
      selected_support_size = selected,
      source = "native final-fit ABESS GIC path", engine_class = class(engine)),
      class = "stormbreaker_algorithm_diagnostic"))
  }
  if (!inherits(engine, c("glmnet", "ncvreg", "ncvsurv"))) sb_abort(
    "This learner does not have a supported native regularization path.",
    "stormbreaker_unsupported_diagnostic")
  beta <- as.matrix(engine$beta)
  lambda <- engine$lambda
  if (!is.numeric(lambda) || ncol(beta) != length(lambda) ||
      !length(lambda) || any(!is.finite(lambda)) || any(lambda <= 0) ||
      any(!is.finite(beta))) sb_abort(
    "Native regularization path is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
  labels <- rownames(beta)
  if (is.null(labels)) sb_abort("Native path feature labels are missing.",
    "stormbreaker_invalid_diagnostic")
  table <- data.frame(lambda = rep(lambda, each = nrow(beta)),
    feature = rep(labels, length(lambda)), coefficient = as.vector(beta))
  structure(list(type = type, table = table, learner = fit$pipeline$learner,
    pipeline_id = fit$pipeline$id, source = "native final-fit regularization path",
    engine_class = class(engine), tuning = fit$tuning),
    class = "stormbreaker_algorithm_diagnostic")
}

#' Plot a saved native algorithm diagnostic
#' @param x A `stormbreaker_algorithm_diagnostic` result, from
#'   [sb_algorithm_diagnostics()] or a family-specific diagnostic extractor such
#'   as [sb_naive_bayes_diagnostics()] or [sb_survival_svm_diagnostics()].
#' @return A ggplot displaying the requested saved native diagnostic, without
#'   fitting or changing the model. Consult the extractor's documentation and
#'   result metadata for the meaning of the values and any native uncertainty;
#'   a development diagnostic is not external validation performance.
#' @export
sb_plot_algorithm_diagnostic <- function(x) {
  sb_require_curve_plotting()
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type == "abess_coefficient_path") {
    support_size <- coefficient <- feature <- NULL
    return(ggplot2::ggplot(x$table, ggplot2::aes(support_size, coefficient, colour = feature)) +
      ggplot2::geom_line() +
      ggplot2::geom_vline(xintercept = x$selected_support_size, linetype = 2, colour = "grey50") +
      ggplot2::labs(x = "Number of selected predictors", y = "Coefficient", title = "ABESS coefficient path") + sb_viz_theme())
  }
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in% c("cancerclass_centroids", "cancerclass_score_mapping")) return(sb_plot_cancerclass_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type == "loclda_neighborhood") return(sb_plot_loclda_neighborhood(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in% paste0("high_dimensional_", c("loadings", "densities", "spectrum", "centers"))) return(sb_plot_high_dimensional_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in% paste0("mclust_da_", c("centers", "correlations", "mixture_weights"))) return(sb_plot_mclust_da_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in% paste0("shrinkage_discriminant_", c("coefficients", "regularization", "loadings", "class_means"))) return(sb_plot_shrinkage_discriminant_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in% c("fda_centers", "fda_query_variates", "fda_mixture_weights")) return(sb_plot_fda_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in% c("gaussian_process_basis_weights", "gaussian_process_query_contributions")) return(sb_plot_gaussian_process_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in% c("ferns_split_usage", "ferns_thresholds")) return(sb_plot_ferns_diagnostic(x))
  if(inherits(x,"stormbreaker_algorithm_diagnostic") && x$type %in% c("rotation_forest_loadings","rotation_forest_tree_size","rotation_forest_tree_predictions")) return(sb_plot_rotation_forest_diagnostic(x))
  if(inherits(x,"stormbreaker_algorithm_diagnostic") && x$type %in% c("c50_usage","c50_splits","c50_size")) return(sb_plot_c50_diagnostic(x))
  if (inherits(x,"stormbreaker_algorithm_diagnostic") && x$type == "classification_tree_leaf_probabilities") return(sb_plot_classification_tree_probabilities(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type == "bart_posterior") return(sb_plot_bart_posterior(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in% c("bagged_knn_bootstrap", "bagged_knn_bag_curves")) return(sb_plot_bagged_knn_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type == "knn_neighborhood") return(sb_plot_knn_neighborhood(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("ahaz_coefficients", "ahaz_bic_path", "ahaz_coefficient_path")) return(sb_plot_ahaz_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("bujar_coefficients", "bujar_residual_survival")) return(sb_plot_bujar_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("penaft_coefficients", "penaft_cv_curve", "penaft_coefficient_path")) return(sb_plot_penaft_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type == "pamr_centroids") return(sb_plot_pamr_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("wsrf_importance", "wsrf_tree_nodes", "wsrf_tree_oob")) return(sb_plot_wsrf_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("rlt_split_usage", "rlt_tree_size")) return(sb_plot_rlt_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("aorsf_tree_size", "aorsf_oob_score")) return(sb_plot_aorsf_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("cforest_split_usage", "cforest_tree_size")) return(sb_plot_cforest_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("bagging_split_usage", "bagging_tree_size", "bagging_importance")) return(sb_plot_bagging_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("timereg_aalen_effects", "timereg_aalen_baseline", "timereg_odds_effects", "timereg_odds_baseline")) return(sb_plot_timereg_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("mtlr_feature_weights", "mtlr_biases")) return(sb_plot_mtlr_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && identical(x$type, "ktsp_pair_scores")) return(sb_plot_ktsp_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("flexsurv_location_coefficients", "flexsurv_spline_baseline")) return(sb_plot_flexsurv_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("aft_time_ratios", "aft_log_time_coefficients")) return(sb_plot_aft_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && identical(x$type, "ctree_structure")) return(sb_plot_ctree_structure(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("earth_basis_coefficients", "earth_pruning_gcv")) return(sb_plot_earth_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("neural_topology", "neural_weights")) return(sb_plot_neural_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("randomglm_inclusion_frequency", "randomglm_bag_coefficients", "randomglm_bag_size")) return(sb_plot_randomglm_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type == "slope_cv_curve") return(sb_plot_slope_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("lvq_prototypes", "lvq_prototype_distances")) return(sb_plot_lvq_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("discriminant_class_means", "discriminant_class_correlation")) return(sb_plot_discriminant_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type == "linear_coefficients") {
    p <- sb_plot_coefficients(x$parameters, include_zero = TRUE) +
      ggplot2::labs(subtitle = NULL, caption = NULL, x = if (x$parameters$score_scale == "class_1_log_odds")
        "Coefficient (log odds/unit)" else "Coefficient (log HR/unit)")
    attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = max(4, 2 + .3 * nrow(x$table)))
    return(p)
  }
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("naive_bayes_conditional_density", "naive_bayes_log_likelihood_ratio")) return(sb_plot_naive_bayes_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") && x$type %in%
      c("survival_svm_score_distribution", "survival_svm_basis_weights", "survival_svm_coefficients")) return(sb_plot_survival_svm_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") &&
      x$type %in% c("svm_support_counts", "svm_decision_distribution", "svm_coefficients")) return(sb_plot_svm_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") &&
      x$type %in% c("pls_loadings", "pls_score_distribution", "pls_cv_curve")) return(sb_plot_pls_diagnostic(x))
  if (inherits(x, "stormbreaker_algorithm_diagnostic") &&
      x$type %in% c("superpc_feature_scores", "superpc_loadings", "superpc_cv_curve")) {
    return(sb_plot_superpc_diagnostic(x))
  }
  if (inherits(x, "stormbreaker_algorithm_diagnostic") &&
      x$type %in% c("hinge_training_loss", "hinge_component_usage", "hinge_coefficient_path")) {
    return(sb_plot_hinge_boosting_diagnostic(x))
  }
  if (inherits(x, "stormbreaker_algorithm_diagnostic") &&
      x$type %in% c("classical_boosting_split_usage", "classical_boosting_tree_weights",
        "classical_boosting_tree_size", "classical_boosting_stump_sequence")) {
    return(sb_plot_classical_boosting_diagnostic(x))
  }
  if (inherits(x, "stormbreaker_algorithm_diagnostic") &&
      x$type %in% c("coxboost_coefficient_path", "coxboost_nonzero_path", "coxboost_cv_curve")) {
    return(sb_plot_coxboost_diagnostic(x))
  }
  if (!inherits(x, "stormbreaker_algorithm_diagnostic") ||
      !x$type %in% c("coefficient_path", "support_criterion", "cv_curve", "nonzero_path", "bic_path", "gam_effects", "tree_structure", "rule_coefficients", "rule_cv_curve", "forest_importance", "forest_tree_size", "forest_oob_error", "boosting_training_loss", "boosting_component_usage", "boosted_tree_importance", "boosted_tree_training_loss", "boosted_tree_cv_curve", "light_cat_importance", "light_cat_tree_size", "light_cat_interaction")) sb_abort(
    "Supply a supported algorithm diagnostic.", "stormbreaker_bad_visualization_argument")
  lambda <- coefficient <- feature <- NULL
  if (x$type %in% c("light_cat_importance", "light_cat_tree_size", "light_cat_interaction")) return(sb_plot_light_cat_diagnostic(x))
  if (x$type %in% c("boosted_tree_importance", "boosted_tree_training_loss", "boosted_tree_cv_curve")) return(sb_plot_gbm_xgboost_diagnostic(x))
  if (x$type %in% c("boosting_training_loss", "boosting_component_usage")) return(sb_plot_boosting_diagnostic(x))
  if (x$type %in% c("forest_importance", "forest_tree_size", "forest_oob_error")) return(sb_plot_forest_diagnostic(x))
  if (x$type %in% c("rule_coefficients", "rule_cv_curve")) return(sb_plot_pre_diagnostic(x))
  if (identical(x$type, "tree_structure")) return(sb_plot_rpart_tree(x))
  if (identical(x$type, "gam_effects")) return(sb_plot_gam_effects(x))
  if (identical(x$type, "bic_path")) {
    bic <- NULL
    return(ggplot2::ggplot(x$table, ggplot2::aes(lambda, bic)) +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::geom_vline(xintercept = x$selected_lambda, linetype = 2, colour = "#D77935") +
      ggplot2::scale_x_log10() +
      ggplot2::labs(title = "MCP/SCAD BIC path", x = "Lambda (log scale)",
        y = "BIC (lower is better)",
        caption = "Dashed line: original selected lambda. Development-fit criterion, not cross-validation.\nNo refitting or replacement by an unpenalized model.") + sb_viz_theme())
  }
  if (x$type %in% c("cv_curve", "nonzero_path")) return(sb_plot_glmnet_tuning_diagnostic(x))
  if (identical(x$type, "support_criterion")) {
    support_size <- gic <- NULL
    return(ggplot2::ggplot(x$table, ggplot2::aes(support_size, gic)) +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::geom_point(colour = "#2878A0") +
      ggplot2::geom_point(data = x$table[x$table$selected, , drop = FALSE],
        shape = 21, fill = "#D77935", size = 3) +
      ggplot2::scale_x_continuous(breaks = x$table$support_size) +
      ggplot2::labs(title = "ABESS support-size selection",
        subtitle = paste("Original selected support size:", x$selected_support_size),
        x = "Number of selected predictors", y = "GIC (lower is better)",
        caption = "Original development-fit criterion, not held-out performance. No refitting.") +
      sb_viz_theme())
  }
  ggplot2::ggplot(x$table, ggplot2::aes(lambda, coefficient, colour = feature)) +
    ggplot2::geom_line() + ggplot2::scale_x_log10() +
    ggplot2::labs(title = "Regularization path", x = "Lambda (log scale)",
      y = "Coefficient", colour = "Feature",
      caption = "Original final-fit path; no refitting. Not a coefficient confidence interval or a cross-validation performance curve.") +
    sb_viz_theme()
}
