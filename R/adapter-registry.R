sb_uniform_learner_adapter_families <- function() {
  list(
    list(sb_akritas_learner_ids(), sb_akritas_learner_adapter),
    list(sb_coxphf_learner_ids(), sb_coxphf_learner_adapter),
    list(sb_mboost_cox_learner_ids(), sb_mboost_cox_learner_adapter),
    list(sb_mboost_aft_learner_ids(), sb_mboost_aft_learner_adapter),
    list(
      sb_mboost_ranking_learner_ids(),
      sb_mboost_ranking_learner_adapter
    ),
    list(sb_xgboost_learner_ids(), sb_xgboost_learner_adapter),
    list(sb_mtlr_survival_learner_ids(), sb_mtlr_survival_learner_adapter),
    list(
      sb_pamr_survival_learner_ids(),
      sb_pamr_survival_learner_adapter
    ),
    list(
      sb_additional_survival_learner_ids(),
      sb_additional_survival_learner_adapter
    ),
    list(
      sb_flexsurv_parametric_learner_ids(),
      sb_flexsurv_parametric_learner_adapter
    ),
    list(
      sb_flexsurvspline_learner_ids(),
      sb_flexsurvspline_learner_adapter
    ),
    list(sb_timereg_survival_learner_ids(), sb_timereg_survival_learner_adapter),
    list(
      sb_classical_binary_learner_ids(),
      sb_classical_binary_learner_adapter
    ),
    list(
      sb_specialized_binary_linear_learner_ids(),
      sb_specialized_binary_linear_learner_adapter
    ),
    list(
      sb_omics_specialist_binary_learner_ids(),
      sb_omics_specialist_binary_learner_adapter
    ),
    list(sb_abess_learner_ids(), sb_abess_learner_adapter),
    list(sb_ncvreg_learner_ids(), sb_ncvreg_learner_adapter),
    list(sb_rlt_learner_ids(), sb_rlt_learner_adapter),
    list(
      sb_bagged_knn_survival_learner_ids(),
      sb_bagged_knn_survival_learner_adapter
    ),
    list(sb_ahaz_survival_learner_ids(), sb_ahaz_survival_learner_adapter),
    list(sb_bujar_linear_learner_ids(), sb_bujar_linear_learner_adapter),
    list(sb_rule_aft_learner_ids(), sb_rule_aft_learner_adapter),
    list(sb_wsrf_binary_learner_ids(), sb_wsrf_binary_learner_adapter),
    list(
      sb_slope_randomglm_binary_learner_ids(),
      sb_slope_randomglm_binary_learner_adapter
    ),
    list(sb_mclust_da_learner_ids(), sb_mclust_da_learner_adapter),
    list(sb_mixomics_da_learner_ids(), sb_mixomics_da_learner_adapter),
    list(sb_hdda_binary_learner_ids(), sb_hdda_binary_learner_adapter)
  )
}

sb_resolve_uniform_learner_adapter <- function(id, params) {
  for (family in sb_uniform_learner_adapter_families()) {
    if (id %in% family[[1L]]) {
      return(sb_constant_sum_learner_adapter(family[[2L]](id, params), family[[2L]]))
    }
  }
  NULL
}

sb_resolve_selector_adapter <- function(id, params) {
  if (id %in% sb_pamr_selector_ids()) return(sb_pamr_selector_adapter(id, params))
  if (id %in% sb_survival_weight_selector_ids()) return(sb_survival_weight_selector_adapter(id, params))
  if (id %in% sb_survival_boosting_selector_ids()) return(sb_survival_boosting_selector_adapter(id, params))
  if (identical(id, "ordered_chain")) return(sb_chain_selector_adapter(params))
  if (id %in% sb_binary_model_selector_ids()) return(sb_binary_model_selector_adapter(id, params))
  if (identical(id, "none")) {
    return(sb_none_selector_adapter())
  }
  if (identical(id, "univ_cox_top_k")) {
    return(sb_univ_cox_top_k_selector_adapter(params$n_features))
  }
  if (identical(id, "univ_binary_auc_top_k")) {
    return(sb_univ_binary_auc_top_k_selector_adapter(params$n_features))
  }
  if (identical(id, "stepcox")) {
    return(sb_stepcox_selector_adapter(params$direction))
  }
  if (identical(id, "rsf_minimal_depth")) {
    return(sb_rsf_minimal_depth_selector_adapter(
      params$nodesize,
      params$seed
    ))
  }
  if (identical(id, "coxboost_nonzero")) {
    return(sb_coxboost_nonzero_selector_adapter(
      params$inner_folds,
      params$seed
    ))
  }
  if (id %in% c("lasso_nonzero", "enet_cox_nonzero")) {
    return(sb_lasso_nonzero_selector_adapter(
      params$inner_folds,
      params$seed, id = id
    ))
  }
  if (identical(id, "rsf_then_stepcox")) {
    return(sb_rsf_then_stepcox_selector_adapter(
      params$nodesize,
      params$seed,
      params$direction
    ))
  }
  if (identical(id, "coxboost_then_stepcox")) {
    return(sb_coxboost_then_stepcox_selector_adapter(
      params$inner_folds,
      params$seed,
      params$direction
    ))
  }
  if (identical(id, "lasso_then_stepcox")) {
    return(sb_lasso_then_stepcox_selector_adapter(
      params$inner_folds,
      params$seed,
      params$direction
    ))
  }
  if (id %in% sb_abess_selector_ids()) {
    return(sb_abess_selector_adapter(id))
  }
  if (id %in% sb_forest_selector_ids()) {
    return(sb_forest_selector_adapter(id, options = params))
  }
  if (id %in% sb_high_dimensional_selector_ids()) {
    return(sb_high_dimensional_selector_adapter(id))
  }
  if (identical(id, "fixed_features")) {
    return(sb_fixed_features_selector_adapter(
      params$features,
      params$feature_set_id
    ))
  }
  sb_abort(
    paste0("No admitted selector adapter is registered for `", id, "`."),
    "stormbreaker_unresolved_selector_adapter"
  )
}

sb_resolve_learner_adapter <- function(id, params = list()) {
  id <- sb_canonical_learner_id(id)
  sb_reject_retired_learner(id)

  if (identical(id, "coxph")) {
    return(sb_coxph_learner_adapter())
  }
  if (identical(id, "glmnet_cox")) {
    return(list(
      id = "glmnet_cox",
      params = params,
      requires_tuning = TRUE
    ))
  }
  if (identical(id, "glmnet_logistic")) {
    return(list(
      id = "glmnet_logistic",
      params = params,
      requires_tuning = TRUE
    ))
  }
  if (identical(id, "logistic")) {
    return(sb_logistic_learner_adapter())
  }
  if (identical(id, "gaussian_naive_bayes")) {
    return(sb_gaussian_naive_bayes_learner_adapter(params$var_smoothing))
  }
  if (id %in% sb_lightgbm_binary_learner_ids()) {
    return(sb_lightgbm_binary_learner_adapter(
      nrounds = params$nrounds,
      learning_rate = params$learning_rate,
      num_leaves = params$num_leaves,
      min_data_in_leaf = params$min_data_in_leaf,
      seed = params$seed,
      learner = id
    ))
  }
  if (identical(id, "catboost_binary")) {
    return(sb_catboost_binary_learner_adapter(
      iterations = params$iterations,
      learning_rate = params$learning_rate,
      depth = params$depth,
      l2_leaf_reg = params$l2_leaf_reg,
      seed = params$random_seed
    ))
  }
  if (identical(id, "catboost_cox")) {
    return(sb_catboost_cox_learner_adapter(
      iterations = params$iterations,
      learning_rate = params$learning_rate,
      depth = params$depth,
      l2_leaf_reg = params$l2_leaf_reg,
      seed = params$random_seed
    ))
  }
  if (identical(id, "catboost_aft_survival")) {
    return(sb_catboost_aft_survival_learner_adapter(
      iterations = params$iterations,
      learning_rate = params$learning_rate,
      depth = params$depth,
      l2_leaf_reg = params$l2_leaf_reg,
      seed = params$random_seed
    ))
  }
  if (identical(id, "cancerclass_binary")) {
    return(sb_cancerclass_binary_learner_adapter())
  }
  if (identical(id, "gbm_cox")) {
    return(sb_gbm_cox_learner_adapter(
      n_trees = params$n_trees,
      interaction_depth = params$interaction_depth,
      min_obs_node = params$min_obs_node,
      shrinkage = params$shrinkage,
      bag_fraction = params$bag_fraction,
      seed = params$seed
    ))
  }
  if (identical(id, "gbm_cox_cv")) {
    return(sb_gbm_cox_cv_learner_adapter(
      inner_folds = params$inner_folds,
      seed = params$seed,
      maximum_trees = params$maximum_trees
    ))
  }
  if (identical(id, "coxboost")) {
    return(sb_coxboost_learner_adapter(
      stepno = params$stepno,
      penalty = params$penalty
    ))
  }
  if (identical(id, "coxboost_cv")) {
    return(sb_coxboost_cv_learner_adapter(
      inner_folds = params$inner_folds,
      seed = params$seed,
      maximum_steps = params$maximum_steps,
      start_penalty = params$start_penalty
    ))
  }
  if (identical(id, "plsrcox")) {
    return(sb_plsrcox_learner_adapter(
      components = params$components
    ))
  }
  if (identical(id, "plsrcox_cv")) {
    return(sb_plsrcox_cv_learner_adapter(
      inner_folds = params$inner_folds,
      seed = params$seed,
      maximum_components = params$maximum_components
    ))
  }
  if (identical(id, "superpc")) {
    return(sb_superpc_learner_adapter(
      threshold = params$threshold,
      components = params$components
    ))
  }
  if (id %in% sb_survivalsvm_learner_ids()) {
    return(sb_survivalsvm_learner_adapter(
      gamma_mu = params$gamma_mu,
      learner = id,
      params = params
    ))
  }
  if (id %in% sb_ranger_binary_learner_ids()) {
    return(sb_ranger_binary_learner_adapter(
      num_trees = params$num_trees,
      mtry = params$mtry,
      min_node_size = params$min_node_size,
      seed = params$seed,
      learner = id
    ))
  }
  if (id %in% sb_classical_binary_boosting_learner_ids()) {
    return(sb_classical_binary_boosting_learner_adapter(
      learner = id,
      iterations = params$iterations
    ))
  }
  if (identical(id, "superpc_cv")) {
    return(sb_superpc_cv_learner_adapter(
      inner_folds = params$inner_folds,
      seed = params$seed
    ))
  }
  if (id %in% sb_ranger_survival_learner_ids()) {
    return(sb_ranger_survival_learner_adapter(
      num_trees = params$num_trees,
      mtry = params$mtry,
      min_node_size = params$min_node_size,
      seed = params$seed,
      learner = id
    ))
  }
  adapter <- sb_resolve_uniform_learner_adapter(id, params)
  if (!is.null(adapter)) {
    return(adapter)
  }
  sb_abort(
    paste0("No admitted learner adapter is registered for `", id, "`."),
    "stormbreaker_unresolved_learner_adapter"
  )
}

sb_resolve_pipeline_runtime <- function(pipeline) {
  pipeline <- sb_upgrade_legacy_names(pipeline)

  if (!inherits(pipeline, "stormbreaker_pipeline")) {
    sb_abort(
      "`pipeline` must be created by `sb_pipeline()`.",
      "stormbreaker_bad_pipeline"
    )
  }
  runtime <- list(
    id = pipeline$id,
    selector = sb_resolve_selector_adapter(
      pipeline$selector,
      if (!is.null(pipeline$selector_params)) {
        pipeline$selector_params
      } else {
        list()
      }
    ),
    learner = sb_resolve_learner_adapter(
      pipeline$learner,
      if (!is.null(pipeline$learner_params)) {
        pipeline$learner_params
      } else {
        list()
      }
    )
  )
  class(runtime) <- "stormbreaker_pipeline_runtime"
  runtime
}

sb_validate_pipeline_runtime <- function(runtime, pipeline) {
  selector_params <- if (is.null(pipeline$selector_params)) {
    list()
  } else {
    pipeline$selector_params
  }
  learner_params <- if (is.null(pipeline$learner_params)) {
    list()
  } else {
    pipeline$learner_params
  }
  valid_runtime <- inherits(runtime, "stormbreaker_pipeline_runtime") &&
    identical(runtime$id, pipeline$id) &&
    is.list(runtime$selector) &&
    identical(runtime$selector$id, pipeline$selector) &&
    identical(runtime$selector$params, selector_params) &&
    is.logical(runtime$selector$needs_outcome) &&
    length(runtime$selector$needs_outcome) == 1L &&
    !is.na(runtime$selector$needs_outcome) &&
    is.function(runtime$selector$fit) &&
    is.function(runtime$selector$apply) &&
    is.list(runtime$learner) &&
    identical(runtime$learner$id, pipeline$learner)

  tuned_placeholder <- pipeline$learner %in% c(
    "glmnet_cox",
    "glmnet_logistic"
  )
  if (valid_runtime && tuned_placeholder) {
    valid_runtime <- identical(runtime$learner$params, learner_params) &&
      identical(runtime$learner$requires_tuning, TRUE)
  }
  if (valid_runtime && !tuned_placeholder) {
    required_functions <- c(
      "fit",
      "predict",
      "extract_state",
      "predict_state"
    )
    valid_runtime <- all(vapply(
      runtime$learner[required_functions],
      is.function,
      logical(1)
    ))
    if (valid_runtime && identical(pipeline$task_type, "survival")) {
      valid_runtime <- is.function(runtime$learner$predict_survival)
    }
    if (valid_runtime && (length(learner_params) > 0L ||
      "params" %in% names(runtime$learner))) {
      valid_runtime <- identical(runtime$learner$params, learner_params)
    }
  }
  if (!valid_runtime) {
    sb_abort(
      "The resolved pipeline runtime is incomplete or inconsistent.",
      "stormbreaker_bad_pipeline_runtime"
    )
  }
  runtime
}
