sb_pipeline_model_id <- function(base_id, learner, learner_params, uses_inner_cv) {
  is_catboost <- learner %in% c(
    "catboost_binary",
    "catboost_cox",
    "catboost_aft_survival"
  )
  is_lightgbm <- learner %in% sb_lightgbm_binary_learner_ids()
  is_gbm_fixed <- identical(learner, "gbm_cox")
  is_analysis_tuned <- learner %in% sb_analysis_tuned_survival_learner_ids()
  is_coxboost <- identical(learner, "coxboost")
  is_plsrcox <- identical(learner, "plsrcox")
  is_superpc <- identical(learner, "superpc")
  is_survivalsvm <- learner %in% sb_survivalsvm_learner_ids()
  is_classical_binary_boosting <- learner %in%
    sb_classical_binary_boosting_learner_ids()
  is_ranger <- learner %in% sb_ranger_learner_ids()
  pipeline_id <- if (uses_inner_cv) {
    standardize_suffix <- if (isTRUE(learner_params$standardize)) {
      ""
    } else {
      "__standardize0"
    }
    lambda_rule_suffix <- if (identical(
      learner_params$selection_rule,
      "lambda.min"
    )) {
      ""
    } else {
      "__lambda1se"
    }
    measure_suffix <- if (identical(learner_params$measure, "deviance")) {
      ""
    } else {
      paste0("__measure", learner_params$measure)
    }
    paste0(
      base_id,
      "__alpha", format(
        learner_params$alpha,
        trim = TRUE,
        scientific = FALSE,
        digits = 15L
      ),
      standardize_suffix,
      lambda_rule_suffix,
      measure_suffix,
      "__inner", learner_params$inner_folds,
      "__seed", learner_params$inner_seed
    )
  } else if (identical(learner, "gaussian_naive_bayes")) {
    paste0(
      base_id,
      "__smooth",
      sb_format_roundtrip_double(learner_params$var_smoothing)
    )
  } else if (is_lightgbm) {
    paste0(
      base_id,
      "__rounds", learner_params$nrounds,
      "__lr", sb_format_roundtrip_double(
        learner_params$learning_rate
      ),
      "__leaves", learner_params$num_leaves,
      "__minleaf", learner_params$min_data_in_leaf,
      "__seed", learner_params$seed
    )
  } else if (is_catboost) {
    paste0(
      base_id,
      "__iter", learner_params$iterations,
      "__lr", sb_format_roundtrip_double(
        learner_params$learning_rate
      ),
      "__depth", learner_params$depth,
      "__l2reg", sb_format_roundtrip_double(
        learner_params$l2_leaf_reg
      ),
      "__seed", learner_params$random_seed
    )
  } else if (identical(learner, "gbm_cox_cv")) {
    paste0(
      base_id,
      "__inner", learner_params$inner_folds,
      "__seed", learner_params$seed,
      "__maxtrees", learner_params$maximum_trees
    )
  } else if (is_analysis_tuned) {
    paste0(
      base_id,
      "__inner", learner_params$inner_folds,
      "__seed", learner_params$seed
    )
  } else if (is_gbm_fixed) {
    paste0(
      base_id,
      "__trees", learner_params$n_trees,
      "__depth", learner_params$interaction_depth,
      "__minnode", learner_params$min_obs_node,
      "__shrink", sb_format_roundtrip_double(
        learner_params$shrinkage
      ),
      "__bag", sb_format_roundtrip_double(
        learner_params$bag_fraction
      ),
      "__seed", learner_params$seed
    )
  } else if (is_ranger) {
    paste0(
      base_id,
      "__trees", learner_params$num_trees,
      "__mtry", if (is.null(learner_params$mtry)) {
        "sqrt"
      } else {
        learner_params$mtry
      },
      "__minnode", learner_params$min_node_size,
      "__seed", learner_params$seed
    )
  } else if (is_coxboost) {
    paste0(
      base_id,
      "__steps", learner_params$stepno,
      "__penalty", if (is.null(learner_params$penalty)) {
        "event9"
      } else {
        sb_format_roundtrip_double(learner_params$penalty)
      }
    )
  } else if (is_plsrcox) {
    paste0(base_id, "__components", learner_params$components)
  } else if (is_superpc) {
    paste0(
      base_id,
      "__threshold", sb_format_roundtrip_double(learner_params$threshold),
      "__components", learner_params$components
    )
  } else if (is_survivalsvm) {
    paste0(
      base_id,
      "__gamma", sb_format_roundtrip_double(learner_params$gamma_mu)
    )
  } else if (is_classical_binary_boosting) {
    paste0(base_id, "__iterations", learner_params$iterations)
  } else if (identical(learner, "stepglm_binary")) {
    paste0(base_id, "__direction", learner_params$direction)
  } else {
    base_id
  }
  pipeline_id
}
