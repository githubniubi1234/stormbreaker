sb_pipeline_tuning_params <- function(learner, pipeline_task_type, uses_inner_cv, options) {
  is_analysis_tuned <- learner %in% sb_analysis_tuned_survival_learner_ids()
  learner_params <- list()
  if (!uses_inner_cv) {
    if (!is.null(options$inner_folds) || !is.null(options$inner_seed)) {
      sb_abort(
        paste0(
          "`inner_folds` and `inner_seed` apply only when ",
          "a glmnet learner is requested."
        ),
        "stormbreaker_irrelevant_inner_parameter"
      )
    }
    if (!is.null(options$glmnet_alpha) || !is.null(options$glmnet_standardize) ||
      !is.null(options$glmnet_lambda_rule) || !is.null(options$glmnet_tuning_measure)) {
      sb_abort(
        paste0(
          "`glmnet_alpha`, `glmnet_standardize`, `glmnet_lambda_rule`, ",
          "and `glmnet_tuning_measure` ",
          "apply only when a glmnet learner is requested."
        ),
        "stormbreaker_irrelevant_glmnet_parameter"
      )
    }
  } else {
    if (is.null(options$inner_folds) || is.null(options$inner_seed)) {
      sb_abort(
        paste0(
          "`inner_folds` and `inner_seed` are required when ",
          "a glmnet learner is requested."
        ),
        "stormbreaker_missing_inner_parameter"
      )
    }
    if (is.null(options$glmnet_alpha)) {
      options$glmnet_alpha <- 1
    }
    if (!is.numeric(options$glmnet_alpha) || length(options$glmnet_alpha) != 1L ||
      is.na(options$glmnet_alpha) || !is.finite(options$glmnet_alpha) ||
      options$glmnet_alpha < 0 || options$glmnet_alpha > 1) {
      sb_abort(
        "`glmnet_alpha` must be one finite number between 0 and 1.",
        "stormbreaker_bad_glmnet_alpha"
      )
    }
    if (is.null(options$glmnet_standardize)) {
      options$glmnet_standardize <- TRUE
    }
    if (!is.logical(options$glmnet_standardize) ||
      length(options$glmnet_standardize) != 1L || is.na(options$glmnet_standardize)) {
      sb_abort(
        "`glmnet_standardize` must be `TRUE` or `FALSE`.",
        "stormbreaker_bad_glmnet_standardize"
      )
    }
    if (is.null(options$glmnet_lambda_rule)) {
      options$glmnet_lambda_rule <- "lambda.min"
    }
    options$glmnet_lambda_rule <- sb_validate_glmnet_lambda_rule(
      options$glmnet_lambda_rule,
      argument = "glmnet_lambda_rule"
    )
    if (is.null(options$glmnet_tuning_measure)) {
      options$glmnet_tuning_measure <- "deviance"
    }
    options$glmnet_tuning_measure <- sb_validate_glmnet_tuning_measure(
      options$glmnet_tuning_measure,
      pipeline_task_type,
      argument = "glmnet_tuning_measure"
    )
    learner_params <- list(
      alpha = as.numeric(options$glmnet_alpha),
      standardize = options$glmnet_standardize,
      inner_folds = sb_check_positive_integer(
        options$inner_folds,
        "inner_folds",
        minimum = 3L
      ),
      inner_seed = sb_check_positive_integer(
        options$inner_seed,
        "inner_seed",
        minimum = 1L
      ),
      selection_rule = options$glmnet_lambda_rule,
      measure = options$glmnet_tuning_measure
    )
  }
  has_analysis_tuning_control <- !is.null(options$learner_tuning_folds) ||
    !is.null(options$learner_seed) || !is.null(options$gbm_cv_max_trees)
  if (!is_analysis_tuned && has_analysis_tuning_control) {
    sb_abort(
      paste0(
        "`learner_tuning_folds`, `learner_seed`, and `gbm_cv_max_trees` ",
        "apply only to an analysis-fold tuned survival learner."
      ),
      "stormbreaker_irrelevant_learner_tuning_parameter"
    )
  }
  if (is_analysis_tuned) {
    if (is.null(options$learner_tuning_folds)) {
      options$learner_tuning_folds <- if (identical(learner, "plsrcox_cv")) {
        5L
      } else {
        10L
      }
    }
    if (is.null(options$learner_seed)) {
      options$learner_seed <- 1L
    }
    if (!identical(learner, "gbm_cox_cv") &&
      !is.null(options$gbm_cv_max_trees)) {
      sb_abort(
        "`gbm_cv_max_trees` applies only to `gbm_cox_cv`.",
        "stormbreaker_irrelevant_learner_tuning_parameter"
      )
    }
    learner_params <- switch(learner,
      coxboost_cv = sb_coxboost_cv_params(
        inner_folds = options$learner_tuning_folds,
        seed = options$learner_seed
      ),
      plsrcox_cv = sb_plsrcox_cv_params(
        inner_folds = options$learner_tuning_folds,
        seed = options$learner_seed
      ),
      superpc_cv = sb_superpc_cv_params(
        inner_folds = options$learner_tuning_folds,
        seed = options$learner_seed
      ),
      gbm_cox_cv = sb_gbm_cox_cv_params(
        inner_folds = options$learner_tuning_folds,
        seed = options$learner_seed,
        maximum_trees = if (is.null(options$gbm_cv_max_trees)) {
          10000L
        } else {
          options$gbm_cv_max_trees
        }
      )
    )
  }
  learner_params
}

sb_pipeline_tree_params <- function(learner, learner_params, options) {
  is_catboost <- learner %in% c(
    "catboost_binary",
    "catboost_cox",
    "catboost_aft_survival"
  )
  is_lightgbm <- learner %in% sb_lightgbm_binary_learner_ids()
  is_gbm_fixed <- identical(learner, "gbm_cox")
  is_ranger_survival <- learner %in% sb_ranger_survival_learner_ids()
  is_ranger <- learner %in% sb_ranger_learner_ids()
  if (identical(learner, "gaussian_naive_bayes")) {
    if (is.null(options$naive_bayes_var_smoothing)) {
      options$naive_bayes_var_smoothing <- 1e-9
    }
    valid_smoothing <- is.numeric(options$naive_bayes_var_smoothing) &&
      length(options$naive_bayes_var_smoothing) == 1L &&
      !is.na(options$naive_bayes_var_smoothing) &&
      is.finite(options$naive_bayes_var_smoothing) &&
      options$naive_bayes_var_smoothing >= 0
    if (!valid_smoothing) {
      sb_abort(
        paste0(
          "`naive_bayes_var_smoothing` must be one finite, ",
          "non-negative number."
        ),
        "stormbreaker_bad_naive_bayes_var_smoothing"
      )
    }
    learner_params <- list(
      var_smoothing = as.numeric(options$naive_bayes_var_smoothing),
      variance_estimator = "population",
      prior = "empirical"
    )
  } else if (!is.null(options$naive_bayes_var_smoothing)) {
    sb_abort(
      paste0(
        "`naive_bayes_var_smoothing` applies only when ",
        "Gaussian Naive Bayes is requested."
      ),
      "stormbreaker_irrelevant_naive_bayes_parameter"
    )
  }
  lightgbm_arguments <- list(
    nrounds = options$lightgbm_nrounds,
    learning_rate = options$lightgbm_learning_rate,
    num_leaves = options$lightgbm_num_leaves,
    min_data_in_leaf = options$lightgbm_min_data_in_leaf,
    seed = options$lightgbm_seed
  )
  if (is_lightgbm) {
    lightgbm_defaults <- list(
      nrounds = 100L,
      learning_rate = 0.1,
      num_leaves = 31L,
      min_data_in_leaf = 20L,
      seed = 1L
    )
    for (name in names(lightgbm_arguments)) {
      if (is.null(lightgbm_arguments[[name]])) {
        lightgbm_arguments[[name]] <- lightgbm_defaults[[name]]
      }
    }
    learner_params <- do.call(
      sb_lightgbm_binary_params,
      c(lightgbm_arguments, list(learner = learner))
    )
  } else if (any(!vapply(lightgbm_arguments, is.null, logical(1)))) {
    sb_abort(
      paste0(
        "LightGBM parameters apply only to fixed LightGBM learners."
      ),
      "stormbreaker_irrelevant_lightgbm_parameter"
    )
  }
  catboost_arguments <- list(
    iterations = options$catboost_iterations,
    learning_rate = options$catboost_learning_rate,
    depth = options$catboost_depth,
    l2_leaf_reg = options$catboost_l2_leaf_reg,
    seed = options$catboost_seed
  )
  if (is_catboost) {
    catboost_defaults <- list(
      iterations = 100L,
      learning_rate = 0.1,
      depth = 6L,
      l2_leaf_reg = 3,
      seed = 1L
    )
    for (name in names(catboost_arguments)) {
      if (is.null(catboost_arguments[[name]])) {
        catboost_arguments[[name]] <- catboost_defaults[[name]]
      }
    }
    parameter_constructor <- if (identical(learner, "catboost_cox")) {
      sb_catboost_cox_params
    } else if (identical(learner, "catboost_aft_survival")) {
      sb_catboost_aft_survival_params
    } else {
      sb_catboost_binary_params
    }
    learner_params <- do.call(parameter_constructor, catboost_arguments)
  } else if (any(!vapply(catboost_arguments, is.null, logical(1)))) {
    sb_abort(
      paste0(
        "CatBoost parameters apply only to a fixed CatBoost learner."
      ),
      "stormbreaker_irrelevant_catboost_parameter"
    )
  }
  gbm_arguments <- list(
    n_trees = options$gbm_n_trees,
    interaction_depth = options$gbm_interaction_depth,
    min_obs_node = options$gbm_min_obs_node,
    shrinkage = options$gbm_shrinkage,
    bag_fraction = options$gbm_bag_fraction,
    seed = options$gbm_seed
  )
  if (is_gbm_fixed) {
    gbm_defaults <- list(
      n_trees = 100L,
      interaction_depth = 1L,
      min_obs_node = 10L,
      shrinkage = 0.1,
      bag_fraction = 1,
      seed = 1L
    )
    for (name in names(gbm_arguments)) {
      if (is.null(gbm_arguments[[name]])) {
        gbm_arguments[[name]] <- gbm_defaults[[name]]
      }
    }
    learner_params <- do.call(sb_gbm_cox_params, gbm_arguments)
  } else if (any(!vapply(gbm_arguments, is.null, logical(1)))) {
    sb_abort(
      "GBM parameters apply only when `learner = 'gbm_cox'`.",
      "stormbreaker_irrelevant_gbm_parameter"
    )
  }
  ranger_arguments <- list(
    num_trees = options$ranger_num_trees,
    mtry = options$ranger_mtry,
    min_node_size = options$ranger_min_node_size,
    seed = options$ranger_seed
  )
  if (is_ranger) {
    ranger_defaults <- list(
      num_trees = 500L,
      min_node_size = if (is_ranger_survival) 3L else 10L,
      seed = 1L
    )
    for (name in names(ranger_defaults)) {
      if (is.null(ranger_arguments[[name]])) {
        ranger_arguments[[name]] <- ranger_defaults[[name]]
      }
    }
    parameter_constructor <- if (is_ranger_survival) {
      sb_ranger_survival_params
    } else {
      sb_ranger_binary_params
    }
    ranger_arguments$learner <- learner
    learner_params <- do.call(parameter_constructor, ranger_arguments)
  } else if (any(!vapply(ranger_arguments, is.null, logical(1)))) {
    sb_abort(
      paste0(
        "ranger parameters apply only to a fixed ranger learner."
      ),
      "stormbreaker_irrelevant_ranger_parameter"
    )
  }
  learner_params
}

sb_pipeline_fixed_params <- function(learner, learner_params, options) {
  is_coxboost <- identical(learner, "coxboost")
  is_plsrcox <- identical(learner, "plsrcox")
  is_superpc <- identical(learner, "superpc")
  is_survivalsvm <- learner %in% sb_survivalsvm_learner_ids()
  is_classical_binary_boosting <- learner %in%
    sb_classical_binary_boosting_learner_ids()
  is_classical_binary_fixed <- learner %in%
    sb_classical_binary_learner_ids()
  is_specialized_binary_linear <- learner %in%
    sb_specialized_binary_linear_learner_ids()
  is_omics_specialist_binary <- learner %in%
    sb_omics_specialist_binary_learner_ids()
  is_cancerclass <- identical(learner, "cancerclass_binary")
  is_xgboost <- learner %in% sb_xgboost_learner_ids()
  is_additional_survival <- learner %in%
    sb_additional_survival_learner_ids()
  is_flexsurv_parametric <- learner %in%
    sb_flexsurv_parametric_learner_ids()
  is_flexsurvspline <- learner %in% sb_flexsurvspline_learner_ids()
  is_mtlr_survival <- learner %in% sb_mtlr_survival_learner_ids()
  is_pamr_survival <- learner %in% sb_pamr_survival_learner_ids()
  is_timereg_survival <- learner %in% sb_timereg_survival_learner_ids()
  is_mboost_cox <- learner %in% sb_mboost_cox_learner_ids()
  is_mboost_aft <- learner %in% sb_mboost_aft_learner_ids()
  is_mboost_ranking <- learner %in% sb_mboost_ranking_learner_ids()
  is_abess <- learner %in% sb_abess_learner_ids()
  is_ncvreg <- learner %in% sb_ncvreg_learner_ids()
  is_rlt <- learner %in% sb_rlt_learner_ids()
  is_coxphf <- learner %in% sb_coxphf_learner_ids()
  is_bagged_knn_survival <- learner %in%
    sb_bagged_knn_survival_learner_ids()
  is_ahaz_survival <- learner %in% sb_ahaz_survival_learner_ids()
  is_bujar_linear <- learner %in% sb_bujar_linear_learner_ids()
  is_rule_aft <- learner %in% sb_rule_aft_learner_ids()
  is_wsrf <- learner %in% sb_wsrf_binary_learner_ids()
  is_slope_randomglm_binary <- learner %in%
    sb_slope_randomglm_binary_learner_ids()
  is_mgcv_gam <- learner %in% sb_mgcv_gam_learner_ids()
  is_mclust_da <- learner %in% sb_mclust_da_learner_ids()
  is_mixomics_da <- learner %in% sb_mixomics_da_learner_ids()
  is_hdda_binary <- learner %in% sb_hdda_binary_learner_ids()
  classical_arguments <- list(
    coxboost_stepno = options$coxboost_stepno,
    coxboost_penalty = options$coxboost_penalty,
    plsrcox_components = options$plsrcox_components,
    superpc_threshold = options$superpc_threshold,
    superpc_components = options$superpc_components,
    survivalsvm_gamma_mu = options$survivalsvm_gamma_mu
  )
  relevant_classical_arguments <- if (is_survivalsvm) {
    "survivalsvm_gamma_mu"
  } else {
    switch(learner,
      coxboost = c("coxboost_stepno", "coxboost_penalty"),
      plsrcox = "plsrcox_components",
      superpc = c("superpc_threshold", "superpc_components"),
      character()
    )
  }
  irrelevant_classical_arguments <- setdiff(
    names(classical_arguments),
    relevant_classical_arguments
  )
  has_irrelevant_classical_argument <- any(!vapply(
    classical_arguments[irrelevant_classical_arguments],
    is.null,
    logical(1)
  ))
  if (has_irrelevant_classical_argument) {
    sb_abort(
      paste0(
        "CoxBoost, plsRcox, SuperPC, and survival-SVM parameters apply ",
        "only to their matching learner."
      ),
      "stormbreaker_irrelevant_classical_survival_parameter"
    )
  }
  if (is_coxboost) {
    learner_params <- sb_coxboost_params(
      stepno = if (is.null(options$coxboost_stepno)) 100L else options$coxboost_stepno,
      penalty = options$coxboost_penalty
    )
  } else if (is_plsrcox) {
    learner_params <- sb_plsrcox_params(
      components = if (is.null(options$plsrcox_components)) {
        2L
      } else {
        options$plsrcox_components
      }
    )
  } else if (is_superpc) {
    learner_params <- sb_superpc_params(
      threshold = if (is.null(options$superpc_threshold)) 1 else options$superpc_threshold,
      components = if (is.null(options$superpc_components)) {
        1L
      } else {
        options$superpc_components
      }
    )
  } else if (is_survivalsvm) {
    learner_params <- sb_survivalsvm_params(
      gamma_mu = if (is.null(options$survivalsvm_gamma_mu)) {
        1
      } else {
        options$survivalsvm_gamma_mu
      },
      learner = learner
    )
  }

  classical_binary_boosting_arguments <- list(
    adaboost_iterations = options$adaboost_iterations,
    logitboost_iterations = options$logitboost_iterations
  )
  relevant_binary_boosting_argument <- switch(learner,
    adaboost_binary = "adaboost_iterations",
    logitboost_binary = "logitboost_iterations",
    character()
  )
  irrelevant_binary_boosting_arguments <- setdiff(
    names(classical_binary_boosting_arguments),
    relevant_binary_boosting_argument
  )
  if (any(!vapply(
    classical_binary_boosting_arguments[
      irrelevant_binary_boosting_arguments
    ],
    is.null,
    logical(1)
  ))) {
    sb_abort(
      paste0(
        "AdaBoost and LogitBoost iteration parameters apply only to their ",
        "matching learner."
      ),
      "stormbreaker_irrelevant_classical_binary_boosting_parameter"
    )
  }
  if (is_classical_binary_boosting) {
    learner_params <- sb_classical_binary_boosting_params(
      learner,
      classical_binary_boosting_arguments[[
        relevant_binary_boosting_argument
      ]]
    )
  }
  if (is_classical_binary_fixed) {
    learner_params <- sb_classical_binary_params(learner)
  }
  if (!identical(learner, "stepglm_binary") &&
    !is.null(options$stepglm_direction)) {
    sb_abort(
      paste0(
        "`stepglm_direction` applies only when ",
        "`learner = 'stepglm_binary'`."
      ),
      "stormbreaker_irrelevant_stepglm_parameter"
    )
  }
  if (is_specialized_binary_linear) {
    learner_params <- sb_specialized_binary_linear_params(
      learner,
      stepglm_direction = options$stepglm_direction
    )
  }
  if (is_omics_specialist_binary) {
    learner_params <- sb_omics_specialist_binary_params(learner)
  }
  if (is_cancerclass) {
    learner_params <- sb_cancerclass_binary_params()
  }
  if (is_xgboost) {
    learner_params <- sb_xgboost_params(learner)
  }
  if (is_additional_survival) {
    learner_params <- sb_additional_survival_params(learner)
  }
  if (is_flexsurv_parametric) {
    learner_params <- sb_flexsurv_parametric_params(learner)
  }
  if (is_flexsurvspline) {
    learner_params <- sb_flexsurvspline_params(learner)
  }
  if (is_mtlr_survival) {
    learner_params <- sb_mtlr_survival_params(learner)
  }
  if (is_pamr_survival) {
    learner_params <- sb_pamr_survival_params(learner)
  }
  if (is_timereg_survival) {
    learner_params <- sb_timereg_survival_params(learner)
  }
  if (is_mboost_cox) {
    learner_params <- sb_mboost_cox_params(learner)
  }
  if (is_mboost_aft) {
    learner_params <- sb_mboost_aft_params(learner)
  }
  if (is_mboost_ranking) {
    learner_params <- sb_mboost_ranking_params(learner)
  }
  if (is_abess) {
    learner_params <- sb_abess_learner_params(learner)
  }
  if (is_ncvreg) {
    learner_params <- sb_ncvreg_learner_params(learner)
  }
  if (is_rlt) {
    learner_params <- sb_rlt_learner_params(learner)
  }
  if (is_coxphf) {
    learner_params <- sb_coxphf_learner_params(learner)
  }
  if (learner %in% sb_akritas_learner_ids()) learner_params <- sb_akritas_params()
  if (is_bagged_knn_survival) {
    learner_params <- sb_bagged_knn_survival_params(learner)
  }
  if (is_ahaz_survival) {
    learner_params <- sb_ahaz_learner_params(learner)
  }
  if (is_bujar_linear) {
    learner_params <- sb_bujar_linear_params(learner)
  }
  if (is_rule_aft) {
    learner_params <- sb_rule_aft_params(learner)
  }
  if (is_wsrf) {
    learner_params <- sb_wsrf_binary_params(learner)
  }
  if (is_slope_randomglm_binary) {
    learner_params <- sb_slope_randomglm_binary_params(learner)
  }
  if (is_mgcv_gam) {
    learner_params <- sb_mgcv_gam_params(learner)
  }
  if (is_mclust_da) {
    learner_params <- sb_mclust_da_params(learner)
  }
  if (is_mixomics_da) {
    learner_params <- sb_mixomics_da_params(learner)
  }
  if (is_hdda_binary) {
    learner_params <- sb_hdda_binary_params(learner)
  }
  learner_params
}
