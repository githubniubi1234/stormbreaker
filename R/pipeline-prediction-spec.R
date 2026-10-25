sb_pipeline_prediction_spec <- function(pipeline_task_type, learner, learner_params) {
  is_gbm <- learner %in% c("gbm_cox", "gbm_cox_cv")
  is_survivalsvm_kernel <- learner %in%
    sb_survivalsvm_kernel_learner_ids()
  is_xgboost_cox <- learner %in% sb_xgboost_cox_learner_ids()
  is_xgboost_aft <- learner %in% sb_xgboost_aft_learner_ids()
  is_flexsurv_parametric <- learner %in%
    sb_flexsurv_parametric_learner_ids()
  is_flexsurvspline <- learner %in% sb_flexsurvspline_learner_ids()
  is_mtlr_survival <- learner %in% sb_mtlr_survival_learner_ids()
  is_pamr_survival <- learner %in% sb_pamr_survival_learner_ids()
  is_mboost_cox <- learner %in% sb_mboost_cox_learner_ids()
  is_mboost_aft <- learner %in% sb_mboost_aft_learner_ids()
  is_mboost_ranking <- learner %in% sb_mboost_ranking_learner_ids()
  is_bagged_knn_survival <- learner %in%
    sb_bagged_knn_survival_learner_ids()
  is_ahaz_survival <- learner %in% sb_ahaz_survival_learner_ids()
  is_bujar_linear <- learner %in% sb_bujar_linear_learner_ids()
  is_rule_aft <- learner %in% sb_rule_aft_learner_ids()
  is_ranger_survival <- learner %in% sb_ranger_survival_learner_ids()
  pipeline <- list()
  if (identical(pipeline_task_type, "binary")) {
    pipeline$prediction_type <- "probability"
    pipeline$positive_class <- 1L
    pipeline$prediction_direction <- "higher_is_more_likely_positive"
  } else if (is_survivalsvm_kernel) {
    pipeline$risk_type <- "calibrated_survival_svm_kernel_score"
    pipeline$risk_reference <- sb_survivalsvm_kernel_reference(
      learner_params$kernel
    )
    pipeline$risk_formula <-
      "negative_native_kernel_score_times_analysis_cox_coefficient"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "catboost_cox")) {
    pipeline$risk_type <- "cox_raw_margin"
    pipeline$risk_reference <- "catboost_native"
    pipeline$risk_formula <- "catboost_raw_formula_value"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "catboost_aft_survival")) {
    pipeline$risk_type <- "aft_log_time"
    pipeline$risk_reference <- "catboost_aft_normal_native"
    pipeline$risk_formula <- "negative_catboost_aft_log_time_margin"
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (is_xgboost_cox) {
    pipeline$risk_type <- "cox_raw_margin"
    pipeline$risk_reference <- "xgboost_native"
    pipeline$risk_formula <- "xgboost_raw_margin"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (is_xgboost_aft) {
    pipeline$risk_type <- "aft_log_time"
    pipeline$risk_reference <- "xgboost_aft_normal_native"
    pipeline$risk_formula <- "negative_xgboost_aft_log_time_margin"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (is_gbm) {
    pipeline$risk_type <- "cox_log_hazard"
    pipeline$risk_reference <- "gbm_native"
    pipeline$risk_formula <- "gbm_link_prediction"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (is_mboost_cox) {
    pipeline$risk_type <- "cox_log_hazard"
    pipeline$risk_reference <- "mboost_native"
    pipeline$risk_formula <- "mboost_coxph_link_prediction"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (is_mboost_aft) {
    pipeline$risk_type <- "aft_log_time"
    pipeline$risk_reference <- paste0(
      "mboost_", learner_params$algorithm, "_",
      learner_params$distribution, "_aft_native"
    )
    pipeline$risk_formula <- "negative_mboost_aft_log_time_link"
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (is_mboost_ranking) {
    pipeline$risk_type <- if (identical(learner_params$family, "Gehan")) {
      "gehan_aft_rank_score"
    } else {
      "uno_ipcw_concordance_score"
    }
    pipeline$risk_reference <- paste0(
      "mboost_", learner_params$algorithm, "_",
      tolower(learner_params$family), "_native"
    )
    pipeline$risk_formula <- "negative_mboost_ranking_link"
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (identical(learner, "gam_cox_reml")) {
    pipeline$risk_type <- "cox_log_hazard"
    pipeline$risk_reference <- "mgcv_native"
    pipeline$risk_formula <- "mgcv_cox_ph_link_prediction"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (learner %in% sb_rlt_survival_learner_ids()) {
    pipeline$risk_type <- "reinforcement_survival_forest_rmst"
    pipeline$risk_reference <- "RLT_native_discrete_hazard"
    pipeline$risk_formula <- "negative_step_RMST_over_native_event_grid"
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (learner %in% sb_akritas_learner_ids()) {
    pipeline$risk_type <- "akritas_conditional_survival_rmst"
    pipeline$risk_reference <- "maximum_analysis_followup"
    pipeline$risk_formula <- "negative_step_RMST_over_analysis_time_grid"
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (is_bagged_knn_survival) {
    pipeline$risk_type <- "bagged_knn_survival_rmst"
    pipeline$risk_reference <- "bnnSurvival_native_training_time_grid"
    pipeline$risk_formula <-
      "negative_step_RMST_over_native_training_time_grid"
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (is_ahaz_survival) {
    pipeline$risk_type <- "additive_hazard_linear_predictor"
    pipeline$risk_reference <- "ahaz_native"
    pipeline$risk_formula <- "x_beta"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (is_bujar_linear) {
    pipeline$risk_type <- "buckley_james_log_time"
    pipeline$risk_reference <- "bujar_native"
    pipeline$risk_formula <- "negative_predicted_log_time"
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (is_rule_aft) {
    risk <- switch(learner,
      bart_aft_survival = list(
        type = "aft_log_time",
        reference = "BART_abart_normal_native",
        formula = "negative_posterior_mean_log_time",
        direction = "higher_is_shorter_predicted_survival"
      ),
      pre_survival_rule_ensemble = list(
        type = "cox_rule_ensemble_link",
        reference = "pre_native_rule_ensemble",
        formula = "pre_cox_link_prediction",
        direction = "higher_is_higher_hazard"
      ),
      penaft_gehan_elastic_net = list(
        type = "semiparametric_gehan_aft_log_time",
        reference = "penAFT_native",
        formula = "negative_centered_elastic_net_log_time_prediction",
        direction = "higher_is_shorter_predicted_survival"
      )
    )
    pipeline$risk_type <- risk$type
    pipeline$risk_reference <- risk$reference
    pipeline$risk_formula <- risk$formula
    pipeline$risk_direction <- risk$direction
  } else if (is_ranger_survival) {
    pipeline$risk_type <- "survival_forest_mortality"
    pipeline$risk_reference <- "ranger_native"
    pipeline$risk_formula <- "row_sums_cumulative_hazard"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "aalen_additive_survival")) {
    pipeline$risk_type <- "aalen_cumulative_hazard"
    pipeline$risk_reference <- "timereg_native_last_analysis_event_time"
    pipeline$risk_formula <- "cumulative_hazard_at_last_analysis_event"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "prop_odds_survival")) {
    pipeline$risk_type <- "proportional_odds_log_failure_odds"
    pipeline$risk_reference <- "zero"
    pipeline$risk_formula <- "x_beta"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (learner %in% sb_survreg_learner_ids()) {
    pipeline$risk_type <- "aft_log_time"
    pipeline$risk_reference <- paste0(
      learner_params$distribution,
      "_aft_native"
    )
    pipeline$risk_formula <- "negative_log_time_location"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (is_flexsurv_parametric) {
    pipeline$risk_type <- if (identical(
      learner_params$location_parameter,
      "mu"
    )) {
      "aft_log_time"
    } else if (identical(learner_params$distribution, "gompertz")) {
      "proportional_hazards_log_rate"
    } else {
      "accelerated_failure_time_log_rate"
    }
    pipeline$risk_reference <- sb_flexsurv_parametric_risk_reference(learner)
    pipeline$risk_formula <- sb_flexsurv_parametric_risk_formula(learner)
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (is_flexsurvspline) {
    pipeline$risk_type <- switch(learner_params$scale,
      hazard = "flexible_parametric_log_cumulative_hazard_location",
      odds = "flexible_parametric_log_cumulative_odds_location",
      normal = "flexible_parametric_normal_survival_location"
    )
    pipeline$risk_reference <- sb_flexsurvspline_risk_reference(learner)
    pipeline$risk_formula <- sb_flexsurvspline_risk_formula(learner)
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (is_mtlr_survival) {
    pipeline$risk_type <- "multitask_logistic_mean_survival"
    pipeline$risk_reference <- "MTLR_native_time_grid"
    pipeline$risk_formula <- "negative_native_mean_predicted_survival_time"
    pipeline$risk_direction <- "higher_is_shorter_predicted_survival"
  } else if (is_pamr_survival) {
    pipeline$risk_type <- "nearest_shrunken_centroid_survival_group_posterior"
    pipeline$risk_reference <- "pamr_native_survival_groups"
    pipeline$risk_formula <- "native_posterior_for_early_survival_group"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "rpart_survival")) {
    pipeline$risk_type <- "survival_tree_log_rate"
    pipeline$risk_reference <- "rpart_native"
    pipeline$risk_formula <- "log_relative_event_rate"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "rsf_survival")) {
    pipeline$risk_type <- "survival_forest_mortality"
    pipeline$risk_reference <- "randomForestSRC_native"
    pipeline$risk_formula <- "mortality_sum_cumulative_hazard"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "ctree_survival")) {
    pipeline$risk_type <- "conditional_survival_cumulative_hazard"
    pipeline$risk_reference <- "partykit_native"
    pipeline$risk_formula <- "terminal_distribution_cumulative_hazard"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "cforest_survival")) {
    pipeline$risk_type <- "conditional_survival_cumulative_hazard"
    pipeline$risk_reference <- "partykit_native"
    pipeline$risk_formula <- "ensemble_distribution_cumulative_hazard"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "bart_survival")) {
    pipeline$risk_type <- "posterior_negative_restricted_mean_survival"
    pipeline$risk_reference <- "BART_analysis_fold_time_grid"
    pipeline$risk_formula <- "negative_step_RMST_from_posterior_mean_survival"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else if (identical(learner, "aorsf_survival")) {
    pipeline$risk_type <- "oblique_survival_forest_mortality"
    pipeline$risk_reference <- "aorsf_native"
    pipeline$risk_formula <- "native_mortality_over_event_times"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  } else {
    pipeline$risk_type <- "linear_predictor"
    pipeline$risk_reference <- "zero"
    pipeline$risk_formula <- "x_beta"
    pipeline$risk_direction <- "higher_is_higher_hazard"
  }
  pipeline
}
