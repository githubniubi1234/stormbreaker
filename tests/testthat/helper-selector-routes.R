expected_fixed_learner_selectors <- function(task) {
  if (task == "binary") {
    c("none", "univ_binary_auc_top_k", "abess_binary_gic_nonzero",
      "pamr_binary_nonzero",
      "lasso_binary_nonzero", "stepglm_binary_aic", "glmboost_binary_nonzero",
      "rf_binary_importance", "enet_binary_nonzero", "gbm_binary_importance",
      "svm_rfe_binary", "gamboost_binary_nonzero", "xgboost_binary_importance",
      "plsrglm_binary_weights", "boruta_binary_roughfix", "varselrf_binary_oob_1se",
      "sis_binary_nlogn", "vsurf_binary_prediction", "mrmr_binary_nlogn",
      "jmi_binary_nlogn", "cmim_binary_nlogn", "disr_binary_nlogn",
      "jmim_binary_nlogn", "lpc_binary_nlogn", "stabsel_binary_lasso",
      "fixed_features")
  } else {
    c("none", "univ_cox_top_k", "stepcox", "abess_cox_gic_nonzero",
      "pamr_survival_nonzero",
      "ridge_cox_weights", "plsrcox_weights", "svm_survival_weights",
      "gbm_cox_importance", "glmboost_cox_nonzero", "superpc_survival_scores",
      "gamboost_cox_nonzero", "xgboost_cox_importance", "boruta_survival_roughfix",
      "lpc_cox_nlogn", "stabsel_cox_lasso", "lasso_nonzero", "enet_cox_nonzero",
      "coxboost_nonzero", "rsf_minimal_depth", "fixed_features")
  }
}

expect_route_dependencies <- function(rows, learner_packages) {
  selector_packages <- c(
    none = "", univ_binary_auc_top_k = "", univ_cox_top_k = "",
    stepcox = "", fixed_features = "", stepglm_binary_aic = "",
    abess_binary_gic_nonzero = "abess", abess_cox_gic_nonzero = "abess",
    pamr_binary_nonzero = "pamr", pamr_survival_nonzero = "pamr",
    lasso_binary_nonzero = "glmnet", enet_binary_nonzero = "glmnet",
    ridge_cox_weights = "glmnet", lasso_nonzero = "glmnet", enet_cox_nonzero = "glmnet",
    glmboost_binary_nonzero = "mboost", gamboost_binary_nonzero = "mboost",
    glmboost_cox_nonzero = "mboost", gamboost_cox_nonzero = "mboost",
    rf_binary_importance = "randomForest", gbm_binary_importance = "gbm",
    gbm_cox_importance = "gbm", svm_rfe_binary = "e1071",
    svm_survival_weights = "survivalsvm", plsrcox_weights = "plsRcox",
    plsrglm_binary_weights = "plsRglm", superpc_survival_scores = "superpc",
    xgboost_binary_importance = "xgboost", xgboost_cox_importance = "xgboost",
    boruta_binary_roughfix = "Boruta", boruta_survival_roughfix = "Boruta",
    varselrf_binary_oob_1se = "varSelRF", sis_binary_nlogn = "SIS",
    vsurf_binary_prediction = "VSURF", mrmr_binary_nlogn = "praznik",
    jmi_binary_nlogn = "praznik", cmim_binary_nlogn = "praznik",
    disr_binary_nlogn = "praznik", jmim_binary_nlogn = "praznik",
    lpc_binary_nlogn = "lpc", lpc_cox_nlogn = "lpc",
    stabsel_binary_lasso = "stabs + glmnet", stabsel_cox_lasso = "stabs + glmnet",
    coxboost_nonzero = "CoxBoost", rsf_minimal_depth = "randomForestSRC"
  )
  for (i in seq_len(nrow(rows))) {
    expected <- paste(selector_packages[[rows$selector[i]]],
      learner_packages[[rows$learner[i]]], sep = " + ")
    packages <- strsplit(expected, " + ", fixed = TRUE)[[1L]]
    expect_setequal(strsplit(rows$optional_dependency[i], " + ", fixed = TRUE)[[1L]],
      unique(packages[nzchar(packages)]))
  }
}
