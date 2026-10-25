sb_display_learner <- function(learner, parameters = list(), model_id = NULL) {
  learner <- sb_canonical_learner_id(learner)

  if (is.null(learner) || !length(learner) || is.na(learner)) return("Model")
  if (learner %in% c("glmnet_logistic", "glmnet_cox")) {
    alpha <- parameters$alpha
    if (is.null(alpha)) alpha <- if (!is.null(model_id) && grepl("^ridge_", model_id)) 0 else
      if (!is.null(model_id) && grepl("^elastic_net_", model_id)) 0.5 else 1
    penalty <- if (identical(as.numeric(alpha), 0)) "Ridge" else
      if (identical(as.numeric(alpha), 1)) "Lasso" else "Elastic net"
    return(paste(penalty, if (learner == "glmnet_cox") "Cox" else "logistic"))
  }
  known <- c(
    logistic = "Logistic regression", coxph = "Cox regression",
    gaussian_naive_bayes = "Gaussian naive Bayes",
    ranger_binary = "Random forest (ranger)",
    ranger_survival = "Random survival forest (ranger)",
    random_forest_binary = "Random forest (randomForest)",
    rsf_survival = "Random survival forest (randomForestSRC)",
    lda_binary = "LDA", qda_binary = "QDA", fda_binary = "FDA",
    pda_binary = "Penalized FDA", mda_binary = "MDA",
    slda_binary = "Stabilized LDA", sda_binary = "Shrinkage DA",
    earth_binary = "MARS", decision_tree_binary = "Decision tree (rpart)",
    tree_binary = "Decision tree (tree)", rpart_survival = "Survival tree (rpart)",
    knn_binary = "k-nearest neighbours", kernel_weighted_knn_binary = "Weighted kNN",
    kernel_naive_bayes_binary = "Kernel naive Bayes",
    nearest_centroid_binary = "Nearest centroid",
    bernoulli_naive_bayes_binary = "Bernoulli naive Bayes",
    gaussian_process_binary = "Gaussian process",
    nnet_binary = "Neural network (nnet)",
    neuralnet_multilayer_binary = "Multilayer neural network",
    aorsf_binary = "Oblique random forest (aorsf)",
    aorsf_survival = "Oblique survival forest (aorsf)",
    brglm2_logistic = "Bias-reduced logistic",
    bayesglm_logistic_cauchy = "Bayesian logistic (Cauchy)",
    robust_logistic_mqle = "Robust logistic (Mqle)",
    coxphf_firth = "Firth Cox regression",
    hda_binary = "High-dimensional discriminant analysis",
    dwd_binary = "Distance-weighted discrimination",
    ktsp_binary = "k-top scoring pairs",
    rotation_forest_binary = "Rotation forest",
    rrf_binary = "Regularized random forest",
    random_ferns_binary = "Random ferns",
    c50_rules_binary = "C5.0 rules", c50_binary = "C5.0",
    bagging_binary = "Bagged trees",
    bagged_knn_survival = "Bagged survival kNN",
    akritas_survival = "Akritas",
    survreg_weibull = "Weibull AFT", survreg_exponential = "Exponential AFT",
    survreg_lognormal = "Log-normal AFT", survreg_loglogistic = "Log-logistic AFT",
    mtlr_survival = "MTLR", ctree_binary = "Conditional inference tree",
    ctree_survival = "Conditional inference survival tree",
    cforest_binary = "Conditional inference forest",
    cforest_survival = "Conditional inference survival forest"
  )
  if (learner %in% names(known)) return(unname(known[[learner]]))
  tokens <- strsplit(learner, "_", fixed = TRUE)[[1L]]
  tokens <- tokens[!tokens %in% c("fixed", "binary", "survival")]
  dictionary <- c(
    xgboost = "XGBoost", lightgbm = "LightGBM", catboost = "CatBoost",
    gbm = "GBM", cox = "Cox", aft = "AFT", dart = "DART", goss = "GOSS",
    rf = "RF", linear = "linear", svm = "SVM", lssvm = "LS-SVM",
    radial = "RBF", rbf = "RBF", polynomial = "polynomial", sigmoid = "sigmoid",
    survivalsvm = "Survival SVM", vanbelle1 = "Van Belle 1", vanbelle2 = "Van Belle 2",
    ranger = "ranger", extratrees = "Extra Trees", hellinger = "Hellinger",
    cindex = "C-index", maxstat = "maxstat", rlt = "RLT",
    adaboost = "AdaBoost", logitboost = "LogitBoost", bst = "BST",
    rda = "RDA", evtree = "Evolutionary tree", coxboost = "CoxBoost",
    plsrcox = "PLS-Cox", plsrglm = "PLS-GLM", superpc = "SuperPC",
    gamboost = "GamBoost", glmboost = "GLMBoost", blackboost = "BlackBoost",
    gehan = "Gehan", unoc = "Uno C", ncvreg = "ncvreg", mcp = "MCP",
    scad = "SCAD", bic = "BIC", gic = "GIC", abess = "ABESS",
    gam = "GAM", reml = "REML", bart = "BART", pamr = "PAM",
    flexsurv = "flexsurv", flexsurvspline = "flexsurvspline",
    liblinear = "LIBLINEAR", l1 = "L1", l2 = "L2", svc = "SVC",
    l2loss = "L2 loss", l1loss = "L1 loss", lvq1 = "LVQ1",
    olvq1 = "OLVQ1", lvq21 = "LVQ2.1", lvq3 = "LVQ3",
    stepglm = "Stepwise logistic", pre = "PRE", ahaz = "ahaz",
    aalen = "Aalen", bujar = "Bujar", penaft = "penAFT",
    cv = "CV-tuned", lasso = "Lasso", sscad = "S-SCAD",
    c50 = "C5.0"
  )
  matched <- match(tokens, names(dictionary))
  tokens[!is.na(matched)] <- unname(dictionary[matched[!is.na(matched)]])
  paste(tokens[nzchar(tokens)], collapse = " ")
}

sb_display_selector <- function(selector, parameters = list()) {
  if (identical(selector, "ordered_chain")) return(paste(vapply(parameters$stages,
    function(p) sb_display_selector(p$selector, p$selector_params), character(1)), collapse = " -> "))
  if (is.null(selector) || !length(selector) || is.na(selector)) return("Not recorded")
  known <- c(
    lasso_binary_nonzero = "LASSO", stepglm_binary_aic = "Stepglm",
    enet_binary_nonzero = "Enet", gbm_binary_importance = "GBM",
    glmboost_binary_nonzero = "glmBoost", rf_binary_importance = "RF importance",
    none = "None", univ_cox_top_k = "Univariate Cox", univ_binary_auc_top_k = "Univariate AUC",
    stepcox = "Stepwise Cox", abess_binary_gic_nonzero = "ABESS",
    abess_cox_gic_nonzero = "ABESS", boruta_binary_roughfix = "Boruta",
    boruta_survival_roughfix = "Boruta", varselrf_binary_oob_1se = "varSelRF",
    sis_binary_nlogn = "SIS", vsurf_binary_prediction = "VSURF",
    mrmr_binary_nlogn = "mRMR", jmi_binary_nlogn = "JMI", cmim_binary_nlogn = "CMIM",
    disr_binary_nlogn = "DISR", jmim_binary_nlogn = "JMIM",
    lpc_binary_nlogn = "LPC", lpc_cox_nlogn = "LPC",
    stabsel_binary_lasso = "Stability selection", stabsel_cox_lasso = "Stability selection",
    rsf_minimal_depth = "RSF minimal depth", coxboost_nonzero = "CoxBoost",
    lasso_nonzero = "Lasso", enet_cox_nonzero = "Enet", rsf_then_stepcox = "RSF + Stepwise Cox",
    coxboost_then_stepcox = "CoxBoost + Stepwise Cox",
    lasso_then_stepcox = "Lasso + Stepwise Cox", fixed_features = "Specified gene set"
  )
  label <- if (selector %in% names(known)) unname(known[[selector]]) else gsub("_", " ", selector)
  k <- if (!is.null(parameters$n_features)) parameters$n_features else parameters$k
  if (!is.null(k)) label <- paste0(label, " (", k, ")")
  if (!is.null(parameters$direction)) label <- paste0(label, " (", parameters$direction, ")")
  label
}

sb_display_parameters <- function(parameters) {
  if (!length(parameters)) return("")
  values <- vapply(parameters, function(value) {
    if (is.atomic(value) && length(value) <= 8L) paste(value, collapse = ", ") else
      paste0("[", length(value), " values]")
  }, character(1))
  paste(paste(names(values), values, sep = " = "), collapse = "; ")
}
