sb_pipeline_registry_data <- function() {
  paired_rows <- function(task_type, learners, selector, dependencies,
                          tuning, base_semantics, selected_semantics) {
    data.frame(
      source = rep("stormbreaker_builtin", length(learners) * 2L),
      task_type = rep(task_type, length(learners) * 2L),
      selector = rep(c("none", selector), times = length(learners)),
      learner = rep(learners, each = 2L),
      admission_status = rep("admitted", length(learners) * 2L),
      selector_parameter = rep(
        c(NA, "selector_n_features"),
        times = length(learners)
      ),
      tuning_semantics = rep(tuning, each = 2L),
      optional_dependency = rep(dependencies, each = 2L),
      parameter_semantics = as.vector(rbind(
        base_semantics,
        selected_semantics
      )),
      stringsAsFactors = FALSE
    )
  }

  survival_learners <- c(
    "coxph", "glmnet_cox", "catboost_cox",
    "catboost_aft_survival", "gbm_cox",
    sb_mboost_cox_learner_ids(),
    sb_mboost_aft_learner_ids(),
    sb_mboost_ranking_learner_ids(),
    sb_ranger_survival_learner_ids(), sb_rlt_survival_learner_ids(),
    "coxboost", "plsrcox",
    "superpc", sb_survivalsvm_learner_ids(),
    sb_xgboost_cox_learner_ids(),
    sb_xgboost_aft_learner_ids(),
    sb_survreg_learner_ids(),
    sb_flexsurv_parametric_learner_ids(),
    sb_flexsurvspline_learner_ids(),
    sb_mtlr_survival_learner_ids(),
    sb_pamr_survival_learner_ids(), "rpart_survival",
    "rsf_survival",
    "ctree_survival", "cforest_survival", "bart_survival",
    "aorsf_survival", sb_timereg_survival_learner_ids(),
    sb_ncvreg_survival_learner_ids(),
    "abess_cox_gic",
    sb_coxphf_learner_ids(), sb_bagged_knn_survival_learner_ids(),
    sb_ahaz_survival_learner_ids(), sb_bujar_linear_learner_ids(), sb_akritas_learner_ids()
  )
  survivalsvm_base <- vapply(
    sb_survivalsvm_learner_ids(),
    function(learner) {
      kernel_type <- sb_survivalsvm_kernel(learner)
      kernel <- switch(kernel_type,
        lin_kernel = "linear",
        add_kernel = "additive",
        rbf_kernel = "radial-basis-function"
      )
      kernel_rule <- if (identical(kernel_type, "rbf_kernel")) {
        paste0(
          "; analysis-median-squared-distance kernel-parameter rule ",
          "fixed before resampling"
        )
      } else {
        ""
      }
      switch(sb_survivalsvm_type(learner),
        regression = paste0(
          kernel,
          " survival-SVM regression gamma fixed before resampling",
          kernel_rule
        ),
        vanbelle1 = paste0(
          kernel,
          " survival-SVM vanbelle1 ranking gamma and makediff3 rule fixed ",
          "before resampling",
          kernel_rule
        ),
        vanbelle2 = paste0(
          kernel,
          " survival-SVM vanbelle2 ranking gamma and makediff3 rule fixed ",
          "before resampling",
          kernel_rule
        ),
        hybrid = paste0(
          kernel,
          " survival-SVM hybrid regression/ranking gamma pair and ",
          "makediff3 rule fixed before resampling",
          kernel_rule
        )
      )
    },
    character(1)
  )
  survival_base <- c(
    "fixed ordinary Cox fit",
    paste0(
      "alpha, standardization setting, lambda rule, and tuning measure fixed; ",
      "lambda selected by inner CV"
    ),
    "all learner parameters fixed before resampling",
    paste0(
      "CatBoost SurvivalAft Normal loss, unit scale, tree controls, and seed ",
      "fixed before resampling"
    ),
    "all learner parameters fixed before resampling",
    paste0(
      "100 CoxPH boosting iterations, shrinkage 0.1, and spline base ",
      "learners fixed before resampling"
    ),
    paste0(
      "100 CoxPH boosting iterations, shrinkage 0.1, centering, and linear ",
      "base learners fixed before resampling"
    ),
    paste0(
      "100 CoxPH boosting iterations, shrinkage 0.1, and tree controls ",
      "fixed before resampling"
    ),
    vapply(
      sb_mboost_aft_learner_ids(),
      sb_mboost_aft_parameter_semantics,
      character(1)
    ),
    vapply(
      sb_mboost_ranking_learner_ids(),
      sb_mboost_ranking_parameter_semantics,
      character(1)
    ),
    vapply(
      sb_ranger_survival_learner_ids(),
      sb_ranger_parameter_semantics,
      character(1)
    ),
    vapply(
      sb_rlt_survival_learner_ids(),
      sb_rlt_parameter_semantics,
      character(1)
    ),
    "boosting steps and penalty rule fixed before resampling",
    "PLS component count fixed before resampling",
    "SuperPC threshold and component count fixed before resampling",
    survivalsvm_base,
    "XGBoost-Cox gbtree controls and seed fixed before resampling",
    paste0(
      "XGBoost-Cox DART tree and dropout controls and seed fixed before ",
      "resampling"
    ),
    paste0(
      "XGBoost-Cox coordinate-descent linear booster penalty, update order, ",
      "and seed fixed before resampling"
    ),
    paste0(
      "XGBoost-AFT normal-loss gbtree controls, scale, and seed fixed ",
      "before resampling"
    ),
    paste0(
      "XGBoost-AFT normal-loss DART tree/dropout controls, scale, and seed ",
      "fixed before resampling"
    ),
    paste0(
      "XGBoost-AFT normal-loss coordinate-descent linear booster penalty, ",
      "scale, update order, and seed fixed before resampling"
    ),
    "Weibull AFT distribution fixed; scale estimated in each analysis fold",
    "exponential AFT distribution and unit scale fixed before resampling",
    "log-normal AFT distribution fixed; scale estimated in each analysis fold",
    "log-logistic AFT distribution fixed; scale estimated in each analysis fold",
    vapply(
      sb_flexsurv_parametric_learner_ids(),
      sb_flexsurv_parametric_parameter_semantics,
      character(1)
    ),
    vapply(
      sb_flexsurvspline_learner_ids(),
      sb_flexsurvspline_parameter_semantics,
      character(1)
    ),
    sb_mtlr_survival_parameter_semantics(),
    sb_pamr_survival_parameter_semantics(),
    "survival-tree controls fixed before resampling",
    "random survival forest controls and seed fixed before resampling",
    "conditional-inference survival-tree controls fixed before resampling",
    "conditional-inference survival-forest controls fixed before resampling",
    paste0(
      "discrete-time probit BART trees, time grid, posterior draws, burn-in, ",
      "and seed fixed before resampling"
    ),
    paste0(
      "oblique survival-forest trees, split controls, mtry rule, sampling, ",
      "seed, and one-thread execution fixed before resampling"
    ),
    paste0(
      "Aalen additive-hazards time-varying coefficients and seeded native ",
      "tie handling fixed before resampling"
    ),
    paste0(
      "semiparametric proportional-odds iterations, profile fitting, ",
      "baseline variance, and seeded native tie handling fixed before ",
      "resampling"
    ),
    vapply(
      sb_ncvreg_survival_learner_ids(),
      sb_ncvreg_parameter_semantics,
      character(1)
    ),
    paste0(
      "ABESS sequence support path 0 through min(10,p), GIC criterion, ",
      "normalization, exact Newton updates, seed, and one-thread execution ",
      "fixed before resampling"
    ),
    sb_coxphf_parameter_semantics(sb_coxphf_learner_ids()),
    sb_bagged_knn_survival_parameter_semantics(
      sb_bagged_knn_survival_learner_ids()
    ),
    vapply(
      sb_ahaz_survival_learner_ids(),
      sb_ahaz_parameter_semantics,
      character(1)
    ),
    sb_bujar_linear_parameter_semantics(),
    "native Akritas conditional product-limit; analysis empirical joint CDF; fixed lambda 0.5; negative step-RMST through maximum analysis follow-up; no search; requires sensitive training predictor basis and paired outcomes"
  )
  survival_selected <- c(
    "fold-local Top-K then fixed ordinary Cox fit",
    paste0(
      "Top-K refit in each inner analysis fold; alpha and ",
      "standardization setting, lambda rule, and tuning measure fixed; ",
      "lambda selected by inner CV"
    ),
    paste0("fold-local Top-K; ", survival_base[-c(1L, 2L)])
  )
  survival <- paired_rows(
    "survival",
    survival_learners,
    "univ_cox_top_k",
    c(
      NA, "glmnet", "catboost", "catboost", "gbm", "mboost", "mboost",
      "mboost + partykit",
      vapply(
        sb_mboost_aft_learner_ids(),
        function(learner) {
          paste(sb_mboost_aft_dependency(learner)$packages, collapse = " + ")
        },
        character(1)
      ),
      vapply(
        sb_mboost_ranking_learner_ids(),
        function(learner) {
          paste(
            sb_mboost_ranking_dependency(learner)$packages,
            collapse = " + "
          )
        },
        character(1)
      ),
      rep("ranger", length(sb_ranger_survival_learner_ids())),
      rep("RLT", length(sb_rlt_survival_learner_ids())),
      "CoxBoost",
      "plsRcox", "superpc",
      rep("survivalsvm", length(sb_survivalsvm_learner_ids())),
      rep("xgboost", 6L),
      rep("survival", 4L),
      rep("flexsurv", length(sb_flexsurv_parametric_learner_ids())),
      rep("flexsurv", length(sb_flexsurvspline_learner_ids())),
      "MTLR",
      "pamr",
      "rpart", "randomForestSRC", "partykit + coin", "partykit + coin",
      "BART", "aorsf", rep("timereg", 2L),
      rep("ncvreg", length(sb_ncvreg_survival_learner_ids())),
      "abess", "coxphf", "bnnSurvival",
      rep("ahaz", length(sb_ahaz_survival_learner_ids())), "bujar", "survivalmodels + distr6"
    ),
    c(
      "none",
      "lambda_selected_by_inner_cv",
      rep("none", length(survival_learners) - 2L)
    ),
    survival_base,
    survival_selected
  )
  survival$tuning_semantics[
    survival$learner == "abess_cox_gic"
  ] <- "analysis_fold_gic"
  survival$tuning_semantics[
    survival$learner %in% sb_ncvreg_survival_learner_ids()
  ] <- "analysis_fold_bic"
  survival$tuning_semantics[
    survival$learner %in% c("ahaz_lasso_bic", "ahaz_sscad_bic")
  ] <- "analysis_fold_bic"
  survival <- rbind(
    survival,
    paired_rows(
      "survival",
      "bart_aft_survival",
      "univ_cox_top_k",
      "BART",
      "none",
      sb_rule_aft_parameter_semantics("bart_aft_survival"),
      paste0(
        "fold-local Top-K; ",
        sb_rule_aft_parameter_semantics("bart_aft_survival")
      )
    )
  )
  tuned_survival_map <- c(
    coxboost = "coxboost_cv",
    plsrcox = "plsrcox_cv",
    superpc = "superpc_cv",
    gbm_cox = "gbm_cox_cv"
  )
  tuned_survival <- do.call(rbind, lapply(
    names(tuned_survival_map),
    function(source_learner) {
      rows <- survival[survival$learner == source_learner, , drop = FALSE]
      rows$learner <- unname(tuned_survival_map[[source_learner]])
      rows$tuning_semantics <- "analysis_fold_internal_cv"
      rows$parameter_semantics <- ifelse(
        rows$selector == "none",
        "reference tuning repeated inside each analysis fold",
        paste0(
          "fold-local Top-K followed by reference tuning inside the same ",
          "analysis fold"
        )
      )
      rows
    }
  ))
  rownames(tuned_survival) <- NULL
  survival <- rbind(survival, tuned_survival)

  binary_learners <- c(
    "logistic", "glmnet_logistic", "gaussian_naive_bayes",
    sb_omics_specialist_binary_learner_ids(),
    sb_lightgbm_binary_learner_ids(),
    "catboost_binary", sb_ranger_binary_learner_ids(),
    sb_rlt_binary_learner_ids(),
    "adaboost_binary", "logitboost_binary",
    sb_classical_binary_learner_ids(),
    sb_specialized_binary_linear_learner_ids(),
    sb_xgboost_binary_learner_ids(), sb_ncvreg_binary_learner_ids(),
    "abess_binary_gic",
    sb_mclust_da_learner_ids(), sb_mixomics_da_learner_ids(),
    sb_hdda_binary_learner_ids(),
    sb_wsrf_binary_learner_ids(), sb_random_glm_binary_learner_ids()
  )
  binary_base <- c(
    "fixed ordinary logistic fit",
    paste0(
      "alpha, standardization setting, lambda rule, and tuning measure fixed; ",
      "lambda selected by inner CV"
    ),
    "variance smoothing fixed before resampling",
    vapply(
      sb_omics_specialist_binary_learner_ids(),
      sb_omics_specialist_binary_parameter_semantics,
      character(1)
    ),
    "LightGBM binary GBDT controls and seed fixed before resampling",
    paste0(
      "LightGBM binary DART tree and dropout controls and seed fixed before ",
      "resampling"
    ),
    paste0(
      "LightGBM binary GOSS top-gradient and random-small-gradient fractions ",
      "and seed fixed before resampling"
    ),
    paste0(
      "LightGBM binary random-forest row and feature fractions, frequency, ",
      "and seed fixed before resampling"
    ),
    "all learner parameters fixed before resampling",
    vapply(
      sb_ranger_binary_learner_ids(),
      sb_ranger_parameter_semantics,
      character(1)
    ),
    vapply(
      sb_rlt_binary_learner_ids(),
      sb_rlt_parameter_semantics,
      character(1)
    ),
    "AdaBoost.M1 maximum iteration count fixed before resampling",
    "LogitBoost iteration count fixed before resampling",
    "Euclidean nearest class mean, uniform priors, no scaling/shrinkage; uncalibrated squared-distance softmax",
    "Bernoulli NB with x > 0, Laplace alpha 1 and training-fold empirical priors",
    "empirical class priors fixed within each analysis fold",
    "empirical class priors fixed within each analysis fold",
    paste0(
      "kernel-density Naive Bayes priors, bandwidth adjustment, and Laplace ",
      "factor fixed before resampling"
    ),
    "RDA gamma and lambda covariance regularization fixed before resampling",
    "standardized Euclidean KNN with k=5",
    paste0(
      "single-hidden-layer neural-network size, decay, iteration limit, ",
      "scaling, and seed fixed before resampling"
    ),
    "classification tree controls fixed before resampling",
    "random-forest controls and seed fixed before resampling",
    "linear SVM cost and scaling fixed before resampling",
    "radial SVM cost, gamma rule, and scaling fixed before resampling",
    paste0(
      "radial SVM cost, gamma rule, scaling, and fold-local inverse-frequency ",
      "class weights fixed before resampling"
    ),
    "polynomial SVM cost, gamma rule, degree, and scaling fixed before resampling",
    "sigmoid SVM cost, gamma rule, offset, and scaling fixed before resampling",
    "GBM controls and seed fixed before resampling",
    paste0(
      "probit BART trees, posterior draws, burn-in, priors, and seed fixed ",
      "before resampling"
    ),
    "boosted C5.0 trials, pruning controls, and seed fixed before resampling",
    "kernel-weighted KNN k, distance, kernel, and scaling fixed before resampling",
    paste0(
      "Gaussian-process RBF kernel, automatic analysis-fold width rule, ",
      "scaling, tolerance, and seed fixed before resampling"
    ),
    paste0(
      "nearest shrunken-centroid threshold, variance scaling, offset, and ",
      "empirical analysis-fold priors fixed before resampling"
    ),
    paste0(
      "two-hidden-layer resilient-backprop architecture, stopping threshold, ",
      "scaling, and seed fixed before resampling"
    ),
    "conditional-inference classification-tree controls fixed before resampling",
    paste0(
      "conditional-inference classification-forest controls, mtry rule, ",
      "sampling, seed, and one-core execution fixed before resampling"
    ),
    "bootstrap Rpart bagging controls and seed fixed before resampling",
    paste0(
      "smooth additive Binomial boosting iterations, shrinkage, and spline ",
      "base learners fixed before resampling"
    ),
    paste0(
      "conditional-tree Binomial boosting iterations, shrinkage, and tree ",
      "controls fixed before resampling"
    ),
    paste0(
      "localized LDA neighbor cap and locally weighted empirical priors ",
      "fixed before resampling"
    ),
    paste0(
      "oblique classification-forest trees, split controls, mtry rule, ",
      "sampling, seed, and one-thread execution fixed before resampling"
    ),
    paste0(
      "Earth hinge-basis degree, basis cap, backward pruning rule, and ",
      "Binomial-logit family fixed before resampling"
    ),
    paste0(
      "flexible discriminant polynomial degree, basis convention, and one ",
      "discriminant dimension fixed before resampling"
    ),
    paste0(
      "mixture discriminant subclasses, iterations, seed, linear polyreg ",
      "map, and discriminant dimension fixed before resampling"
    ),
    vapply(
      sb_liblinear_binary_learner_ids(),
      sb_liblinear_binary_parameter_semantics,
      character(1)
    ),
    paste0(
      "evolutionary-tree population, iteration budget, depth controls, ",
      "complexity weight, and seed fixed before resampling"
    ),
    paste0(
      "shrinkage discriminant covariance, variance, and class-frequency ",
      "intensities estimated analytically inside each analysis fold"
    ),
    paste0(
      "rotation-forest tree count, variable-subset rule, tree controls, and ",
      "seed fixed before resampling"
    ),
    paste0(
      "regularized random-forest tree count, mtry rule, node size, ",
      "regularization coefficient, and seed fixed before resampling"
    ),
    paste0(
      "penalized discriminant ridge lambda 1, one discriminant dimension, ",
      "and empirical analysis-fold priors fixed before resampling"
    ),
    paste0(
      "stabilised LDA native PC_q score-dimension rule and empirical class ",
      "priors applied inside each analysis fold"
    ),
    paste0(
      "classic tree deviance split, minimum node sizes, and minimum deviance ",
      "reduction fixed before resampling"
    ),
    paste0(
      "C5.0 rule-set trials, pruning controls, and seed fixed before ",
      "resampling"
    ),
    paste0(
      "random-fern depth, ensemble size, seed, one-thread execution, and ",
      "fixed softmax mapping of native class scores"
    ),
    paste0(
      "robust Mqle logistic tuning constant, convergence controls, and ",
      "predictor-weighting rule fixed before resampling"
    ),
    paste0(
      "full radial least-squares SVM tau and scaling fixed; predictor-only ",
      "kernel width estimated by native sigest inside each analysis fold ",
      "with local seed 1; fixed softmax mapping of native decision scores"
    ),
    paste0(
      "LVQ1 equal-prior-target five-neighbor native initialization with the ",
      "native size rule and minimum two, analysis-fold scaling, 100 updates ",
      "per prototype, alpha 0.03, seed 1, and fixed nearest-prototype score mapping"
    ),
    paste0(
      "optimized LVQ1 equal-prior-target five-neighbor native initialization ",
      "with the native size rule and minimum two, analysis-fold scaling, 40 ",
      "updates per prototype, alpha 0.3, seed 1, and fixed nearest-prototype ",
      "score mapping"
    ),
    paste0(
      "LVQ2.1 equal-prior-target five-neighbor native initialization with the ",
      "native size rule and minimum two, analysis-fold scaling, 100 updates per ",
      "prototype, alpha 0.03, window 0.3, seed 1, and fixed nearest-prototype ",
      "score mapping"
    ),
    paste0(
      "LVQ3 equal-prior-target five-neighbor native initialization with the ",
      "native size rule and minimum two, analysis-fold scaling, 100 updates per ",
      "prototype, alpha 0.03, window 0.3, epsilon 0.1, seed 1, and fixed ",
      "nearest-prototype score mapping"
    ),
    "direction and AIC stopping rule fixed before resampling",
    paste0(
      "iteration count, shrinkage, family, and centering fixed before ",
      "resampling"
    ),
    "component cap, scaling, and logistic PLS model fixed before resampling",
    "XGBoost binary gbtree controls and seed fixed before resampling",
    paste0(
      "XGBoost binary DART tree and dropout controls and seed fixed before ",
      "resampling"
    ),
    paste0(
      "XGBoost binary coordinate-descent linear booster penalty, update ",
      "order, and seed fixed before resampling"
    ),
    vapply(
      sb_ncvreg_binary_learner_ids(),
      sb_ncvreg_parameter_semantics,
      character(1)
    ),
    paste0(
      "ABESS sequence support path 0 through min(10,p), GIC criterion, ",
      "normalization, exact Newton updates, seed, and one-thread execution ",
      "fixed before resampling"
    ),
    paste0(
      "one Gaussian component per class with a common covariance structure ",
      "selected by BIC inside each analysis fold; known class labels require ",
      "no random initialization"
    ),
    paste0(
      "one through three Gaussian components and covariance structure ",
      "selected separately by class using BIC inside each analysis fold; ",
      "hierarchical subset cap 2000, VVV/SVD initialization, and seed 1 fixed"
    ),
    paste0(
      "two-component cap, analysis-fold scaling, and fixed softmax mapping ",
      "of native PLS-DA dummy-response scores"
    ),
    paste0(
      "two-component cap, up to ten nonzero loadings per component, ",
      "analysis-fold scaling, and fixed softmax mapping of native sparse ",
      "PLS-DA dummy-response scores"
    ),
    sb_hdda_binary_parameter_semantics(sb_hdda_binary_learner_ids()),
    sb_wsrf_binary_parameter_semantics(sb_wsrf_binary_learner_ids()),
    sb_slope_randomglm_binary_parameter_semantics(
      sb_random_glm_binary_learner_ids()
    )
  )
  binary_selected <- c(
    "fold-local Top-K then fixed ordinary logistic fit",
    paste0(
      "Top-K refit in each inner analysis fold; alpha and ",
      "standardization setting, lambda rule, and tuning measure fixed; ",
      "lambda selected by inner CV"
    ),
    paste0("fold-local Top-K; ", binary_base[-c(1L, 2L)])
  )
  binary <- paired_rows(
    "binary",
    binary_learners,
    "univ_binary_auc_top_k",
    c(
      NA, "glmnet", NA,
      vapply(
        sb_omics_specialist_binary_learner_ids(),
        function(learner) {
          sb_omics_specialist_binary_dependency(learner)$package
        },
        character(1)
      ),
      rep("lightgbm", 4L), "catboost",
      rep("ranger", length(sb_ranger_binary_learner_ids())),
      rep("RLT", length(sb_rlt_binary_learner_ids())),
      "fastAdaboost", "caTools", NA, NA, "MASS", "MASS", "klaR", "klaR", NA,
      "nnet", "rpart", "randomForest", "e1071", "e1071", "e1071",
      "e1071", "e1071", "gbm", "BART", "C50", "kknn", "kernlab",
      "pamr", "neuralnet",
      "partykit", "partykit", "ipred", "mboost", "mboost", "klaR",
      "aorsf", "earth", "mda", "mda",
      rep("LiblineaR", length(sb_liblinear_binary_learner_ids())),
      "evtree", "sda",
      "rotationForest", "RRF", "mda", "ipred", "tree", "C50",
      "rFerns", "robustbase", "kernlab", rep("class", 4L),
      NA, "mboost", "plsRglm",
      rep("xgboost", 3L),
      rep("ncvreg", length(sb_ncvreg_binary_learner_ids())),
      "abess",
      "mclust", "mclust",
      "mixOmics", "mixOmics", "HDclassif", "wsrf", "randomGLM"
    ),
    c(
      "none",
      "lambda_selected_by_inner_cv",
      rep("none", length(binary_learners) - 2L)
    ),
    binary_base,
    binary_selected
  )
  binary$tuning_semantics[
    binary$learner == "abess_binary_gic"
  ] <- "analysis_fold_gic"
  binary$tuning_semantics[
    binary$learner %in% sb_mclust_da_learner_ids()
  ] <- "analysis_fold_bic"
  binary$tuning_semantics[
    binary$learner %in% sb_ncvreg_binary_learner_ids()
  ] <- "analysis_fold_bic"
  binary$tuning_semantics[
    binary$learner == "lssvm_radial_binary"
  ] <- "analysis_fold_predictor_only_sigest"
  binary$tuning_semantics[
    binary$learner %in% sb_svm_probability_learner_ids()
  ] <- "analysis_fold_oof_probability_calibration"
  registry <- rbind(survival, binary)
  stepcox <- registry[
    registry$task_type == "survival" & registry$selector == "none", ,
    drop = FALSE
  ]
  stepcox$selector <- "stepcox"
  stepcox$selector_parameter <- "stepcox_direction"
  stepcox$parameter_semantics <- ifelse(
    stepcox$tuning_semantics == "lambda_selected_by_inner_cv",
    paste0(
      "StepCox AIC selection refit in each inner analysis fold; ",
      "direction, alpha, standardization setting, lambda rule, and ",
      "tuning measure fixed; lambda selected by inner CV"
    ),
    paste0(
      "fold-local StepCox AIC selection with direction fixed before ",
      "resampling; ",
      stepcox$parameter_semantics
    )
  )
  registry <- rbind(registry, stepcox)

  combine_dependencies <- function(selector_dependency, learner_dependency) {
    if (is.na(selector_dependency) || !nzchar(selector_dependency)) {
      return(learner_dependency)
    }
    if (is.na(learner_dependency) || !nzchar(learner_dependency)) {
      return(selector_dependency)
    }
    dependencies <- unique(c(
      strsplit(selector_dependency, " + ", fixed = TRUE)[[1L]],
      strsplit(learner_dependency, " + ", fixed = TRUE)[[1L]]
    ))
    paste(dependencies, collapse = " + ")
  }
  abess_selector_rows <- function(task_type, selector, excluded_learners) {
    rows <- registry[
      registry$task_type == task_type & registry$selector == "none" &
        !registry$learner %in% excluded_learners, ,
      drop = FALSE
    ]
    rows$source <- "stormbreaker_abess_runtime"
    rows$selector <- selector
    rows$selector_parameter <- NA_character_
    rows$optional_dependency <- mapply(
      combine_dependencies,
      "abess",
      rows$optional_dependency,
      USE.NAMES = FALSE
    )
    rows$parameter_semantics <- paste0(
      "fold-local ABESS sequence path with GIC-selected nonzero support ",
      "from 2 through min(10,p); ",
      rows$parameter_semantics
    )
    rows
  }
  registry <- rbind(
    registry,
    abess_selector_rows(
      "binary",
      "abess_binary_gic_nonzero",
      "abess_binary_gic"
    ),
    abess_selector_rows(
      "survival",
      "abess_cox_gic_nonzero",
      c("abess_cox_gic", sb_analysis_tuned_survival_learner_ids())
    )
  )
  learned_selector_rows <- function(task_type, selector, selector_dependency,
                                    excluded_learners, selection_semantics,
                                    source) {
    rows <- registry[
      registry$task_type == task_type & registry$selector == "none" &
        !registry$learner %in% excluded_learners, ,
      drop = FALSE
    ]
    rows$source <- source
    rows$selector <- selector
    rows$selector_parameter <- NA_character_
    rows$optional_dependency <- mapply(
      combine_dependencies,
      selector_dependency,
      rows$optional_dependency,
      USE.NAMES = FALSE
    )
    rows$parameter_semantics <- paste0(
      selection_semantics,
      "; ",
      rows$parameter_semantics
    )
    rows
  }
  information_criteria <- c(
    jmi_binary_nlogn = "joint mutual-information (JMI)",
    cmim_binary_nlogn = "conditional mutual-information maximization (CMIM)",
    disr_binary_nlogn = "double-input symmetrical relevance (DISR)",
    jmim_binary_nlogn = "minimal joint mutual-information maximization (JMIM)"
  )
  model_selectors <- sb_binary_model_selector_ids()
  for (id in sb_pamr_selector_ids()) registry <- rbind(registry,
    learned_selector_rows(if (id == "pamr_survival_nonzero") "survival" else "binary",
      id, "pamr", sb_analysis_tuned_survival_learner_ids(),
      paste0("fold-local native PAMR nonzero shrunken-centroid features; fixed threshold 1; ",
        if (id == "pamr_survival_nonzero") "two native censoring-aware survival groups" else "empirical binary class priors",
        "; scale.sd TRUE, offset 50, no threshold search or feature padding"),
      "stormbreaker_pamr_selector_runtime"))
  for (i in seq_along(sb_survival_weight_selector_ids())) registry <- rbind(registry,
    learned_selector_rows("survival", sb_survival_weight_selector_ids()[i],
      c("glmnet", "plsRcox", "survivalsvm")[i], sb_analysis_tuned_survival_learner_ids(),
      paste0("fold-local ", c("Ridge Cox", "PLS-Cox", "linear survival-SVM")[i],
        " absolute coefficient times training SD; retain top half, minimum 2 positive weights when available"),
      "stormbreaker_survival_weight_selector_runtime"))
  for (i in seq_along(sb_survival_boosting_selector_ids())) registry <- rbind(registry,
    learned_selector_rows("survival", sb_survival_boosting_selector_ids()[i],
      c("gbm", "mboost", "superpc", "mboost", "xgboost")[i], sb_analysis_tuned_survival_learner_ids(),
      c("fold-local Cox GBM positive relative influence; 100 trees, depth 1, shrinkage 0.1",
        "fold-local linear Cox glmBoost nonzero coefficients; 100 steps, shrinkage 0.1",
        "fold-local SuperPC native survival scores; absolute threshold 1, s0 percentile 0.5",
        "fold-local Cox GamBoost active bbs spline terms; dfbase 4, 100 steps, shrinkage 0.1; score is spline coefficient L2 norm",
        "fold-local native XGBoost Cox positive gain; 100 rounds, eta 0.1, depth 3, full sampling, one thread, seed 1")[i],
      "stormbreaker_survival_boosting_selector_runtime"))
  for (i in seq_along(model_selectors)) registry <- rbind(registry,
    learned_selector_rows("binary", model_selectors[i],
      c("glmnet", NA, "mboost", "randomForest", "glmnet", "gbm", "e1071", "mboost", "xgboost", "plsRglm")[i], character(),
      c("fold-local binomial LASSO CV nonzero coefficients",
        "fold-local stepwise logistic AIC variables",
        "fold-local glmBoost nonzero coefficients",
        "fold-local RF importance strictly above declared threshold",
        "fold-local binomial elastic-net CV nonzero coefficients",
        "fold-local GBM positive relative influence",
        "linear SVM recursive elimination with grouped inner-CV feature count selection",
        "fold-local binomial GamBoost active bbs spline terms; dfbase 4, 100 steps, shrinkage 0.1; score is spline coefficient L2 norm",
        "fold-local native XGBoost binary positive gain; 100 rounds, eta 0.1, depth 3, full sampling, one thread, seed 1",
        "fold-local native logistic PLS; at most two scaled components; absolute coefficient times training SD; retain top half, minimum two positive weights when available; no search")[i],
      "stormbreaker_binary_model_selector_runtime"))
  information_selector_rows <- do.call(rbind, lapply(
    names(information_criteria),
    function(selector) {
      learned_selector_rows(
        "binary",
        selector,
        "praznik",
        character(),
        paste0(
          "fold-local ", information_criteria[[selector]],
          " selection using the native praznik ranking, ",
          "min(p, max(2, floor(n/log(n)))) retained predictors, all ",
          "finite native scores, and one thread"
        ),
        "stormbreaker_high_dimensional_runtime"
      )
    }
  ))
  registry <- rbind(
    registry,
    learned_selector_rows(
      "binary",
      "boruta_binary_roughfix",
      "Boruta",
      character(),
      paste0(
        "fold-local Boruta all-relevant selection using ranger normalized ",
        "permutation importance, multiplicity-adjusted p=0.01 decisions, ",
        "at most 100 runs, 500 trees, one thread, seed 1, and ",
        "TentativeRoughFix"
      ),
      "stormbreaker_forest_selector_runtime"
    ),
    learned_selector_rows(
      "binary",
      "varselrf_binary_oob_1se",
      "varSelRF",
      character(),
      paste0(
        "fold-local varSelRF backward elimination using fixed initial ",
        "MeanDecreaseAccuracy, the one-standard-error OOB subset rule, ",
        "5,000 initial trees, 2,000 iteration trees, 20% variable drops, ",
        "whole-range search, and seed 1"
      ),
      "stormbreaker_forest_selector_runtime"
    ),
    learned_selector_rows(
      "survival",
      "boruta_survival_roughfix",
      "Boruta",
      sb_analysis_tuned_survival_learner_ids(),
      paste0(
        "fold-local survival Boruta all-relevant selection using ranger ",
        "normalized permutation importance, multiplicity-adjusted p=0.01 ",
        "decisions, at most 100 runs, 500 trees, one thread, seed 1, and ",
        "TentativeRoughFix"
      ),
      "stormbreaker_forest_selector_runtime"
    ),
    learned_selector_rows(
      "binary",
      "sis_binary_nlogn",
      "SIS",
      character(),
      paste0(
        "fold-local binomial SIS marginal-utility ranking with ",
        "min(p, max(2, floor(n/log(n)))) retained predictors, native ",
        "standardization, non-iterative screening, and seed 1"
      ),
      "stormbreaker_high_dimensional_runtime"
    ),
    learned_selector_rows(
      "binary",
      "vsurf_binary_prediction",
      "VSURF",
      character(),
      paste0(
        "fold-local VSURF threshold, interpretation, and prediction stages ",
        "using its randomForest backend, native default forest counts, ",
        "serial execution, and seed 1"
      ),
      "stormbreaker_high_dimensional_runtime"
    ),
    learned_selector_rows(
      "binary",
      "mrmr_binary_nlogn",
      "praznik",
      character(),
      paste0(
        "fold-local minimum-redundancy maximum-relevance selection using ",
        "the native praznik ranking, min(p, max(2, floor(n/log(n)))) ",
        "retained predictors, all finite native scores, and one thread"
      ),
      "stormbreaker_high_dimensional_runtime"
    ),
    information_selector_rows,
    learned_selector_rows(
      "binary",
      "lpc_binary_nlogn",
      "lpc",
      character(),
      paste0(
        "fold-local Lassoed Principal Components using native two-class ",
        "t scores, the native 10 repeated half-split predictive-advantage ",
        "soft threshold, decreasing absolute LPC-score ranking, ",
        "min(p, max(2, floor(n/log(n)))) retained predictors, and seed 1"
      ),
      "stormbreaker_high_dimensional_runtime"
    ),
    learned_selector_rows(
      "survival",
      "lpc_cox_nlogn",
      "lpc",
      sb_analysis_tuned_survival_learner_ids(),
      paste0(
        "fold-local Lassoed Principal Components using native survival Cox ",
        "scores, the native 10 repeated half-split predictive-advantage ",
        "soft threshold, decreasing absolute LPC-score ranking, ",
        "min(p, max(2, floor(n/log(n)))) retained predictors, and seed 1"
      ),
      "stormbreaker_high_dimensional_runtime"
    ),
    learned_selector_rows(
      "binary",
      "stabsel_binary_lasso",
      "stabs + glmnet",
      "glmnet_logistic",
      paste0(
        "fold-local complementary-pairs stability selection with ",
        "conservative standardized binomial Lasso paths, 50 pairs, ",
        "min(max(2, floor(n/log(n))), floor(p/2)) variables per ",
        "subsample with at least four candidate predictors, ",
        "selection-probability cutoff 0.75, unimodality assumption, ",
        "serial execution, and seed 1"
      ),
      "stormbreaker_high_dimensional_runtime"
    ),
    learned_selector_rows(
      "survival",
      "stabsel_cox_lasso",
      "stabs + glmnet",
      c(sb_analysis_tuned_survival_learner_ids(), "glmnet_cox"),
      paste0(
        "fold-local complementary-pairs stability selection with ",
        "conservative standardized Cox Lasso paths with Breslow ties, ",
        "50 pairs, ",
        "min(max(2, floor(n/log(n))), floor(p/2)) variables per ",
        "subsample with at least four candidate predictors, ",
        "selection-probability cutoff 0.75, unimodality assumption, ",
        "serial execution, and seed 1"
      ),
      "stormbreaker_high_dimensional_runtime"
    )
  )
  mime_selector_rows <- function(selector, learners, selector_parameter,
                                 selector_dependency, selection_semantics) {
    rows <- registry[
      registry$task_type == "survival" &
        registry$selector == "none" &
        registry$learner %in% learners, ,
      drop = FALSE
    ]
    rows <- rows[match(learners, rows$learner), , drop = FALSE]
    rows$source <- "stormbreaker_mime117_runtime"
    rows$selector <- selector
    rows$selector_parameter <- selector_parameter
    rows$optional_dependency <- mapply(
      combine_dependencies,
      selector_dependency,
      rows$optional_dependency,
      USE.NAMES = FALSE
    )
    rows$parameter_semantics <- paste0(
      selection_semantics,
      "; ",
      rows$parameter_semantics
    )
    rows
  }

  registry <- rbind(
    registry,
    mime_selector_rows(
      "rsf_minimal_depth",
      c(
        "coxboost_cv", "glmnet_cox", "gbm_cox_cv",
        "plsrcox_cv", "superpc_cv", sb_survivalsvm_learner_ids()
      ),
      "selector_rsf_nodesize_and_seed",
      "randomForestSRC",
      paste0(
        "fold-local 1,000-tree RSF minimal-depth selection with ",
        "conservative='high'"
      )
    ),
    mime_selector_rows(
      "coxboost_nonzero",
      c(
        "glmnet_cox", "gbm_cox_cv", "plsrcox_cv",
        "superpc_cv", sb_survivalsvm_learner_ids()
      ),
      "selector_inner_folds_and_seed",
      "CoxBoost",
      paste0(
        "fold-local CoxBoost penalty and step tuning followed by nonzero ",
        "coefficient selection"
      )
    ),
    mime_selector_rows(
      "enet_cox_nonzero",
      unique(c(registry$learner[registry$selector == "lasso_nonzero"],
        "coxph", "glmnet_cox", "ranger_survival",
        "coxboost_cv", "gbm_cox_cv", "plsrcox_cv",
        "rsf_survival", "superpc_cv", sb_survivalsvm_learner_ids())),
      "selector_inner_folds_and_seed", "glmnet",
      "fold-local Cox Enet alpha=0.5 lambda.min nonzero coefficients"
    ),
    mime_selector_rows(
      "lasso_nonzero",
      c(
        "coxboost_cv", "gbm_cox_cv", "plsrcox_cv",
        "rsf_survival", "superpc_cv", sb_survivalsvm_learner_ids()
      ),
      "selector_inner_folds_and_seed",
      "glmnet",
      paste0(
        "fold-local Cox Lasso lambda.min tuning followed by nonzero ",
        "coefficient selection"
      )
    ),
    mime_selector_rows(
      "rsf_then_stepcox",
      "coxph",
      "selector_rsf_nodesize_seed_and_stepcox_direction",
      "randomForestSRC",
      paste0(
        "fold-local 1,000-tree RSF minimal-depth selection followed by ",
        "fold-local AIC StepCox"
      )
    ),
    mime_selector_rows(
      "coxboost_then_stepcox",
      "coxph",
      "selector_inner_folds_seed_and_stepcox_direction",
      "CoxBoost",
      paste0(
        "fold-local tuned CoxBoost nonzero selection followed by fold-local ",
        "AIC StepCox"
      )
    ),
    mime_selector_rows(
      "lasso_then_stepcox",
      "coxph",
      "selector_inner_folds_seed_and_stepcox_direction",
      "glmnet",
      paste0(
        "fold-local tuned Cox Lasso nonzero selection followed by fold-local ",
        "AIC StepCox"
      )
    )
  )
  for (selector in c("lasso_nonzero", "enet_cox_nonzero", "coxboost_nonzero", "rsf_minimal_depth")) {
    dependency <- switch(selector, lasso_nonzero = "glmnet", enet_cox_nonzero = "glmnet", coxboost_nonzero = "CoxBoost",
      rsf_minimal_depth = "randomForestSRC")
    extra <- learned_selector_rows("survival", selector, dependency,
      sb_analysis_tuned_survival_learner_ids(),
      "fold-local native feature selection; existing selector defaults retained",
      "stormbreaker_mime117_runtime")
    extra$selector_parameter <- if (selector == "rsf_minimal_depth")
      "selector_rsf_nodesize_and_seed" else "selector_inner_folds_and_seed"
    extra <- extra[!extra$learner %in% registry$learner[registry$selector == selector], , drop = FALSE]
    registry <- rbind(registry, extra)
  }
  restricted_internal_cv <- data.frame(
    source = rep("stormbreaker_builtin", 4L),
    task_type = c("binary", "survival", "survival", "binary"),
    selector = rep("none", 4L),
    learner = c(
      "pre_binary_rule_ensemble",
      "pre_survival_rule_ensemble",
      "penaft_gehan_elastic_net",
      "cancerclass_binary"
    ),
    admission_status = rep("admitted", 4L),
    selector_parameter = rep(NA_character_, 4L),
    tuning_semantics = c(
      "analysis_fold_native_coefficient_cv_after_rule_induction",
      "analysis_fold_native_coefficient_cv_after_rule_induction",
      "analysis_fold_internal_cv_fixed",
      "analysis_fold_oof_probability_calibration"
    ),
    optional_dependency = c(
      "pre", "pre + mboost", "penAFT",
      "cancerclass + Biobase"
    ),
    parameter_semantics = c(
      vapply(
        c(
          "pre_binary_rule_ensemble",
          "pre_survival_rule_ensemble",
          "penaft_gehan_elastic_net"
        ),
        sb_rule_aft_parameter_semantics,
        character(1)
      ),
      paste0(
        "native Welch ranking and correlation-distance z score over all ",
        "supplied features; five-fold OOF logistic calibration inside each ",
        "analysis fold"
      )
    ),
    stringsAsFactors = FALSE
  )
  registry <- rbind(registry, restricted_internal_cv)
  first_learner_row <- !duplicated(registry$learner)
  fixed <- registry[first_learner_row, , drop = FALSE]
  fixed$selector <- "fixed_features"
  fixed$selector_parameter <- "selector_features_and_feature_set_id"
  fixed$parameter_semantics <- paste0(
    "predeclared outcome-independent feature set; ",
    fixed$parameter_semantics
  )
  rownames(fixed) <- NULL
  rbind(registry, fixed)
}

sb_pipeline_registry_match <- function(selector, learner,
                                       registry = sb_pipeline_registry_data()) {
  registry[
    registry$selector == selector & registry$learner == learner, ,
    drop = FALSE
  ]
}

sb_pipeline_uses_inner_cv <- function(learner,
                                      registry = sb_pipeline_registry_data()) {
  semantics <- unique(registry$tuning_semantics[
    registry$learner == learner
  ])
  identical(semantics, "lambda_selected_by_inner_cv")
}

sb_pipeline_registry_learner_view <- function(registry) {
  learners <- registry[registry$selector == "none", , drop = FALSE]
  route_counts <- table(registry$learner)
  selectors <- vapply(
    learners$learner,
    function(learner) {
      paste(registry$selector[registry$learner == learner], collapse = ", ")
    },
    character(1)
  )
  data.frame(
    task_type = learners$task_type,
    learner = learners$learner,
    pipeline_count = unname(as.integer(route_counts[learners$learner])),
    selectors = selectors,
    tuning_semantics = learners$tuning_semantics,
    optional_dependency = learners$optional_dependency,
    parameter_semantics = learners$parameter_semantics,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' List pipeline combinations admitted by the built-in runtime
#'
#' Returns the selector and learner combinations that the installed package
#' knows how to construct and resolve. `admitted` means that STORMBREAKER has a
#' matching built-in specification and adapter route. It does not mean that an
#' optional package is installed, that a model will fit a particular dataset,
#' or that the statistical design is appropriate. `source` identifies package
#' code, not literature provenance. At learner level, dependency columns
#' describe the learner itself (the `selector = "none"` route); use route level
#' when dependencies of a learner plus a feature selector must both be checked.
#'
#' @param task_type Optional task filter: `"survival"` or `"binary"`.
#' @param level Return one row per executable `"route"` or one row per
#'   `"learner"`.
#' @param check_installed Whether to add a side-effect-free check of installed
#'   R packages. This uses package locations only: it does not load optional R
#'   namespaces, fit a model, or claim that a method suits the supplied data.
#'
#' @return A data frame. At route level, each row is an admitted selector and
#'   learner combination. At learner level, each row summarizes one learner,
#'   its available selectors, pipeline count, tuning semantics, and optional
#'   dependency. With `check_installed = TRUE`, four additional columns list R
#'   packages, a reserved module field, missing R packages, and one dependency
#'   status. The current runtime is R-only, so `python_modules` is empty for
#'   every registered route. Admission does not imply that optional software is
#'   installed or operational.
#' @examples
#' registry <- sb_pipeline_registry()
#' registry[c("task_type", "selector", "learner")]
#' sb_pipeline_registry(task_type = "binary", level = "learner")
#' sb_pipeline_registry(level = "learner", check_installed = TRUE)
#' @export
sb_pipeline_registry <- function(task_type = NULL,
                                 level = c("route", "learner"),
                                 check_installed = FALSE) {
  level <- match.arg(level)
  if (!is.null(task_type)) {
    task_type <- match.arg(task_type, c("survival", "binary"))
  }
  if (!is.logical(check_installed) || length(check_installed) != 1L ||
    is.na(check_installed)) {
    sb_abort(
      "`check_installed` must be `TRUE` or `FALSE`.",
      "stormbreaker_bad_registry_check_installed"
    )
  }
  registry <- sb_pipeline_registry_data()
  if (!is.null(task_type)) {
    registry <- registry[registry$task_type == task_type, , drop = FALSE]
    rownames(registry) <- NULL
  }
  if (identical(level, "learner")) {
    registry <- sb_pipeline_registry_learner_view(registry)
  }
  if (check_installed) {
    registry <- sb_add_dependency_status(registry)
  }
  registry
}

sb_pipeline_lint_report <- function(status, pipeline = NULL, spec = NULL,
                                    task_type = NULL, error = NULL) {
  valid_named_spec <- !inherits(spec, "stormbreaker_pipeline") &&
    is.list(spec) && (length(spec) == 0L ||
    (!is.null(names(spec)) && !any(is.na(names(spec))) &&
      !any(!nzchar(names(spec))) && anyDuplicated(names(spec)) == 0L))
  raw_value <- function(name, default = NA_character_) {
    if (!valid_named_spec) {
      return(NA_character_)
    }
    if (!name %in% names(spec)) {
      return(default)
    }
    value <- spec[[name]]
    if (!is.character(value) || length(value) != 1L || is.na(value)) {
      return(NA_character_)
    }
    value
  }
  pipeline_value <- function(name) {
    if (!is.list(pipeline) || !name %in% names(pipeline)) {
      return(NA_character_)
    }
    value <- pipeline[[name]]
    if (!is.character(value) || length(value) != 1L || is.na(value)) {
      return(NA_character_)
    }
    value
  }
  has_pipeline <- is.list(pipeline)
  pipeline_id <- pipeline_value("id")
  reported_task <- if (!is.null(task_type)) {
    task_type
  } else {
    pipeline_value("task_type")
  }
  selector <- if (has_pipeline) {
    pipeline_value("selector")
  } else {
    raw_value("selector", "none")
  }
  learner <- if (has_pipeline) {
    pipeline_value("learner")
  } else {
    raw_value("learner", "coxph")
  }
  issue_code <- if (is.null(error)) NA_character_ else class(error)[[1L]]
  message <- if (is.null(error)) NA_character_ else conditionMessage(error)

  data.frame(
    status = status,
    pipeline_id = pipeline_id,
    task_type = reported_task,
    selector = selector,
    learner = learner,
    issue_code = issue_code,
    message = message,
    stringsAsFactors = FALSE
  )
}

sb_is_expected_pipeline_lint_error <- function(error) {
  sb_is_expected_pipeline_constructor_error(error) ||
    inherits(error, c(
      "stormbreaker_bad_pipeline_spec",
      "stormbreaker_unknown_pipeline_argument",
      "stormbreaker_bad_pipeline_task_type",
      "stormbreaker_pipeline_task_mismatch",
      "stormbreaker_unresolved_pipeline"
    ))
}

sb_pipeline_from_lint_spec <- function(spec) {
  if (!is.list(spec)) {
    sb_abort(
      paste0(
        "`spec` must be one `stormbreaker_pipeline` or one named list ",
        "of arguments for `sb_pipeline()`."
      ),
      "stormbreaker_bad_pipeline_spec"
    )
  }
  if (length(spec) > 0L && (is.null(names(spec)) ||
    any(is.na(names(spec))) || any(!nzchar(names(spec))) ||
    anyDuplicated(names(spec)) > 0L)) {
    sb_abort(
      "A pipeline argument list must have unique, non-empty names.",
      "stormbreaker_bad_pipeline_spec"
    )
  }
  unknown <- setdiff(names(spec), names(formals(sb_pipeline)))
  if (length(unknown) > 0L) {
    sb_abort(
      paste0(
        "Unknown `sb_pipeline()` argument", if (length(unknown) > 1L) {
          "s"
        } else {
          ""
        }, ": ", paste0("`", unknown, "`", collapse = ", "), "."
      ),
      "stormbreaker_unknown_pipeline_argument"
    )
  }
  do.call(sb_pipeline, spec)
}

#' Check whether one pipeline specification is admitted
#'
#' Checks one existing pipeline object or one named list of arguments for
#' [sb_pipeline()]. The check reuses the constructor, the pipeline contract,
#' and the built-in adapter resolver; it does not fit a model or inspect data.
#' Expected STORMBREAKER specification errors are returned as `rejected`.
#' Unexpected internal errors are not converted and remain visible.
#'
#' `admitted` means that the specification has a registered adapter route. It
#' does not check optional-package availability, study design quality, data
#' suitability, literature provenance, or whether fitting will succeed.
#'
#' @param spec One `stormbreaker_pipeline` object or one named list of
#'   arguments accepted by [sb_pipeline()].
#' @param task_type Optional target task, either `"survival"` or `"binary"`.
#'   When supplied, the pipeline must match it.
#'
#' @return A one-row data frame containing status, pipeline identity, task,
#'   selector, learner, and an issue code and message for a rejected
#'   specification.
#' @examples
#' sb_lint_pipeline(list(learner = "logistic"))
#' sb_lint_pipeline(list(learner = "glmnet_cox"))
#' @export
sb_lint_pipeline <- function(spec, task_type = NULL) {
  if (!is.null(task_type) && (!is.character(task_type) ||
    length(task_type) != 1L || is.na(task_type) ||
    !task_type %in% c("survival", "binary"))) {
    sb_abort(
      "`task_type` must be `NULL`, `\"survival\"`, or `\"binary\"`.",
      "stormbreaker_bad_pipeline_lint_task"
    )
  }

  pipeline <- NULL
  result <- tryCatch(
    {
      pipeline <- if (inherits(spec, "stormbreaker_pipeline")) {
        spec
      } else {
        sb_pipeline_from_lint_spec(spec)
      }
      target_task <- if (is.null(task_type)) {
        pipeline[["task_type"]]
      } else {
        task_type
      }
      if (!is.character(target_task) || length(target_task) != 1L ||
        is.na(target_task) ||
        !target_task %in% c("survival", "binary")) {
        sb_abort(
          "The pipeline does not declare a supported task type.",
          "stormbreaker_bad_pipeline_task_type"
        )
      }
      pipeline <- sb_validate_pipeline_contract(pipeline, target_task)
      pipeline
    },
    stormbreaker_error = function(error) {
      if (!sb_is_expected_pipeline_lint_error(error)) {
        stop(error)
      }
      error
    }
  )

  if (inherits(result, "stormbreaker_error")) {
    return(sb_pipeline_lint_report(
      status = "rejected",
      pipeline = pipeline,
      spec = spec,
      task_type = task_type,
      error = result
    ))
  }
  runtime <- sb_resolve_pipeline_runtime(result)
  sb_validate_pipeline_runtime(runtime, result)
  sb_pipeline_lint_report(
    status = "admitted",
    pipeline = result,
    spec = spec,
    task_type = task_type
  )
}
