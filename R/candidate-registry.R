sb_candidate_catalog_text <- function(value) {
  if (length(value) == 0L || is.na(value)) {
    return("")
  }
  as.character(value)
}

sb_candidate_stepcox_direction <- function(label) {
  matched <- regexec(
    "^StepCox\\[(forward|backward|both)\\]$",
    label
  )
  parts <- regmatches(label, matched)[[1L]]
  if (length(parts) != 2L) {
    return(NULL)
  }
  parts[[2L]]
}

sb_candidate_stepglm_direction <- function(label) {
  matched <- regexec(
    "^Stepglm\\[(forward|backward|both)\\]$",
    label
  )
  parts <- regmatches(label, matched)[[1L]]
  if (length(parts) != 2L) {
    return(NULL)
  }
  parts[[2L]]
}

sb_candidate_glmnet_alpha <- function(label) {
  if (identical(label, "Ridge")) {
    return(0)
  }
  if (identical(label, "Lasso")) {
    return(1)
  }
  if (!grepl("^Enet\\[alpha=[0-9.]+\\]$", label)) {
    return(NULL)
  }
  alpha <- as.numeric(sub(
    "^Enet\\[alpha=([0-9.]+)\\]$",
    "\\1",
    label
  ))
  if (is.na(alpha) || !is.finite(alpha) || alpha < 0 || alpha > 1) {
    return(NULL)
  }
  alpha
}

sb_candidate_runtime_parameter <- function(name, value) {
  paste0(
    name,
    "=",
    format(value, trim = TRUE, scientific = FALSE)
  )
}

sb_candidate_source_spec_audit <- function(task_type, model_type,
                                           learner_label, selector) {
  if (identical(selector, "stepcox")) {
    return(list(
      status = "audited_source_identity_ambiguous",
      reference = paste0(
        "Mime1 0.0.0.9000 ML.Dev.Prog.Sig.R and the historical Fig.5A ",
        "reconstruction ML296_Fig5_Backend.R; audited 2026-08-17 and ",
        "2026-08-24."
      ),
      difference = paste0(
        "The shared catalog label does not identify one operational StepCox ",
        "specification. Mime1 applies a common full-training univariable-Cox ",
        "prefilter, starts every direction from a full Cox model, and omits a ",
        "two-stage result when fewer than two variables remain. The historical ",
        "Fig.5A reconstruction uses different forward-start, minimum-feature, ",
        "truncation, warning, and fallback behavior. STORMBREAKER has no ",
        "implicit common prefilter, refits inside analysis resampling, starts ",
        "forward from the null model and backward/both from the full model, ",
        "and reports failures instead of silently substituting a result."
      )
    ))
  }

  direct_survival_glmnet <- identical(task_type, "survival") &&
    identical(model_type, "single") && identical(selector, "none") &&
    !is.null(sb_candidate_glmnet_alpha(learner_label))
  if (direct_survival_glmnet) {
    difference <- if (identical(learner_label, "Ridge")) {
      paste0(
        "Official Mime1 fits alpha=0 but selects lambda.min from a separate ",
        "cv.glmnet call that omits alpha and therefore uses glmnet's default ",
        "alpha=1; the historical Fig.5A reconstruction cross-validates ",
        "alpha=0 with a sample-size-dependent fold count."
      )
    } else if (identical(learner_label, "Lasso")) {
      paste0(
        "Official Mime1 cross-validates alpha=1 with fixed 10 folds and ",
        "lambda.min; the historical Fig.5A reconstruction uses alpha=1 with ",
        "a sample-size-dependent fold count."
      )
    } else {
      paste0(
        "Official Mime1 cross-validates the label alpha with fixed 10 folds ",
        "and lambda.min; the historical Fig.5A reconstruction uses the same ",
        "label alpha with a sample-size-dependent fold count."
      )
    }
    difference <- paste0(
      difference,
      " Mime1's wrapper also applies a common full-training univariable-Cox ",
      "prefilter, intersects predictors across supplied cohorts, and ",
      "mean-imputes validation cohorts. The catalog supplies only ",
      "historical_fig5a_label_inferred and no DOI or figure, so it does not ",
      "establish which implementation generated the row. STORMBREAKER ",
      "cross-validates the declared alpha with declared inner folds, seed, ",
      "measure, lambda rule, and standardization inside each outer analysis ",
      "fold, without an implicit common prefilter or fallback."
    )
    return(list(
      status = "audited_source_identity_ambiguous",
      reference = paste0(
        "Mime1 0.0.0.9000 ML.Dev.Prog.Sig.R lines 764-786, 2013-2057, ",
        "2424-2447, and 2638-2697; historical Fig.5A reconstruction ",
        "ML296_Fig5_Backend.R lines 323-350. Catalog evidence is ",
        "historical_fig5a_label_inferred with no source DOI or figure, so ",
        "source identity is not uniquely established; audited 2026-08-24."
      ),
      difference = difference
    ))
  }

  direct_survival_fixed <- identical(task_type, "survival") &&
    identical(model_type, "single") && identical(selector, "none") &&
    learner_label %in% c(
      "survival-SVM", "CoxBoost", "GBM", "SuperPC", "plsRcox", "Ranger",
      "GamBoost", "GlmBoost", "BlackBoost"
    )
  if (direct_survival_fixed) {
    ambiguous <- learner_label %in% c(
      "survival-SVM", "CoxBoost", "GBM", "SuperPC"
    )
    difference <- switch(learner_label,
      CoxBoost = paste0(
        "Official Mime1 optimizes a penalty from start.penalty=500 and then ",
        "selects the boosting step by fixed 10-fold cross-validation with ",
        "maxstepno=500. The historical Fig.5A reconstruction instead uses ",
        "maxstepno=200 and a sample-size-dependent fold count. STORMBREAKER ",
        "uses a declared fixed step count and either a declared penalty or ",
        "its documented analysis-event-count penalty rule, without an ",
        "automatic penalty or step search."
      ),
      GBM = paste0(
        "Official Mime1 searches 10,000 trees by fixed 10-fold ",
        "cross-validation and refits at the selected tree count. The ",
        "historical Fig.5A reconstruction searches 2,000 trees with a ",
        "sample-size-dependent fold count and different shrinkage, without ",
        "Mime1's second fit. STORMBREAKER uses declared fixed tree, depth, ",
        "node-size, shrinkage, bag-fraction, and seed parameters with no ",
        "early stopping or automatic tree-count search."
      ),
      SuperPC = paste0(
        "Official Mime1 selects a threshold by fixed 10-fold superpc.cv, ",
        "retries indefinitely after an error, and supplies validation ",
        "outcomes in each superpc.predict test object. The historical Fig.5A ",
        "reconstruction uses a sample-size-dependent fold count and no retry. ",
        "STORMBREAKER uses a declared fixed threshold and component count, ",
        "does not receive validation outcomes for assessment prediction, and ",
        "does not retry or substitute a fallback result."
      ),
      `survival-SVM` = paste0(
        "Official Mime1 calls survivalsvm with gamma.mu=1 and the package ",
        "defaults, then reports the raw predicted score. The historical ",
        "Fig.5A reconstruction makes the regression, linear-kernel, and ",
        "quadprog choices explicit but retains raw package scores and its ",
        "optional outer fallback. STORMBREAKER uses an explicit linear ",
        "regression SVM, fold-local standardization, and analysis-fold Cox ",
        "calibration to a log-hazard score, with failures reported instead ",
        "of fallback."
      ),
      plsRcox = paste0(
        "The 101/117 implementations select the requested component count ",
        "with cv.plsRcox before fitting plsRcox; the historical Fig.5A ",
        "reconstruction exposes that learner only behind a selector and also ",
        "adjusts the candidate component limit to the fitted feature count. ",
        "STORMBREAKER uses a declared fixed component count and performs no ",
        "automatic component search. The compatibility label alone does not ",
        "identify which operational definition the direct row should reproduce."
      ),
      Ranger = paste0(
        "The cited Table 1 identifies the Ranger label and ranger package but ",
        "does not specify tree count, mtry, node size, sampling, split rule, ",
        "seed, or risk-score semantics. Mime1's RSF uses randomForestSRC and ",
        "is not ranger. STORMBREAKER uses its own declared fixed ranger ",
        "survival parameters and native survival-curve prediction; neither ",
        "source establishes that this is the catalog row's implementation."
      ),
      GamBoost = paste0(
        "Official Mime1 does not implement GamBoost. The historical Fig.5A ",
        "backend uses mboost::gamboost with CoxPH, mstop=100, and nu=0.1, ",
        "whereas a later local all-in-one reconstruction incorrectly dispatches ",
        "its GamBoost label to glmboost. STORMBREAKER uses the distinct native ",
        "gamboost spline learner with those fixed boosting controls."
      ),
      GlmBoost = paste0(
        "Official Mime1 does not implement survival GlmBoost. Local historical ",
        "reconstructions use mboost::glmboost with CoxPH, mstop=100, and nu=0.1. ",
        "STORMBREAKER uses that distinct native linear boosting family with ",
        "explicit predictors and no automatic iteration search."
      ),
      BlackBoost = paste0(
        "Official Mime1 and the historical Fig.5A backend do not implement ",
        "BlackBoost. A later local all-in-one reconstruction uses ",
        "mboost::blackboost with CoxPH, mstop=100, and nu=0.1. STORMBREAKER ",
        "uses the native tree-boosting learner with those fixed controls and ",
        "explicit default tree settings."
      )
    )
    if (ambiguous) {
      difference <- paste0(
        difference,
        " Mime1 also applies a common full-training univariable-Cox prefilter, ",
        "intersects predictors across supplied cohorts, and mean-imputes each ",
        "validation cohort. The reconstructed wrapper and STORMBREAKER do not ",
        "share that complete operational procedure."
      )
    }
    reference <- if (ambiguous) {
      paste0(
        "Mime1 0.0.0.9000 ML.Dev.Prog.Sig.R all-mode lines 1521-2008 and ",
        "single-mode lines 2474-2634; historical Fig.5A reconstruction ",
        "ML296_Fig5_Backend.R lines 279-305 and 366-409. The bundled row ",
        "records historical_fig5a_label_inferred but no source DOI or figure, ",
        "so it does not identify one operational implementation; audited ",
        "2026-08-24."
      )
    } else if (identical(learner_label, "plsRcox")) {
      paste0(
        "Catalog-construction audit: direct plsRcox is the sole 101/117 ",
        "compatibility addition absent from the reconstructed 296 single-model ",
        "set. Comparison references are the Mime paper DOI ",
        "10.1016/j.csbj.2024.06.035, Mime1 0.0.0.9000 ",
        "ML.Dev.Prog.Sig.R lines 1907-1923 and 2501-2524, and historical ",
        "ML296_Fig5_Backend.R lines 279-305 and 373-377. The bundled row ",
        "records published_101_117_compatibility_candidate but no source DOI ",
        "or figure, so its exact source identity remains unverified; audited ",
        "2026-08-24."
      )
    } else if (learner_label %in% c(
      "GamBoost", "GlmBoost", "BlackBoost"
    )) {
      paste0(
        "Comparison references: historical ML296_Fig5_Backend.R lines ",
        "234-248 and 386-393, plus the later local ml296_all_in_one.R ",
        "mboost dispatcher and fixed parameter table. Official Mime1 ",
        "0.0.0.9000 ML.Dev.Prog.Sig contains no mboost learner. The bundled ",
        "row supplies no source DOI or figure, so its exact source identity ",
        "remains unverified; audited 2026-08-26."
      )
    } else {
      paste0(
        "Comparison reference only: npj Precision Oncology 2026 Table 1 ",
        "lists Random Forest (Ranger implementation), abbreviation Ranger, ",
        "and package ranger. The bundled row records ",
        "table1_catalog_candidate_endpoint_unverified but no source DOI or ",
        "figure; Table 1 supplies no operational settings, and official Mime1 ",
        "uses RSF/randomForestSRC rather than ranger. Source identity remains ",
        "unverified; audited 2026-08-24."
      )
    }
    return(list(
      status = if (ambiguous) {
        "audited_source_identity_ambiguous"
      } else {
        "audited_source_identity_unverified"
      },
      reference = reference,
      difference = difference
    ))
  }

  direct_binary_linear <- identical(task_type, "binary") &&
    identical(model_type, "single") && identical(selector, "none") &&
    (identical(learner_label, "LogisticRegression") ||
      !is.null(sb_candidate_glmnet_alpha(learner_label)))
  if (direct_binary_linear) {
    difference <- if (identical(learner_label, "LogisticRegression")) {
      paste0(
        "The catalog generator contributes only a normalized ",
        "LogisticRegression label; it does not specify an engine, formula ",
        "options, preprocessing, resampling, convergence handling, or ",
        "probability convention. STORMBREAKER fits an unpenalized ",
        "binomial-logit GLM on each analysis fold, requires a converged ",
        "full-rank fit with finite coefficients, and reports any warning or ",
        "invalid fit as a failure."
      )
    } else {
      paste0(
        "The catalog generator contributes only the Ridge, Lasso, or Enet ",
        "label and alpha; it does not specify an engine, lambda selection, ",
        "inner folds, seed, tuning measure, lambda rule, standardization, ",
        "preprocessing, or failure handling. STORMBREAKER cross-validates ",
        "the declared alpha inside each outer analysis fold using declared ",
        "inner folds, seed, measure, lambda rule, and standardization; any ",
        "supervised selector is refitted in the corresponding inner analysis ",
        "folds, and failures are reported without fallback."
      )
    }
    difference <- paste0(
      difference,
      " The catalog does not establish that this operational implementation ",
      "generated the candidate."
    )
    return(list(
      status = "audited_source_identity_unverified",
      reference = paste0(
        "Catalog construction reference: build_classification_candidate_space() ",
        "lines 609-681 and summary lines 1035-1044 in the 2026-07-18 ",
        "union-space generator, plus its algorithm-label provenance CSV. The ",
        "rows record union_generated_candidate_endpoint_unverified, ",
        "source_mime_binary=no, and no source DOI or figure. Official Mime1's ",
        "category module exposes neither LogisticRegression nor glmnet; source ",
        "identity remains unverified; audited 2026-08-24."
      ),
      difference = difference
    ))
  }

  mime_reference_binary <- identical(task_type, "binary") &&
    identical(model_type, "single") && identical(selector, "none") &&
    learner_label %in% c("NaiveBayes", "AdaBoost", "LogitBoost")
  if (!mime_reference_binary) {
    return(list(
      status = "not_yet_audited",
      reference = "",
      difference = ""
    ))
  }

  difference <- switch(learner_label,
    NaiveBayes = paste0(
      "Mime1 calls caret's nb method with usekernel=TRUE and tunes 25 ",
      "fL/adjust combinations by 5-fold, 10-repeat ROC. STORMBREAKER uses ",
      "fixed Gaussian Naive Bayes with empirical priors and a declared ",
      "variance-smoothing value, not the kernel model or Mime1 grid."
    ),
    AdaBoost = paste0(
      "Mime1 tunes 10 iteration/method combinations, including AdaBoost.M1 ",
      "and Real AdaBoost, by 5-fold, 10-repeat ROC. STORMBREAKER uses fixed ",
      "AdaBoost.M1 only and does not search the Mime1 grid."
    ),
    LogitBoost = paste0(
      "Mime1 tunes ten iteration counts from 11 through 101 by 5-fold, ",
      "10-repeat ROC. STORMBREAKER uses one caller-declared fixed iteration ",
      "count and does not search the Mime1 grid."
    )
  )
  difference <- paste0(
    difference,
    " Mime1's surrounding wrapper also intersects predictors across supplied ",
    "training and validation cohorts and mean-imputes each cohort; ",
    "STORMBREAKER does neither implicitly. The catalog does not establish ",
    "that Mime1 generated this candidate."
  )
  list(
    status = "audited_source_identity_unverified",
    reference = paste0(
      "Comparison reference only: Mime1 0.0.0.9000 ",
      "ML.Dev.Pred.Category.Sig.R lines 48-115 and 265-379. The catalog ",
      "rows provide no source DOI or figure and record ",
      "union_generated_candidate_endpoint_unverified, so source identity ",
      "remains unverified; audited 2026-08-24."
    ),
    difference = difference
  )
}

sb_candidate_learner_map <- function(label, task_type) {
  alpha <- sb_candidate_glmnet_alpha(label)
  if (!is.null(alpha)) {
    return(list(
      learner = if (identical(task_type, "survival")) {
        "glmnet_cox"
      } else {
        "glmnet_logistic"
      },
      parameters = sb_candidate_runtime_parameter("glmnet_alpha", alpha),
      configuration_required = TRUE,
      unsupported_shape = FALSE,
      missing = character()
    ))
  }

  if (identical(task_type, "survival")) {
    direction <- sb_candidate_stepcox_direction(label)
    if (!is.null(direction)) {
      return(list(
        learner = NA_character_,
        parameters = character(),
        configuration_required = FALSE,
        unsupported_shape = TRUE,
        missing = paste0("learner:", label)
      ))
    }
    map <- c(
      Coxph = "coxph",
      CoxBoost = "coxboost",
      plsRcox = "plsrcox",
      SuperPC = "superpc",
      `survival-SVM` = "survivalsvm",
      GBM = "gbm_cox",
      GamBoost = "gamboost_cox",
      GlmBoost = "glmboost_cox",
      BlackBoost = "blackboost_cox",
      Ranger = "ranger_survival",
      XGBoost = "xgboost_cox",
      SurvReg = "survreg_weibull",
      Rpart = "rpart_survival",
      RSF = "rsf_survival",
      Ctree = "ctree_survival",
      CForest = "cforest_survival"
    )
    configuration_required <- !identical(label, "Coxph")
  } else {
    direction <- sb_candidate_stepglm_direction(label)
    if (!is.null(direction)) {
      return(list(
        learner = "stepglm_binary",
        parameters = sb_candidate_runtime_parameter(
          "stepglm_direction",
          direction
        ),
        configuration_required = TRUE,
        unsupported_shape = FALSE,
        missing = character()
      ))
    }
    map <- c(
      LogisticRegression = "logistic",
      NaiveBayes = "gaussian_naive_bayes",
      AdaBoost = "adaboost_binary",
      LogitBoost = "logitboost_binary",
      LDA = "lda_binary",
      QDA = "qda_binary",
      KNN = "knn_binary",
      DecisionTree = "decision_tree_binary",
      RF = "random_forest_binary",
      SVM = "svm_radial_binary",
      svmRadialWeights = "svm_radial_weighted_binary",
      GBM = "gbm_binary",
      glmBoost = "glmboost_binary",
      plsRglm = "plsrglm_binary",
      XGBoost = "xgboost_binary"
    )
    configuration_required <- !identical(label, "LogisticRegression")
  }

  learner <- unname(map[label])
  if (length(learner) == 0L || is.na(learner)) {
    return(list(
      learner = NA_character_,
      parameters = character(),
      configuration_required = FALSE,
      unsupported_shape = FALSE,
      missing = paste0("learner:", label)
    ))
  }
  list(
    learner = learner,
    parameters = character(),
    configuration_required = configuration_required,
    unsupported_shape = FALSE,
    missing = character()
  )
}

sb_candidate_selector_map <- function(label, task_type) {
  if (!nzchar(label)) {
    return(list(
      selector = "none",
      parameters = character(),
      missing = character()
    ))
  }
  direction <- if (identical(task_type, "survival")) {
    sb_candidate_stepcox_direction(label)
  } else {
    NULL
  }
  if (!is.null(direction)) {
    return(list(
      selector = "stepcox",
      parameters = sb_candidate_runtime_parameter(
        "stepcox_direction",
        direction
      ),
      missing = character()
    ))
  }
  list(
    selector = NA_character_,
    parameters = character(),
    missing = paste0("selector:", label)
  )
}

sb_candidate_mapping_result <- function(task_type, route_status,
                                        admission_status, reason,
                                        selector = NA_character_,
                                        learner = NA_character_,
                                        parameters = NA_character_,
                                        route_available = FALSE,
                                        configuration_status = "not_applicable",
                                        source_spec_status = "not_applicable",
                                        source_spec_reference = "",
                                        source_spec_difference = "",
                                        missing = NA_character_) {
  data.frame(
    runtime_task_type = task_type,
    runtime_selector = selector,
    runtime_learner = learner,
    runtime_parameters = parameters,
    runtime_route_status = route_status,
    runtime_route_available = route_available,
    runtime_configuration_status = configuration_status,
    runtime_source_spec_status = source_spec_status,
    runtime_source_spec_reference = source_spec_reference,
    runtime_source_spec_difference = source_spec_difference,
    missing_runtime_components = missing,
    runtime_admission_status = admission_status,
    runtime_mapping_reason = reason,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

sb_candidate_mapping_row <- function(row, route_keys) {
  task_type <- if (identical(
    sb_candidate_catalog_text(row$task_type),
    "survival_right_censored"
  )) {
    "survival"
  } else {
    "binary"
  }
  admission_status <- sb_candidate_catalog_text(row$admission_status)
  if (grepl("^blocked", admission_status)) {
    return(sb_candidate_mapping_result(
      task_type,
      "methodologically_blocked",
      "not_admitted_methodologically_blocked",
      paste0(
        "The source catalog blocks this candidate until its stated ",
        "methodological requirement is implemented."
      )
    ))
  }
  if (identical(
    sb_candidate_catalog_text(row$role_semantics_status),
    "stage_roles_unresolved"
  )) {
    return(sb_candidate_mapping_result(
      task_type,
      "stage_roles_unresolved",
      "not_admitted_stage_roles_unresolved",
      paste0(
        "The source catalog does not resolve stage roles; STORMBREAKER ",
        "does not infer selector or learner roles from stage order."
      )
    ))
  }

  model_type <- sb_candidate_catalog_text(row$model_type)
  learner_label <- sb_candidate_catalog_text(row$learner)
  if (identical(model_type, "single") &&
    identical(task_type, "survival")) {
    direction <- sb_candidate_stepcox_direction(learner_label)
  } else {
    direction <- NULL
  }
  if (!is.null(direction)) {
    selector <- list(
      selector = "stepcox",
      parameters = sb_candidate_runtime_parameter(
        "stepcox_direction",
        direction
      ),
      missing = character()
    )
    learner <- list(
      learner = "coxph",
      parameters = character(),
      configuration_required = FALSE,
      unsupported_shape = FALSE,
      missing = character()
    )
  } else if (identical(model_type, "single")) {
    selector <- sb_candidate_selector_map("", task_type)
    learner <- sb_candidate_learner_map(learner_label, task_type)
  } else if (identical(model_type, "selector_learner")) {
    selector <- sb_candidate_selector_map(
      sb_candidate_catalog_text(row$selector_1),
      task_type
    )
    learner <- sb_candidate_learner_map(learner_label, task_type)
  } else {
    return(sb_candidate_mapping_result(
      task_type,
      "unsupported_pipeline_shape",
      "not_admitted_runtime_route_unavailable",
      "The current runtime does not represent this catalog pipeline shape."
    ))
  }

  parameters <- c(selector$parameters, learner$parameters)
  parameters <- if (length(parameters) == 0L) {
    ""
  } else {
    paste(parameters, collapse = "; ")
  }
  if (learner$unsupported_shape) {
    return(sb_candidate_mapping_result(
      task_type = task_type,
      route_status = "unsupported_pipeline_shape",
      admission_status = "not_admitted_runtime_route_unavailable",
      reason = paste0(
        "The catalog learner itself performs StepCox selection; the current ",
        "one-selector/one-learner grammar does not silently collapse that ",
        "additional stage."
      ),
      selector = selector$selector,
      parameters = parameters,
      missing = learner$missing
    ))
  }

  missing <- c(selector$missing, learner$missing)
  if (length(missing) > 0L) {
    return(sb_candidate_mapping_result(
      task_type = task_type,
      route_status = "component_unavailable",
      admission_status = "not_admitted_runtime_route_unavailable",
      reason = paste0(
        "No current runtime adapter is assigned to: ",
        paste(missing, collapse = "; "),
        "."
      ),
      selector = selector$selector,
      learner = learner$learner,
      parameters = parameters,
      missing = paste(missing, collapse = "; ")
    ))
  }

  route_key <- paste(
    task_type,
    selector$selector,
    learner$learner,
    sep = "::"
  )
  route_available <- route_key %in% route_keys
  if (!route_available) {
    return(sb_candidate_mapping_result(
      task_type = task_type,
      route_status = "route_not_admitted",
      admission_status = "not_admitted_runtime_route_unavailable",
      reason = paste0(
        "Mapped components exist, but their task-aware selector-learner ",
        "route is not admitted by the current runtime."
      ),
      selector = selector$selector,
      learner = learner$learner,
      parameters = parameters
    ))
  }

  source_spec <- sb_candidate_source_spec_audit(
    task_type,
    model_type,
    learner_label,
    selector$selector
  )
  source_identity_ambiguous <- identical(
    source_spec$status,
    "audited_source_identity_ambiguous"
  )
  source_identity_unverified <- identical(
    source_spec$status,
    "audited_source_identity_unverified"
  )

  sb_candidate_mapping_result(
    task_type = task_type,
    route_status = "route_available",
    admission_status = if (source_identity_ambiguous) {
      "not_admitted_source_identity_ambiguous"
    } else if (source_identity_unverified) {
      "not_admitted_source_identity_unverified"
    } else {
      "not_admitted_source_to_spec_unverified"
    },
    reason = if (source_identity_ambiguous) {
      paste0(
        "A current runtime route matches the task and mapped stage families, ",
        "but the shared historical label does not identify one ",
        "operational source specification."
      )
    } else if (source_identity_unverified) {
      paste0(
        "A current runtime route matches the task and named learner family, ",
        "but the candidate row does not identify one operational source ",
        "specification; any documented implementation remains only a ",
        "comparison reference."
      )
    } else {
      paste0(
        "A current runtime route matches this candidate's task and mapped ",
        "stage families, but source-to-spec semantics and any required ",
        "configuration are not yet verified."
      )
    },
    selector = selector$selector,
    learner = learner$learner,
    parameters = parameters,
    route_available = TRUE,
    configuration_status = if (learner$configuration_required) {
      "additional_parameters_required"
    } else {
      "no_additional_runtime_parameters"
    },
    source_spec_status = source_spec$status,
    source_spec_reference = source_spec$reference,
    source_spec_difference = source_spec$difference
  )
}

sb_candidate_registry_data <- function() {
  catalog <- utils::read.csv(
    system.file(
      "extdata",
      "unified_pipeline_registry_1213_20260806.csv",
      package = "stormbreaker"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  routes <- sb_pipeline_registry_data()
  route_keys <- paste(
    routes$task_type,
    routes$selector,
    routes$learner,
    sep = "::"
  )
  mapping <- lapply(seq_len(nrow(catalog)), function(index) {
    sb_candidate_mapping_row(catalog[index, , drop = FALSE], route_keys)
  })
  mapping <- do.call(rbind, mapping)
  rownames(mapping) <- NULL
  result <- data.frame(
    catalog,
    mapping,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(result) <- NULL
  result
}

sb_candidate_route_statuses <- function() {
  c(
    "route_available",
    "component_unavailable",
    "route_not_admitted",
    "unsupported_pipeline_shape",
    "stage_roles_unresolved",
    "methodologically_blocked"
  )
}

#' Inspect the historical candidate catalog and current runtime mapping
#'
#' `sb_candidate_registry()` returns the bundled task-aware catalog of 1,175 retained
#' historical candidate specifications and appends a descriptive mapping to
#' STORMBREAKER's current runtime. The source columns are preserved; in
#' particular, historical adapter, execution, and admission fields are not
#' upgraded because a similarly named runtime route now exists.
#'
#' `runtime_route_status = "route_available"` means only that the current
#' package has a task-matched selector and learner family route capable of
#' representing the catalog shape. It does not verify that the source used the
#' same operational definition, package, defaults, tuning procedure, or
#' parameters. Every currently mapped row therefore remains non-admitted. The
#' appended `runtime_source_spec_status` separates routes that have not yet
#' been audited from audits that found a name-level ambiguity or operational
#' semantics that differ; it never promotes a family match to an executable
#' candidate. The function never constructs, fits, ranks, tunes, or selects a
#' pipeline.
#'
#' The first source audit covers available routes that contain the StepCox
#' selector. The historical catalog merges the label across sources with
#' different operational definitions: official Mime1 and the historical Fig.5A
#' reconstruction differ in their common prefilter, forward start, minimum-
#' feature, warning, truncation, and fallback behavior. STORMBREAKER's current
#' fold-local AIC selector is therefore not silently called an exact historical
#' implementation. These rows are marked
#' `audited_source_identity_ambiguous` and remain
#' `not_admitted_source_identity_ambiguous`. This records provenance; it does
#' not change the admitted StepCox runtime.
#'
#' The second source audit uses official Mime1 only as an operational comparison
#' for the direct binary `NaiveBayes`, `AdaBoost`, and `LogitBoost` labels. Mime1
#' tunes each method by 5-fold, 10-repeat ROC. Its Naive Bayes route uses a
#' kernel model, its AdaBoost grid includes both AdaBoost.M1 and Real AdaBoost,
#' and its LogitBoost route searches iteration counts. STORMBREAKER instead
#' exposes fixed Gaussian Naive Bayes, fixed AdaBoost.M1, and fixed LogitBoost
#' methods without those searches. Because the catalog rows provide no source
#' DOI or figure and identify themselves as union-generated, this comparison
#' does not attribute them to Mime1. They are marked
#' `audited_source_identity_unverified` and remain
#' `not_admitted_source_identity_unverified`. The three current STORMBREAKER
#' methods remain available under their own explicit specifications.
#'
#' The third source audit covers the 11 direct survival `Ridge`, `Lasso`, and
#' elastic-net labels. Their catalog evidence is only
#' `historical_fig5a_label_inferred`, and the two available historical
#' implementations differ: official Mime1 uses fixed 10-fold cross-validation,
#' while the Fig.5A reconstruction uses a sample-size-dependent fold count.
#' Mime1's Ridge code also fits `alpha = 0` but obtains `lambda.min` from a
#' separate cross-validation call that omits `alpha`; the reconstruction
#' cross-validates `alpha = 0`. STORMBREAKER requires explicit alpha, inner
#' folds, seed, measure, lambda rule, and standardization inside outer analysis
#' resampling. The rows are therefore marked
#' `audited_source_identity_ambiguous`, remain non-admitted, and do not change
#' the current glmnet Cox route.
#'
#' The fourth source audit covers the 12 direct binary `LogisticRegression`,
#' `Ridge`, `Lasso`, and elastic-net labels. The catalog generator created them
#' as a label union with `union_generated_candidate_endpoint_unverified`; its
#' provenance marks these families outside the Mime1 binary method set, and the
#' rows provide no DOI or figure. The labels do not specify ordinary-logistic
#' convergence semantics or glmnet lambda selection, folds, seed, objective,
#' rule, standardization, preprocessing, or failures. They are marked
#' `audited_source_identity_unverified` and remain non-admitted. The current
#' explicit logistic and nested glmnet-logistic routes are unchanged.
#'
#' The fifth source audit covers nine direct fixed survival labels.
#' `survival-SVM`, `CoxBoost`, `GBM`, and `SuperPC` carry only
#' `historical_fig5a_label_inferred`; official Mime1 and the historical Fig.5A
#' reconstruction differ in tuning, retry, preprocessing, or score semantics,
#' so these four rows are `audited_source_identity_ambiguous`. Direct `plsRcox`
#' is a 101/117 compatibility addition outside the reconstructed 296
#' single-model set, and `Ranger` is supported only by a Table 1 family/package
#' label while Mime1 uses RSF/randomForestSRC. Those two rows are
#' `audited_source_identity_unverified`. All six remain non-admitted. The
#' current fixed learners keep their own explicit parameters and are not
#' presented as exact reproductions of the catalog rows. The same audit records
#' that official Mime1 contains no GamBoost, GlmBoost, or BlackBoost survival
#' learner; local reconstructions provide comparison implementations, including
#' one incorrect GamBoost-to-glmboost dispatch. STORMBREAKER therefore exposes
#' three distinct native mboost algorithms while keeping their catalog rows
#' source-unverified.
#'
#' Seven direct binary labels (`LDA`, `QDA`, `KNN`, `DecisionTree`, `RF`, `SVM`,
#' and `GBM`) currently have family-level runtime routes but have not yet
#' received source-to-specification audits. They therefore carry
#' `runtime_source_spec_status = "not_yet_audited"` together with
#' `runtime_admission_status = "not_admitted_source_to_spec_unverified"`.
#' Route availability for these labels is not a claim of exact historical
#' reproduction.
#'
#' Fixed runtime implementations now also create descriptive family matches for
#' binary and survival `XGBoost`, `RSF`, `Rpart`, `Ctree`, `CForest`, and the
#' compatible `SurvReg` selector routes. Direct matches without a completed
#' source audit remain `not_yet_audited`; StepCox-containing matches retain the
#' StepCox source-identity decision. No such row is promoted to an executable
#' historical specification.
#'
#' Unresolved multi-stage labels remain unresolved rather than having roles
#' inferred from their order. Methodologically blocked rows remain blocked.
#' Names such as RSF and RF are not silently equated with ranger, and binary GBM
#' is not equated with LightGBM or CatBoost.
#'
#' @param task_type Optional task filter: `NULL`, `"survival"`, or `"binary"`.
#' @param route_status Optional exact runtime-route status returned by the
#'   function. `NULL` retains every status.
#'
#' @return A data frame retaining every source catalog column and adding the
#'   current runtime task, mapped selector and learner, parameters declared by
#'   the label, route status and availability, configuration status, source-
#'   specification audit status/reference/difference, missing components,
#'   current non-admission status, and mapping reason.
#' @export
#'
#' @examples
#' candidates <- sb_candidate_registry()
#' table(candidates$runtime_task_type, candidates$runtime_route_status)
#' sb_candidate_registry(
#'   task_type = "survival",
#'   route_status = "route_available"
#' )
sb_candidate_registry <- function(task_type = NULL, route_status = NULL) {
  if (!is.null(task_type) && (!is.character(task_type) ||
    length(task_type) != 1L || is.na(task_type) ||
    !task_type %in% c("survival", "binary"))) {
    sb_abort(
      "`task_type` must be `NULL`, `\"survival\"`, or `\"binary\"`.",
      "stormbreaker_bad_candidate_registry_task"
    )
  }
  if (!is.null(route_status) && (!is.character(route_status) ||
    length(route_status) != 1L || is.na(route_status) ||
    !route_status %in% sb_candidate_route_statuses())) {
    sb_abort(
      "`route_status` must be `NULL` or one documented route status.",
      "stormbreaker_bad_candidate_registry_route_status"
    )
  }

  registry <- sb_candidate_registry_data()
  if (!is.null(task_type)) {
    registry <- registry[
      registry$runtime_task_type == task_type, ,
      drop = FALSE
    ]
  }
  if (!is.null(route_status)) {
    registry <- registry[
      registry$runtime_route_status == route_status, ,
      drop = FALSE
    ]
  }
  rownames(registry) <- NULL
  registry
}
