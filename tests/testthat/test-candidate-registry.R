test_that("candidate registry preserves the retained historical catalog", {
  registry <- sb_candidate_registry()

  expect_s3_class(registry, "data.frame")
  expect_identical(nrow(registry), 1175L)
  task_counts <- table(registry$runtime_task_type)
  expect_identical(names(task_counts), c("binary", "survival"))
  expect_identical(as.integer(task_counts), c(584L, 591L))
  expect_identical(anyDuplicated(registry$pipeline_id), 0L)
  expect_identical(registry$pipeline_id[[1L]], "surv_0001")
  expect_identical(registry$pipeline_id[[1175L]], "class_0584")
  expect_identical(unique(registry$adapter_status), "not_verified")
  expect_identical(any(registry$execution_status == "executable"), FALSE)
  expect_identical(
    all(grepl("^(not_admitted|blocked)", registry$admission_status)),
    TRUE
  )
  expect_identical(
    any(registry$runtime_admission_status == "admitted"),
    FALSE
  )
  expect_identical(sum(registry$runtime_route_available), 138L)
})

test_that("candidate registry maps only explicit current runtime roles", {
  registry <- sb_candidate_registry()
  row_by_id <- function(id) {
    registry[registry$pipeline_id == id, , drop = FALSE]
  }

  logistic <- row_by_id("class_0001")
  adaboost <- row_by_id("class_0029")
  stepcox <- row_by_id("surv_0013")
  xgboost <- row_by_id("class_0024")
  decision_tree <- registry[
    registry$task_type == "binary_classification" &
      registry$pipeline_label == "DecisionTree" &
      registry$model_type == "single", ,
    drop = FALSE
  ]
  stepcox_ridge <- registry[
    registry$task_type == "survival_right_censored" &
      registry$pipeline_label == "StepCox[forward] + Ridge", ,
    drop = FALSE
  ]
  blocked <- registry[grepl("^blocked", registry$admission_status), ,
    drop = FALSE
  ][1L, , drop = FALSE]
  unresolved <- registry[
    registry$role_semantics_status == "stage_roles_unresolved", ,
    drop = FALSE
  ][1L, , drop = FALSE]
  nested_stepcox <- registry[
    registry$selector_1 == "StepCox[forward]" &
      registry$learner == "StepCox[both]", ,
    drop = FALSE
  ]

  expect_identical(logistic$runtime_selector, "none")
  expect_identical(logistic$runtime_learner, "logistic")
  expect_identical(logistic$runtime_route_status, "route_available")
  expect_identical(
    logistic$runtime_configuration_status,
    "no_additional_runtime_parameters"
  )
  expect_identical(adaboost$runtime_learner, "adaboost_binary")
  expect_identical(
    adaboost$runtime_configuration_status,
    "additional_parameters_required"
  )
  expect_identical(stepcox$runtime_selector, "stepcox")
  expect_identical(stepcox$runtime_learner, "coxph")
  expect_identical(
    stepcox$runtime_parameters,
    "stepcox_direction=forward"
  )
  expect_identical(nrow(stepcox_ridge), 1L)
  expect_identical(stepcox_ridge$runtime_route_status, "route_available")
  expect_identical(stepcox_ridge$runtime_learner, "glmnet_cox")
  expect_identical(
    stepcox_ridge$runtime_parameters,
    "stepcox_direction=forward; glmnet_alpha=0"
  )
  expect_identical(xgboost$runtime_route_status, "route_available")
  expect_identical(xgboost$runtime_learner, "xgboost_binary")
  expect_true(is.na(xgboost$missing_runtime_components))
  expect_identical(nrow(decision_tree), 1L)
  expect_identical(
    decision_tree$runtime_learner,
    "decision_tree_binary"
  )
  expect_identical(decision_tree$runtime_route_status, "route_available")
  expect_identical(blocked$runtime_route_status, "methodologically_blocked")
  expect_identical(blocked$runtime_route_available, FALSE)
  expect_identical(unresolved$runtime_route_status, "stage_roles_unresolved")
  expect_identical(unresolved$runtime_route_available, FALSE)
  expect_identical(nrow(nested_stepcox), 1L)
  expect_identical(
    nested_stepcox$runtime_route_status,
    "unsupported_pipeline_shape"
  )
})

test_that("new specialized binary labels map to real runtime learners", {
  registry <- sb_candidate_registry()
  labels <- c(
    "Stepglm[forward]", "Stepglm[backward]", "Stepglm[both]",
    "svmRadialWeights", "glmBoost", "plsRglm"
  )
  direct <- registry[
    registry$runtime_task_type == "binary" &
      registry$model_type == "single" &
      registry$pipeline_label %in% labels, ,
    drop = FALSE
  ]
  expected_learners <- c(
    `Stepglm[forward]` = "stepglm_binary",
    `Stepglm[backward]` = "stepglm_binary",
    `Stepglm[both]` = "stepglm_binary",
    svmRadialWeights = "svm_radial_weighted_binary",
    glmBoost = "glmboost_binary",
    plsRglm = "plsrglm_binary"
  )

  expect_identical(nrow(direct), 6L)
  expect_identical(sort(direct$pipeline_label), sort(labels))
  expect_identical(
    direct$runtime_learner,
    unname(expected_learners[direct$pipeline_label])
  )
  expect_identical(unique(direct$runtime_selector), "none")
  expect_identical(unique(direct$runtime_route_status), "route_available")
  expect_identical(unique(direct$runtime_route_available), TRUE)
  expect_identical(
    direct$runtime_parameters[grepl("^Stepglm", direct$pipeline_label)],
    paste0(
      "stepglm_direction=",
      sub("^Stepglm\\[(.*)\\]$", "\\1", direct$pipeline_label[
        grepl("^Stepglm", direct$pipeline_label)
      ])
    )
  )
  expect_identical(
    unique(direct$runtime_admission_status),
    "not_admitted_source_to_spec_unverified"
  )
  expect_identical(any(direct$execution_status == "executable"), FALSE)
})

test_that("three mboost survival labels map to distinct real learners", {
  registry <- sb_candidate_registry()
  labels <- c("GamBoost", "GlmBoost", "BlackBoost")
  direct <- registry[
    registry$runtime_task_type == "survival" &
      registry$model_type == "single" &
      registry$pipeline_label %in% labels, ,
    drop = FALSE
  ]
  expected <- c(
    GamBoost = "gamboost_cox",
    GlmBoost = "glmboost_cox",
    BlackBoost = "blackboost_cox"
  )

  expect_identical(nrow(direct), 3L)
  expect_identical(
    direct$runtime_learner,
    unname(expected[direct$pipeline_label])
  )
  expect_identical(unique(direct$runtime_route_status), "route_available")
  expect_identical(unique(direct$runtime_route_available), TRUE)
  expect_identical(
    unique(direct$runtime_source_spec_status),
    "audited_source_identity_unverified"
  )
  expect_identical(
    unique(direct$runtime_admission_status),
    "not_admitted_source_identity_unverified"
  )
  expect_true(grepl(
    "distinct native gamboost",
    direct$runtime_source_spec_difference[direct$pipeline_label == "GamBoost"],
    fixed = TRUE
  ))
})

test_that("available candidate routes refer only to admitted runtime templates", {
  candidates <- sb_candidate_registry()
  routes <- sb_pipeline_registry()
  available <- candidates[candidates$runtime_route_available, , drop = FALSE]
  candidate_keys <- paste(
    available$runtime_task_type,
    available$runtime_selector,
    available$runtime_learner,
    sep = "::"
  )
  route_keys <- paste(
    routes$task_type,
    routes$selector,
    routes$learner,
    sep = "::"
  )

  expect_gt(nrow(available), 0L)
  expect_identical(all(candidate_keys %in% route_keys), TRUE)
  expect_identical(
    sort(unique(available$runtime_admission_status)),
    sort(c(
      "not_admitted_source_identity_ambiguous",
      "not_admitted_source_identity_unverified",
      "not_admitted_source_to_spec_unverified"
    ))
  )
  expect_identical(
    any(available$execution_status == "executable"),
    FALSE
  )
})

test_that("StepCox source identity ambiguity is explicit without admission", {
  registry <- sb_candidate_registry()
  available <- registry$runtime_route_available
  stepcox <- available & registry$runtime_selector == "stepcox"
  mime_binary <- available &
    registry$runtime_task_type == "binary" &
    registry$model_type == "single" &
    registry$runtime_selector == "none" &
    registry$pipeline_label %in% c("NaiveBayes", "AdaBoost", "LogitBoost")
  survival_glmnet <- available &
    registry$runtime_task_type == "survival" &
    registry$model_type == "single" &
    registry$runtime_selector == "none" &
    registry$runtime_learner == "glmnet_cox"
  binary_linear <- available &
    registry$runtime_task_type == "binary" &
    registry$model_type == "single" &
    registry$runtime_selector == "none" &
    registry$runtime_learner %in% c("logistic", "glmnet_logistic")
  direct_survival_fixed <- available &
    registry$runtime_task_type == "survival" &
    registry$model_type == "single" &
    registry$runtime_selector == "none" &
    registry$runtime_learner %in% c(
      "survivalsvm",
      "coxboost",
      "gbm_cox",
      "superpc",
      "plsrcox",
      "ranger_survival",
      sb_mboost_cox_learner_ids()
    )
  new_binary <- available &
    registry$runtime_task_type == "binary" &
    registry$model_type == "single" &
    registry$runtime_selector == "none" &
    registry$runtime_learner %in% sb_classical_binary_learner_ids()
  new_specialized_binary <- available &
    registry$runtime_task_type == "binary" &
    registry$model_type == "single" &
    registry$runtime_selector == "none" &
    registry$runtime_learner %in%
      sb_specialized_binary_linear_learner_ids()
  new_expansion <- available &
    registry$model_type == "single" &
    registry$runtime_selector == "none" &
    registry$runtime_learner %in% c(
      "xgboost_binary", "xgboost_cox",
      "survreg_weibull", "rpart_survival",
      "rsf_survival", "ctree_survival",
      "cforest_survival"
    )
  former_python <- registry$learner %in% c(
      "CoxTime", "DeepHit", "DeepSurv", "Logistic-Hazard", "PC-Hazard"
    )
  other_available <- available & !stepcox & !mime_binary &
    !survival_glmnet & !binary_linear & !direct_survival_fixed & !new_binary &
    !new_specialized_binary & !new_expansion
  unavailable <- !available

  expect_gt(sum(stepcox), 0L)
  expect_identical(sum(other_available), 0L)
  expect_identical(sum(former_python), 60L)
  expect_true(all(!registry$runtime_route_available[former_python]))
  expect_true(all(
    registry$runtime_route_status[former_python] %in%
      c("component_unavailable", "methodologically_blocked")
  ))
  former_unblocked <- former_python &
    registry$runtime_route_status == "component_unavailable"
  expect_true(all(
    registry$runtime_admission_status[former_unblocked] ==
      "not_admitted_runtime_route_unavailable"
  ))
  expect_true(all(grepl(
    "learner:",
    registry$missing_runtime_components[former_unblocked]
  )))
  expect_identical(sum(new_binary), 8L)
  expect_identical(
    sort(registry$pipeline_label[new_binary]),
    sort(c(
      "DecisionTree", "GBM", "KNN", "LDA", "QDA", "RF", "SVM",
      "svmRadialWeights"
    ))
  )
  expect_identical(sum(new_specialized_binary), 5L)
  expect_identical(
    sort(registry$pipeline_label[new_specialized_binary]),
    sort(c(
      "Stepglm[forward]", "Stepglm[backward]", "Stepglm[both]",
      "glmBoost", "plsRglm"
    ))
  )
  expect_identical(sum(new_expansion), 6L)
  expect_identical(
    unique(registry$runtime_source_spec_status[new_expansion]),
    "not_yet_audited"
  )
  expect_identical(
    unique(registry$runtime_source_spec_status[stepcox]),
    "audited_source_identity_ambiguous"
  )
  expect_identical(
    unique(registry$runtime_source_spec_status[new_binary]),
    "not_yet_audited"
  )
  expect_identical(
    unique(registry$runtime_source_spec_status[new_specialized_binary]),
    "not_yet_audited"
  )
  expect_identical(
    unique(registry$runtime_source_spec_status[unavailable]),
    "not_applicable"
  )
  expect_identical(
    unique(registry$runtime_admission_status[stepcox]),
    "not_admitted_source_identity_ambiguous"
  )
  expect_identical(
    all(nzchar(registry$runtime_source_spec_reference[stepcox])),
    TRUE
  )
  expect_identical(
    all(grepl(
      "Mime1.*Fig.5A",
      registry$runtime_source_spec_reference[stepcox]
    )),
    TRUE
  )
  expect_identical(
    all(grepl(
      "forward.*null",
      registry$runtime_source_spec_difference[stepcox]
    )),
    TRUE
  )
  expect_identical(
    all(registry$adapter_status[stepcox] == "not_verified"),
    TRUE
  )
  expect_identical(
    any(registry$execution_status[stepcox] == "executable"),
    FALSE
  )
})

test_that("direct survival glmnet source identity ambiguity is explicit", {
  registry <- sb_candidate_registry()
  labels <- c(
    "Ridge",
    sprintf("Enet[alpha=%.1f]", seq(0.1, 0.9, by = 0.1)),
    "Lasso"
  )
  direct <- registry[
    registry$runtime_route_available &
      registry$runtime_task_type == "survival" &
      registry$model_type == "single" &
      registry$runtime_selector == "none" &
      registry$runtime_learner == "glmnet_cox", ,
    drop = FALSE
  ]

  expect_identical(nrow(direct), 11L)
  expect_identical(sort(direct$pipeline_label), sort(labels))
  expect_identical(anyDuplicated(direct$pipeline_label), 0L)
  expect_identical(
    unique(direct$evidence_status),
    "historical_fig5a_label_inferred"
  )
  expect_identical(
    unique(direct$runtime_source_spec_status),
    "audited_source_identity_ambiguous"
  )
  expect_identical(
    unique(direct$runtime_admission_status),
    "not_admitted_source_identity_ambiguous"
  )
  expect_identical(
    all(grepl(
      "Mime1.*historical Fig.5A.*historical_fig5a_label_inferred",
      direct$runtime_source_spec_reference
    )),
    TRUE
  )
  expect_identical(
    all(grepl(
      "does not establish which implementation generated",
      direct$runtime_source_spec_difference
    )),
    TRUE
  )

  ridge <- direct[direct$pipeline_label == "Ridge", , drop = FALSE]
  lasso <- direct[direct$pipeline_label == "Lasso", , drop = FALSE]
  enet <- direct[grepl("^Enet", direct$pipeline_label), , drop = FALSE]
  expect_match(ridge$runtime_source_spec_difference, "omits alpha")
  expect_match(ridge$runtime_source_spec_difference, "default alpha=1")
  expect_match(ridge$runtime_source_spec_difference, "cross-validates alpha=0")
  expect_match(lasso$runtime_source_spec_difference, "fixed 10 folds")
  expect_match(lasso$runtime_source_spec_difference, "sample-size-dependent")
  expect_identical(
    all(grepl(
      "fixed 10 folds.*sample-size-dependent",
      enet$runtime_source_spec_difference
    )),
    TRUE
  )
  expect_identical(
    all(grepl(
      paste0(
        "declared inner folds, seed, measure, lambda rule, and ",
        "standardization"
      ),
      direct$runtime_source_spec_difference
    )),
    TRUE
  )
  expect_identical(unique(direct$adapter_status), "not_verified")
  expect_identical(any(direct$execution_status == "executable"), FALSE)
})

test_that("direct binary linear source identity remains unverified", {
  registry <- sb_candidate_registry()
  labels <- c(
    "LogisticRegression",
    "Ridge",
    "Lasso",
    sprintf("Enet[alpha=%.1f]", seq(0.1, 0.9, by = 0.1))
  )
  direct <- registry[
    registry$runtime_route_available &
      registry$runtime_task_type == "binary" &
      registry$model_type == "single" &
      registry$runtime_selector == "none" &
      registry$runtime_learner %in% c("logistic", "glmnet_logistic"), ,
    drop = FALSE
  ]

  expect_identical(nrow(direct), 12L)
  expect_identical(sort(direct$pipeline_label), sort(labels))
  expect_identical(anyDuplicated(direct$pipeline_label), 0L)
  expect_identical(
    unique(direct$evidence_status),
    "union_generated_candidate_endpoint_unverified"
  )
  expect_identical(all(direct$source_doi == ""), TRUE)
  expect_identical(all(direct$source_figure == ""), TRUE)
  expect_identical(
    unique(direct$runtime_source_spec_status),
    "audited_source_identity_unverified"
  )
  expect_identical(
    unique(direct$runtime_admission_status),
    "not_admitted_source_identity_unverified"
  )
  expect_identical(
    all(grepl(
      paste0(
        "Catalog construction reference:.*",
        "union_generated_candidate_endpoint_unverified.*",
        "source_mime_binary=no"
      ),
      direct$runtime_source_spec_reference
    )),
    TRUE
  )
  expect_identical(
    all(grepl(
      "Official Mime1.*neither LogisticRegression nor glmnet",
      direct$runtime_source_spec_reference
    )),
    TRUE
  )
  expect_identical(
    all(grepl(
      "does not establish that this operational implementation generated",
      direct$runtime_source_spec_difference
    )),
    TRUE
  )

  logistic <- direct[direct$runtime_learner == "logistic", , drop = FALSE]
  glmnet <- direct[direct$runtime_learner == "glmnet_logistic", , drop = FALSE]
  expect_identical(nrow(logistic), 1L)
  expect_identical(nrow(glmnet), 11L)
  expect_match(
    logistic$runtime_source_spec_difference,
    "unpenalized binomial-logit"
  )
  expect_match(logistic$runtime_source_spec_difference, "converged full-rank")
  expect_match(logistic$runtime_source_spec_difference, "warning or invalid fit")
  expect_identical(
    all(grepl(
      "label and alpha.*does not specify.*lambda selection",
      glmnet$runtime_source_spec_difference
    )),
    TRUE
  )
  expect_identical(
    all(grepl(
      paste0(
        "declared inner folds, seed, measure, lambda rule, and ",
        "standardization"
      ),
      glmnet$runtime_source_spec_difference
    )),
    TRUE
  )
  expect_identical(
    all(grepl(
      "inside each outer analysis fold.*without fallback",
      glmnet$runtime_source_spec_difference
    )),
    TRUE
  )
  expect_identical(unique(direct$adapter_status), "not_verified")
  expect_identical(any(direct$execution_status == "executable"), FALSE)
})

test_that("remaining direct survival source identities are explicit", {
  registry <- sb_candidate_registry()
  learner_map <- c(
    `survival-SVM` = "survivalsvm",
    CoxBoost = "coxboost",
    GBM = "gbm_cox",
    SuperPC = "superpc",
    plsRcox = "plsrcox",
    Ranger = "ranger_survival"
  )
  direct <- registry[
    registry$runtime_route_available &
      registry$runtime_task_type == "survival" &
      registry$model_type == "single" &
      registry$runtime_selector == "none" &
      registry$runtime_learner %in% unname(learner_map), ,
    drop = FALSE
  ]

  expect_identical(nrow(direct), 6L)
  expect_identical(sort(direct$pipeline_label), sort(names(learner_map)))
  expect_identical(anyDuplicated(direct$pipeline_label), 0L)
  expect_identical(
    unname(learner_map[direct$pipeline_label]),
    direct$runtime_learner
  )

  inferred_labels <- c("survival-SVM", "CoxBoost", "GBM", "SuperPC")
  inferred <- direct[direct$pipeline_label %in% inferred_labels, , drop = FALSE]
  unverified <- direct[
    direct$pipeline_label %in% c("plsRcox", "Ranger"), ,
    drop = FALSE
  ]
  expect_identical(nrow(inferred), 4L)
  expect_identical(
    unique(inferred$evidence_status),
    "historical_fig5a_label_inferred"
  )
  expect_identical(
    unique(inferred$runtime_source_spec_status),
    "audited_source_identity_ambiguous"
  )
  expect_identical(
    unique(inferred$runtime_admission_status),
    "not_admitted_source_identity_ambiguous"
  )
  expect_identical(
    unverified$evidence_status[match(
      c("plsRcox", "Ranger"),
      unverified$pipeline_label
    )],
    c(
      "published_101_117_compatibility_candidate",
      "table1_catalog_candidate_endpoint_unverified"
    )
  )
  expect_identical(
    unique(unverified$runtime_source_spec_status),
    "audited_source_identity_unverified"
  )
  expect_identical(
    unique(unverified$runtime_admission_status),
    "not_admitted_source_identity_unverified"
  )

  by_label <- function(label) {
    direct[direct$pipeline_label == label, , drop = FALSE]
  }
  coxboost <- by_label("CoxBoost")
  gbm <- by_label("GBM")
  superpc <- by_label("SuperPC")
  survival_svm <- by_label("survival-SVM")
  plsrcox <- by_label("plsRcox")
  ranger <- by_label("Ranger")

  expect_match(coxboost$runtime_source_spec_difference, "maxstepno=500")
  expect_match(coxboost$runtime_source_spec_difference, "maxstepno=200")
  expect_match(coxboost$runtime_source_spec_difference, "fixed step count")
  expect_match(gbm$runtime_source_spec_difference, "10,000")
  expect_match(gbm$runtime_source_spec_difference, "2,000")
  expect_match(gbm$runtime_source_spec_difference, "no early stopping")
  expect_match(superpc$runtime_source_spec_difference, "retries indefinitely")
  expect_match(superpc$runtime_source_spec_difference, "validation outcomes")
  expect_match(superpc$runtime_source_spec_difference, "fixed threshold")
  expect_match(survival_svm$runtime_source_spec_difference, "package defaults")
  expect_match(survival_svm$runtime_source_spec_difference, "raw predicted")
  expect_match(survival_svm$runtime_source_spec_difference, "Cox calibration")
  expect_match(plsrcox$runtime_source_spec_reference, "101/117 compatibility")
  expect_match(plsrcox$runtime_source_spec_difference, "cv.plsRcox")
  expect_match(plsrcox$runtime_source_spec_difference, "fixed component count")
  expect_match(ranger$runtime_source_spec_reference, "Table 1")
  expect_match(ranger$runtime_source_spec_difference, "randomForestSRC")
  expect_match(ranger$runtime_source_spec_difference, "not ranger")

  expect_identical(unique(direct$adapter_status), "not_verified")
  expect_identical(any(direct$execution_status == "executable"), FALSE)
})

test_that("Mime1 is only a reference when binary source identity is unverified", {
  registry <- sb_candidate_registry()
  labels <- c("NaiveBayes", "AdaBoost", "LogitBoost")
  direct <- registry[
    registry$runtime_route_available &
      registry$runtime_task_type == "binary" &
      registry$model_type == "single" &
      registry$runtime_selector == "none" &
      registry$pipeline_label %in% labels, ,
    drop = FALSE
  ]

  expect_identical(sort(direct$pipeline_label), sort(labels))
  expect_identical(anyDuplicated(direct$pipeline_label), 0L)
  expect_identical(
    unique(direct$runtime_source_spec_status),
    "audited_source_identity_unverified"
  )
  expect_identical(
    unique(direct$runtime_admission_status),
    "not_admitted_source_identity_unverified"
  )
  expect_identical(
    all(grepl(
      "Comparison reference only: Mime1.*source identity remains unverified",
      direct$runtime_source_spec_reference
    )),
    TRUE
  )

  naive_bayes <- direct[direct$pipeline_label == "NaiveBayes", , drop = FALSE]
  adaboost <- direct[direct$pipeline_label == "AdaBoost", , drop = FALSE]
  logitboost <- direct[direct$pipeline_label == "LogitBoost", , drop = FALSE]
  expect_match(naive_bayes$runtime_source_spec_difference, "usekernel=TRUE")
  expect_match(naive_bayes$runtime_source_spec_difference, "Gaussian")
  expect_match(adaboost$runtime_source_spec_difference, "Real AdaBoost")
  expect_match(adaboost$runtime_source_spec_difference, "fixed AdaBoost.M1")
  expect_match(logitboost$runtime_source_spec_difference, "11 through 101")
  expect_match(logitboost$runtime_source_spec_difference, "fixed iteration")
  expect_identical(
    all(grepl(
      "does not establish that Mime1 generated this candidate",
      direct$runtime_source_spec_difference
    )),
    TRUE
  )
  expect_identical(unique(direct$adapter_status), "not_verified")
  expect_identical(
    any(direct$execution_status == "executable"),
    FALSE
  )
})

test_that("candidate registry filters one task or route status", {
  all_candidates <- sb_candidate_registry()
  binary <- sb_candidate_registry(task_type = "binary")
  available <- sb_candidate_registry(route_status = "route_available")
  bad_task <- capture_sb_error(sb_candidate_registry(
    task_type = "diagnostic"
  ))
  bad_status <- capture_sb_error(sb_candidate_registry(
    route_status = "executable"
  ))

  expect_identical(nrow(binary), 584L)
  expect_identical(unique(binary$runtime_task_type), "binary")
  expect_identical(
    nrow(available),
    sum(all_candidates$runtime_route_status == "route_available")
  )
  expect_identical(unique(available$runtime_route_status), "route_available")
  expect_s3_class(bad_task, "stormbreaker_bad_candidate_registry_task")
  expect_s3_class(
    bad_status,
    "stormbreaker_bad_candidate_registry_route_status"
  )
})
