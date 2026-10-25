test_that("MIME117 generator returns 117 real admitted pipelines", {
  suite <- sb_mime117_pipelines(
    inner_folds = 3L,
    inner_seed = 8701L,
    selector_inner_folds = 3L,
    selector_seed = 8702L,
    rsf_nodesize = 5L
  )

  expect_s3_class(suite, "data.frame")
  expect_identical(nrow(suite), 117L)
  expect_identical(suite$official_order, seq_len(117L))
  expect_identical(anyDuplicated(suite$config_id), 0L)
  expect_identical(anyDuplicated(suite$official_model_label), 0L)
  expect_identical(
    anyDuplicated(vapply(suite$pipeline, `[[`, character(1), "id")),
    0L
  )
  expect_identical(
    as.integer(table(factor(
      suite$block,
      levels = c(
        "RSF", "Enet", "StepCox", "CoxBoost", "plsRcox", "SuperPC",
        "GBM", "survival-SVM", "Ridge", "Lasso"
      )
    ))),
    c(20L, 9L, 54L, 19L, 1L, 1L, 1L, 1L, 1L, 10L)
  )
  expect_identical(
    as.integer(table(suite$pipeline_type)),
    c(20L, 97L)
  )
  expect_true(all(vapply(
    suite$pipeline,
    inherits,
    logical(1),
    what = "stormbreaker_pipeline"
  )))

  lint_status <- vapply(suite$pipeline, function(pipeline) {
    sb_lint_pipeline(pipeline)$status
  }, character(1))
  expect_identical(unique(lint_status), "admitted")
})

test_that("MIME117 generator preserves variants and stage roles", {
  suite <- sb_mime117_pipelines(
    inner_folds = 3L,
    inner_seed = 8711L,
    selector_inner_folds = 3L,
    selector_seed = 8712L,
    rsf_nodesize = 5L
  )

  enet <- suite[grepl("Enet", suite$official_model_label, fixed = TRUE), ]
  stepcox <- suite[grepl("StepCox", suite$official_model_label, fixed = TRUE), ]
  expect_setequal(
    unique(enet$enet_alpha),
    seq(0.1, 0.9, by = 0.1)
  )
  expect_setequal(
    unique(stats::na.omit(stepcox$stepcox_direction)),
    c("both", "backward", "forward")
  )
  expect_identical(suite$official_model_label[[1L]], "RSF")
  expect_identical(suite$pipeline[[1L]]$selector, "none")
  expect_identical(suite$pipeline[[1L]]$learner, "rsf_survival")

  rsf_coxboost <- match("RSF + CoxBoost", suite$official_model_label)
  expect_identical(
    suite$pipeline[[rsf_coxboost]]$selector,
    "rsf_minimal_depth"
  )
  expect_identical(suite$pipeline[[rsf_coxboost]]$learner, "coxboost_cv")
  coxboost <- suite$pipeline[[match("CoxBoost", suite$official_model_label)]]
  plsrcox <- suite$pipeline[[match("plsRcox", suite$official_model_label)]]
  superpc <- suite$pipeline[[match("SuperPC", suite$official_model_label)]]
  gbm <- suite$pipeline[[match("GBM", suite$official_model_label)]]
  expect_identical(coxboost$learner, "coxboost_cv")
  expect_identical(coxboost$learner_params$start_penalty, 500)
  expect_identical(coxboost$learner_params$maximum_steps, 500L)
  expect_identical(plsrcox$learner, "plsrcox_cv")
  expect_identical(plsrcox$learner_params$maximum_components, 10L)
  expect_identical(plsrcox$learner_params$inner_folds, 5L)
  expect_identical(superpc$learner, "superpc_cv")
  expect_identical(superpc$learner_params$threshold_count, 20L)
  expect_identical(superpc$learner_params$tuning_components, 3L)
  expect_identical(superpc$learner_params$prediction_components, 1L)
  expect_identical(gbm$learner, "gbm_cox_cv")
  expect_identical(gbm$learner_params$maximum_trees, 10000L)
  expect_identical(gbm$learner_params$interaction_depth, 3L)
  expect_identical(gbm$learner_params$shrinkage, 0.001)

  lasso_stepcox <- match(
    "Lasso + StepCox[forward]",
    suite$official_model_label
  )
  expect_identical(
    suite$pipeline[[lasso_stepcox]]$selector,
    "lasso_then_stepcox"
  )
  expect_identical(suite$pipeline[[lasso_stepcox]]$learner, "coxph")
  expect_identical(
    suite$pipeline[[lasso_stepcox]]$selector_params$direction,
    "forward"
  )
  expect_identical(
    unique(suite$operational_fidelity),
    "fold_local_mime_training_reproduction_not_bug_for_bug"
  )
})
