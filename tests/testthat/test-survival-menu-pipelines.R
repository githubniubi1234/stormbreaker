test_that("source-backed survival menu preserves registry pairs and MIME settings", {
  features <- paste0("gene_", seq_len(60L))
  menu <- sb_survival_menu_pipelines(
    selector_features = features,
    feature_set_id = "candidate60_test"
  )
  registry <- sb_pipeline_registry(task_type = "survival")
  expect_setequal(
    paste(menu$selector, menu$learner),
    paste(registry$selector, registry$learner)
  )
  expect_identical(anyDuplicated(menu$config_id), 0L)
  expect_identical(anyDuplicated(vapply(menu$pipeline, `[[`, "", "id")), 0L)
  reference <- sb_mime117_pipelines(
    inner_seed = 1L, selector_seed = 1L, learner_seed = 1L
  )
  mapping <- attr(menu, "mime117_mapping")
  expect_identical(nrow(mapping), 117L)
  expect_identical(mapping$official_model_label, reference$official_model_label)
  indices <- match(mapping$comparison_config_id, menu$config_id)
  expect_identical(anyNA(indices), FALSE)
  expect_identical(unname(menu$pipeline[indices]), unname(reference$pipeline))
  expect_identical(sum(menu$comparison_role == "MIME_structure_matched"), 117L)
  expect_identical(attr(menu, "registry_pairs"), nrow(registry))
  for (p in menu$pipeline) {
    if (p$selector == "fixed_features") {
      expect_identical(p$selector_params$features, features)
      expect_identical(p$selector_params$feature_set_id, "candidate60_test")
    }
    if (p$selector %in% c("lasso_nonzero", "coxboost_nonzero",
      "lasso_then_stepcox", "coxboost_then_stepcox", "enet_cox_nonzero")) {
      expect_identical(p$selector_params$inner_folds, 10L)
      expect_identical(p$selector_params$seed, 1L)
    }
    if (p$learner == "glmnet_cox") {
      expect_identical(p$learner_params$inner_folds, 10L)
      expect_identical(p$learner_params$inner_seed, 1L)
      expect_identical(p$learner_params$selection_rule, "lambda.min")
      expect_identical(p$learner_params$standardize, TRUE)
    }
    if (p$learner %in% c("coxboost_cv", "gbm_cox_cv", "superpc_cv")) {
      expect_identical(p$learner_params$inner_folds, 10L)
      expect_identical(p$learner_params$seed, 1L)
    }
    if (p$learner == "plsrcox_cv") {
      expect_identical(p$learner_params$inner_folds, 5L)
      expect_identical(p$learner_params$maximum_components, 10L)
    }
    if (p$learner == "gbm_cox_cv") {
      expect_identical(p$learner_params$maximum_trees, 10000L)
    }
    if (p$selector == "boruta_survival_roughfix") {
      expect_null(p$selector_params)
    }
  }
})

test_that("survival menu validates caller-declared features without data access", {
  invalid <- tryCatch(sb_survival_menu_pipelines(
    selector_features = c("g", "g"), feature_set_id = "duplicate"
  ), error = identity)
  expect_s3_class(invalid, "stormbreaker_bad_fixed_features")
  invalid <- tryCatch(sb_survival_menu_pipelines(
    selector_features = c("a", "b"), feature_set_id = "valid",
    selector_n_features = 0L
  ), error = identity)
  expect_s3_class(invalid, "error")
})
