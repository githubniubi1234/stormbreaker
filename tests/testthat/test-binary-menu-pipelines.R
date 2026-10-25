test_that("binary menu covers the registry without fitting or changing defaults", {
  genes <- paste0("gene_", seq_len(55L))
  menu <- sb_binary_menu_pipelines(genes, "fixed55")
  registry <- sb_pipeline_registry("binary")
  expect_setequal(paste(menu$selector, menu$learner),
    paste(registry$selector, registry$learner))
  expect_identical(nrow(menu), nrow(registry))
  expect_identical(anyDuplicated(menu$config_id), 0L)
  expect_identical(anyDuplicated(vapply(menu$pipeline, `[[`, "", "id")), 0L)
  expect_identical(attr(menu, "registry_pairs"), nrow(registry))
  for (p in menu$pipeline) {
    expect_identical(p$task_type, "binary")
    if (p$selector == "fixed_features") {
      expect_identical(p$selector_params$features, genes)
      expect_identical(p$selector_params$feature_set_id, "fixed55")
    }
    if (p$selector == "univ_binary_auc_top_k") {
      expect_identical(p$selector_params$n_features, 10L)
    }
    if (p$learner == "glmnet_logistic") {
      expect_identical(p$learner_params$inner_folds, 5L)
      expect_identical(p$learner_params$inner_seed, 1L)
    }
  }
  reference <- attr(menu, "mime7_family_mapping")
  expect_setequal(reference$official_method,
    c("nb", "svmRadialWeights", "rf", "kknn", "adaboost", "LogitBoost", "cancerclass"))
  expect_identical(anyNA(match(reference$config_id, menu$config_id)), FALSE)
  expect_identical(unique(reference$relationship), "algorithm_family_only")
  expect_identical(sum(nzchar(menu$mime_family)), 7L)
})

test_that("binary menu controls are explicit and invalid inputs are rejected", {
  invalid <- tryCatch(sb_binary_menu_pipelines(c("a", "a"), "duplicate"), error = identity)
  expect_s3_class(invalid, "stormbreaker_bad_fixed_features")
  invalid <- tryCatch(sb_binary_menu_pipelines(c("a", "b"), "valid", inner_folds = 2L), error = identity)
  expect_s3_class(invalid, "error")
  row <- sb_pipeline_registry("binary")
  row <- row[row$selector == "none" & row$learner == "glmnet_logistic", ]
  args <- sb_binary_menu_arguments(row, c("a", "b"), "fixed", 2L, 3L, 42L)
  p <- do.call(sb_pipeline, args)
  expect_identical(p$learner_params$inner_folds, 3L)
  expect_identical(p$learner_params$inner_seed, 42L)
  direct <- sb_pipeline(learner = "glmnet_logistic", inner_folds = 3L, inner_seed = 42L)
  expect_identical(p, direct)
})
