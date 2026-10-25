test_that("fixed feature pipelines retain an explicit feature block", {
  binary <- sb_pipeline(
    selector = "fixed_features",
    learner = "logistic",
    selector_features = c("age", "gene_a"),
    feature_set_id = "clinical_omics"
  )
  survival <- sb_pipeline(
    selector = "fixed_features",
    learner = "coxph",
    selector_features = c("gene_b", "gene_a"),
    feature_set_id = "two_gene"
  )

  expect_identical(binary$id, "fixed_features__logistic__setclinical_omics")
  expect_identical(binary$task_type, "binary")
  expect_identical(binary$selector_params, list(
    features = c("age", "gene_a"),
    feature_set_id = "clinical_omics",
    n_features = 2L
  ))
  expect_identical(survival$task_type, "survival")
  expect_identical(
    sb_lint_pipeline(binary, task_type = "binary")$status,
    "admitted"
  )

  registry <- sb_pipeline_registry()
  fixed <- registry[registry$selector == "fixed_features", , drop = FALSE]
  expect_identical(nrow(fixed), length(unique(registry$learner)))
  for (task_type in c("binary", "survival")) {
    expect_identical(
      sum(fixed$task_type == task_type),
      length(unique(registry$learner[registry$task_type == task_type]))
    )
  }
})

test_that("fixed feature arguments preserve the existing positional API", {
  positional <- sb_pipeline("none", "glmnet_cox", NULL, 3L, 101L)
  named <- sb_pipeline(
    selector = "none",
    learner = "glmnet_cox",
    selector_n_features = NULL,
    inner_folds = 3L,
    inner_seed = 101L
  )

  expect_identical(positional, named)
  expect_identical(
    positional$id,
    "none__glmnet_cox__alpha1__inner3__seed101"
  )

  ranger_positional <- do.call(
    sb_pipeline,
    c(
      list("none", "ranger_survival"),
      rep(list(NULL), 15L),
      list(100L, 2L, 3L, 7L)
    )
  )
  ranger_named <- sb_pipeline(
    learner = "ranger_survival",
    ranger_num_trees = 100L,
    ranger_mtry = 2L,
    ranger_min_node_size = 3L,
    ranger_seed = 7L
  )

  expect_identical(ranger_positional, ranger_named)
  expect_identical(
    ranger_positional$id,
    "none__ranger_survival__trees100__mtry2__minnode3__seed7"
  )

  gbm_positional <- do.call(
    sb_pipeline,
    c(
      list("none", "gbm_cox"),
      rep(list(NULL), 25L),
      list(100L, 2L, 10L, 0.1, 1, 1L)
    )
  )
  gbm_named <- sb_pipeline(
    learner = "gbm_cox",
    gbm_n_trees = 100L,
    gbm_interaction_depth = 2L,
    gbm_min_obs_node = 10L,
    gbm_shrinkage = 0.1,
    gbm_bag_fraction = 1,
    gbm_seed = 1L
  )

  expect_identical(gbm_positional, gbm_named)
  expect_identical(
    gbm_positional$id,
    paste0(
      "none__gbm_cox__trees100__depth2__minnode10__",
      "shrink0.1__bag1__seed1"
    )
  )
})

test_that("fixed feature specifications reject incomplete or ambiguous input", {
  missing_features <- capture_sb_error(sb_pipeline(
    selector = "fixed_features",
    learner = "logistic",
    feature_set_id = "set_a"
  ))
  missing_id <- capture_sb_error(sb_pipeline(
    selector = "fixed_features",
    learner = "logistic",
    selector_features = c("gene_a", "gene_b")
  ))
  duplicate_features <- capture_sb_error(sb_pipeline(
    selector = "fixed_features",
    learner = "logistic",
    selector_features = c("gene_a", "gene_a"),
    feature_set_id = "set_a"
  ))
  blank_id <- capture_sb_error(sb_pipeline(
    selector = "fixed_features",
    learner = "logistic",
    selector_features = c("gene_a", "gene_b"),
    feature_set_id = " "
  ))
  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "logistic",
    selector_features = c("gene_a", "gene_b"),
    feature_set_id = "set_a"
  ))
  one_feature_glmnet <- capture_sb_error(sb_pipeline(
    selector = "fixed_features",
    learner = "glmnet_logistic",
    selector_features = "gene_a",
    feature_set_id = "one_gene",
    inner_folds = 3L,
    inner_seed = 1L
  ))

  expect_s3_class(
    missing_features,
    "stormbreaker_missing_selector_parameter"
  )
  expect_s3_class(missing_id, "stormbreaker_missing_selector_parameter")
  expect_s3_class(
    duplicate_features,
    "stormbreaker_bad_fixed_features"
  )
  expect_s3_class(blank_id, "stormbreaker_bad_feature_set_id")
  expect_s3_class(
    irrelevant,
    "stormbreaker_irrelevant_selector_parameter"
  )
  expect_s3_class(
    one_feature_glmnet,
    "stormbreaker_glmnet_too_few_selected_features"
  )
})

test_that("fixed feature selection is outcome-blind and order-preserving", {
  x <- data.frame(
    gene_a = c(1, 2, 3),
    gene_b = c(4, 5, 6),
    gene_c = c(7, 8, 9)
  )
  adapter <- sb_resolve_selector_adapter(
    "fixed_features",
    list(
      features = c("gene_c", "gene_a"),
      feature_set_id = "ordered",
      n_features = 2L
    )
  )
  first <- sb_fit_runtime_selector(
    adapter,
    x,
    list(task_type = "binary", outcome = c(0L, 1L, 0L))
  )
  changed <- sb_fit_runtime_selector(
    adapter,
    x,
    list(task_type = "binary", outcome = c(1L, 0L, 1L))
  )
  selected <- adapter$apply(first, x)

  expect_identical(adapter$needs_outcome, FALSE)
  expect_identical(changed, first)
  expect_identical(first$feature_set_id, "ordered")
  expect_identical(first$selected_features, c("gene_c", "gene_a"))
  expect_identical(names(selected), c("gene_c", "gene_a"))
  expect_identical(selected$gene_c, x$gene_c)
  expect_s3_class(
    capture_sb_error(adapter$fit(x[c("gene_a", "gene_b")])),
    "stormbreaker_missing_fixed_feature"
  )
  expect_s3_class(
    capture_sb_error(adapter$apply(first, x[c("gene_b", "gene_a", "gene_c")])),
    "stormbreaker_selector_schema_mismatch"
  )
})

test_that("fixed feature blocks compare on identical outer folds", {
  data <- make_sb_binary_fixture(n = 60L, seed = 9961L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(
      selector = "fixed_features",
      learner = "logistic",
      selector_features = c("gene_a", "gene_b"),
      feature_set_id = "omics_pair"
    ),
    sb_pipeline(
      selector = "fixed_features",
      learner = "logistic",
      selector_features = "gene_c",
      feature_set_id = "single_numeric"
    )
  )

  observed <- sb_compare_cv(task, pipelines, folds = 3L, seed = 9962L)
  assignments <- lapply(observed$results, function(result) {
    result$predictions[c("patient_id", "fold")]
  })

  expect_identical(observed$status, "success")
  expect_identical(assignments[[1L]], assignments[[2L]])
  expect_identical(
    observed$comparison$selection_stability_status,
    rep("not_applicable", 2L)
  )
  expect_identical(
    observed$comparison$n_successful_selection_folds,
    rep(0L, 2L)
  )
  expect_identical(
    lapply(observed$results[[1L]]$selections, function(selection) {
      selection$selected_features
    }),
    rep(list(c("gene_a", "gene_b")), 3L)
  )
  expect_identical(
    lapply(observed$results[[2L]]$selections, function(selection) {
      selection$selected_features
    }),
    rep(list("gene_c"), 3L)
  )
})

test_that("fixed subsets feed final binary and survival prediction", {
  binary_data <- make_sb_binary_fixture(n = 60L, seed = 9971L)
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  binary_pipeline <- sb_pipeline(
    selector = "fixed_features",
    learner = "logistic",
    selector_features = c("gene_c", "gene_a"),
    feature_set_id = "binary_pair"
  )
  binary_fit <- sb_fit_final(binary_task, binary_pipeline)
  binary_new <- binary_data[1:6, c("patient_id", "gene_c", "gene_a")]
  binary_prediction <- sb_predict_final(
    binary_fit,
    binary_new,
    "patient_id"
  )
  binary_direct_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    c("gene_c", "gene_a")
  )
  binary_direct_fit <- sb_fit_final(
    binary_direct_task,
    sb_pipeline(learner = "logistic")
  )
  binary_direct_prediction <- sb_predict_final(
    binary_direct_fit,
    binary_new,
    "patient_id"
  )

  survival_data <- make_sb_fixture(n = 60L, seed = 9972L)
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  survival_pipeline <- sb_pipeline(
    selector = "fixed_features",
    learner = "coxph",
    selector_features = c("gene_b", "gene_a"),
    feature_set_id = "survival_pair"
  )
  survival_fit <- sb_fit_final(survival_task, survival_pipeline)
  survival_new <- survival_data[1:6, c("patient_id", "gene_b", "gene_a")]
  survival_prediction <- sb_predict_final(
    survival_fit,
    survival_new,
    "patient_id"
  )
  survival_direct_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    c("gene_b", "gene_a"),
    "month"
  )
  survival_direct_fit <- sb_fit_final(
    survival_direct_task,
    sb_pipeline(learner = "coxph")
  )
  survival_direct_prediction <- sb_predict_final(
    survival_direct_fit,
    survival_new,
    "patient_id"
  )
  survival_cv <- sb_run_cv(
    survival_task,
    survival_pipeline,
    folds = 3L,
    seed = 9973L
  )

  expect_identical(binary_fit$status, "success")
  expect_identical(binary_fit$required_features, c("gene_c", "gene_a"))
  expect_identical(binary_prediction$status, "success")
  expect_equal(
    binary_prediction$predictions$probability,
    binary_direct_prediction$predictions$probability,
    tolerance = 0
  )
  expect_identical(
    binary_prediction$predictions$patient_id,
    binary_new$patient_id
  )
  expect_identical(survival_fit$status, "success")
  expect_identical(
    survival_fit$required_features,
    c("gene_b", "gene_a")
  )
  expect_identical(survival_prediction$status, "success")
  expect_equal(
    survival_prediction$predictions$risk,
    survival_direct_prediction$predictions$risk,
    tolerance = 0
  )
  expect_identical(
    survival_prediction$predictions$patient_id,
    survival_new$patient_id
  )
  expect_identical(survival_cv$status, "success")
  expect_identical(
    survival_cv$selection_stability$status,
    "not_applicable"
  )
  expect_identical(
    lapply(survival_cv$selections, function(selection) {
      selection$selected_features
    }),
    rep(list(c("gene_b", "gene_a")), 3L)
  )
})

test_that("fixed glmnet subsets do not masquerade as learned inner selection", {
  skip_if_not_installed("glmnet")

  binary_data <- make_sb_binary_fixture(n = 60L, seed = 9981L)
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  binary_pipeline <- sb_pipeline(
    selector = "fixed_features",
    learner = "glmnet_logistic",
    selector_features = c("gene_b", "gene_a"),
    feature_set_id = "binary_glmnet_pair",
    glmnet_alpha = 0,
    inner_folds = 3L,
    inner_seed = 9982L
  )
  binary_fit <- sb_fit_final(binary_task, binary_pipeline)

  survival_data <- make_sb_fixture(n = 60L, seed = 9983L)
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  survival_pipeline <- sb_pipeline(
    selector = "fixed_features",
    learner = "glmnet_cox",
    selector_features = c("gene_c", "gene_a"),
    feature_set_id = "survival_glmnet_pair",
    glmnet_alpha = 0,
    inner_folds = 3L,
    inner_seed = 9984L
  )
  survival_fit <- sb_fit_final(survival_task, survival_pipeline)

  expect_identical(binary_fit$status, "success")
  expect_identical(
    binary_fit$tuning$inner_selection_method,
    "predeclared_feature_subset"
  )
  expect_length(binary_fit$tuning$inner_selections, 0L)
  expect_identical(survival_fit$status, "success")
  expect_identical(
    survival_fit$tuning$inner_selection_method,
    "predeclared_feature_subset"
  )
  expect_length(survival_fit$tuning$inner_selections, 0L)
})
