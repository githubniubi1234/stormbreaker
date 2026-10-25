test_that("glmnet active sets are summarized from outer-fold coefficients", {
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 1,
    inner_folds = 3L,
    inner_seed = 1L
  )
  coefficients <- list(
    list(
      fold = 2L,
      status = "success",
      nonzero_features = c("gene_b", "gene_c")
    ),
    list(
      fold = 1L,
      status = "success",
      nonzero_features = c("gene_a", "gene_b")
    ),
    list(
      fold = 3L,
      status = "success",
      nonzero_features = "gene_b"
    )
  )

  observed <- sb_glmnet_active_set_stability(
    pipeline,
    coefficients,
    c("gene_a", "gene_b", "gene_c")
  )

  expect_identical(observed$learner, "glmnet_logistic")
  expect_identical(observed$alpha, 1)
  expect_identical(observed$status, "success")
  expect_identical(observed$per_fold$fold, 1:3)
  expect_identical(observed$per_fold$n_selected, c(2L, 2L, 1L))
  expect_equal(observed$mean_pairwise_jaccard, 4 / 9, tolerance = 0)
  expect_identical(
    observed$feature_frequency$selected_count,
    c(1L, 3L, 1L)
  )
  expect_equal(
    observed$feature_frequency$selection_frequency,
    c(1 / 3, 1, 1 / 3),
    tolerance = 0
  )
})

test_that("ridge and non-glmnet fits have no sparse active-set stability", {
  coefficients <- list(
    list(
      fold = 1L,
      status = "success",
      nonzero_features = c("gene_a", "gene_b")
    ),
    list(
      fold = 2L,
      status = "success",
      nonzero_features = c("gene_a", "gene_b")
    )
  )
  ridge <- sb_pipeline(
    learner = "glmnet_cox",
    glmnet_alpha = 0,
    inner_folds = 3L,
    inner_seed = 1L
  )
  ordinary <- sb_pipeline()

  ridge_audit <- sb_glmnet_active_set_stability(
    ridge,
    coefficients,
    c("gene_a", "gene_b")
  )
  ordinary_audit <- sb_glmnet_active_set_stability(
    ordinary,
    coefficients,
    c("gene_a", "gene_b")
  )

  expect_identical(ridge_audit$status, "not_applicable")
  expect_identical(ridge_audit$alpha, 0)
  expect_identical(ridge_audit$n_successful_folds, 0L)
  expect_true(all(is.na(
    ridge_audit$feature_frequency$selection_frequency
  )))
  expect_identical(ordinary_audit$status, "not_applicable")
  expect_true(is.na(ordinary_audit$alpha))
})

test_that("empty successful glmnet active sets remain visibly unevaluable", {
  pipeline <- sb_pipeline(
    selector = "fixed_features",
    learner = "glmnet_logistic",
    glmnet_alpha = 1,
    inner_folds = 3L,
    inner_seed = 1L,
    selector_features = c("gene_b", "gene_a"),
    feature_set_id = "pair"
  )
  coefficients <- lapply(1:3, function(fold) {
    list(
      fold = fold,
      status = "success",
      nonzero_features = character()
    )
  })

  observed <- sb_glmnet_active_set_stability(
    pipeline,
    coefficients,
    c("gene_a", "gene_b", "gene_c")
  )

  expect_identical(observed$status, "not_evaluable")
  expect_identical(
    observed$feature_frequency$feature,
    c("gene_b", "gene_a")
  )
  expect_identical(observed$per_fold$n_selected, rep(0L, 3L))
  expect_true(is.na(observed$mean_pairwise_jaccard))
})

test_that("binary and survival CV expose glmnet active-set audits", {
  skip_if_not_installed("glmnet")

  binary_data <- make_sb_binary_fixture(n = 60L, seed = 2404L)
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  binary <- sb_run_cv(
    binary_task,
    sb_pipeline(
      learner = "glmnet_logistic",
      glmnet_alpha = 1,
      inner_folds = 3L,
      inner_seed = 2405L
    ),
    folds = 3L,
    seed = 2406L
  )

  survival_data <- make_sb_fixture(n = 60L, seed = 2407L)
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  survival <- sb_run_cv(
    survival_task,
    sb_pipeline(
      learner = "glmnet_cox",
      glmnet_alpha = 1,
      inner_folds = 3L,
      inner_seed = 2408L
    ),
    folds = 3L,
    seed = 2409L
  )

  for (result in list(binary, survival)) {
    audit <- result$glmnet_active_set_stability
    expect_identical(audit$learner, result$pipeline$learner)
    expect_identical(audit$alpha, 1)
    expect_identical(
      audit$per_fold$n_selected,
      vapply(result$coefficients, function(x) x$n_nonzero, integer(1))
    )
    expect_identical(
      audit$n_successful_folds,
      as.integer(sum(vapply(
        result$coefficients,
        function(x) identical(x$status, "success"),
        logical(1)
      )))
    )
  }
})
