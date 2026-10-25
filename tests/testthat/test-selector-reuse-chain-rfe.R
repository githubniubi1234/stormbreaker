test_that("SVM-RFE shares the same learned prefix while retaining inner refits", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(120L, 1208L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  pipelines <- lapply(c("logistic", "gaussian_naive_bayes"), function(learner) {
    sb_new_chain_pipeline(sb_pipeline("none", learner), list(
      sb_pipeline("univ_binary_auc_top_k", learner, selector_n_features = 2L),
      sb_pipeline("svm_rfe_binary", learner, selector_options = list(inner_folds = 3L))))
  })
  native <- sb_fit_svm_rfe
  calls <- 0L
  local_mocked_bindings(sb_fit_svm_rfe = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L, reuse_selectors = FALSE)
  expect_identical(calls, 6L)
  calls <- 0L
  shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L)
  expect_identical(calls, 3L)
  expect_identical(shared, plain)
  expect_identical(unique(shared$comparison$status), "success")

  # Changing an upstream input must not hit merely because selected columns match.
  base <- sb_pipeline("none", "logistic")
  chain <- sb_new_chain_pipeline(base, list(
    sb_pipeline("fixed_features", "logistic", selector_features = features[1:2],
      feature_set_id = "first_two"),
    sb_pipeline("svm_rfe_binary", "logistic", selector_options = list(inner_folds = 3L))))
  adapter <- sb_resolve_pipeline_runtime(chain)$selector
  endpoint <- list(task_type = "binary", outcome = d$outcome, group_ids = d$patient_id)
  x <- d[features]; changed <- x; changed$gene_c <- -changed$gene_c
  calls <- 0L
  sb_with_selector_reuse({
    a <- adapter$fit(x, endpoint)
    expect_identical(adapter$fit(x, endpoint), a)
    adapter$fit(changed, endpoint)
  })
  expect_identical(calls, 2L)
})
