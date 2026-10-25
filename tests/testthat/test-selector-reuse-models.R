test_that("model selectors reuse native fits without changing results or RNG", {
  for (package in c("glmnet", "mboost", "randomForest", "gbm", "e1071",
                    "xgboost", "plsRglm", "abess")) {
    skip_if_not_installed(package)
    requireNamespace(package, quietly = TRUE)
  }
  native <- sb_fit_runtime_selector_uncached
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    native(selector, x, endpoint)
  })
  d <- sb_example_binary(160L, 1208L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "binary", outcome = d$outcome,
    group_ids = d$patient_id)
  withr::local_seed(37L)
  before <- .Random.seed
  for (id in c(sb_binary_model_selector_ids(), sb_abess_selector_ids())) {
    if (id == "abess_cox_gic_nonzero") {
      s <- sb_example_survival(160L, 11032L)
      endpoint <- list(task_type = "survival", time = s$time, event = s$event,
        group_ids = s$patient_id)
      x <- s[c("gene_a", "gene_b", "gene_c")]
    }
    adapter <- if (id %in% sb_binary_model_selector_ids()) {
      sb_binary_model_selector_adapter(id, sb_binary_model_selector_options(id))
    } else sb_abess_selector_adapter(id)
    plain <- native(adapter, x, endpoint)
    calls <- 0L
    shared <- sb_with_selector_reuse({
      a <- sb_fit_runtime_selector(adapter, x, endpoint)
      b <- sb_fit_runtime_selector(adapter, x, endpoint)
      list(a, b)
    }, diagnostics = TRUE)
    expect_identical(shared$result[[1]], plain, info = id)
    expect_identical(shared$result[[2]], plain, info = id)
    expect_identical(calls, 1L, info = id)
    expect_identical(shared$reuse$hits, 1L, info = id)
    expect_identical(.Random.seed, before, info = id)
  }
})

test_that("shared selectors keep diagnostic history requests separate", {
  for (package in c("glmnet", "randomForest", "abess")) skip_if_not_installed(package)
  d <- sb_example_binary(160L, 1208L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "binary", outcome = d$outcome)
  for (id in c("lasso_binary_nonzero", "rf_binary_importance", "abess_binary_gic_nonzero")) {
    adapters <- lapply(c(FALSE, TRUE), function(history) {
      if (id == "abess_binary_gic_nonzero") sb_abess_selector_adapter(id, history)
      else sb_binary_model_selector_adapter(id, sb_binary_model_selector_options(id), history)
    })
    plain <- lapply(adapters, sb_fit_runtime_selector_uncached, x = x, endpoint = endpoint)
    shared <- sb_with_selector_reuse({
      lapply(c(adapters, adapters), sb_fit_runtime_selector, x = x, endpoint = endpoint)
    }, diagnostics = TRUE)
    expect_identical(shared$result, c(plain, plain), info = id)
    expect_identical(shared$reuse$misses, 2L, info = id)
    expect_identical(shared$reuse$hits, 2L, info = id)
  }
})

test_that("ABESS shares distinct outer and inner fits without sharing across folds", {
  skip_if_not_installed("abess")
  skip_if_not_installed("glmnet")
  native <- sb_fit_runtime_selector_uncached
  calls <- list()
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls[[length(calls) + 1L]] <<- rownames(x)
    native(selector, x, endpoint)
  })
  for (type in c("binary", "survival")) {
    features <- c("gene_a", "gene_b", "gene_c")
    if (type == "binary") {
      d <- sb_example_binary(160L, 1208L)
      task <- sb_binary_task(d, "patient_id", "outcome", features)
      selector <- "abess_binary_gic_nonzero"; learner <- "glmnet_logistic"
    } else {
      d <- sb_example_survival(160L, 11032L)
      task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
      selector <- "abess_cox_gic_nonzero"; learner <- "glmnet_cox"
    }
    pipelines <- lapply(c(1, .5), function(alpha) sb_pipeline(selector, learner,
      glmnet_alpha = alpha, inner_folds = 3L, inner_seed = 101L))
    calls <- list()
    plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L, reuse_selectors = FALSE)
    expect_length(calls, 24L)
    calls <- list()
    shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L)
    expect_identical(shared, plain)
    expect_identical(unique(shared$comparison$status), "success")
    expect_length(calls, 12L)
    expected <- c(lapply(shared$results[[1]]$folds, function(f) rownames(task$data)[f$analysis_rows]),
      unlist(lapply(shared$results[[1]]$tuning, function(t)
        lapply(t$inner_selections, function(s) s$analysis_patient_ids)), recursive = FALSE))
    expect_identical(vapply(expected, function(ids) sum(vapply(calls,
      function(rows) identical(rows, ids), logical(1))), integer(1)), rep(1L, 12L))
  }
})

test_that("binary model selectors preserve complete comparisons with sharing", {
  for (package in c("glmnet", "mboost", "randomForest", "gbm", "e1071",
                    "xgboost", "plsRglm")) skip_if_not_installed(package)
  d <- sb_example_binary(160L, 1208L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  native <- sb_fit_runtime_selector_uncached
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    native(selector, x, endpoint)
  })
  for (id in sb_binary_model_selector_ids()) {
    pipelines <- lapply(c("logistic", "gaussian_naive_bayes"),
      function(learner) sb_pipeline(id, learner))
    calls <- 0L
    plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L, reuse_selectors = FALSE)
    expect_identical(calls, 6L, info = id)
    calls <- 0L
    shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L)
    expect_identical(shared, plain, info = id)
    failed <- sum(vapply(shared$results[[1]]$selections,
      function(s) !identical(s$status, "success"), logical(1)))
    expect_identical(calls, 3L + failed, info = id)
    # Empty native selections are retried for each consumer, not cached or filled.
    if (failed) expect_identical(unique(shared$comparison$status), "failed", info = id)
  }
})
