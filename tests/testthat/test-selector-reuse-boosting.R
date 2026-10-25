test_that("survival boosting selectors share only identical successful fits", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    list(selected_features = names(x), params = selector$params)
  })
  x <- data.frame(a = 1:6, b = 6:1)
  endpoint <- list(task_type = "survival", time = 1:6, event = rep(0:1, 3))
  out <- sb_with_selector_reuse({
    for (id in sb_survival_boosting_selector_ids()) {
      p <- sb_pipeline(id, "coxph")
      s <- sb_resolve_selector_adapter(id, p$selector_params)
      first <- sb_fit_runtime_selector(s, x, endpoint)
      expect_identical(sb_fit_runtime_selector(s, x, endpoint), first)
      s$params$seed <- 9L
      sb_fit_runtime_selector(s, x, endpoint)
      sb_fit_runtime_selector(s, x[1:4, ], list(task_type = "survival", time = 1:4, event = rep(0:1, 2)))
    }
  }, diagnostics = TRUE)
  expect_identical(calls, 15L)
  expect_identical(out$reuse$hits, 5L)
  expect_null(.sb_selector_reuse_context$active)
})

test_that("boosting warnings and failures are still evaluated on every request", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    if (selector$params$mode == "warning") warning("original warning")
    if (selector$params$mode == "error") sb_abort("original error", "stormbreaker_probe_failure")
    list(selected_features = names(x))
  })
  s <- sb_resolve_selector_adapter("gbm_cox_importance", list())
  sb_with_selector_reuse({
    s$params <- list(mode = "warning")
    for (i in 1:2) expect_warning(sb_fit_runtime_selector(s, data.frame(x = 1:4), list()), "original warning")
    s$params <- list(mode = "error")
    for (i in 1:2) expect_error(sb_fit_runtime_selector(s, data.frame(x = 1:4), list()), class = "stormbreaker_probe_failure")
  })
  expect_identical(calls, 4L)
})

test_that("native boosting selection, final predictions and RNG survive reuse", {
  for (pkg in c("gbm", "mboost", "superpc", "xgboost")) skip_if_not_installed(pkg)
  d <- sb_example_survival(150, 82)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
  pipelines <- unlist(lapply(sb_survival_boosting_selector_ids(), function(id) {
    lapply(c("coxph", "survreg_weibull"), function(learner) sb_pipeline(id, learner))
  }), recursive = FALSE)
  run <- function(enabled) sb_with_selector_reuse(lapply(pipelines, function(p) {
    f <- sb_fit_final(task, p, capture_diagnostics = TRUE)
    expect_identical(f$status, "success")
    prediction <- sb_predict_final(unserialize(serialize(f, NULL)), d[c("patient_id", task$predictors)], "patient_id")
    list(selection = f$selection, prediction = prediction)
  }), enabled = enabled, diagnostics = TRUE)
  withr::local_seed(847L)
  before <- .Random.seed
  plain <- run(FALSE)
  expect_identical(.Random.seed, before)
  shared <- run(TRUE)
  expect_identical(.Random.seed, before)
  expect_identical(shared$result, plain$result)
  expect_identical(shared$reuse$hits, 5L)
})

test_that("boosting reuse retains separate outer and glmnet inner analysis folds", {
  skip_if_not_installed("mboost")
  skip_if_not_installed("glmnet")
  d <- sb_example_survival(150, 82)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
  calls <- list()
  original <- sb_fit_runtime_selector_uncached
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls[[length(calls) + 1L]] <<- rownames(x)
    original(selector, x, endpoint)
  })
  pipelines <- lapply(c(0, .5), function(alpha) sb_pipeline("gamboost_cox_nonzero", "glmnet_cox",
    glmnet_alpha = alpha, inner_folds = 3L, inner_seed = 7L))
  plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 82L, reuse_selectors = FALSE)
  plain_calls <- calls; calls <- list()
  shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 82L)
  expect_identical(shared, plain)
  expect_length(plain_calls, 24L)
  expect_length(calls, 12L)
  expected <- c(lapply(shared$results[[1L]]$folds, function(f) rownames(task$data)[f$analysis_rows]),
    unlist(lapply(shared$results[[1L]]$tuning, function(t) lapply(t$inner_selections, function(s) s$analysis_patient_ids)), recursive = FALSE))
  expect_true(all(lengths(calls) < nrow(task$data)))
  expect_true(all(vapply(expected, function(ids) sum(vapply(calls, identical, logical(1), y = ids)) == 1L, logical(1))))
})
