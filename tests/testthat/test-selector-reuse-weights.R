weight_reuse_pipeline <- function(id, learner = "coxph") {
  if (id %in% c("lasso_nonzero", "enet_cox_nonzero")) {
    sb_pipeline(id, learner, selector_inner_folds = 10L, selector_seed = 1L)
  } else sb_pipeline(id, learner)
}

test_that("weight and Elastic Net selectors reuse only the same analysis inputs", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    list(selected_features = names(x), params = selector$params)
  })
  x <- data.frame(a = 1:6, b = 6:1, row.names = paste0("p", 1:6))
  endpoint <- list(task_type = "survival", time = 1:6, event = rep(0:1, 3),
    group_ids = c(1, 1, 2, 2, 3, 3))
  ids <- c(sb_survival_weight_selector_ids(), "enet_cox_nonzero")
  out <- sb_with_selector_reuse({
    for (id in ids) {
      p <- weight_reuse_pipeline(id)
      s <- sb_resolve_selector_adapter(id, p$selector_params)
      first <- sb_fit_runtime_selector(s, x, endpoint)
      expect_identical(sb_fit_runtime_selector(s, x, endpoint), first)
      changed <- endpoint; changed$group_ids <- 1:6
      sb_fit_runtime_selector(s, x, changed)
      changed <- endpoint; changed$event[1] <- 1L
      sb_fit_runtime_selector(s, x, changed)
      sb_fit_runtime_selector(s, x[c("b", "a")], endpoint)
      changed_s <- s; changed_s$params$probe_parameter <- 2L
      sb_fit_runtime_selector(changed_s, x, endpoint)
    }
  }, diagnostics = TRUE)
  expect_identical(calls, 5L * length(ids))
  expect_identical(out$reuse$hits, length(ids))
  expect_null(.sb_selector_reuse_context$active)
})

test_that("weight and Elastic Net warning or failed fits are not reused", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    if (selector$params$probe == "warning") warning("native warning")
    if (selector$params$probe == "error") sb_abort("native error", "stormbreaker_probe_failure")
    list(selected_features = names(x))
  })
  sb_with_selector_reuse({
    for (id in c(sb_survival_weight_selector_ids(), "enet_cox_nonzero")) {
      p <- weight_reuse_pipeline(id)
      s <- sb_resolve_selector_adapter(id, p$selector_params)
      s$params$probe <- "warning"
      for (i in 1:2) expect_warning(sb_fit_runtime_selector(s, data.frame(x = 1:4), list()), "native warning")
      s$params$probe <- "error"
      for (i in 1:2) expect_error(sb_fit_runtime_selector(s, data.frame(x = 1:4), list()), class = "stormbreaker_probe_failure")
    }
  })
  expect_identical(calls, 16L)
})

test_that("sharing keeps captured coefficient paths and loadings distinct", {
  for (pkg in c("glmnet", "plsRcox")) skip_if_not_installed(pkg)
  d <- sb_example_survival(150, 82)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "survival", time = d$time, event = d$event)
  for (id in c("ridge_cox_weights", "plsrcox_weights", "lasso_nonzero", "enet_cox_nonzero")) {
    p <- weight_reuse_pipeline(id)
    ordinary <- sb_resolve_selector_adapter(id, p$selector_params)
    captured <- sb_captured_selector_adapter(id, p$selector_params)
    run <- function(enabled) sb_with_selector_reuse({
      a <- sb_fit_runtime_selector(ordinary, x, endpoint)
      b <- sb_fit_runtime_selector(captured, x, endpoint)
      c <- sb_fit_runtime_selector(captured, x, endpoint)
      list(a, b, c)
    }, enabled = enabled, diagnostics = TRUE)
    plain <- run(FALSE); shared <- run(TRUE)
    expect_identical(shared$result, plain$result)
    expect_identical(shared$reuse$hits, 1L)
    history <- if (id == "plsrcox_weights") shared$result[[2]]$evidence$loadings else
      if (id == "ridge_cox_weights") shared$result[[2]]$evidence$coefficient_path else shared$result[[2]]$coefficient_path
    expect_false(is.null(history))
  }
})

test_that("CoxBoost reuse also distinguishes requested diagnostic history", {
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    list(history = environment(selector$fit)$capture_history)
  })
  out <- sb_with_selector_reuse({
    x <- data.frame(a = 1:6, b = 6:1)
    ordinary <- sb_coxboost_nonzero_selector_adapter(3L, 1L)
    captured <- sb_coxboost_nonzero_selector_adapter(3L, 1L, TRUE)
    list(sb_fit_runtime_selector(ordinary, x, list()),
      sb_fit_runtime_selector(captured, x, list()),
      sb_fit_runtime_selector(captured, x, list()))
  }, diagnostics = TRUE)
  expect_identical(out$result, list(list(history = FALSE), list(history = TRUE), list(history = TRUE)))
  expect_identical(out$reuse$hits, 1L)
})

test_that("native weight reuse preserves RNG and outcome-blind final predictions", {
  for (pkg in c("glmnet", "plsRcox", "survivalsvm")) skip_if_not_installed(pkg)
  d <- sb_example_survival(150, 82)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
  pipelines <- unlist(lapply(c(sb_survival_weight_selector_ids(), "enet_cox_nonzero"), function(id)
    lapply(c("coxph", "survreg_weibull"), function(learner) weight_reuse_pipeline(id, learner))), recursive = FALSE)
  run <- function(enabled) sb_with_selector_reuse(lapply(pipelines, function(p) {
    f <- sb_fit_final(task, p, capture_diagnostics = TRUE)
    expect_identical(f$status, "success")
    pred <- sb_predict_final(unserialize(serialize(f, NULL)), d[c("patient_id", task$predictors)], "patient_id")
    expect_identical(pred$status, "success")
    list(selection = f$selection, prediction = pred)
  }), enabled = enabled, diagnostics = TRUE)
  withr::local_seed(847L)
  before <- .Random.seed
  plain <- run(FALSE); expect_identical(.Random.seed, before)
  shared <- run(TRUE); expect_identical(.Random.seed, before)
  expect_identical(shared$result, plain$result)
  expect_identical(shared$reuse$hits, 4L)
})

test_that("weight reuse does not substitute outer selection for glmnet inner selection", {
  for (pkg in c("glmnet", "plsRcox")) skip_if_not_installed(pkg)
  d <- sb_example_survival(150, 82)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
  calls <- list()
  original <- sb_fit_runtime_selector_uncached
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls[[length(calls) + 1L]] <<- rownames(x)
    original(selector, x, endpoint)
  })
  pipelines <- lapply(c(0, .5), function(alpha) sb_pipeline("plsrcox_weights", "glmnet_cox",
    glmnet_alpha = alpha, inner_folds = 3L, inner_seed = 7L))
  plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 82L, reuse_selectors = FALSE)
  plain_calls <- calls; calls <- list()
  shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 82L)
  expect_identical(shared, plain)
  expect_length(plain_calls, 24L)
  expect_length(calls, 12L)
  expected <- c(lapply(shared$results[[1]]$folds, function(f) rownames(task$data)[f$analysis_rows]),
    unlist(lapply(shared$results[[1]]$tuning, function(t) lapply(t$inner_selections, function(s) s$analysis_patient_ids)), recursive = FALSE))
  expect_true(all(lengths(calls) < nrow(task$data)))
  expect_true(all(vapply(expected, function(ids) sum(vapply(calls, identical, logical(1), y = ids)) == 1L, logical(1))))
})
