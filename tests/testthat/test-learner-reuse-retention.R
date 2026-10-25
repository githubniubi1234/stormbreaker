test_that("final fits reuse a nonconsecutive identical analysis input", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls <<- calls + 1L
    list(features = names(x), means = colMeans(x))
  })
  learner <- list(id = "logistic", params = list(), fit = function(...) NULL)
  x <- data.frame(a = 1:6, b = 6:1)
  endpoint <- list(task_type = "binary", outcome = rep(0:1, 3L))
  ids <- paste0("p", 1:6)
  out <- sb_with_selector_reuse({
    first <- sb_fit_final_runtime_learner(learner, x, endpoint, ids, "none")
    different <- sb_fit_final_runtime_learner(learner, x["a"], endpoint, ids, "fixed_features")
    again <- sb_fit_final_runtime_learner(learner, x, endpoint, ids, "none")
    expect_identical(again, first)
    expect_identical(different$features, "a")
  }, diagnostics = TRUE)
  expect_identical(calls, 2L)
  expect_identical(out$reuse$final_learner_hits, 1L)
  expect_identical(out$reuse$final_learner_entries, 2L)
})

test_that("retained learner payloads obey a shared CV and final memory budget", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls <<- calls + 1L
    list(value = sum(x), payload = raw(20000L))
  })
  learner <- list(id = "logistic", params = list(), fit = function(...) NULL)
  endpoint <- list(task_type = "binary", outcome = rep(0:1, 3L))
  ids <- paste0("p", 1:6)
  out <- sb_with_selector_reuse({
    context <- .sb_selector_reuse_context$active
    context$learner_byte_limit <- 70000
    fit <- function(v, phase = "final") sb_fit_reusable_runtime_learner(
      learner, data.frame(a = seq_len(6L) + v), endpoint, ids, "none", phase)
    a <- fit(0); fit(10, "cv"); expect_identical(fit(0), a)
    fit(20) # A was recently touched: evict the older CV entry, not A.
    expect_identical(fit(0), a)
    expect_identical(calls, 3L)
    fit(10, "cv")
    expect_identical(calls, 4L)
    expect_lte(context$learner_bytes, context$learner_byte_limit)
    expect_gt(context$learner_evictions, 0L)
  }, diagnostics = TRUE)
  expect_lte(out$reuse$learner_bytes, 70000)
  expect_null(.sb_selector_reuse_context$active)
})

test_that("oversized learner states are returned but not retained", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(...) {
    calls <<- calls + 1L
    list(payload = raw(20000L))
  })
  learner <- list(id = "logistic", params = list(), fit = function(...) NULL)
  out <- sb_with_selector_reuse({
    .sb_selector_reuse_context$active$learner_byte_limit <- 1000
    fit <- function() sb_fit_final_runtime_learner(learner, data.frame(a = 1:6),
      list(outcome = rep(0:1, 3L)), paste0("p", 1:6), "none")
    expect_identical(fit(), fit())
  }, diagnostics = TRUE)
  expect_identical(calls, 2L)
  expect_identical(out$reuse$final_learner_entries, 0L)
  expect_identical(out$reuse$learner_oversized, 2L)
})

test_that("eviction preserves native final model predictions", {
  d <- sb_example_binary(100L, 1208L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  a <- sb_pipeline("none", "logistic")
  b <- sb_pipeline("fixed_features", "logistic",
    selector_features = features[1:2], feature_set_id = "two")
  plain <- lapply(list(a, b, a), function(p) sb_fit_final(task, p))
  shared <- sb_with_selector_reuse({
    first <- sb_fit_final(task, a)
    context <- .sb_selector_reuse_context$active
    context$learner_byte_limit <- context$learner_bytes
    list(first, sb_fit_final(task, b), sb_fit_final(task, a))
  }, diagnostics = TRUE)
  expect_gt(shared$reuse$learner_evictions, 0L)
  expect_lte(shared$reuse$learner_bytes, shared$reuse$learner_byte_limit)
  expect_equal(shared$result, plain)
  for (i in 1:3) expect_identical(
    sb_predict_final(unserialize(serialize(shared$result[[i]], NULL)),
      d[c("patient_id", features)], "patient_id"),
    sb_predict_final(plain[[i]], d[c("patient_id", features)], "patient_id"))
})
