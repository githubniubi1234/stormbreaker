test_that("Boruta worker handoff preserves matching and analysis boundaries", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    list(selected_features = names(x), value = calls)
  })
  selector <- sb_resolve_selector_adapter("boruta_binary_roughfix", list())
  x <- data.frame(a = 1:8, b = 8:1)
  endpoint <- list(task_type = "binary", outcome = rep(0:1, 4))
  saved <- sb_with_selector_reuse({
    fit <- sb_fit_runtime_selector(selector, x, endpoint)
    list(fit = fit, checkpoint = sb_boruta_worker_checkpoint(selector$id))
  })
  saved <- unserialize(serialize(saved, NULL))
  result <- sb_with_selector_reuse({
    sb_restore_boruta_worker_checkpoint(saved$checkpoint)
    expect_identical(sb_fit_runtime_selector(selector, x, endpoint), saved$fit)
    sb_fit_runtime_selector(selector, x[-1, ],
      list(task_type = "binary", outcome = endpoint$outcome[-1]))
  }, diagnostics = TRUE)
  expect_identical(calls, 2L)
  expect_identical(result$reuse$hits, 1L)
  expect_identical(result$reuse$misses, 1L)
  expect_null(.sb_selector_reuse_context$active)
})

test_that("failed consumers do not alter saved Boruta entries", {
  selector <- sb_resolve_selector_adapter("boruta_binary_roughfix", list())
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(...) list(selected_features = "a"))
  x <- data.frame(a = 1:6)
  endpoint <- list(task_type = "binary", outcome = rep(0:1, 3))
  saved <- sb_with_selector_reuse({
    sb_fit_runtime_selector(selector, x, endpoint)
    sb_boruta_worker_checkpoint(selector$id)
  })
  expect_error(sb_with_selector_reuse({
    sb_restore_boruta_worker_checkpoint(saved)
    stop("downstream failure")
  }), "downstream failure")
  result <- sb_with_selector_reuse({
    sb_restore_boruta_worker_checkpoint(saved)
    sb_fit_runtime_selector(selector, x, endpoint)
  }, diagnostics = TRUE)
  expect_identical(result$reuse$hits, 1L)
  expect_null(.sb_selector_reuse_context$active)
})
