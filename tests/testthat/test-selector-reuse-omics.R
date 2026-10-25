test_that("seeded omics and PAMR selectors share complete native selections", {
  for (package in c("SIS", "VSURF", "lpc", "stabs", "pamr")) {
    skip_if_not_installed(package)
    requireNamespace(package, quietly = TRUE)
  }
  fixture <- withr::with_seed(9201L, {
    n <- 96L; p <- 60L
    y <- rep(0:1, length.out = n)
    x <- as.data.frame(matrix(rnorm(n * p), n, p))
    x[[1]] <- x[[1]] + 2 * y
    x[[2]] <- x[[2]] - 1.5 * y
    t <- rexp(n, exp(x[[1]] - x[[2]]))
    c <- rexp(n, .3)
    list(x = x, y = y, time = pmin(t, c) + .01, event = as.integer(t <= c))
  })
  native <- sb_fit_runtime_selector_uncached
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    native(selector, x, endpoint)
  })
  ids <- c("sis_binary_nlogn", "vsurf_binary_prediction", "lpc_binary_nlogn",
    "lpc_cox_nlogn", "stabsel_binary_lasso", "stabsel_cox_lasso", sb_pamr_selector_ids())
  for (id in ids) {
    adapter <- if (id %in% sb_pamr_selector_ids()) {
      sb_pamr_selector_adapter(id, sb_pamr_selector_params(id))
    } else sb_high_dimensional_selector_adapter(id)
    endpoint <- if (grepl("cox|survival", id)) {
      list(task_type = "survival", time = fixture$time, event = fixture$event)
    } else list(task_type = "binary", outcome = fixture$y)
    withr::with_seed(37L, {
      before <- .Random.seed
      plain <- native(adapter, fixture$x, endpoint)
      calls <- 0L
      shared <- sb_with_selector_reuse({
        list(sb_fit_runtime_selector(adapter, fixture$x, endpoint),
          sb_fit_runtime_selector(adapter, fixture$x, endpoint))
      }, diagnostics = TRUE)
      expect_identical(shared$result, list(plain, plain), info = id)
      expect_identical(calls, 1L, info = id)
      expect_identical(shared$reuse$hits, 1L, info = id)
      expect_identical(.Random.seed, before, info = id)
    })
  }
})
