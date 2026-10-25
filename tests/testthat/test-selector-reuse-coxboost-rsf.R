test_that("native CoxBoost and RSF selections replay without RNG side effects", {
  skip_if_not_installed("CoxBoost")
  skip_if_not_installed("randomForestSRC")
  withr::local_seed(8601L)
  x <- as.data.frame(matrix(rnorm(90 * 6), nrow = 90))
  names(x) <- paste0("gene", 1:6)
  rownames(x) <- paste0("patient", 1:90)
  time <- rexp(90, exp(1.2 * x$gene1 - 0.9 * x$gene2 + 0.6 * x$gene3) / 15)
  censor <- rexp(90, 1 / 20)
  endpoint <- list(task_type = "survival", time = pmax(pmin(time, censor), 0.05),
    event = as.integer(time <= censor))
  for (selector in list(sb_coxboost_nonzero_selector_adapter(3L, 1L),
    sb_rsf_minimal_depth_selector_adapter(5L, 1L))) {
    before <- .Random.seed
    plain <- sb_with_selector_reuse(sb_fit_runtime_selector(selector, x, endpoint),
      enabled = FALSE)
    out <- sb_with_selector_reuse({
      first <- sb_fit_runtime_selector(selector, x, endpoint)
      second <- sb_fit_runtime_selector(selector, x, endpoint)
      expect_identical(first, second)
      first
    }, diagnostics = TRUE)
    expect_s3_class(plain, "stormbreaker_selector_fit")
    expect_identical(out$result, plain)
    expect_identical(out$reuse$hits, 1L)
    expect_identical(out$reuse$misses, 1L)
    expect_identical(.Random.seed, before)
  }
})

test_that("CoxBoost and RSF reuse distinguishes settings and never stores failures", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    if (identical(endpoint$mode, "warning")) warning("native warning")
    if (identical(endpoint$mode, "error")) stop("native error")
    list(selected_features = names(x))
  })
  x <- data.frame(a = 1:12, b = 12:1, row.names = paste0("p", 1:12))
  endpoint <- list(task_type = "survival", time = as.numeric(1:12), event = rep(0:1, 6))
  for (selector in list(sb_coxboost_nonzero_selector_adapter(3L, 1L),
    sb_rsf_minimal_depth_selector_adapter(5L, 1L))) {
    out <- sb_with_selector_reuse({
      sb_fit_runtime_selector(selector, x, endpoint)
      sb_fit_runtime_selector(selector, x, endpoint)
      changed <- selector
      changed$params$seed <- 2L
      sb_fit_runtime_selector(changed, x, endpoint)
      changed <- selector
      if (selector$id == "coxboost_nonzero") changed$params$inner_folds <- 4L
      else changed$params$nodesize <- 6L
      sb_fit_runtime_selector(changed, x, endpoint)
      changed_endpoint <- endpoint
      changed_endpoint$time[1] <- 0.5
      sb_fit_runtime_selector(selector, x, changed_endpoint)
      sb_fit_runtime_selector(selector, x[12:1, ], endpoint)
      for (i in 1:2) {
        expect_warning(sb_fit_runtime_selector(selector, x,
          c(endpoint, list(mode = "warning"))), "native warning")
        expect_error(sb_fit_runtime_selector(selector, x,
          c(endpoint, list(mode = "error"))), "native error")
      }
    }, diagnostics = TRUE)
    expect_identical(out$reuse$hits, 1L)
    expect_identical(out$reuse$entries, 5L)
    expect_identical(out$reuse$misses, 9L)
  }
  expect_identical(calls, 18L)
})
