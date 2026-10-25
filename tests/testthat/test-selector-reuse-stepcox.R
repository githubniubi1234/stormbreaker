test_that("real StepCox fits share only identical analysis and direction", {
  withr::local_seed(241L)
  x <- data.frame(a = rnorm(120), b = rnorm(120), c = rnorm(120))
  rownames(x) <- paste0("patient_", seq_len(nrow(x)))
  endpoint <- list(task_type = "survival",
    time = rexp(120, exp(1.5 * x$a - x$b)), event = rep(1L, 120))
  before <- .Random.seed
  out <- sb_with_selector_reuse({
    for (direction in c("forward", "backward", "both")) {
      selector <- sb_stepcox_selector_adapter(direction)
      first <- sb_fit_runtime_selector(selector, x, endpoint)
      expect_identical(sb_fit_runtime_selector(selector, x, endpoint), first)
      expect_identical(first, sb_fit_runtime_selector_uncached(selector, x, endpoint))
    }
    sb_fit_runtime_selector(selector, x[120:1, ],
      list(task_type = "survival", time = rev(endpoint$time), event = endpoint$event))
  }, diagnostics = TRUE)
  expect_identical(out$reuse$hits, 3L)
  expect_identical(out$reuse$misses, 4L)
  expect_identical(.Random.seed, before)
  expect_null(.sb_selector_reuse_context$active)
})

test_that("StepCox warning and error executions are not reused", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls <<- calls + 1L
    if (selector$params$direction == "forward") warning("native warning")
    else stop("native failure")
    list(selected_features = names(x))
  })
  out <- sb_with_selector_reuse({
    for (i in 1:2) {
      expect_warning(sb_fit_runtime_selector(sb_stepcox_selector_adapter("forward"),
        data.frame(a = 1:4), list(task_type = "survival")), "native warning")
      expect_error(sb_fit_runtime_selector(sb_stepcox_selector_adapter("both"),
        data.frame(a = 1:4), list(task_type = "survival")), "native failure")
    }
  }, diagnostics = TRUE)
  expect_identical(calls, 4L)
  expect_identical(out$reuse$misses, 4L)
  expect_identical(out$reuse$entries, 0L)
})

test_that("StepCox reuse preserves complete nested Cox comparisons", {
  skip_if_not_installed("glmnet")
  withr::local_seed(241L)
  d <- data.frame(patient_id = paste0("p", 1:120), a = rnorm(120), b = rnorm(120))
  d$time <- rexp(120, exp(1.5 * d$a - d$b))
  censor <- rexp(120, 0.12)
  d$event <- as.integer(d$time <= censor)
  d$time <- pmin(d$time, censor)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("a", "b"), "day")
  pipelines <- lapply(c(0, 1), function(alpha) sb_pipeline("stepcox", "glmnet_cox",
    stepcox_direction = "forward", glmnet_alpha = alpha, inner_folds = 3L,
    inner_seed = 17L))
  original <- sb_fit_runtime_selector_uncached
  calls <- list()
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls[[length(calls) + 1L]] <<- rownames(x)
    original(selector, x, endpoint)
  })
  plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 81L, reuse_selectors = FALSE)
  plain_calls <- length(calls)
  calls <- list()
  shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 81L)
  expect_true(all(shared$comparison$status == "success"))
  expect_identical(shared, plain)
  expect_lt(length(calls), plain_calls)
  expect_true(all(lengths(calls) < nrow(d)))
  expected <- c(lapply(shared$results[[1]]$folds,
    function(f) rownames(task$data)[f$analysis_rows]),
    unlist(lapply(shared$results[[1]]$tuning, function(t)
      lapply(t$inner_selections, function(s) s$analysis_patient_ids)), recursive = FALSE))
  expect_true(all(vapply(expected, function(ids)
    sum(vapply(calls, function(rows) identical(rows, ids), logical(1))) == 1L,
    logical(1))))
})
