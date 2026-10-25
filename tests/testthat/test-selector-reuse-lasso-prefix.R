test_that("Lasso selection and its StepCox prefix share actual fitted states", {
  skip_if_not_installed("glmnet")
  withr::local_seed(241L)
  x <- data.frame(a = rnorm(120), b = rnorm(120), noise = rnorm(120))
  rownames(x) <- paste0("patient_", seq_len(nrow(x)))
  time <- rexp(120, exp(1.5 * x$a - x$b))
  censor <- rexp(120, 0.12)
  e <- list(task_type = "survival", time = pmin(time, censor),
    event = as.integer(time <= censor))
  first <- sb_lasso_nonzero_selector_adapter(3L, 17L)
  chain <- sb_lasso_then_stepcox_selector_adapter(3L, 17L, "forward")
  before <- .Random.seed
  expected <- sb_with_selector_reuse({
    list(first = sb_fit_runtime_selector(first, x, e),
      chain = sb_fit_runtime_selector(chain, x, e))
  }, enabled = FALSE)
  actual <- sb_with_selector_reuse({
    list(first = sb_fit_runtime_selector(first, x, e),
      chain = sb_fit_runtime_selector(chain, x, e))
  }, diagnostics = TRUE)
  expect_identical(actual$result, expected)
  expect_identical(actual$result$first, actual$result$chain$stages$first)
  expect_identical(actual$reuse$hits, 1L)
  expect_identical(actual$reuse$entries, 2L)
  expect_identical(.Random.seed, before)

  settings <- sb_with_selector_reuse({
    sb_fit_runtime_selector(first, x, e)
    sb_fit_runtime_selector(sb_lasso_nonzero_selector_adapter(3L, 18L), x, e)
    sb_fit_runtime_selector(sb_lasso_nonzero_selector_adapter(4L, 17L), x, e)
  }, diagnostics = TRUE)
  expect_identical(settings$reuse$hits, 0L)
  expect_identical(settings$reuse$misses, 3L)
})

test_that("Lasso reuse preserves outer selection and complete predictions", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("survivalsvm")
  withr::local_seed(241L)
  d <- data.frame(patient_id = paste0("p", 1:120), a = rnorm(120), b = rnorm(120))
  time <- rexp(120, exp(1.5 * d$a - d$b))
  censor <- rexp(120, 0.12)
  d$time <- pmin(time, censor)
  d$event <- as.integer(time <= censor)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("a", "b"), "day")
  pipelines <- lapply(c(1, 2), function(gamma)
    sb_pipeline("lasso_nonzero", "survivalsvm",
      selector_inner_folds = 3L, selector_seed = 17L,
      survivalsvm_gamma_mu = gamma))
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
  expected <- lapply(shared$results[[1]]$folds,
    function(f) rownames(task$data)[f$analysis_rows])
  expect_true(all(vapply(expected, function(ids)
    sum(vapply(calls, function(rows) identical(rows, ids), logical(1))) == 1L,
    logical(1))))
  expect_true(all(lengths(calls) < nrow(d)))
})
