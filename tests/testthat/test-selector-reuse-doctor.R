test_that("doctor comparison reuses selectors without additional user settings", {
  d <- sb_example_survival(120, 37)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month")
  pipelines <- list(
    cox = surv_pipeline("Cox", surv_selector("StepCox")),
    chain = surv_pipeline("Cox", surv_selector(c("StepCox", "StepCox"))))
  original <- sb_fit_runtime_selector_uncached
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    if (identical(selector$id, "stepcox")) calls <<- calls + 1L
    original(selector, x, endpoint)
  })
  plain <- sb_compare_cv(task, unname(pipelines), folds = 3, seed = 19,
    reuse_selectors = FALSE, progress = FALSE)
  plain_calls <- calls
  calls <- 0L
  run <- surv_run(task, pipelines = pipelines, folds = 3, seed = 19,
    plots = FALSE)
  expect_identical(run$comparison, plain)
  expect_identical(run$comparison$comparison$status, rep("success", 2))
  expect_lt(calls, plain_calls)
  expect_length(run$details, 0L)
})
