test_that("survival forest predictions agree across assessment batch sizes", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("coin")
  withr::local_seed(123)
  x <- data.frame(a = rnorm(66), b = rnorm(66), c = rnorm(66))
  time <- rexp(66, exp(x$a) / 10)
  event <- rep(c(1L, 1L, 0L), 22)
  adapter <- sb_resolve_pipeline_runtime(sb_pipeline("none", "cforest_survival"))$learner
  model <- adapter$fit(x, time, event)
  large <- x[rep(seq_len(66), length.out = 518), , drop = FALSE]
  actual <- sb_partykit_survival_distributions(model, large)
  expected <- lapply(seq_len(nrow(large)), function(i) {
    predict(model$engine, newdata = large[i, , drop = FALSE], type = "prob")[[1L]]
  })
  expect_length(actual, 518L)
  expect_equal(lapply(actual, `[[`, "surv"), lapply(expected, `[[`, "surv"))
  expect_equal(lapply(actual, `[[`, "time"), lapply(expected, `[[`, "time"))
  small <- sb_partykit_survival_distributions(model, x[1:10, ])
  expect_equal(small, predict(model$engine, newdata = x[1:10, ], type = "prob"))
})
