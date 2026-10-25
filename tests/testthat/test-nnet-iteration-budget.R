test_that("nnet gains optimization time without changing saved budgets", {
  p <- sb_pipeline(learner = "nnet_binary")
  expect_identical(p$learner_params$maxit, 2000L)
  expect_identical(sb_validate_pipeline_contract(p, "binary"), p)
  old <- p
  old$learner_params$maxit <- 200L
  expect_identical(sb_validate_pipeline_contract(old, "binary"), old)
  old$learner_params$decay <- 10
  error <- tryCatch(sb_validate_pipeline_contract(old, "binary"), error = identity)
  expect_s3_class(error, "stormbreaker_unresolved_pipeline")
})
