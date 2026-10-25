test_that("one pipeline construction prepares its registry once", {
  native <- sb_pipeline_registry_data
  count <- 0L
  expected <- sb_pipeline("none", "coxph")
  actual <- with_mocked_bindings(
    sb_pipeline("none", "coxph"),
    sb_pipeline_registry_data = function() {
      count <<- count + 1L
      native()
    }
  )
  expect_identical(actual, expected)
  expect_identical(count, 1L)
})

test_that("registry helpers agree for local and standalone lookup", {
  registry <- sb_pipeline_registry_data()
  for (learner in c("coxph", "logistic", "glmnet_cox", "glmnet_logistic")) {
    expect_identical(
      sb_pipeline_registry_match("none", learner, registry),
      sb_pipeline_registry_match("none", learner)
    )
    expect_identical(
      sb_pipeline_uses_inner_cv(learner, registry),
      sb_pipeline_uses_inner_cv(learner)
    )
  }
  expect_identical(
    sb_pipeline_registry_match("unregistered", "coxph", registry),
    sb_pipeline_registry_match("unregistered", "coxph")
  )
})
