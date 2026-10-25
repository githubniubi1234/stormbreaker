test_that("RSF sampling callback does not retain the native training frame", {
  skip_if_not_installed("randomForestSRC")
  d <- sb_example_survival(80L, 11032L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_additional_survival_learner_adapter("rsf_survival")
  model <- adapter$fit(x, d$time, d$event)
  expect_false(exists("nativeOutput", environment(model$engine$forest$sampsize), inherits = FALSE))
  expect_identical(model$engine$forest$sampsize(80L), 80L * 0.632)
  prepared <- sb_additional_survival_formula_data(x, d$time, d$event,
    random_forest_src = TRUE)
  native <- withr::with_options(list(rf.cores = 1L, mc.cores = 1L),
    randomForestSRC::rfsrc(prepared$formula,
      data = prepared$data, ntree = 1000L, nodesize = 5L, splitrule = "logrank",
      importance = FALSE, proximity = FALSE, forest = TRUE, seed = 1L))
  reference <- model
  reference$engine <- native
  expect_identical(withr::with_seed(22L, adapter$predict(model, x)),
    withr::with_seed(22L, adapter$predict(reference, x)))
  expect_identical(withr::with_seed(23L,
    adapter$predict_survival(model, x, d$time, d$event, x, 7)),
    withr::with_seed(23L,
      adapter$predict_survival(reference, x, d$time, d$event, x, 7)))
})
