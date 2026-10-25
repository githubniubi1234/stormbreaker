test_that("retired learners have no registry, constructor or adapter route", {
  retired <- sb_retired_learner_ids()
  registry <- sb_pipeline_registry()
  expect_identical(intersect(registry$learner, retired), character())
  expect_length(unique(registry$learner), 176L)
  for (id in retired) {
    expect_error(sb_pipeline(learner = id), class = "stormbreaker_retired_learner")
    expect_error(sb_pipeline("ordered_chain", id), class = "stormbreaker_retired_learner")
    expect_error(sb_resolve_learner_adapter(id), class = "stormbreaker_retired_learner")
    p <- sb_pipeline(learner = "logistic")
    p$learner <- id
    expect_error(sb_validate_pipeline_contract(p, "binary"),
      class = "stormbreaker_retired_learner")
  }
  expect_identical(sb_mgcv_gam_learner_ids(), character())
  expect_identical(sb_slope_binary_learner_ids(), character())
  expect_identical(sb_slope_randomglm_binary_learner_ids(), "random_glm_binary")
})

test_that("removed native fitting implementations cannot be invoked", {
  expect_error(sb_fit_mgcv_gam(NULL, NULL, "gam_cox_reml", NULL),
    class = "stormbreaker_retired_learner")
  expect_error(sb_fit_mgcv_gam(NULL, NULL, "gam_binary_reml", NULL),
    class = "stormbreaker_retired_learner")
  expect_error(sb_fit_slope_binary(NULL, NULL), class = "stormbreaker_retired_learner")
})

test_that("randomGLM and GAM boosting remain constructible", {
  for (id in c("random_glm_binary", "gamboost_binary", "gamboost_cox",
    "gamboost_gehan_aft", "gamboost_loglogistic_aft", "gamboost_lognormal_aft",
    "gamboost_weibull_aft", "gamboost_unoc")) {
    expect_identical(sb_pipeline(learner = id)$learner, id)
  }
  expect_identical(sb_slope_randomglm_binary_params("random_glm_binary")$n_bags, 100L)
})
