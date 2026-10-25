test_that("a single allowed GBM tree is fitted without inventing a CV search", {
  skip_if_not_installed("gbm")
  d <- make_sb_fixture(n = 100L, seed = 8612L)
  features <- c("gene_a", "gene_b", "gene_c")
  x <- d[features]
  params <- sb_gbm_cox_cv_params(inner_folds = 3L, maximum_trees = 1L)
  model <- sb_fit_gbm_cox_cv(x, d$time, d$event, params)
  expect_identical(model$n_trees, 1L)
  expect_identical(model$tuning$cv_error, numeric())
  expect_identical(model$tuning$search_status, "single_candidate_no_search")
  expect_null(model$booster$cv.error)
  formula_data <- sb_additional_survival_formula_data(x, d$time, d$event)
  native <- withr::with_seed(params$seed, gbm::gbm(
    formula = formula_data$formula, data = formula_data$data,
    distribution = params$distribution, n.trees = 1L,
    interaction.depth = params$interaction_depth, n.minobsinnode = params$min_obs_node,
    shrinkage = params$shrinkage, bag.fraction = params$bag_fraction,
    cv.folds = 0L, n.cores = 1L, keep.data = FALSE, verbose = FALSE))
  expect_equal(sb_predict_gbm_cox(model, x),
    unname(as.numeric(stats::predict(native, newdata = x, n.trees = 1L, type = "link"))), tolerance = 0)
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
  fit <- sb_fit_final(task, sb_pipeline(learner = "gbm_cox_cv", gbm_cv_max_trees = 1L,
    learner_tuning_folds = 3L), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  expect_identical(fit$learner_state$tuning$search_status, "single_candidate_no_search")
  err <- tryCatch(sb_boosting_diagnostics(fit, "cv_curve"), error = identity)
  expect_s3_class(err, "stormbreaker_unavailable_diagnostic")
  loss <- sb_boosting_diagnostics(fit, "training_loss")
  expect_equal(loss$table$loss, fit$diagnostic_artifacts$learner$booster$train.error)
  expect_equal(nrow(loss$table), 1L)
  expect_identical(serialize(fit, NULL, version = 2L), before)
  pipeline <- fit$pipeline
  compact <- sb_fit_final(task, pipeline)
  expect_identical(compact$status, "success")
  expect_identical(compact$learner_state, fit$learner_state)
  cv <- sb_run_cv(task, pipeline, folds = 3L, seed = 7733L)
  expect_identical(cv$status, "success")
  expect_true(cv$complete_oof)
  expect_true(all(is.finite(cv$predictions$risk)))
})
