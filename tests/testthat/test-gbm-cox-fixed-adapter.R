test_that("fixed GBM-Cox state reconstructs risk and horizon probability", {
  skip_if_not_installed("gbm", minimum_version = "2.2.0")

  d <- make_sb_fixture(n = 54L, seed = 7701L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(4L, 19L, 37L), , drop = FALSE]
  adapter <- sb_gbm_cox_learner_adapter(
    n_trees = 12L,
    interaction_depth = 2L,
    min_obs_node = 3L,
    shrinkage = 0.05,
    bag_fraction = 1,
    seed = 7702L
  )
  model <- adapter$fit(x, d$time, d$event)
  observed_risk <- adapter$predict(model, new_x)
  probability <- adapter$predict_survival(
    model,
    x,
    d$time,
    d$event,
    new_x,
    prediction_horizon = 7
  )
  audit <- adapter$extract_tree_model(model)

  restored <- unserialize(audit$model_raw)
  restored_risk <- as.numeric(stats::predict(
    restored,
    newdata = as.matrix(new_x),
    n.trees = 12L,
    type = "link"
  ))
  analysis_risk <- adapter$predict(model, x)
  expected_hazard <- sb_breslow_baseline_hazard(
    d$time,
    d$event,
    analysis_risk,
    prediction_horizon = 7
  )
  expected_probability <- sb_survival_probability_from_hazard(
    observed_risk,
    expected_hazard
  )
  state <- adapter$extract_state(model, list(
    x = x,
    time = d$time,
    event = d$event,
    prediction_horizon = 7,
    time_unit = "month"
  ))
  state_prediction <- adapter$predict_state(state, new_x)

  expect_s3_class(model, "stormbreaker_gbm_cox")
  expect_identical(as.integer(model$booster$n.trees), 12L)
  expect_null(model$booster$data)
  expect_equal(observed_risk, restored_risk, tolerance = 0)
  expect_identical(probability$prediction_horizon, 7)
  expect_identical(probability$cumulative_baseline_hazard, expected_hazard)
  expect_identical(
    probability$source,
    "gbm_cox_breslow_log_hazard"
  )
  expect_identical(probability$ties, "breslow")
  expect_equal(probability$probability, expected_probability, tolerance = 0)
  expect_identical(audit$n_trees, 12L)
  expect_identical(audit$features, names(x))
  expect_identical(audit$params, adapter$params)
  expect_identical(audit$importance_type, "relative_influence")
  expect_identical(audit$model_serialization, "R_xdr_v3_gbm")
  expect_type(audit$model_raw, "raw")
  expect_gt(length(audit$model_raw), 0L)
  expect_identical(audit$feature_importance$feature, names(x))
  expect_true(all(is.finite(audit$feature_importance$importance)))
  expect_identical(state$learner, "gbm_cox")
  expect_false("booster" %in% names(state))
  expect_equal(state_prediction$risk, observed_risk, tolerance = 0)
  expect_equal(
    state_prediction$survival_probability,
    expected_probability,
    tolerance = 0
  )
  expect_identical(state_prediction$prediction_horizon, 7)
  expect_identical(state_prediction$time_unit, "month")
})

test_that("fixed GBM-Cox settings are explicit and reproducible", {
  skip_if_not_installed("gbm", minimum_version = "2.2.0")

  pipeline <- sb_pipeline(learner = "gbm_cox")
  expect_identical(
    pipeline$id,
    paste0(
      "none__gbm_cox__trees100__depth1__minnode10",
      "__shrink0.1__bag1__seed1"
    )
  )
  expect_identical(pipeline$task_type, "survival")
  expect_identical(pipeline$risk_type, "cox_log_hazard")
  expect_identical(pipeline$risk_reference, "gbm_native")
  expect_identical(pipeline$risk_formula, "gbm_link_prediction")
  expect_identical(pipeline$learner_params$distribution, "coxph")
  expect_identical(pipeline$learner_params$keep_data, FALSE)

  d <- make_sb_fixture(n = 50L, seed = 7711L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_gbm_cox_learner_adapter(
    n_trees = 10L,
    interaction_depth = 1L,
    min_obs_node = 3L,
    shrinkage = 0.1,
    bag_fraction = 0.8,
    seed = 7712L
  )
  first <- adapter$fit(x[1:40, , drop = FALSE], d$time[1:40], d$event[1:40])
  second <- adapter$fit(x[1:40, , drop = FALSE], d$time[1:40], d$event[1:40])
  new_x <- x[41:50, , drop = FALSE]

  expect_equal(
    adapter$predict(first, new_x),
    adapter$predict(second, new_x),
    tolerance = 0
  )
  expect_error(
    adapter$predict(first, new_x[c("gene_b", "gene_a", "gene_c")]),
    class = "stormbreaker_learner_schema_mismatch"
  )

  expect_error(
    sb_pipeline(learner = "gbm_cox", gbm_shrinkage = 0),
    class = "stormbreaker_bad_gbm_shrinkage"
  )
  expect_error(
    sb_pipeline(learner = "gbm_cox", gbm_bag_fraction = 1.1),
    class = "stormbreaker_bad_gbm_bag_fraction"
  )
  expect_error(
    sb_pipeline(learner = "coxph", gbm_n_trees = 10L),
    class = "stormbreaker_irrelevant_gbm_parameter"
  )
})

test_that("fixed GBM-Cox failures remain explicit", {
  d <- make_sb_fixture(n = 24L, seed = 7721L)
  x <- d[c("gene_a", "gene_b")]
  adapter <- sb_gbm_cox_learner_adapter(
    n_trees = 5L,
    min_obs_node = 3L,
    seed = 7722L
  )

  local_mocked_bindings(sb_gbm_available = function() FALSE)
  expect_error(
    adapter$fit(x, d$time, d$event),
    class = "stormbreaker_missing_gbm"
  )
  expect_error(
    adapter$fit(x, d$time, rep(0L, nrow(d))),
    class = "stormbreaker_gbm_cox_no_event"
  )
})

test_that("GBM-Cox runs selected survival CV and final prediction", {
  skip_if_not_installed("gbm", minimum_version = "2.2.0")

  d <- make_sb_fixture(n = 72L, seed = 7731L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "gbm_cox",
    selector_n_features = 2L,
    gbm_n_trees = 10L,
    gbm_min_obs_node = 3L,
    gbm_seed = 7732L
  )
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 7733L)

  expect_identical(result$status, "success")
  expect_true(result$complete_oof)
  expect_true(result$complete_survival_probability_oof)
  expect_true(all(is.finite(result$predictions$risk)))
  expect_true(all(is.finite(result$predictions$survival_probability)))
  expect_identical(unique(result$predictions$learner), "gbm_cox")
  expect_length(result$selections, 3L)
  expect_true(all(vapply(
    result$selections,
    function(selection) length(selection$selected_features) == 2L,
    logical(1)
  )))
  expect_length(result$tree_models, 3L)
  expect_true(all(vapply(
    result$tree_models,
    function(audit) identical(audit$status, "success") &&
      is.raw(audit$model_raw) && length(audit$model_raw) > 0L,
    logical(1)
  )))
  expect_true(all(vapply(
    result$baseline_hazards,
    function(audit) identical(
      audit$source,
      "gbm_cox_breslow_log_hazard"
    ),
    logical(1)
  )))
  expect_identical(
    result$provenance$gbm_version,
    as.character(utils::packageVersion("gbm"))
  )

  fit <- sb_fit_final(task, pipeline)
  rows <- c(5L, 21L, 49L)
  new_data <- d[rows, c("patient_id", rev(predictors)), drop = FALSE]
  prediction <- sb_predict_final(fit, new_data, "patient_id")
  direct <- sb_predict_gbm_cox_learner_state(
    fit$learner_state,
    new_data[fit$required_features]
  )

  expect_identical(fit$status, "success")
  expect_identical(fit$learner_state$learner, "gbm_cox")
  expect_type(fit$learner_state$model_raw, "raw")
  expect_identical(prediction$status, "success")
  expect_equal(prediction$predictions$risk, direct$risk, tolerance = 0)
  expect_equal(
    prediction$predictions$survival_probability,
    direct$survival_probability,
    tolerance = 0
  )
  expect_identical(
    fit$provenance$gbm_version,
    as.character(utils::packageVersion("gbm"))
  )
})
