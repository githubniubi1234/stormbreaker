test_that("fixed CatBoost Cox state reconstructs risk and horizon probability", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_fixture(n = 48L, seed = 7501L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(3L, 17L, 29L), , drop = FALSE]
  adapter <- sb_catboost_cox_learner_adapter(
    iterations = 12L,
    learning_rate = 0.1,
    depth = 3L,
    l2_leaf_reg = 2.5,
    seed = 7502L
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
  rebuilt_pool <- catboost::catboost.load_pool(as.matrix(new_x))
  rebuilt_risk <- as.numeric(catboost::catboost.predict(
    restored,
    rebuilt_pool,
    prediction_type = "RawFormulaVal",
    ntree_start = 0L,
    ntree_end = 12L,
    thread_count = 1L
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

  expect_s3_class(model, "stormbreaker_catboost_cox")
  expect_equal(observed_risk, rebuilt_risk, tolerance = 0)
  expect_identical(model$booster$tree_count, 12L)
  expect_identical(adapter$params$loss_function, "Cox")
  expect_identical(adapter$params$eval_metric, "Cox")
  expect_identical(probability$prediction_horizon, 7)
  expect_identical(probability$cumulative_baseline_hazard, expected_hazard)
  expect_identical(probability$source, "catboost_cox_breslow_raw_margin")
  expect_identical(probability$ties, "breslow")
  expect_equal(probability$probability, expected_probability, tolerance = 0)
  expect_identical(audit$features, names(x))
  expect_identical(audit$params, adapter$params)
  expect_identical(audit$model_serialization, "R_xdr_v3_catboost.Model")
  expect_type(audit$model_raw, "raw")
  expect_gt(length(audit$model_raw), 0L)
  expect_identical(state$learner, "catboost_cox")
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

test_that("fixed CatBoost Cox is deterministic and batching independent", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_fixture(n = 50L, seed = 7511L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_catboost_cox_learner_adapter(
    iterations = 10L,
    learning_rate = 0.1,
    depth = 3L,
    l2_leaf_reg = 3,
    seed = 7512L
  )
  first <- adapter$fit(x[1:40, , drop = FALSE], d$time[1:40], d$event[1:40])
  second <- adapter$fit(x[1:40, , drop = FALSE], d$time[1:40], d$event[1:40])
  new_x <- x[41:50, , drop = FALSE]

  batch <- adapter$predict(first, new_x)
  repeated <- adapter$predict(second, new_x)
  one_at_a_time <- vapply(seq_len(nrow(new_x)), function(row) {
    adapter$predict(first, new_x[row, , drop = FALSE])
  }, numeric(1))

  expect_equal(batch, repeated, tolerance = 0)
  expect_equal(batch, one_at_a_time, tolerance = 0)
  expect_error(
    adapter$predict(first, new_x[c("gene_b", "gene_a", "gene_c")]),
    class = "stormbreaker_learner_schema_mismatch"
  )
})

test_that("fixed CatBoost Cox failures remain explicit", {
  d <- make_sb_fixture(n = 24L, seed = 7521L)
  x <- d[c("gene_a", "gene_b")]
  adapter <- sb_catboost_cox_learner_adapter(
    iterations = 5L,
    learning_rate = 0.1,
    depth = 2L,
    l2_leaf_reg = 3,
    seed = 7522L
  )

  local_mocked_bindings(sb_catboost_available = function() FALSE)
  expect_error(
    adapter$fit(x, d$time, d$event),
    class = "stormbreaker_missing_catboost"
  )

  local_mocked_bindings(sb_catboost_available = function() TRUE)
  nonpositive <- d$time
  nonpositive[[1L]] <- 0
  expect_error(
    adapter$fit(x, nonpositive, d$event),
    class = "stormbreaker_catboost_cox_nonpositive_time"
  )
  expect_error(
    adapter$fit(x, d$time, rep(0L, nrow(d))),
    class = "stormbreaker_catboost_cox_no_event"
  )
})

test_that("CatBoost Cox runs audited survival CV and final prediction", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_fixture(n = 72L, seed = 7531L)
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
    learner = "catboost_cox",
    catboost_iterations = 8L,
    catboost_learning_rate = 0.1,
    catboost_depth = 3L,
    catboost_l2_leaf_reg = 3,
    catboost_seed = 7532L
  )
  screened <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "catboost_cox",
    selector_n_features = 2L,
    catboost_iterations = 8L,
    catboost_depth = 3L,
    catboost_seed = 7533L
  )
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 7534L)
  screened_result <- sb_run_cv(task, screened, folds = 3L, seed = 7534L)

  expect_identical(pipeline$task_type, "survival")
  expect_identical(pipeline$risk_type, "cox_raw_margin")
  expect_identical(pipeline$risk_reference, "catboost_native")
  expect_identical(pipeline$risk_formula, "catboost_raw_formula_value")
  expect_identical(result$status, "success")
  expect_true(result$complete_oof)
  expect_true(result$complete_survival_probability_oof)
  expect_true(all(is.finite(result$predictions$risk)))
  expect_true(all(is.finite(result$predictions$survival_probability)))
  expect_identical(unique(result$predictions$learner), "catboost_cox")
  expect_identical(result$metrics$harrell_c$status, "success")
  expect_identical(result$metrics$ipcw_brier$status, "success")
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
      "catboost_cox_breslow_raw_margin"
    ),
    logical(1)
  )))
  expect_identical(
    result$provenance$catboost_version,
    as.character(utils::packageVersion("catboost"))
  )
  expect_identical(screened_result$status, "success")
  expect_length(screened_result$selections, 3L)
  expect_true(all(vapply(
    screened_result$selections,
    function(selection) length(selection$selected_features) == 2L,
    logical(1)
  )))

  fit <- sb_fit_final(task, pipeline)
  rows <- c(5L, 21L, 49L)
  new_data <- d[rows, c("patient_id", rev(predictors)), drop = FALSE]
  prediction <- sb_predict_final(fit, new_data, "patient_id")
  selected_x <- new_data[fit$required_features]
  direct <- sb_predict_catboost_cox_learner_state(
    fit$learner_state,
    selected_x
  )

  expect_identical(fit$status, "success")
  expect_s3_class(fit$learner_state, "stormbreaker_learner_state")
  expect_identical(fit$learner_state$learner, "catboost_cox")
  expect_type(fit$learner_state$model_raw, "raw")
  expect_identical(prediction$status, "success")
  expect_equal(prediction$predictions$risk, direct$risk, tolerance = 0)
  expect_equal(
    prediction$predictions$survival_probability,
    direct$survival_probability,
    tolerance = 0
  )
  expect_identical(
    fit$provenance$catboost_version,
    as.character(utils::packageVersion("catboost"))
  )
})
