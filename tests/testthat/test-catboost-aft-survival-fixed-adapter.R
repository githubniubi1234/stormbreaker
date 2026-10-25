test_that("fixed CatBoost AFT reproduces native risk and Normal survival", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  withr::local_seed(8801L)
  x <- data.frame(
    gene_a = stats::rnorm(48L),
    gene_b = stats::rnorm(48L),
    gene_c = stats::rnorm(48L),
    check.names = FALSE
  )
  time <- exp(1 - 0.7 * x$gene_a + 0.25 * x$gene_b +
    stats::rnorm(48L, sd = 0.35))
  event <- rep(c(1L, 1L, 0L), length.out = nrow(x))
  new_x <- x[c(4L, 19L, 37L), , drop = FALSE]
  adapter <- sb_catboost_aft_survival_learner_adapter(
    iterations = 10L,
    learning_rate = 0.1,
    depth = 3L,
    l2_leaf_reg = 2.5,
    seed = 8802L
  )
  model <- adapter$fit(x, time, event)

  native_pool <- catboost::catboost.load_pool(
    data = as.matrix(x),
    label = cbind(time, ifelse(event == 1L, time, -1))
  )
  native <- catboost::catboost.train(
    learn_pool = native_pool,
    test_pool = NULL,
    params = adapter$params
  )
  native_margin <- as.numeric(catboost::catboost.predict(
    native,
    catboost::catboost.load_pool(as.matrix(new_x)),
    prediction_type = "RawFormulaVal",
    ntree_start = 0L,
    ntree_end = 10L,
    thread_count = 1L
  ))
  risk <- adapter$predict(model, new_x)
  survival <- adapter$predict_survival(
    model,
    x,
    time,
    event,
    new_x,
    prediction_horizon = 3
  )
  audit <- adapter$extract_tree_model(model)
  state <- unserialize(serialize(adapter$extract_state(model, list(
    x = x,
    time = time,
    event = event,
    prediction_horizon = 3,
    time_unit = "month"
  )), NULL))
  replay <- adapter$predict_state(state, new_x)

  expect_s3_class(model, "stormbreaker_catboost_aft_survival")
  expect_identical(adapter$params$loss_function,
    "SurvivalAft:dist=Normal;scale=1"
  )
  expect_equal(risk, -native_margin, tolerance = 0)
  expect_equal(
    survival$probability,
    stats::pnorm((native_margin - log(3)) / 1),
    tolerance = 0
  )
  expect_identical(survival$source, "catboost_aft_normal_native")
  expect_identical(survival$ties, "not_applicable")
  expect_identical(audit$tree_count, 10L)
  expect_identical(audit$features, names(x))
  expect_equal(replay$risk, risk, tolerance = 0)
  expect_equal(
    replay$survival_probability,
    survival$probability,
    tolerance = 0
  )
  expect_identical(replay$prediction_horizon, 3)
  expect_identical(replay$time_unit, "month")
})

test_that("fixed CatBoost AFT keeps dependency and schema failures explicit", {
  d <- make_sb_fixture(n = 24L, seed = 8811L)
  x <- d[c("gene_a", "gene_b")]
  adapter <- sb_catboost_aft_survival_learner_adapter(
    iterations = 5L,
    learning_rate = 0.1,
    depth = 2L,
    l2_leaf_reg = 3,
    seed = 8812L
  )

  local_mocked_bindings(sb_catboost_available = function() FALSE)
  missing <- capture_sb_error(adapter$fit(x, d$time, d$event))
  expect_s3_class(missing, "stormbreaker_missing_catboost")

  local_mocked_bindings(sb_catboost_available = function() TRUE)
  nonpositive <- d$time
  nonpositive[[1L]] <- 0
  bad_time <- capture_sb_error(adapter$fit(x, nonpositive, d$event))
  no_event <- capture_sb_error(adapter$fit(x, d$time, rep(0L, nrow(d))))
  expect_s3_class(bad_time, "stormbreaker_catboost_aft_nonpositive_time")
  expect_s3_class(no_event, "stormbreaker_catboost_aft_no_event")
})

test_that("CatBoost AFT runs audited survival CV and final prediction", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_fixture(n = 60L, seed = 8821L)
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
    learner = "catboost_aft_survival",
    catboost_iterations = 6L,
    catboost_learning_rate = 0.1,
    catboost_depth = 3L,
    catboost_l2_leaf_reg = 3,
    catboost_seed = 8822L
  )
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 8823L)

  expect_identical(pipeline$task_type, "survival")
  expect_identical(pipeline$risk_type, "aft_log_time")
  expect_identical(pipeline$risk_reference, "catboost_aft_normal_native")
  expect_identical(
    pipeline$risk_formula,
    "negative_catboost_aft_log_time_margin"
  )
  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(result$complete_survival_probability_oof, TRUE)
  expect_length(result$tree_models, 3L)
  expect_identical(
    vapply(result$tree_models, `[[`, "", "status"),
    rep("success", 3L)
  )
  expect_identical(
    vapply(result$baseline_hazards, `[[`, "", "status"),
    rep("not_applicable", 3L)
  )
  expect_identical(
    vapply(result$baseline_hazards, `[[`, "", "source"),
    rep("catboost_aft_normal_native", 3L)
  )
  expect_identical(
    result$provenance$catboost_version,
    as.character(utils::packageVersion("catboost"))
  )

  fit <- sb_fit_final(task, pipeline)
  rows <- c(5L, 21L, 49L)
  new_data <- d[rows, c("patient_id", rev(predictors)), drop = FALSE]
  prediction <- sb_predict_final(fit, new_data, "patient_id")
  direct <- sb_predict_catboost_aft_survival_learner_state(
    fit$learner_state,
    new_data[fit$required_features]
  )

  expect_identical(fit$status, "success")
  expect_identical(
    fit$learner_state$learner,
    "catboost_aft_survival"
  )
  expect_identical(prediction$status, "success")
  expect_equal(prediction$predictions$risk, direct$risk, tolerance = 0)
  expect_equal(
    prediction$predictions$survival_probability,
    direct$survival_probability,
    tolerance = 0
  )
})
