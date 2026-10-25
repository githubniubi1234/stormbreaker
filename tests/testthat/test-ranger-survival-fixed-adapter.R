test_that("C-index execution threads preserve forest and portable predictions", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")
  d <- make_sb_fixture(n = 54L, seed = 7601L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_ranger_survival_learner_adapter(
    num_trees = 18L, learner = "ranger_cindex_survival")
  one <- withr::with_options(list(stormbreaker.ranger_cindex_threads = 1L),
    adapter$fit(x, d$time, d$event))
  four <- withr::with_options(list(stormbreaker.ranger_cindex_threads = 4L),
    adapter$fit(x, d$time, d$event))
  expect_identical(one$params$num_threads, 1L)
  expect_identical(four$params$num_threads, 4L)
  expect_identical(one$forest, four$forest)
  expect_identical(adapter$predict(one, x), adapter$predict(four, x))
  context <- list(x = x, time = d$time, event = d$event,
    prediction_horizon = 7, time_unit = "month")
  s1 <- adapter$extract_state(one, context)
  s4 <- unserialize(serialize(adapter$extract_state(four, context), NULL))
  expect_identical(adapter$predict_state(s1, x), adapter$predict_state(s4, x))
  expect_equal(one$feature_importance, four$feature_importance, tolerance = 1e-12)
})

test_that("fixed ranger survival state reconstructs risk and direct survival", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  d <- make_sb_fixture(n = 54L, seed = 7601L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(4L, 19L, 37L), , drop = FALSE]
  adapter <- sb_ranger_survival_learner_adapter(
    num_trees = 18L,
    mtry = NULL,
    min_node_size = 3L,
    seed = 7602L
  )
  model <- adapter$fit(x, d$time, d$event)
  observed_risk <- adapter$predict(model, new_x)
  observed_survival <- adapter$predict_survival(
    model,
    x,
    d$time,
    d$event,
    new_x,
    prediction_horizon = 7
  )
  audit <- adapter$extract_tree_model(model)

  restored <- unserialize(audit$model_raw)
  oracle <- stats::predict(
    restored,
    data = new_x,
    num.trees = 18L,
    type = "response",
    seed = 0L,
    num.threads = 1L,
    verbose = FALSE
  )
  expected_risk <- unname(rowSums(oracle$chf))
  eligible <- which(oracle$unique.death.times <= 7)
  expected_survival <- if (length(eligible) == 0L) {
    rep(1, nrow(new_x))
  } else {
    unname(oracle$survival[, max(eligible)])
  }
  state <- adapter$extract_state(model, list(
    x = x,
    time = d$time,
    event = d$event,
    prediction_horizon = 7,
    time_unit = "month"
  ))
  state_prediction <- adapter$predict_state(state, new_x)

  expect_s3_class(model, "stormbreaker_ranger_survival")
  expect_identical(model$effective_mtry, 1L)
  expect_equal(observed_risk, expected_risk, tolerance = 0)
  expect_equal(observed_survival$probability, expected_survival, tolerance = 0)
  expect_identical(observed_survival$prediction_horizon, 7)
  expect_identical(
    observed_survival$source,
    "ranger_survival_forest_direct"
  )
  expect_identical(observed_survival$ties, "not_applicable")
  expect_false("cumulative_baseline_hazard" %in% names(observed_survival))
  expect_identical(audit$num_trees, 18L)
  expect_identical(audit$effective_mtry, 1L)
  expect_identical(audit$features, names(x))
  expect_identical(audit$params, adapter$params)
  expect_identical(audit$importance_type, "impurity")
  expect_identical(audit$model_serialization, "R_xdr_v3_ranger.forest")
  expect_type(audit$model_raw, "raw")
  expect_gt(length(audit$model_raw), 0L)
  expect_identical(audit$feature_importance$feature, names(x))
  expect_true(all(is.finite(audit$feature_importance$importance)))
  expect_identical(state$learner, "ranger_survival")
  expect_false("forest" %in% names(state))
  expect_identical(
    state$survival_probability$source,
    "ranger_survival_forest_direct"
  )
  expect_equal(state_prediction$risk, expected_risk, tolerance = 0)
  expect_equal(
    state_prediction$survival_probability,
    expected_survival,
    tolerance = 0
  )
  expect_identical(state_prediction$prediction_horizon, 7)
  expect_identical(state_prediction$time_unit, "month")

  risk_only_state <- adapter$extract_state(model, list())
  expect_null(risk_only_state$survival_probability)
  expect_identical(
    names(adapter$predict_state(risk_only_state, new_x)),
    "risk"
  )
})

test_that("fixed ranger survival is deterministic and batching independent", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  d <- make_sb_fixture(n = 50L, seed = 7611L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_ranger_survival_learner_adapter(
    num_trees = 15L,
    mtry = 2L,
    min_node_size = 3L,
    seed = 7612L
  )
  first <- adapter$fit(x[1:40, , drop = FALSE], d$time[1:40], d$event[1:40])
  second <- adapter$fit(x[1:40, , drop = FALSE], d$time[1:40], d$event[1:40])
  new_x <- x[41:50, , drop = FALSE]

  batch <- adapter$predict(first, new_x)
  repeated <- adapter$predict(second, new_x)
  one_at_a_time <- vapply(seq_len(nrow(new_x)), function(row) {
    adapter$predict(first, new_x[row, , drop = FALSE])
  }, numeric(1))
  first_survival <- adapter$predict_survival(
    first,
    x[1:40, , drop = FALSE],
    d$time[1:40],
    d$event[1:40],
    new_x,
    7
  )$probability
  second_survival <- adapter$predict_survival(
    second,
    x[1:40, , drop = FALSE],
    d$time[1:40],
    d$event[1:40],
    new_x,
    7
  )$probability

  expect_equal(batch, repeated, tolerance = 0)
  expect_equal(batch, one_at_a_time, tolerance = 0)
  expect_equal(first_survival, second_survival, tolerance = 0)
  expect_error(
    adapter$predict(first, new_x[c("gene_b", "gene_a", "gene_c")]),
    class = "stormbreaker_learner_schema_mismatch"
  )
})

test_that("fixed ranger survival failures remain explicit", {
  d <- make_sb_fixture(n = 24L, seed = 7621L)
  x <- d[c("gene_a", "gene_b")]
  adapter <- sb_ranger_survival_learner_adapter(
    num_trees = 5L,
    min_node_size = 3L,
    seed = 7622L
  )

  local_mocked_bindings(sb_ranger_available = function() FALSE)
  expect_error(
    adapter$fit(x, d$time, d$event),
    class = "stormbreaker_missing_ranger"
  )

  local_mocked_bindings(sb_ranger_available = function() TRUE)
  expect_error(
    adapter$fit(x, d$time, rep(0L, nrow(d))),
    class = "stormbreaker_ranger_survival_no_event"
  )
  too_wide <- sb_ranger_survival_learner_adapter(
    num_trees = 5L,
    mtry = 3L,
    min_node_size = 3L,
    seed = 7623L
  )
  expect_error(
    too_wide$fit(x, d$time, d$event),
    class = "stormbreaker_ranger_mtry_exceeds_features"
  )
})

test_that("ranger survival runs audited CV and final direct prediction", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  d <- make_sb_fixture(n = 72L, seed = 7631L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  horizon_task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )
  risk_task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month"
  )
  pipeline <- sb_pipeline(
    learner = "ranger_survival",
    ranger_num_trees = 10L,
    ranger_min_node_size = 3L,
    ranger_seed = 7632L
  )
  screened <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "ranger_survival",
    selector_n_features = 2L,
    ranger_num_trees = 10L,
    ranger_mtry = 1L,
    ranger_min_node_size = 3L,
    ranger_seed = 7633L
  )
  result <- sb_run_cv(horizon_task, pipeline, folds = 3L, seed = 7634L)
  screened_result <- sb_run_cv(
    horizon_task,
    screened,
    folds = 3L,
    seed = 7634L
  )
  risk_result <- sb_run_cv(risk_task, pipeline, folds = 3L, seed = 7634L)

  expect_identical(pipeline$task_type, "survival")
  expect_identical(pipeline$risk_type, "survival_forest_mortality")
  expect_identical(pipeline$risk_reference, "ranger_native")
  expect_identical(pipeline$risk_formula, "row_sums_cumulative_hazard")
  expect_identical(pipeline$learner_params$mtry_rule, "floor_sqrt_p")
  expect_identical(result$status, "success")
  expect_true(result$complete_oof)
  expect_true(result$complete_survival_probability_oof)
  expect_true(all(is.finite(result$predictions$risk)))
  expect_true(all(is.finite(result$predictions$survival_probability)))
  expect_identical(
    unique(result$predictions$learner),
    "ranger_survival"
  )
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
    function(audit) identical(audit$status, "not_applicable") &&
      identical(audit$source, "ranger_survival_forest_direct") &&
      is.na(audit$cumulative_baseline_hazard),
    logical(1)
  )))
  expect_identical(
    result$provenance$ranger_version,
    as.character(utils::packageVersion("ranger"))
  )
  expect_identical(screened_result$status, "success")
  expect_length(screened_result$selections, 3L)
  expect_true(all(vapply(
    screened_result$selections,
    function(selection) length(selection$selected_features) == 2L,
    logical(1)
  )))
  expect_identical(risk_result$status, "success")
  expect_true(risk_result$complete_oof)
  expect_false("survival_probability" %in% names(risk_result$predictions))
  expect_false("baseline_hazards" %in% names(risk_result))

  fit <- sb_fit_final(horizon_task, pipeline)
  rows <- c(5L, 21L, 49L)
  new_data <- d[rows, c("patient_id", rev(predictors)), drop = FALSE]
  prediction <- sb_predict_final(fit, new_data, "patient_id")
  selected_x <- new_data[fit$required_features]
  direct <- sb_predict_ranger_survival_learner_state(
    fit$learner_state,
    selected_x
  )

  expect_identical(fit$status, "success")
  expect_s3_class(fit$learner_state, "stormbreaker_learner_state")
  expect_identical(fit$learner_state$learner, "ranger_survival")
  expect_type(fit$learner_state$model_raw, "raw")
  expect_null(fit$learner_state$baseline_hazard)
  expect_identical(prediction$status, "success")
  expect_equal(prediction$predictions$risk, direct$risk, tolerance = 0)
  expect_equal(
    prediction$predictions$survival_probability,
    direct$survival_probability,
    tolerance = 0
  )
  expect_identical(
    fit$provenance$ranger_version,
    as.character(utils::packageVersion("ranger"))
  )
})

test_that("ranger survival informational warnings retain native results", {
  skip_if_not_installed("ranger")
  x <- data.frame(a = sin(1:30), b = cos(1:30))
  adapter <- sb_ranger_survival_learner_adapter(num_trees = 5L)
  reference <- adapter$fit(x, seq_len(30), rep(0:1, 15L))
  native <- ranger::ranger
  local_mocked_bindings(ranger = function(...) {
    warning("informational native warning")
    native(...)
  }, .package = "ranger")
  seen <- character()
  observed <- withCallingHandlers(adapter$fit(x, seq_len(30), rep(0:1, 15L)), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(seen, "informational native warning")
  expect_identical(observed$forest, reference$forest)
  expect_identical(adapter$predict(observed, x), adapter$predict(reference, x))
})

test_that("ranger survival still rejects ignored model arguments", {
  skip_if_not_installed("ranger")
  x <- data.frame(a = sin(1:30), b = cos(1:30))
  adapter <- sb_ranger_survival_learner_adapter(num_trees = 5L)
  message <- "Unused arguments: misspelled_parameter"
  local_mocked_bindings(ranger = function(...) {
    warning(message)
    stop("ignored argument warning should stop fitting first")
  }, .package = "ranger")
  expect_s3_class(tryCatch(adapter$fit(x, seq_len(30), rep(0:1, 15L)), error = identity), "stormbreaker_ranger_survival_warning")
  message <- "Argument 'num.random.splits' ignored if splitrule is not 'extratrees'."
  expect_s3_class(tryCatch(adapter$fit(x, seq_len(30), rep(0:1, 15L)), error = identity), "stormbreaker_ranger_survival_warning")
})

test_that("ranger survival still rejects incomplete forests", {
  skip_if_not_installed("ranger")
  x <- data.frame(a = sin(1:30), b = cos(1:30))
  adapter <- sb_ranger_survival_learner_adapter(num_trees = 5L)
  native <- ranger::ranger
  local_mocked_bindings(ranger = function(...) {
    fitted <- native(...)
    fitted$forest <- NULL
    fitted
  }, .package = "ranger")
  expect_s3_class(tryCatch(adapter$fit(x, seq_len(30), rep(0:1, 15L)), error = identity), "stormbreaker_incomplete_ranger_fit")
})
