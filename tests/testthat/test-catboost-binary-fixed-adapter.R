test_that("fixed CatBoost raw model reconstructs class-one probabilities", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  x <- data.frame(
    gene_a = c(-2, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
    gene_b = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
    unused = 1,
    check.names = FALSE
  )
  outcome <- c(0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 1L)
  new_x <- data.frame(
    gene_a = c(-1.5, 0.75, 3.25),
    gene_b = c(1, 0, 1),
    unused = 1,
    check.names = FALSE
  )
  adapter <- sb_catboost_binary_learner_adapter(
    iterations = 12L,
    learning_rate = 0.1,
    depth = 3L,
    l2_leaf_reg = 2.5,
    seed = 101L
  )
  model <- adapter$fit(x, outcome)
  observed <- adapter$predict(model, new_x)
  audit <- adapter$extract_tree_model(model)
  engine_params <- catboost::catboost.get_model_params(model$booster)
  expected_flat_params <- list(
    iterations = 12L,
    learning_rate = 0.1,
    depth = 3L,
    l2_leaf_reg = 2.5,
    random_seed = 101L,
    loss_function = "Logloss",
    eval_metric = "Logloss",
    task_type = "CPU",
    thread_count = 1L,
    bootstrap_type = "Bayesian",
    bagging_temperature = 1L,
    random_strength = 1L,
    rsm = 1L,
    grow_policy = "SymmetricTree",
    boosting_type = "Plain",
    use_best_model = FALSE,
    boost_from_average = TRUE,
    allow_writing_files = FALSE,
    save_snapshot = FALSE,
    logging_level = "Silent"
  )
  restored <- unserialize(audit$model_raw)
  rebuilt_pool <- catboost::catboost.load_pool(as.matrix(new_x))
  rebuilt <- as.numeric(catboost::catboost.predict(
    restored,
    rebuilt_pool,
    prediction_type = "Probability",
    thread_count = 1L
  ))

  expect_s3_class(model, "stormbreaker_catboost_binary")
  expect_identical(
    engine_params$flat_params[names(expected_flat_params)],
    expected_flat_params
  )
  expect_identical(engine_params$boosting_options$iterations, 12L)
  expect_identical(engine_params$boosting_options$od_config$type, "None")
  expect_identical(model$booster$tree_count, 12L)
  expect_equal(observed, rebuilt, tolerance = 0)
  expect_true(all(is.finite(observed)))
  expect_true(all(observed >= 0 & observed <= 1))
  expect_identical(audit$iterations, 12L)
  expect_identical(audit$tree_count, 12L)
  expect_identical(audit$features, names(x))
  expect_identical(audit$params, adapter$params)
  expect_identical(audit$model_serialization, "R_xdr_v3_catboost.Model")
  expect_type(audit$model_raw, "raw")
  expect_gt(length(audit$model_raw), 0L)
  expect_named(
    audit$feature_importance,
    c("feature", "importance"),
    ignore.order = FALSE
  )
  expect_identical(audit$feature_importance$feature, names(x))
  expect_identical(
    audit$feature_importance$importance[
      audit$feature_importance$feature == "unused"
    ],
    0
  )
})

test_that("fixed CatBoost predictions do not depend on assessment batching", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_binary_fixture(n = 50L, seed = 4201L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_catboost_binary_learner_adapter(
    iterations = 10L,
    learning_rate = 0.1,
    depth = 3L,
    l2_leaf_reg = 3,
    seed = 4301L
  )
  model <- adapter$fit(x[1:40, , drop = FALSE], d$outcome[1:40])
  new_x <- x[41:50, , drop = FALSE]

  batch <- adapter$predict(model, new_x)
  one_at_a_time <- vapply(seq_len(nrow(new_x)), function(row) {
    adapter$predict(model, new_x[row, , drop = FALSE])
  }, numeric(1))
  order <- c(10L, 2L, 8L, 4L, 6L, 1L, 9L, 3L, 7L, 5L)
  reordered <- adapter$predict(model, new_x[order, , drop = FALSE])
  restored_order <- reordered[match(seq_along(order), order)]

  expect_equal(batch, one_at_a_time, tolerance = 0)
  expect_equal(batch, restored_order, tolerance = 0)
})

test_that("fixed CatBoost is deterministic and locks predictor schema", {
  skip_if_not_installed("catboost", minimum_version = "1.2.0")

  d <- make_sb_binary_fixture(n = 40L, seed = 4401L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_catboost_binary_learner_adapter(
    iterations = 8L,
    learning_rate = 0.1,
    depth = 3L,
    l2_leaf_reg = 3,
    seed = 4501L
  )
  first <- adapter$fit(x, d$outcome)
  second <- adapter$fit(x, d$outcome)

  expect_equal(
    adapter$predict(first, x),
    adapter$predict(second, x),
    tolerance = 0
  )
  expect_equal(
    adapter$extract_tree_model(first)$feature_importance,
    adapter$extract_tree_model(second)$feature_importance,
    tolerance = 0
  )

  schema_error <- capture_sb_error(
    adapter$predict(first, x[c("gene_b", "gene_a", "gene_c")])
  )
  expect_s3_class(schema_error, "stormbreaker_learner_schema_mismatch")
})

test_that("fixed CatBoost dependency and class failures remain explicit", {
  x <- data.frame(
    gene_a = seq_len(8),
    gene_b = rev(seq_len(8))
  )
  adapter <- sb_catboost_binary_learner_adapter(
    iterations = 5L,
    learning_rate = 0.1,
    depth = 2L,
    l2_leaf_reg = 3,
    seed = 1L
  )

  local_mocked_bindings(sb_catboost_available = function() FALSE)
  missing_error <- capture_sb_error(
    adapter$fit(x, rep(c(0L, 1L), each = 4L))
  )
  expect_s3_class(missing_error, "stormbreaker_missing_catboost")

  local_mocked_bindings(sb_catboost_available = function() TRUE)
  class_error <- capture_sb_error(adapter$fit(x, rep(0L, 8L)))
  expect_s3_class(class_error, "stormbreaker_catboost_single_class")
})
