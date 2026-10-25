test_that("fixed LightGBM model text reconstructs class-one probabilities", {
  skip_if_not_installed("lightgbm", minimum_version = "4.0.0")

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
  adapter <- sb_lightgbm_binary_learner_adapter(
    nrounds = 12L,
    learning_rate = 0.1,
    num_leaves = 7L,
    min_data_in_leaf = 2L,
    seed = 101L
  )
  model <- adapter$fit(x, outcome)
  observed <- adapter$predict(model, new_x)
  audit <- adapter$extract_tree_model(model)
  engine_params <- model$booster$params
  expected_engine_params <- list(
    learning_rate = 0.1,
    num_leaves = 7L,
    min_data_in_leaf = 2L,
    seed = 101L,
    objective = "binary",
    boosting = "gbdt",
    device_type = "cpu",
    num_threads = 1L,
    deterministic = TRUE,
    force_col_wise = TRUE,
    feature_fraction = 1,
    bagging_fraction = 1,
    bagging_freq = 0L
  )
  restored <- lightgbm::lgb.load(model_str = audit$model_string)
  rebuilt <- as.numeric(stats::predict(
    restored,
    as.matrix(new_x),
    type = "response",
    num_iteration = audit$nrounds,
    params = list(num_threads = 1L)
  ))

  expect_s3_class(model, "stormbreaker_lightgbm_binary")
  expect_identical(
    engine_params[names(expected_engine_params)],
    expected_engine_params
  )
  expect_identical(engine_params$num_iterations, 12L)
  expect_identical(as.integer(model$booster$num_iter()), 12L)
  expect_identical(as.integer(model$booster$num_trees()), 12L)
  expect_identical(model$booster$best_iter, -1L)
  expect_identical(model$booster$record_evals, list())
  expect_equal(observed, rebuilt, tolerance = 0)
  expect_true(all(is.finite(observed)))
  expect_true(all(observed >= 0 & observed <= 1))
  expect_identical(audit$nrounds, 12L)
  expect_identical(audit$features, names(x))
  expect_identical(audit$params, adapter$params)
  expect_type(audit$model_string, "character")
  expect_length(audit$model_string, 1L)
  expect_named(
    audit$feature_importance,
    c("feature", "gain", "cover", "frequency"),
    ignore.order = FALSE
  )
  expect_identical(audit$feature_importance$feature, names(x))
  expect_identical(
    audit$feature_importance$frequency[audit$feature_importance$feature ==
      "unused"],
    0L
  )
})

test_that("fixed LightGBM predictions do not depend on assessment batching", {
  skip_if_not_installed("lightgbm", minimum_version = "4.0.0")

  d <- make_sb_binary_fixture(n = 50L, seed = 2801L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_lightgbm_binary_learner_adapter(
    nrounds = 10L,
    learning_rate = 0.1,
    num_leaves = 7L,
    min_data_in_leaf = 3L,
    seed = 2901L
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

test_that("fixed LightGBM is deterministic and locks predictor schema", {
  skip_if_not_installed("lightgbm", minimum_version = "4.0.0")

  d <- make_sb_binary_fixture(n = 40L, seed = 3001L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_lightgbm_binary_learner_adapter(
    nrounds = 8L,
    learning_rate = 0.1,
    num_leaves = 7L,
    min_data_in_leaf = 3L,
    seed = 3101L
  )
  first <- adapter$fit(x, d$outcome)
  second <- adapter$fit(x, d$outcome)

  expect_identical(
    adapter$extract_tree_model(first)$model_string,
    adapter$extract_tree_model(second)$model_string
  )
  expect_equal(
    adapter$predict(first, x),
    adapter$predict(second, x),
    tolerance = 0
  )

  schema_error <- capture_sb_error(
    adapter$predict(first, x[c("gene_b", "gene_a", "gene_c")])
  )
  expect_s3_class(schema_error, "stormbreaker_learner_schema_mismatch")
})

test_that("fixed LightGBM dependency and class failures remain explicit", {
  x <- data.frame(
    gene_a = seq_len(8),
    gene_b = rev(seq_len(8))
  )
  adapter <- sb_lightgbm_binary_learner_adapter(
    nrounds = 5L,
    learning_rate = 0.1,
    num_leaves = 3L,
    min_data_in_leaf = 1L,
    seed = 1L
  )

  local_mocked_bindings(sb_lightgbm_available = function() FALSE)
  missing_error <- capture_sb_error(
    adapter$fit(x, rep(c(0L, 1L), each = 4L))
  )
  expect_s3_class(missing_error, "stormbreaker_missing_lightgbm")

  local_mocked_bindings(sb_lightgbm_available = function() TRUE)
  class_error <- capture_sb_error(adapter$fit(x, rep(0L, 8L)))
  expect_s3_class(class_error, "stormbreaker_lightgbm_single_class")
})
