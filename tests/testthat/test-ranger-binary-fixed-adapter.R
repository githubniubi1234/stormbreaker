test_that("fixed ranger binary state reconstructs class-one probabilities", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  d <- make_sb_binary_fixture(n = 54L, seed = 7701L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(4L, 19L, 37L), , drop = FALSE]
  adapter <- sb_ranger_binary_learner_adapter(
    num_trees = 18L,
    mtry = NULL,
    min_node_size = 10L,
    seed = 7702L
  )
  model <- adapter$fit(x, d$outcome)
  observed <- adapter$predict(model, new_x)
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
  )$predictions
  expected <- unname(oracle[, "1"])
  state <- adapter$extract_state(model)

  expect_s3_class(model, "stormbreaker_ranger_binary")
  expect_identical(model$effective_mtry, 1L)
  expect_equal(observed, expected, tolerance = 0)
  expect_identical(adapter$params$probability, TRUE)
  expect_identical(adapter$params$splitrule, "gini")
  expect_identical(adapter$params$min_node_size, 10L)
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
  expect_identical(state$learner, "ranger_binary")
  expect_s3_class(state, "stormbreaker_learner_state")
  expect_false("forest" %in% names(state))
  expect_equal(
    adapter$predict_state(state, new_x)$probability,
    expected,
    tolerance = 0
  )
})

test_that("fixed ranger binary is deterministic and fails explicitly", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  d <- make_sb_binary_fixture(n = 50L, seed = 7711L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_ranger_binary_learner_adapter(
    num_trees = 15L,
    mtry = 2L,
    min_node_size = 10L,
    seed = 7712L
  )
  first <- adapter$fit(x[1:40, , drop = FALSE], d$outcome[1:40])
  second <- adapter$fit(x[1:40, , drop = FALSE], d$outcome[1:40])
  new_x <- x[41:50, , drop = FALSE]
  batch <- adapter$predict(first, new_x)
  one_at_a_time <- vapply(seq_len(nrow(new_x)), function(row) {
    adapter$predict(first, new_x[row, , drop = FALSE])
  }, numeric(1))

  expect_equal(batch, adapter$predict(second, new_x), tolerance = 0)
  expect_equal(batch, one_at_a_time, tolerance = 0)
  expect_error(
    adapter$predict(first, new_x[c("gene_b", "gene_a", "gene_c")]),
    class = "stormbreaker_learner_schema_mismatch"
  )

  local_mocked_bindings(sb_ranger_available = function() FALSE)
  expect_error(
    adapter$fit(x, d$outcome),
    class = "stormbreaker_missing_ranger"
  )
  local_mocked_bindings(sb_ranger_available = function() TRUE)
  expect_error(
    adapter$fit(x, rep(0L, nrow(x))),
    class = "stormbreaker_ranger_binary_single_class"
  )
  too_wide <- sb_ranger_binary_learner_adapter(
    num_trees = 5L,
    mtry = 4L,
    seed = 7713L
  )
  expect_error(
    too_wide$fit(x, d$outcome),
    class = "stormbreaker_ranger_mtry_exceeds_features"
  )
})

test_that("ranger binary runs all-feature and fold-local Top-K workflows", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  d <- make_sb_binary_fixture(n = 72L, seed = 7721L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  all_features <- sb_pipeline(
    learner = "ranger_binary",
    ranger_num_trees = 12L,
    ranger_seed = 7722L
  )
  screened <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "ranger_binary",
    selector_n_features = 2L,
    ranger_num_trees = 12L,
    ranger_mtry = 1L,
    ranger_seed = 7723L
  )
  result <- sb_run_cv(task, all_features, folds = 3L, seed = 7724L)
  screened_result <- sb_run_cv(task, screened, folds = 3L, seed = 7724L)

  expect_identical(all_features$task_type, "binary")
  expect_identical(
    all_features$id,
    "none__ranger_binary__trees12__mtrysqrt__minnode10__seed7722"
  )
  expect_identical(
    screened$id,
    paste0(
      "univ_binary_auc_top_k__ranger_binary__k2",
      "__trees12__mtry1__minnode10__seed7723"
    )
  )
  expect_identical(all_features$learner_params$min_node_size, 10L)
  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_true(all(is.finite(result$predictions$probability)))
  expect_length(result$tree_models, 3L)
  expect_identical(
    vapply(result$tree_models, `[[`, "", "status"),
    rep("success", 3L)
  )
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("not_applicable", 3L)
  )
  expect_identical(
    result$provenance$ranger_version,
    as.character(utils::packageVersion("ranger"))
  )
  for (fold in seq_along(result$folds)) {
    split <- result$folds[[fold]]
    audit <- result$tree_models[[fold]]
    rows <- match(split$assessment_ids, d$patient_id)
    rebuilt <- stats::predict(
      unserialize(audit$model_raw),
      data = d[rows, audit$features, drop = FALSE],
      num.trees = audit$num_trees,
      type = "response",
      seed = 0L,
      num.threads = 1L,
      verbose = FALSE
    )$predictions[, "1"]
    observed <- result$predictions$probability[
      match(split$assessment_ids, result$predictions$patient_id)
    ]
    expect_equal(observed, unname(rebuilt), tolerance = 0)
  }
  expect_identical(screened_result$status, "success")
  expect_true(all(vapply(
    screened_result$selections,
    function(selection) length(selection$selected_features) == 2L,
    logical(1)
  )))

  fit <- sb_fit_final(task, screened)
  rows <- c(5L, 21L, 49L)
  new_data <- d[rows, c("patient_id", rev(predictors)), drop = FALSE]
  prediction <- sb_predict_final(fit, new_data, "patient_id")
  direct <- sb_predict_ranger_binary_learner_state(
    fit$learner_state,
    new_data[fit$required_features]
  )

  expect_identical(fit$status, "success")
  expect_identical(fit$learner_state$learner, "ranger_binary")
  expect_type(fit$learner_state$model_raw, "raw")
  expect_identical(prediction$status, "success")
  expect_equal(prediction$predictions$probability, direct$probability, tolerance = 0)
  expect_identical(
    fit$provenance$ranger_version,
    as.character(utils::packageVersion("ranger"))
  )
})

test_that("ranger binary informational warnings retain native results", {
  skip_if_not_installed("ranger")
  x <- data.frame(a = sin(1:30), b = cos(1:30))
  adapter <- sb_ranger_binary_learner_adapter(num_trees = 5L)
  reference <- adapter$fit(x, rep(0:1, 15L))
  native <- ranger::ranger
  local_mocked_bindings(ranger = function(...) {
    warning("informational native warning")
    native(...)
  }, .package = "ranger")
  seen <- character()
  observed <- withCallingHandlers(adapter$fit(x, rep(0:1, 15L)), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(seen, "informational native warning")
  expect_identical(observed$forest, reference$forest)
  expect_identical(adapter$predict(observed, x), adapter$predict(reference, x))
})

test_that("ranger binary still rejects ignored model arguments", {
  skip_if_not_installed("ranger")
  x <- data.frame(a = sin(1:30), b = cos(1:30))
  adapter <- sb_ranger_binary_learner_adapter(num_trees = 5L)
  message <- "Unused arguments: misspelled_parameter"
  local_mocked_bindings(ranger = function(...) {
    warning(message)
    stop("ignored argument warning should stop fitting first")
  }, .package = "ranger")
  expect_s3_class(tryCatch(adapter$fit(x, rep(0:1, 15L)), error = identity), "stormbreaker_ranger_binary_warning")
  message <- "Argument 'num.random.splits' ignored if splitrule is not 'extratrees'."
  expect_s3_class(tryCatch(adapter$fit(x, rep(0:1, 15L)), error = identity), "stormbreaker_ranger_binary_warning")
})

test_that("ranger binary still rejects incomplete forests", {
  skip_if_not_installed("ranger")
  x <- data.frame(a = sin(1:30), b = cos(1:30))
  adapter <- sb_ranger_binary_learner_adapter(num_trees = 5L)
  native <- ranger::ranger
  local_mocked_bindings(ranger = function(...) {
    fitted <- native(...)
    fitted$forest <- NULL
    fitted
  }, .package = "ranger")
  expect_s3_class(tryCatch(adapter$fit(x, rep(0:1, 15L)), error = identity), "stormbreaker_incomplete_ranger_fit")
})
