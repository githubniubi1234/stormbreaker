test_that("LightGBM DART, GOSS, and RF are distinct training regimes", {
  learners <- c(
    "lightgbm_binary_dart",
    "lightgbm_binary_goss",
    "lightgbm_binary_rf"
  )
  expect_setequal(
    sb_lightgbm_binary_learner_ids(),
    c("lightgbm_binary", learners)
  )

  dart <- sb_lightgbm_binary_params(
    learner = "lightgbm_binary_dart"
  )
  expect_identical(dart$boosting, "dart")
  expect_identical(dart$drop_rate, 0.1)
  expect_identical(dart$max_drop, 50L)
  expect_identical(dart$skip_drop, 0)
  expect_identical(dart$uniform_drop, FALSE)
  expect_identical(dart$xgboost_dart_mode, FALSE)

  goss <- sb_lightgbm_binary_params(
    learner = "lightgbm_binary_goss"
  )
  expect_identical(goss$boosting, "gbdt")
  expect_identical(goss$data_sample_strategy, "goss")
  expect_identical(goss$top_rate, 0.2)
  expect_identical(goss$other_rate, 0.1)
  expect_identical(goss$bagging_fraction, 1)
  expect_identical(goss$bagging_freq, 0L)

  forest <- sb_lightgbm_binary_params(
    learner = "lightgbm_binary_rf"
  )
  expect_identical(forest$boosting, "rf")
  expect_identical(forest$feature_fraction, 0.8)
  expect_identical(forest$bagging_fraction, 0.8)
  expect_identical(forest$bagging_freq, 1L)

  registry <- sb_pipeline_registry()
  expected_selectors <- expected_fixed_learner_selectors("binary")
  for (learner in learners) {
    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_identical(unique(rows$task_type), "binary")
    expect_setequal(rows$selector, expected_selectors)
    expect_true(all(grepl(
      "lightgbm", rows$optional_dependency,
      fixed = TRUE
    )))
    expect_identical(unique(rows$tuning_semantics), "none")

    pipeline <- sb_pipeline(
      learner = learner,
      lightgbm_nrounds = 12L,
      lightgbm_learning_rate = 0.05,
      lightgbm_num_leaves = 7L,
      lightgbm_min_data_in_leaf = 3L,
      lightgbm_seed = 29L
    )
    expect_identical(pipeline$learner, learner)
    expect_identical(pipeline$learner_params$nrounds, 12L)
    expect_identical(pipeline$learner_params$seed, 29L)
    expect_match(pipeline$id, paste0("none__", learner), fixed = TRUE)
  }
})

test_that("new LightGBM regimes preserve native probabilities and state", {
  skip_if_not_installed("lightgbm", minimum_version = "4.0.0")

  data <- make_sb_binary_fixture(n = 72L, seed = 8911L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("gene A", "gene-B", "gene_c")
  assessment <- x[c(5L, 21L, 46L, 67L), , drop = FALSE]
  learners <- c(
    "lightgbm_binary_dart",
    "lightgbm_binary_goss",
    "lightgbm_binary_rf"
  )
  withr::local_seed(8912L)

  for (learner in learners) {
    adapter <- sb_lightgbm_binary_learner_adapter(
      nrounds = 12L,
      learning_rate = 0.1,
      num_leaves = 7L,
      min_data_in_leaf = 3L,
      seed = 1L,
      learner = learner
    )
    rng_before <- .Random.seed
    model <- adapter$fit(x, data$outcome)
    expect_identical(.Random.seed, rng_before)
    repeated <- adapter$fit(x, data$outcome)
    expect_identical(.Random.seed, rng_before)
    expect_identical(
      model$booster$save_model_to_string(num_iteration = 12L),
      repeated$booster$save_model_to_string(num_iteration = 12L)
    )

    observed <- adapter$predict(model, assessment)
    native <- as.numeric(stats::predict(
      model$booster,
      as.matrix(assessment),
      type = "response",
      num_iteration = 12L,
      params = list(num_threads = 1L)
    ))
    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    replay <- adapter$predict_state(state, assessment)$probability

    expect_equal(observed, native, tolerance = 0)
    expect_equal(replay, observed, tolerance = 0)
    expect_identical(state$learner, learner)
    expect_type(state$model_string, "character")
    expect_identical("booster" %in% names(state), FALSE)
    expect_true(all(observed >= 0 & observed <= 1))
  }
})

test_that("new LightGBM regimes execute through fold-local CV and final fit", {
  skip_if_not_installed("lightgbm", minimum_version = "4.0.0")

  data <- make_sb_binary_fixture(n = 72L, seed = 8921L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  learners <- c(
    "lightgbm_binary_dart",
    "lightgbm_binary_goss",
    "lightgbm_binary_rf"
  )

  for (learner in learners) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 2L,
      lightgbm_nrounds = 12L,
      lightgbm_num_leaves = 7L,
      lightgbm_min_data_in_leaf = 3L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 8922L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      data[1:6, c("patient_id", "gene_c", "gene_a", "gene_b")],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(sum(comparison$runs$status == "success"), 3L)
    expect_identical(
      comparison$provenance$lightgbm_version,
      as.character(utils::packageVersion("lightgbm"))
    )
    expect_identical(fit$status, "success")
    expect_identical(fit$learner_state$learner, learner)
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$probability, 6L)
  }
})
