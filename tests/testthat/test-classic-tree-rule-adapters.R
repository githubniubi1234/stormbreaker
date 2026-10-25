test_that("classic tree and C5.0 rules declare distinct native specifications", {
  expected_dependencies <- c(
    tree_binary = "tree",
    c50_rules_binary = "C50"
  )

  expect_setequal(
    intersect(
      names(expected_dependencies),
      sb_classical_binary_learner_ids()
    ),
    names(expected_dependencies)
  )
  for (learner in names(expected_dependencies)) {
    dependency <- sb_classical_binary_dependency(learner)
    expect_identical(
      dependency$package,
      unname(expected_dependencies[[learner]])
    )
    expect_identical(
      dependency$available,
      requireNamespace(dependency$package, quietly = TRUE)
    )
  }

  expect_identical(
    sb_classical_binary_params("tree_binary"),
    list(
      split = "deviance",
      mincut = 5L,
      minsize = 10L,
      mindev = 0.01,
      keep_model = FALSE,
      keep_x = FALSE,
      keep_y = FALSE,
      keep_weights = FALSE
    )
  )
  expect_identical(
    sb_classical_binary_params("c50_rules_binary"),
    list(
      trials = 10L,
      rules = TRUE,
      winnow = FALSE,
      early_stopping = FALSE,
      seed = 1L
    )
  )

  selectors <- expected_fixed_learner_selectors("binary")
  registry <- sb_pipeline_registry()
  for (learner in names(expected_dependencies)) {
    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_identical(nrow(rows), length(selectors))
    expect_setequal(rows$selector, selectors)
    expect_identical(unique(rows$task_type), "binary")
    expect_identical(unique(rows$tuning_semantics), "none")
  }
})

test_that("classic tree and C5.0 rules reproduce native probabilities", {
  skip_if_not_installed("tree", minimum_version = "1.0-40")
  skip_if_not_installed("C50", minimum_version = "0.2.0")

  data <- make_sb_binary_fixture(n = 96L, seed = 8631L)
  names(data)[names(data) == "gene_a"] <- "gene A"
  names(data)[names(data) == "gene_b"] <- "gene-B"
  predictors <- setdiff(names(data), c("patient_id", "outcome"))
  x <- data[predictors]
  grouping <- factor(data$outcome, levels = c(0L, 1L))
  assessment <- x[1:9, , drop = FALSE]

  tree_adapter <- sb_classical_binary_learner_adapter(
    "tree_binary"
  )
  tree_model <- tree_adapter$fit(x, data$outcome)
  formula_data <- sb_tree_binary_formula_data(x, data$outcome)
  native_tree <- tree::tree(
    formula_data$formula,
    data = formula_data$data,
    control = tree::tree.control(
      nobs = nrow(x),
      mincut = 5L,
      minsize = 10L,
      mindev = 0.01
    ),
    split = "deviance",
    model = FALSE,
    x = FALSE,
    y = FALSE,
    wts = FALSE
  )
  tree_assessment <- assessment
  names(tree_assessment) <- sb_tree_binary_internal_names(ncol(assessment))
  expected_tree <- stats::predict(
    native_tree,
    newdata = tree_assessment,
    type = "vector"
  )[, "1"]

  expect_s3_class(tree_model$engine, "tree")
  expect_equal(
    tree_adapter$predict(tree_model, assessment),
    unname(expected_tree),
    tolerance = 0
  )

  rules_adapter <- sb_classical_binary_learner_adapter(
    "c50_rules_binary"
  )
  rules_model <- rules_adapter$fit(x, data$outcome)
  native_rules <- C50::C5.0(
    x = x,
    y = grouping,
    trials = 10L,
    rules = TRUE,
    control = C50::C5.0Control(
      winnow = FALSE,
      seed = 1L,
      earlyStopping = FALSE
    )
  )
  expected_rules <- get(
    "predict.C5.0",
    envir = asNamespace("C50"),
    inherits = FALSE
  )(native_rules, newdata = assessment, type = "prob")[, "1"]

  expect_s3_class(rules_model$engine, "C5.0")
  expect_type(rules_model$engine$rules, "character")
  expect_true(nzchar(rules_model$engine$rules))
  expect_equal(
    rules_adapter$predict(rules_model, assessment),
    unname(expected_rules),
    tolerance = 0
  )
})

test_that("classic tree and C5.0 rule states replay without fitted rows", {
  skip_if_not_installed("tree", minimum_version = "1.0-40")
  skip_if_not_installed("C50", minimum_version = "0.2.0")

  data <- make_sb_binary_fixture(n = 84L, seed = 8632L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  x <- data[predictors]
  assessment <- x[1:7, , drop = FALSE]

  for (learner in c("tree_binary", "c50_rules_binary")) {
    adapter <- sb_classical_binary_learner_adapter(learner)
    model <- adapter$fit(x, data$outcome)
    expected <- adapter$predict(model, assessment)
    state <- adapter$extract_state(model)
    engine <- unserialize(state$engine_raw)
    replay <- adapter$predict_state(
      unserialize(serialize(state, NULL)),
      assessment
    )$probability

    expect_equal(replay, expected, tolerance = 0)
    expect_null(state$training_data_retention)
    if (identical(learner, "tree_binary")) {
      expect_identical(names(engine), c("frame", "terms"))
      expect_identical(environment(engine$terms), baseenv())
    } else {
      expect_identical(
        intersect(c("call", "output", "boostResults", "size"), names(engine)),
        character()
      )
    }
  }
})

test_that("classic tree and C5.0 rules execute through CV and final prediction", {
  skip_if_not_installed("tree", minimum_version = "1.0-40")
  skip_if_not_installed("C50", minimum_version = "0.2.0")

  data <- make_sb_binary_fixture(n = 72L, seed = 8633L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    predictors = c("gene_a", "gene_b", "gene_c")
  )

  for (learner in c("tree_binary", "c50_rules_binary")) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 2L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 8634L)
    final <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      final,
      data[1:8, c("patient_id", "gene_a", "gene_b", "gene_c")],
      id_col = "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(final$status, "success")
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$probability, 8L)
  }
})
