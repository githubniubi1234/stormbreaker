test_that("mixOmics PLS-DA learners declare two distinct fixed methods", {
  expected_dependencies <- c(
    plsda_binary = "mixOmics",
    splsda_binary = "mixOmics"
  )

  expect_setequal(
    sb_mixomics_da_learner_ids(),
    names(expected_dependencies)
  )
  expect_identical(
    sb_mixomics_da_params("plsda_binary"),
    list(
      method = "PLS-DA",
      requested_components = 2L,
      component_rule = "min_requested_features_analysis_rows_minus_one",
      scale = TRUE,
      tolerance = 1e-6,
      maximum_iterations = 100L,
      near_zero_variance_filter = FALSE,
      logratio = "none",
      keep_x_rule = "all_features",
      probability_mapping = "binary_softmax_of_native_dummy_scores"
    )
  )
  expect_identical(
    sb_mixomics_da_params("splsda_binary"),
    list(
      method = "sPLS-DA",
      requested_components = 2L,
      component_rule = "min_requested_features_analysis_rows_minus_one",
      scale = TRUE,
      tolerance = 1e-6,
      maximum_iterations = 100L,
      near_zero_variance_filter = FALSE,
      logratio = "none",
      keep_x_rule = "min_10_and_p_per_component",
      probability_mapping = "binary_softmax_of_native_dummy_scores"
    )
  )

  registry <- sb_pipeline_registry()
  for (learner in names(expected_dependencies)) {
    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_setequal(
      rows$selector,
      expected_fixed_learner_selectors("binary")
    )
    expect_true(all(grepl("mixOmics", rows$optional_dependency, fixed = TRUE)))
    expect_length(unique(rows$task_type), 1L)
    expect_identical(unique(rows$task_type), "binary")

    pipelines <- list(
      sb_pipeline(learner = learner),
      sb_pipeline(
        selector = "univ_binary_auc_top_k",
        learner = learner,
        selector_n_features = 3L
      ),
      sb_pipeline(
        selector = "abess_binary_gic_nonzero",
        learner = learner
      ),
      sb_pipeline(
        selector = "boruta_binary_roughfix",
        learner = learner
      ),
      sb_pipeline(
        selector = "varselrf_binary_oob_1se",
        learner = learner
      ),
      sb_pipeline(
        selector = "fixed_features",
        learner = learner,
        selector_features = c("gene_a", "gene_b", "gene_c"),
        feature_set_id = "mixomics_fixed"
      )
    )
    expect_length(pipelines, 6L)
    expect_true(all(vapply(
      pipelines,
      inherits,
      logical(1),
      what = "stormbreaker_pipeline"
    )))
  }
})

test_that("mixOmics adapters reproduce native dummy scores and softmax", {
  skip_if_not_installed("mixOmics", minimum_version = "6.32.0")

  data <- make_sb_high_dim_binary_fixture(n = 90L, p = 12L, seed = 8641L)
  x <- data[, -(1:2), drop = FALSE]
  names(x) <- c("gene A", "gene-B", paste0("gene_", 3:12))
  rownames(x) <- data$patient_id
  assessment <- x[1:9, , drop = FALSE]
  grouping <- factor(data$outcome, levels = c(0L, 1L))

  for (learner in sb_mixomics_da_learner_ids()) {
    adapter <- sb_mixomics_da_learner_adapter(learner)
    model <- adapter$fit(x, data$outcome)
    native <- if (identical(learner, "plsda_binary")) {
      mixOmics::plsda(
        data.matrix(x), grouping,
        ncomp = 2L, scale = TRUE, tol = 1e-6, max.iter = 100L,
        near.zero.var = FALSE, logratio = "none", all.outputs = TRUE
      )
    } else {
      mixOmics::splsda(
        data.matrix(x), grouping,
        ncomp = 2L, keepX = c(10L, 10L), scale = TRUE,
        tol = 1e-6, max.iter = 100L, near.zero.var = FALSE,
        logratio = "none", all.outputs = TRUE
      )
    }
    native_scores <- stats::predict(
      native,
      data.matrix(assessment)
    )$predict[, , "dim2"]
    expected <- stats::plogis(
      native_scores[, "1"] - native_scores[, "0"]
    )

    expect_equal(
      adapter$predict(model, assessment),
      unname(expected),
      tolerance = 1e-12
    )
    expect_identical(
      adapter$predict(model, assessment) > 0.5,
      unname(native_scores[, "1"] > native_scores[, "0"])
    )
    expect_identical(model$audit$effective_components, 2L)
    expect_identical(
      model$audit$probability_mapping,
      "binary_softmax_of_native_dummy_scores"
    )
    if (identical(learner, "splsda_binary")) {
      expect_identical(
        vapply(model$audit$component_features, length, integer(1)),
        c(10L, 10L)
      )
    } else {
      expect_identical(model$audit$component_features, list())
    }
  }
})

test_that("mixOmics component count follows available dimensions", {
  skip_if_not_installed("mixOmics", minimum_version = "6.32.0")

  data <- make_sb_binary_fixture(n = 70L, seed = 8642L)
  x <- data.frame("gene A" = data$gene_a, check.names = FALSE)
  rownames(x) <- data$patient_id

  for (learner in sb_mixomics_da_learner_ids()) {
    adapter <- sb_mixomics_da_learner_adapter(learner)
    model <- adapter$fit(x, data$outcome)
    probability <- adapter$predict(model, x[1:5, , drop = FALSE])

    expect_identical(model$audit$effective_components, 1L)
    expect_identical(dim(model$engine$coefficients), c(1L, 2L))
    expect_length(probability, 5L)
    expect_true(all(is.finite(probability)))
  }
})

test_that("mixOmics states replay without analysis rows", {
  skip_if_not_installed("mixOmics", minimum_version = "6.32.0")

  data <- make_sb_high_dim_binary_fixture(n = 84L, p = 4L, seed = 8643L)
  x <- data[, -(1:2), drop = FALSE]
  names(x) <- c("gene A", "gene-B", "gene_c", "gene_d")
  rownames(x) <- data$patient_id
  assessment <- x[1:7, , drop = FALSE]

  for (learner in sb_mixomics_da_learner_ids()) {
    adapter <- sb_mixomics_da_learner_adapter(learner)
    model <- adapter$fit(x, data$outcome)
    expected <- adapter$predict(model, assessment)
    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    replay <- adapter$predict_state(state, assessment)$probability

    expect_equal(replay, expected, tolerance = 0)
    expect_identical(
      names(state$engine),
      c("coefficients", "x_center", "x_scale", "y_center", "y_scale")
    )
    expect_identical(
      dim(state$engine$coefficients),
      c(ncol(x), 2L)
    )
    expect_length(state$engine$x_center, ncol(x))
    expect_length(state$engine$x_scale, ncol(x))
    expect_length(state$engine$y_center, 2L)
    expect_length(state$engine$y_scale, 2L)
  }
})

test_that("mixOmics learners execute through CV and final prediction", {
  skip_if_not_installed("mixOmics", minimum_version = "6.32.0")

  data <- make_sb_high_dim_binary_fixture(n = 72L, p = 4L, seed = 8644L)
  task <- sb_binary_task(
    data = data,
    id_col = "patient_id",
    outcome_col = "outcome",
    predictors = c("gene_001", "gene_002", "gene_003", "gene_004")
  )

  for (learner in sb_mixomics_da_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 3L
    )
    cv <- sb_run_binary_cv(task, pipeline, folds = 3L, seed = 8645L)
    final <- sb_fit_final(task, pipeline)
    predicted <- sb_predict_final(
      final,
      data[1:6, c(
        "patient_id", "gene_001", "gene_002", "gene_003", "gene_004"
      )],
      id_col = "patient_id"
    )

    expect_identical(cv$status, "success")
    expect_identical(cv$complete_oof, TRUE)
    expect_identical(final$status, "success")
    expect_identical(predicted$status, "success")
    expect_length(predicted$predictions$probability, 6L)
    expect_true(all(is.finite(predicted$predictions$probability)))
    expect_true(all(
      predicted$predictions$probability >= 0 &
        predicted$predictions$probability <= 1
    ))
  }
})
