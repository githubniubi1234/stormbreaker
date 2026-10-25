test_that("PDA and stabilised LDA declare distinct native specifications", {
  expected_dependencies <- c(
    pda_binary = "mda",
    slda_binary = "ipred"
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
    sb_classical_binary_params("pda_binary"),
    list(
      method = "gen.ridge",
      dimension = 1L,
      lambda = 1,
      prior = "empirical",
      keep_fitted = FALSE
    )
  )
  expect_identical(
    sb_classical_binary_params("slda_binary"),
    list(
      score_rule = "native_PC_q",
      q = NULL,
      prior = "empirical"
    )
  )

  registry <- sb_pipeline_registry()
  selectors <- expected_fixed_learner_selectors("binary")
  for (learner in names(expected_dependencies)) {
    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_identical(nrow(rows), length(selectors))
    expect_setequal(rows$selector, selectors)
    expect_identical(unique(rows$task_type), "binary")
    expect_identical(unique(rows$tuning_semantics), "none")
  }
})

test_that("PDA and stabilised LDA reproduce native high-dimensional probabilities", {
  skip_if_not_installed("mda", minimum_version = "0.5.0")
  skip_if_not_installed("ipred", minimum_version = "0.9.0")

  data <- make_sb_high_dim_binary_fixture(n = 64L, p = 90L, seed = 8611L)
  names(data)[names(data) == "gene_001"] <- "gene A"
  names(data)[names(data) == "gene_002"] <- "gene-B"
  predictors <- setdiff(names(data), c("patient_id", "outcome"))
  x <- data[predictors]
  assessment <- x[1:9, , drop = FALSE]
  grouping <- factor(data$outcome, levels = c(0L, 1L))
  priors <- as.numeric(table(grouping) / length(grouping))

  pda_adapter <- sb_classical_binary_learner_adapter(
    "pda_binary"
  )
  pda <- pda_adapter$fit(x, data$outcome)
  formula_data <- sb_classical_binary_formula_data(x, data$outcome)
  native_pda <- mda::fda(
    formula_data$formula,
    data = formula_data$data,
    dimension = 1L,
    method = mda::gen.ridge,
    keep.fitted = FALSE,
    lambda = 1
  )
  expected_pda <- stats::predict(
    native_pda,
    newdata = assessment,
    type = "posterior"
  )[, "1"]

  expect_s3_class(pda$engine, "fda")
  expect_s3_class(pda$engine$fit, "gen.ridge")
  expect_equal(as.numeric(pda$engine$prior), priors, tolerance = 0)
  expect_identical(pda$engine$fit$lambda, 1)
  expect_gt(ncol(x), nrow(x))
  expect_equal(
    pda_adapter$predict(pda, assessment),
    unname(expected_pda),
    tolerance = 1e-13
  )

  slda_adapter <- sb_classical_binary_learner_adapter(
    "slda_binary"
  )
  slda <- slda_adapter$fit(x, data$outcome)
  native_slda <- ipred::slda(
    grouping,
    X = as.matrix(x),
    q = NULL,
    prior = priors
  )
  projected <- as.matrix(assessment) %*% native_slda$scores
  colnames(projected) <- paste0("V", seq_len(ncol(projected)))
  expected_slda <- stats::predict(
    native_slda$mylda,
    newdata = as.data.frame(projected)
  )$posterior[, "1"]

  expect_s3_class(slda$engine, "slda")
  expect_s3_class(slda$engine$mylda, "lda")
  expect_identical(
    ncol(slda$engine$scores),
    ncol(native_slda$scores)
  )
  expect_equal(unname(slda$engine$mylda$prior), priors, tolerance = 0)
  expect_equal(
    slda_adapter$predict(slda, assessment),
    unname(expected_slda),
    tolerance = 1e-13
  )
})

test_that("PDA and stabilised LDA states replay without fitted patient rows", {
  skip_if_not_installed("mda", minimum_version = "0.5.0")
  skip_if_not_installed("ipred", minimum_version = "0.9.0")

  data <- make_sb_high_dim_binary_fixture(n = 64L, p = 90L, seed = 8612L)
  names(data)[names(data) == "gene_001"] <- "gene A"
  names(data)[names(data) == "gene_002"] <- "gene-B"
  predictors <- setdiff(names(data), c("patient_id", "outcome"))
  x <- data[predictors]
  assessment <- x[1:7, , drop = FALSE]

  for (learner in c("pda_binary", "slda_binary")) {
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
    if (identical(learner, "pda_binary")) {
      expect_identical(
        intersect(c("call", "confusion"), names(engine)),
        character()
      )
      expect_identical(
        intersect(
          c("fitted.values", "x", "residuals", "call"),
          names(engine$fit)
        ),
        character()
      )
      expect_identical(environment(engine$terms), baseenv())
    } else {
      expect_identical(names(engine), c("scores", "mylda"))
      expect_identical("call" %in% names(engine$mylda), FALSE)
    }
  }
})

test_that("PDA and stabilised LDA execute through CV and final prediction", {
  skip_if_not_installed("mda", minimum_version = "0.5.0")
  skip_if_not_installed("ipred", minimum_version = "0.9.0")

  data <- make_sb_high_dim_binary_fixture(n = 72L, p = 12L, seed = 8613L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    predictors = sprintf("gene_%03d", seq_len(12L))
  )

  for (learner in c("pda_binary", "slda_binary")) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 4L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 8614L)
    final <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      final,
      data[1:8, c("patient_id", sprintf("gene_%03d", seq_len(12L)))],
      id_col = "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(final$status, "success")
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$probability, 8L)
  }
})

test_that("PDA no-discrimination fits fail before prediction", {
  skip_if_not_installed("mda", minimum_version = "0.5.0")

  data <- data.frame(
    patient_id = sprintf("PD%02d", seq_len(30L)),
    outcome = rep(c(0L, 1L), each = 15L),
    constant_gene = 0,
    check.names = FALSE
  )
  adapter <- sb_classical_binary_learner_adapter("pda_binary")
  error <- capture_sb_error(suppressWarnings(adapter$fit(
    data["constant_gene"],
    data$outcome
  )))
  task <- sb_binary_task(
    data,
    id_col = "patient_id",
    outcome_col = "outcome",
    predictors = "constant_gene"
  )
  pipeline <- sb_pipeline(learner = "pda_binary")
  comparison <- suppressWarnings(sb_run_cv(
    task,
    pipeline,
    folds = 3L,
    seed = 8615L
  ))
  final <- suppressWarnings(sb_fit_final(task, pipeline))

  expect_s3_class(error, "stormbreaker_pda_no_discrimination")
  expect_s3_class(error$parent, "warning")
  expect_identical(comparison$status, "failed")
  expect_identical(comparison$complete_oof, FALSE)
  expect_identical(
    unique(comparison$runs$error_class),
    "stormbreaker_pda_no_discrimination"
  )
  expect_identical(final$status, "failed")
  expect_identical(final$error_stage, "learner_fit")
  expect_identical(
    final$error_class,
    "stormbreaker_pda_no_discrimination"
  )
})
