make_liblinear_formulation_fixture <- function(n = 72L,
                                               p = 5L,
                                               seed = 9311L) {
  withr::with_seed(seed, {
    outcome <- sample(rep(c(0L, 1L), length.out = n))
    x <- matrix(stats::rnorm(n * p), nrow = n)
    x[, 1L] <- x[, 1L] + 1.5 * outcome
    if (p >= 2L) {
      x[, 2L] <- x[, 2L] - 1.1 * outcome
    }
    names <- c(
      "gene one",
      if (p >= 2L) "2gene" else NULL,
      if (p >= 3L) paste0("gene-", seq_len(p - 2L)) else NULL
    )
    colnames(x) <- names
    data.frame(
      patient_id = sprintf("LL%03d", seq_len(n)),
      outcome = outcome,
      x,
      check.names = FALSE
    )
  })
}

test_that("LIBLINEAR formulations declare four distinct new objectives", {
  expected <- c(
    "liblinear_l2_logistic_binary",
    "liblinear_l1_logistic_binary",
    "liblinear_l2_l2loss_svc_binary",
    "liblinear_l2_l1loss_svc_binary",
    "liblinear_l1_l2loss_svc_binary"
  )
  expect_identical(sb_liblinear_binary_learner_ids(), expected)

  expected_types <- c(0L, 6L, 1L, 3L, 5L)
  expected_tolerances <- c(0.01, 0.01, 0.1, 0.1, 0.01)
  expected_probabilities <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
  for (i in seq_along(expected)) {
    params <- sb_classical_binary_params(expected[[i]])
    expect_identical(params$solver_type, expected_types[[i]])
    expect_identical(params$tolerance, expected_tolerances[[i]])
    expect_identical(params$probability, expected_probabilities[[i]])
    expect_identical(params$cost, 1)
    expect_identical(params$bias, 1)
    expect_identical(params$standardize, TRUE)
  }
  for (learner in expected[-1L]) {
    expect_identical(
      sb_classical_binary_params(learner)$seed,
      1L
    )
  }

  registry <- sb_pipeline_registry()
  expected_selectors <- expected_fixed_learner_selectors("binary")
  rows <- registry[
    registry$learner %in% expected[-1L], ,
    drop = FALSE
  ]
  expect_identical(
    nrow(rows),
    length(expected[-1L]) * length(expected_selectors)
  )
  expect_identical(
    sort(unique(rows$selector)),
    sort(expected_selectors)
  )
  expect_true(all(grepl("LiblineaR", rows$optional_dependency, fixed = TRUE)))
})

test_that("LIBLINEAR formulations reproduce native probability or margin", {
  skip_if_not_installed("LiblineaR", minimum_version = "2.10.0")

  data <- make_liblinear_formulation_fixture()
  features <- setdiff(names(data), c("patient_id", "outcome"))
  x <- data[features]
  assessment <- x[c(2L, 19L, 51L), , drop = FALSE]

  for (learner in sb_liblinear_binary_learner_ids()[-1L]) {
    adapter <- sb_classical_binary_learner_adapter(learner)

    set.seed(9312L)
    caller_rng <- .Random.seed
    model <- adapter$fit(x, data$outcome)
    expect_identical(.Random.seed, caller_rng)

    set.seed(9313L)
    second_caller_rng <- .Random.seed
    repeated <- adapter$fit(x, data$outcome)
    expect_identical(.Random.seed, second_caller_rng)
    expect_identical(
      serialize(model$engine, NULL),
      serialize(repeated$engine, NULL)
    )

    standardized <- sweep(
      sweep(as.matrix(assessment), 2L, model$engine$center, "-"),
      2L,
      model$engine$scale,
      "/"
    )
    native <- if (isTRUE(adapter$params$probability)) {
      stats::predict(
        model$engine$model,
        standardized,
        proba = TRUE,
        decisionValues = FALSE
      )$probabilities[, "1"]
    } else {
      prediction <- stats::predict(
        model$engine$model,
        standardized,
        proba = FALSE,
        decisionValues = TRUE
      )
      decision <- prediction$decisionValues
      mapped <- stats::plogis(decision[, "1"] - decision[, "0"])
      expect_identical(
        as.integer(mapped >= 0.5),
        as.integer(prediction$predictions)
      )
      mapped
    }
    observed <- adapter$predict(model, assessment)
    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    replay <- adapter$predict_state(state, assessment)$probability

    expect_identical(model$engine$model$Type, adapter$params$solver_type)
    expect_equal(observed, unname(native), tolerance = 0)
    expect_equal(replay, observed, tolerance = 0)
    expect_type(state$engine_raw, "raw")
    expect_null(state$training_data_retention)
  }
})

test_that("LIBLINEAR formulations support p greater than n and one feature", {
  skip_if_not_installed("LiblineaR", minimum_version = "2.10.0")

  cases <- list(
    high_dimensional = make_liblinear_formulation_fixture(
      n = 40L,
      p = 48L,
      seed = 9314L
    ),
    one_feature = make_liblinear_formulation_fixture(
      n = 48L,
      p = 1L,
      seed = 9315L
    )
  )
  for (data in cases) {
    features <- setdiff(names(data), c("patient_id", "outcome"))
    x <- data[features]
    for (learner in sb_liblinear_binary_learner_ids()[-1L]) {
      adapter <- sb_classical_binary_learner_adapter(learner)
      model <- adapter$fit(x, data$outcome)
      probability <- adapter$predict(model, x[1:5, , drop = FALSE])

      expect_length(probability, 5L)
      expect_true(all(is.finite(probability)))
      expect_true(all(probability >= 0 & probability <= 1))
    }
  }
})

test_that("LIBLINEAR formulations execute through CV and final prediction", {
  skip_if_not_installed("LiblineaR", minimum_version = "2.10.0")

  data <- make_liblinear_formulation_fixture(n = 72L, seed = 9316L)
  features <- setdiff(names(data), c("patient_id", "outcome"))
  task <- sb_binary_task(data, "patient_id", "outcome", features)

  for (learner in sb_liblinear_binary_learner_ids()[-1L]) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 3L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 9317L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      data[1:6, c("patient_id", rev(features)), drop = FALSE],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_identical(fit$learner_state$learner, learner)
  }
})
