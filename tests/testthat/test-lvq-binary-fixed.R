test_that("four LVQ update algorithms are admitted as binary learners", {
  learners <- c(
    "lvq1_binary",
    "olvq1_binary",
    "lvq21_binary",
    "lvq3_binary"
  )
  expect_setequal(sb_lvq_binary_learner_ids(), learners)

  expected_methods <- c("lvq1", "olvq1", "lvq2.1", "lvq3")
  expected_iterations <- c(100L, 40L, 100L, 100L)
  expected_alpha <- c(0.03, 0.3, 0.03, 0.03)
  expected_selectors <- expected_fixed_learner_selectors("binary")
  registry <- sb_pipeline_registry()

  for (index in seq_along(learners)) {
    learner <- learners[[index]]
    params <- sb_classical_binary_params(learner)
    expect_identical(params$method, expected_methods[[index]])
    expect_identical(
      params$iterations_per_prototype,
      expected_iterations[[index]]
    )
    expect_identical(params$alpha, expected_alpha[[index]])
    expect_identical(params$codebook_size_rule, "max_two_native_default")
    expect_identical(params$prototype_allocation, "equal_prior_target")
    expect_identical(params$initialization_k, 5L)
    expect_identical(params$standardize, TRUE)
    expect_identical(params$seed, 1L)
    expect_identical(
      params$score_mapping,
      "logistic_nearest_prototype_squared_distance_margin"
    )
    if (learner %in% c("lvq21_binary", "lvq3_binary")) {
      expect_identical(params$window, 0.3)
    } else {
      expect_null(params$window)
    }
    if (identical(learner, "lvq3_binary")) {
      expect_identical(params$epsilon, 0.1)
    } else {
      expect_null(params$epsilon)
    }

    dependency <- sb_classical_binary_dependency(learner)
    expect_identical(dependency$package, "class")
    expect_identical(dependency$available, TRUE)

    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_identical(unique(rows$task_type), "binary")
    expect_setequal(rows$selector, expected_selectors)
    expect_identical(unique(rows$tuning_semantics), "none")
    expect_identical(
      all(grepl("class", rows$optional_dependency, fixed = TRUE)),
      TRUE
    )
  }
})

test_that("portable LVQ states preserve the native nearest-prototype boundary", {
  skip_if_not_installed("class", minimum_version = "7.3")

  data <- make_sb_binary_fixture(n = 84L, seed = 11101L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("gene one", "gene-two", "gene_c")
  assessment <- x[c(3L, 19L, 47L, 72L), , drop = FALSE]
  learners <- sb_lvq_binary_learner_ids()

  withr::local_seed(11102L)
  rng_before <- .Random.seed
  for (learner in learners) {
    adapter <- sb_classical_binary_learner_adapter(learner)
    model <- adapter$fit(x, data$outcome)
    expect_identical(.Random.seed, rng_before)
    repeated <- adapter$fit(x, data$outcome)
    expect_identical(.Random.seed, rng_before)
    expect_equal(repeated$engine$prototypes, model$engine$prototypes)

    params <- sb_classical_binary_params(learner)
    standardized_training <- sweep(
      sweep(as.matrix(x), 2L, model$engine$center, "-"),
      2L,
      model$engine$scale,
      "/"
    )
    native <- withr::with_seed(params$seed, {
      size <- max(
        2L,
        as.integer(min(
          round(0.4 * 2 * (2 - 1 + ncol(x) / 2), 0),
          nrow(x)
        ))
      )
      initial <- class::lvqinit(
        standardized_training,
        factor(data$outcome, levels = c(0L, 1L)),
        size = size,
        prior = c(0.5, 0.5),
        k = params$initialization_k
      )
      iterations <- params$iterations_per_prototype * nrow(initial$x)
      switch(params$method,
        lvq1 = class::lvq1(
          standardized_training,
          factor(data$outcome, levels = c(0L, 1L)),
          initial,
          niter = iterations,
          alpha = params$alpha
        ),
        olvq1 = class::olvq1(
          standardized_training,
          factor(data$outcome, levels = c(0L, 1L)),
          initial,
          niter = iterations,
          alpha = params$alpha
        ),
        `lvq2.1` = class::lvq2(
          standardized_training,
          factor(data$outcome, levels = c(0L, 1L)),
          initial,
          niter = iterations,
          alpha = params$alpha,
          win = params$window
        ),
        lvq3 = class::lvq3(
          standardized_training,
          factor(data$outcome, levels = c(0L, 1L)),
          initial,
          niter = iterations,
          alpha = params$alpha,
          win = params$window,
          epsilon = params$epsilon
        )
      )
    })
    expect_equal(
      model$engine$prototypes,
      unname(as.matrix(native$x)),
      tolerance = 0
    )
    expect_identical(
      model$engine$prototype_class,
      unname(as.character(native$cl))
    )

    observed <- adapter$predict(model, assessment)
    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    replay <- adapter$predict_state(state, assessment)$probability
    engine <- unserialize(state$engine_raw)
    standardized <- sweep(
      sweep(as.matrix(assessment), 2L, engine$center, "-"),
      2L,
      engine$scale,
      "/"
    )
    codebook <- list(
      x = engine$prototypes,
      cl = factor(engine$prototype_class, levels = c("0", "1"))
    )
    native_class <- class::lvqtest(codebook, standardized)
    squared_distance <- outer(
      rowSums(standardized^2),
      rowSums(engine$prototypes^2),
      "+"
    ) - 2 * tcrossprod(standardized, engine$prototypes)
    squared_distance <- pmax(squared_distance, 0)
    distance_zero <- apply(
      squared_distance[, engine$prototype_class == "0", drop = FALSE],
      1L,
      min
    )
    distance_one <- apply(
      squared_distance[, engine$prototype_class == "1", drop = FALSE],
      1L,
      min
    )
    expected <- stats::plogis(distance_zero - distance_one)

    expect_equal(observed, unname(expected), tolerance = 1e-12)
    expect_equal(replay, observed, tolerance = 0)
    expect_identical(
      as.integer(observed > 0.5),
      as.integer(as.character(native_class))
    )
    expect_setequal(
      names(engine),
      c(
        "center", "scale", "prototypes", "prototype_class",
        "codebook_size", "iterations"
      )
    )
    expect_null(rownames(engine$prototypes))
    expect_identical(
      any(c("outcome", "training", "native_fit") %in% names(engine)),
      FALSE
    )
    expect_identical("training_data_retention" %in% names(state), FALSE)
  }
})

test_that("the LVQ initialization rule supports a one-feature analysis fit", {
  skip_if_not_installed("class", minimum_version = "7.3")

  outcome <- rep(c(0L, 1L), each = 30L)
  x <- data.frame(signal = c(
    seq(-3, -1, length.out = 30L),
    seq(1, 3, length.out = 30L)
  ))
  for (learner in sb_lvq_binary_learner_ids()) {
    adapter <- sb_classical_binary_learner_adapter(learner)
    model <- adapter$fit(x, outcome)
    probability <- adapter$predict(model, x[c(5L, 55L), , drop = FALSE])

    expect_identical(model$engine$codebook_size, 2L)
    expect_setequal(model$engine$prototype_class, c("0", "1"))
    expect_length(probability, 2L)
    expect_identical(all(is.finite(probability)), TRUE)
  }
})

test_that("LVQ learners execute through fold-local CV and final fit", {
  skip_if_not_installed("class", minimum_version = "7.3")

  data <- make_sb_binary_fixture(n = 72L, seed = 11111L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )

  for (learner in sb_lvq_binary_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 2L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 11112L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      data[1:6, c("patient_id", "gene_c", "gene_a", "gene_b")],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(sum(comparison$runs$status == "success"), 3L)
    expect_identical("class_version" %in% names(comparison$provenance), TRUE)
    expect_identical(fit$status, "success")
    expect_identical(fit$learner_state$learner, learner)
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$probability, 6L)
  }
})
