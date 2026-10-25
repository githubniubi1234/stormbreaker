test_that("batch-nine binary learners declare distinct native methods", {
  expected <- c(
    earth_binary = "earth",
    fda_binary = "mda",
    mda_binary = "mda",
    liblinear_l2_logistic_binary = "LiblineaR",
    evtree_binary = "evtree",
    sda_binary = "sda",
    rotation_forest_binary = "rotationForest",
    rrf_binary = "RRF"
  )

  expect_true(all(names(expected) %in% sb_classical_binary_learner_ids()))
  for (learner in names(expected)) {
    dependency <- sb_classical_binary_dependency(learner)
    expect_identical(dependency$package, unname(expected[[learner]]))
    expect_identical(
      dependency$available,
      requireNamespace(dependency$package, quietly = TRUE)
    )
  }

  expect_identical(
    sb_classical_binary_params("earth_binary"),
    list(
      degree = 1L,
      nk_rule = "max_21_2p_plus1",
      nprune = 15L,
      pmethod = "backward",
      glm_family = "binomial_logit",
      nfold = 0L,
      keep_xy = FALSE
    )
  )
  expect_identical(
    sb_classical_binary_params("fda_binary"),
    list(
      method = "polyreg",
      dimension = 1L,
      degree = 1L,
      monomial = FALSE,
      keep_fitted = FALSE
    )
  )
  expect_identical(
    sb_classical_binary_params("mda_binary"),
    list(
      subclasses_per_class = 2L,
      dimension = 3L,
      iterations = 10L,
      method = "polyreg",
      degree = 1L,
      keep_fitted = FALSE,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params(
      "liblinear_l2_logistic_binary"
    ),
    list(
      solver_type = 0L,
      cost = 1,
      tolerance = 0.01,
      bias = 1,
      standardize = TRUE,
      probability = TRUE
    )
  )
  expect_identical(
    sb_classical_binary_params("evtree_binary"),
    list(
      minbucket = 7L,
      minsplit = 20L,
      maxdepth = 5L,
      iterations = 1000L,
      population_size = 50L,
      complexity_weight = 1,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("sda_binary"),
    list(
      covariance_shrinkage = "analytic",
      variance_shrinkage = "analytic",
      class_frequency_shrinkage = "analytic",
      diagonal = FALSE
    )
  )
  expect_identical(
    sb_classical_binary_params("rotation_forest_binary"),
    list(
      num_trees = 10L,
      subset_rule = "round_p_over_3_bounded_by_1_and_floor_p_over_2",
      cp = 0.01,
      minsplit = 10L,
      maxdepth = 10L,
      xval = 0L,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("rrf_binary"),
    list(
      num_trees = 500L,
      mtry_rule = "floor_sqrt_p",
      min_node_size = 1L,
      replace = TRUE,
      regularization_mode = 1L,
      regularization_coefficient = 0.5,
      seed = 1L
    )
  )
})

test_that("batch-nine binary predictions equal native engine predictions", {
  data <- sb_example_binary(n = 84L, seed = 8511L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:8, , drop = FALSE]
  learners <- c(
    "earth_binary",
    "fda_binary",
    "mda_binary",
    "liblinear_l2_logistic_binary",
    "evtree_binary",
    "sda_binary",
    "rotation_forest_binary",
    "rrf_binary"
  )

  learners <- learners[vapply(learners, function(learner) {
    sb_classical_binary_dependency(learner)$available
  }, logical(1))]
  if (!length(learners)) skip("No optional engines for this batch are installed")

  for (learner in learners) {
    params <- sb_classical_binary_params(learner)
    if (identical(learner, "evtree_binary")) {
      params$iterations <- 200L
      params$population_size <- 20L
    }
    if (identical(learner, "rotation_forest_binary")) {
      params$num_trees <- 3L
    }
    if (identical(learner, "rrf_binary")) {
      params$num_trees <- 40L
    }
    adapter <- sb_classical_binary_learner_adapter(learner, params)
    model <- adapter$fit(x, data$outcome)
    observed <- adapter$predict(model, assessment)

    native <- switch(learner,
      earth_binary = {
        expect_s3_class(model$engine, "earth")
        stats::predict(model$engine, assessment, type = "response")
      },
      fda_binary = {
        expect_s3_class(model$engine, "fda")
        stats::predict(model$engine, assessment, type = "posterior")[, "1"]
      },
      mda_binary = {
        expect_s3_class(model$engine, "mda")
        stats::predict(model$engine, assessment, type = "posterior")[, "1"]
      },
      liblinear_l2_logistic_binary = {
        expect_s3_class(model$engine$model, "LiblineaR")
        matrix_x <- as.matrix(assessment)
        standardized <- sweep(
          sweep(matrix_x, 2L, model$engine$center, "-"),
          2L,
          model$engine$scale,
          "/"
        )
        stats::predict(
          model$engine$model,
          standardized,
          proba = TRUE,
          decisionValues = FALSE
        )$probabilities[, "1"]
      },
      evtree_binary = {
        expect_s3_class(model$engine, "constparty")
        stats::predict(model$engine, assessment, type = "prob")[, "1"]
      },
      sda_binary = {
        expect_s3_class(model$engine, "sda")
        sda::predict.sda(
          model$engine,
          as.matrix(assessment),
          verbose = FALSE
        )$posterior[, "1"]
      },
      rotation_forest_binary = {
        expect_s3_class(model$engine, "rotationForest")
        stats::predict(model$engine, assessment)
      },
      rrf_binary = {
        expect_s3_class(model$engine, "RRF")
        stats::predict(
          model$engine,
          as.matrix(assessment),
          type = "prob"
        )[, "1"]
      }
    )
    state <- adapter$extract_state(model)
    replay <- adapter$predict_state(
      unserialize(serialize(state, NULL)),
      assessment
    )$probability

    expect_equal(observed, unname(as.numeric(native)), tolerance = 1e-13)
    expect_equal(replay, observed, tolerance = 0)
  }
})

test_that("batch-nine prediction states omit avoidable fitted-row payloads", {
  for (package in c("earth", "mda", "rotationForest", "RRF", "evtree")) {
    skip_if_not_installed(package)
  }
  data <- sb_example_binary(n = 72L, seed = 8512L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  state_engine <- function(learner, params = NULL) {
    if (is.null(params)) {
      params <- sb_classical_binary_params(learner)
    }
    adapter <- sb_classical_binary_learner_adapter(learner, params)
    model <- adapter$fit(x, data$outcome)
    list(
      state = adapter$extract_state(model),
      engine = unserialize(adapter$extract_state(model)$engine_raw)
    )
  }

  earth <- state_engine("earth_binary")
  expect_identical(
    intersect(
      c("bx", "fitted.values", "residuals", "glm.list", "weights", "call"),
      names(earth$engine)
    ),
    character()
  )

  fda <- state_engine("fda_binary")
  expect_identical(
    intersect(c("x", "residuals", "call"), names(fda$engine$fit)),
    character()
  )

  mda <- state_engine("mda_binary")
  expect_false("weights" %in% names(mda$engine))

  params <- sb_classical_binary_params("rotation_forest_binary")
  params$num_trees <- 3L
  rotation <- state_engine("rotation_forest_binary", params)
  expect_true(all(vapply(rotation$engine$models, function(tree) {
    identical(intersect(c("model", "x", "y", "where", "call"), names(tree)),
      character()) && identical(environment(tree$terms), baseenv())
  }, logical(1))))

  params <- sb_classical_binary_params("rrf_binary")
  params$num_trees <- 30L
  rrf <- state_engine("rrf_binary", params)
  expect_identical(
    intersect(
      c("predicted", "votes", "oob.times", "y", "inbag", "call"),
      names(rrf$engine)
    ),
    character()
  )

  evtree <- state_engine("evtree_binary")
  expect_identical(
    evtree$state$training_data_retention,
    "required_by_evtree_terminal_class_probability_state"
  )
})

test_that("MDA fixes local initialization without changing caller RNG", {
  skip_if_not_installed("mda")
  data <- sb_example_binary(n = 84L, seed = 8513L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:8, , drop = FALSE]
  adapter <- sb_classical_binary_learner_adapter("mda_binary")

  set.seed(19L)
  caller_state <- .Random.seed
  first <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, caller_state)

  set.seed(991L)
  second <- adapter$fit(x, data$outcome)
  expect_equal(
    adapter$predict(first, assessment),
    adapter$predict(second, assessment),
    tolerance = 0
  )
  expect_identical(first$params$seed, 1L)
})
