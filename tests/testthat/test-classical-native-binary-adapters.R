test_that("batch-eight binary learners declare distinct native methods", {
  expected <- c(
    ctree_binary = "partykit",
    cforest_binary = "partykit",
    bagging_binary = "ipred",
    gamboost_binary = "mboost",
    blackboost_binary = "mboost",
    loclda_binary = "klaR",
    aorsf_binary = "aorsf"
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
    sb_classical_binary_params("ctree_binary"),
    list(
      teststat = "quad",
      testtype = "Bonferroni",
      mincriterion = 0.95,
      minsplit = 20L,
      minbucket = 7L,
      maxdepth = Inf,
      saveinfo = FALSE
    )
  )
  expect_identical(
    sb_classical_binary_params("cforest_binary"),
    list(
      num_trees = 500L,
      mtry_rule = "ceiling_sqrt_p",
      replace = TRUE,
      sample_fraction = 1,
      teststat = "quad",
      testtype = "Univ",
      mincriterion = 0,
      minsplit = 20L,
      minbucket = 7L,
      maxdepth = Inf,
      saveinfo = FALSE,
      seed = 1L,
      cores = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("bagging_binary"),
    list(
      num_trees = 100L,
      cp = 0,
      minsplit = 20L,
      minbucket = 7L,
      maxdepth = 30L,
      xval = 0L,
      seed = 1L,
      coob = FALSE,
      keep_x = FALSE
    )
  )
  expect_identical(
    sb_classical_binary_params("gamboost_binary"),
    list(
      mstop = 100L,
      nu = 0.1,
      risk = "inbag",
      family = "Binomial_logit",
      baselearner = "bbs",
      dfbase = 4,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("blackboost_binary"),
    list(
      mstop = 100L,
      nu = 0.1,
      risk = "inbag",
      family = "Binomial_logit",
      teststat = "quad",
      testtype = "Teststatistic",
      mincriterion = 0,
      minsplit = 20L,
      minbucket = 7L,
      maxdepth = 2L,
      saveinfo = FALSE,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("loclda_binary"),
    list(
      k_max = 30L,
      k_rule = "dimension_aware",
      weighted_apriori = TRUE,
      prior = "local_empirical"
    )
  )
  expect_identical(
    sb_classical_binary_params("aorsf_binary"),
    list(
      num_trees = 500L,
      num_split = 5L,
      num_retry = 3L,
      mtry_rule = "ceiling_sqrt_p",
      sample_with_replacement = TRUE,
      sample_fraction = 0.632,
      control_method = "glm",
      control_scale_x = TRUE,
      control_max_iter = 1L,
      seed = 1L,
      threads = 1L,
      importance = "none",
      attach_data = FALSE
    )
  )
})

test_that("batch-eight binary predictions equal native engine predictions", {
  data <- sb_example_binary(n = 72L, seed = 8421L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:8, , drop = FALSE]
  learners <- c(
    "ctree_binary",
    "cforest_binary",
    "bagging_binary",
    "gamboost_binary",
    "blackboost_binary",
    "loclda_binary",
    "aorsf_binary"
  )

  learners <- learners[vapply(learners, function(learner) {
    sb_classical_binary_dependency(learner)$available
  }, logical(1))]
  if (!length(learners)) skip("No optional engines for this batch are installed")

  for (learner in learners) {
    params <- sb_classical_binary_params(learner)
    if (identical(learner, "cforest_binary")) {
      params$num_trees <- 30L
    }
    if (identical(learner, "bagging_binary")) {
      params$num_trees <- 20L
    }
    if (learner %in% c("gamboost_binary", "blackboost_binary")) {
      params$mstop <- 20L
    }
    if (identical(learner, "aorsf_binary")) {
      params$num_trees <- 30L
    }
    adapter <- sb_classical_binary_learner_adapter(learner, params)
    model <- adapter$fit(x, data$outcome)
    observed <- adapter$predict(model, assessment)

    native <- switch(learner,
      ctree_binary = {
        expect_s3_class(model$engine, "constparty")
        stats::predict(model$engine, assessment, type = "prob")[, "1"]
      },
      cforest_binary = {
        expect_s3_class(model$engine, "cforest")
        stats::predict(model$engine, assessment, type = "prob")[, "1"]
      },
      bagging_binary = {
        expect_s3_class(model$engine, "classbagg")
        stats::predict(
          model$engine,
          assessment,
          type = "prob",
          aggregation = "average"
        )[, "1"]
      },
      gamboost_binary = {
        expect_s3_class(model$engine, "gamboost")
        stats::predict(model$engine, assessment, type = "response")
      },
      blackboost_binary = {
        expect_s3_class(model$engine, "blackboost")
        stats::predict(model$engine, assessment, type = "response")
      },
      loclda_binary = {
        expect_s3_class(model$engine, "loclda")
        get(
          "predict.loclda",
          envir = asNamespace("klaR"),
          inherits = FALSE
        )(model$engine, assessment)$posterior[, "1"]
      },
      aorsf_binary = {
        expect_s3_class(model$engine, "ObliqueForestClassification")
        stats::predict(
          model$engine,
          new_data = assessment,
          pred_type = "prob",
          n_thread = 1L
        )[, "1"]
      }
    )
    state <- adapter$extract_state(model)
    replay <- adapter$predict_state(
      unserialize(serialize(state, NULL)),
      assessment
    )$probability

    expect_equal(observed, unname(as.numeric(native)), tolerance = 1e-14)
    expect_equal(replay, observed, tolerance = 0)
  }
})

test_that("lazy native binary states disclose required analysis data", {
  for (package in c("partykit", "mboost", "klaR", "aorsf")) {
    skip_if_not_installed(package)
  }
  data <- sb_example_binary(n = 60L, seed = 8422L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  expected <- c(
    ctree_binary =
      "required_by_partykit_class_probability_state",
    cforest_binary =
      "required_by_partykit_class_probability_state",
    gamboost_binary =
      "required_by_mboost_prediction_state",
    blackboost_binary =
      "required_by_mboost_prediction_state",
    loclda_binary =
      "required_by_klar_localized_lda_prediction_state"
  )

  for (learner in names(expected)) {
    params <- sb_classical_binary_params(learner)
    if (identical(learner, "cforest_binary")) {
      params$num_trees <- 20L
    }
    if (learner %in% c("gamboost_binary", "blackboost_binary")) {
      params$mstop <- 15L
    }
    adapter <- sb_classical_binary_learner_adapter(learner, params)
    state <- adapter$extract_state(adapter$fit(x, data$outcome))
    expect_identical(state$training_data_retention, unname(expected[[learner]]))
  }

  params <- sb_classical_binary_params("aorsf_binary")
  params$num_trees <- 20L
  adapter <- sb_classical_binary_learner_adapter(
    "aorsf_binary",
    params
  )
  model <- adapter$fit(x, data$outcome)
  state <- adapter$extract_state(model)
  expect_null(model$engine$data)
  expect_false("training_data_retention" %in% names(state))
})

test_that("bagging state drops fitted-row outcomes and replays exactly", {
  skip_if_not_installed("ipred")
  data <- sb_example_binary(n = 60L, seed = 8423L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:8, , drop = FALSE]
  params <- sb_classical_binary_params("bagging_binary")
  params$num_trees <- 20L
  adapter <- sb_classical_binary_learner_adapter(
    "bagging_binary",
    params
  )
  model <- adapter$fit(x, data$outcome)
  expected <- adapter$predict(model, assessment)
  state <- adapter$extract_state(model)
  engine <- unserialize(state$engine_raw)
  replay <- adapter$predict_state(state, assessment)$probability
  has_full_row_payload <- function(value, n) {
    if (is.environment(value) || is.function(value)) {
      return(FALSE)
    }
    if ((is.data.frame(value) || is.matrix(value)) && nrow(value) == n) {
      return(TRUE)
    }
    if (is.atomic(value) && is.null(dim(value)) && length(value) == n) {
      return(TRUE)
    }
    if (is.list(value) || is.pairlist(value) || is.call(value)) {
      return(any(vapply(as.list(value), has_full_row_payload, logical(1), n)))
    }
    FALSE
  }

  expect_identical(levels(engine$y), c("0", "1"))
  expect_length(engine$y, 0L)
  expect_null(engine$X)
  expect_null(engine$OOB)
  expect_null(engine$call)
  expect_true(all(vapply(engine$mtrees, function(tree) {
    is.null(tree$bindx) && is.null(tree$btree$y) &&
      is.null(tree$btree$where) && is.null(tree$btree$call) &&
      identical(environment(tree$btree$terms), baseenv())
  }, logical(1))))
  expect_false(has_full_row_payload(engine, nrow(x)))
  expect_false("training_data_retention" %in% names(state))
  expect_equal(replay, expected, tolerance = 0)
})
