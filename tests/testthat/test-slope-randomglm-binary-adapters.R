make_slope_randomglm_fixture <- function(n = 90L, p = 100L, seed = 2901L) {
  set.seed(seed)
  x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("gene_", seq_len(p))
  linear_predictor <- 1.3 * x[, 1L] - x[, 2L] + 0.6 * x[, 3L]
  outcome <- stats::rbinom(n, 1L, stats::plogis(linear_predictor))
  data <- data.frame(
    patient_id = paste0("patient_", seq_len(n)),
    outcome = outcome,
    x,
    check.names = FALSE
  )
  list(data = data, x = data[colnames(x)], outcome = outcome)
}


test_that("missing SLOPE and randomGLM engines fail explicitly", {
  fixture <- make_slope_randomglm_fixture(n = 30L, p = 5L)

  for (learner in sb_slope_randomglm_binary_learner_ids()) {
    package <- if (identical(learner, "slope_binary")) "SLOPE" else "randomGLM"
    expected_class <- if (identical(learner, "slope_binary")) {
      "stormbreaker_missing_slope"
    } else {
      "stormbreaker_missing_randomglm"
    }
    local_mocked_bindings(
      sb_slope_randomglm_binary_dependency = function(learner) {
        list(package = package, available = FALSE)
      }
    )
    expect_error(
      sb_slope_randomglm_binary_learner_adapter(learner)$fit(
        fixture$x,
        fixture$outcome
      ),
      class = expected_class
    )
  }
})


test_that("compact randomGLM probabilities equal its native ensemble", {
  skip_if_not_installed("randomGLM")
  fixture <- make_slope_randomglm_fixture(n = 70L, p = 100L, seed = 2903L)
  params <- sb_slope_randomglm_binary_params("random_glm_binary")
  params$n_bags <- 10L
  native <- randomGLM::randomGLM(
    x = as.matrix(fixture$x),
    y = factor(fixture$outcome, levels = c(0L, 1L)),
    type = "binary",
    maxInteractionOrder = params$max_interaction_order,
    nBags = params$n_bags,
    replace = params$replace,
    nCandidateCovariates = params$n_candidate_covariates,
    keepModels = TRUE,
    randomSeed = params$seed,
    nThreads = params$threads,
    verbose = 0
  )
  compact <- sb_compact_random_glm_binary(native, names(fixture$x), params)

  expected <- unname(stats::predict(
    native,
    newdata = as.matrix(fixture$x),
    type = "response"
  )[, "1"])
  observed <- sb_predict_slope_randomglm_binary(compact, fixture$x)
  state <- sb_extract_slope_randomglm_binary_state(compact)

  expect_equal(observed, expected, tolerance = 1e-14)
  expect_equal(
    sb_predict_slope_randomglm_binary_state(state, fixture$x)$probability,
    observed,
    tolerance = 0
  )
  expect_false(any(c("models", "training_data", "outcome") %in% names(state)))
})

test_that("randomGLM fitting preserves an absent caller RNG state", {
  skip_if_not_installed("randomGLM")
  fixture <- make_slope_randomglm_fixture(n = 50L, p = 30L, seed = 2906L)
  params <- sb_slope_randomglm_binary_params("random_glm_binary")
  params$n_bags <- 3L
  previous_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit(
    {
      if (is.null(previous_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  sb_fit_random_glm_binary(
    list(x = fixture$x, outcome = fixture$outcome),
    params
  )

  expect_false(exists(".Random.seed", envir = .GlobalEnv))
})

test_that("both learners run through binary CV and final prediction", {
  skip_if_not_installed("SLOPE")
  skip_if_not_installed("randomGLM")
  production_params <- sb_slope_randomglm_binary_params
  local_mocked_bindings(
    sb_slope_randomglm_binary_params = function(learner) {
      params <- production_params(learner)
      if (identical(learner, "slope_binary")) {
        params$q <- c(0.1, 0.2)
        params$inner_folds <- 3L
        params$path_length <- 10L
      } else {
        params$n_bags <- 10L
      }
      params
    }
  )
  fixture <- make_slope_randomglm_fixture(n = 90L, p = 30L, seed = 2904L)
  task <- sb_binary_task(
    fixture$data,
    "patient_id",
    "outcome",
    names(fixture$x)
  )

  for (learner in sb_slope_randomglm_binary_learner_ids()) {
    pipeline <- sb_pipeline(learner = learner)
    cv <- sb_run_cv(task, pipeline = pipeline, folds = 3L, seed = 2905L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      fixture$data[1:6, c("patient_id", names(fixture$x))],
      "patient_id"
    )
    version_field <- if (identical(learner, "slope_binary")) {
      "slope_version"
    } else {
      "randomglm_version"
    }

    expect_identical(cv$status, "success", info = learner)
    expect_identical(cv$complete_oof, TRUE, info = learner)
    expect_true(all(is.finite(cv$predictions$probability)), info = learner)
    expect_true(version_field %in% names(cv$provenance), info = learner)
    expect_identical(fit$status, "success", info = learner)
    expect_identical(prediction$status, "success", info = learner)
    expect_true(all(is.finite(prediction$predictions$probability)), info = learner)
    expect_true(version_field %in% names(fit$provenance), info = learner)
  }
})
