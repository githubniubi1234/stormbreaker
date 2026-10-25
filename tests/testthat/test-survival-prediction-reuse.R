test_that("bagged kNN inverts each bootstrap covariance once per prediction batch", {
  skip_if_not_installed("bnnSurvival")
  d <- sb_example_survival(72L, 42L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_bagged_knn_survival_learner_adapter("bagged_knn_survival")
  model <- adapter$fit(x, d$time, d$event)
  native_solve <- base::solve.default
  count <- 0L
  local_mocked_bindings(solve.default = function(...) {
    count <<- count + 1L
    native_solve(...)
  }, .package = "base")
  curves <- sb_predict_bagged_knn_survival_curves(model, x[1:4, ])
  expect_identical(count, model$params$num_base_learners)
  expect_identical(nrow(curves$survival), 4L)
})

test_that("Akritas multivariate prediction avoids rebuilding an R6 distribution", {
  skip_if_not_installed("distr6"); skip_if_not_installed("survivalmodels")
  d <- sb_example_survival(72L, 42L); x <- d[c("gene_a", "gene_b", "gene_c")]
  model <- sb_fit_akritas_survival(x, d$time, d$event, "akritas_survival", sb_akritas_params())
  old <- sb_akritas_distribution
  calls <- 0L
  local_mocked_bindings(sb_akritas_distribution = function(basis) {
    calls <<- calls + 1L
    old(basis)
  })
  sb_akritas_curves(model, x[1:4, ])
  expect_identical(calls, 0L)
})

test_that("Akritas state prediction evaluates its risk grid and horizon in one call", {
  skip_if_not_installed("distr6"); skip_if_not_installed("survivalmodels")
  d <- sb_example_survival(72L, 42L); x <- d[c("gene_a", "gene_b")]
  model <- sb_fit_akritas_survival(x, d$time, d$event, "akritas_survival", sb_akritas_params())
  h <- as.numeric(median(d$time))
  state <- sb_extract_akritas_state(model, list(time = d$time, prediction_horizon = h, time_unit = "month"))
  original <- sb_akritas_curves
  calls <- 0L
  local_mocked_bindings(sb_akritas_curves = function(...) {
    calls <<- calls + 1L
    original(...)
  })
  result <- sb_predict_akritas_state(state, x[1:4, ])
  expect_identical(calls, 1L)
  expect_identical(result$survival_probability, as.numeric(original(state, x[1:4, ], h)))
  expect_identical(result$risk, sb_predict_akritas_risk(model, x[1:4, ]))
})

test_that("bagged kNN preserves native curves, ties, one-feature inputs and failures", {
  skip_if_not_installed("bnnSurvival")
  d <- sb_example_survival(72L, 42L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("gene A", "gene-B", "gene C")
  x[2L, ] <- x[1L, ]
  adapter <- sb_bagged_knn_survival_learner_adapter("bagged_knn_survival")
  native <- get("predict", asNamespace("bnnSurvival"))
  old_method <- methods::getMethod("predict", "bnnSurvivalEnsemble", where = asNamespace("bnnSurvival"))
  withr::local_seed(19)
  withr::local_options(mc.cores = 7L)
  before <- .Random.seed
  for (p in c(1L, 3L)) {
    train <- x[seq_len(p)]
    model <- adapter$fit(train, d$time, d$event)
    for (rows in list(1L, c(1L, 2L, 1L, 4L, 8L))) {
      query <- train[rows, , drop = FALSE]
      internal <- sb_bagged_knn_internal_data(query)$data
      expected <- sb_with_bagged_knn_single_core(native(model$engine, internal))
      observed <- sb_predict_bagged_knn_survival_curves(model, query)
      expect_identical(observed$survival, as.matrix(bnnSurvival::predictions(expected)))
      expect_identical(observed$time, as.numeric(bnnSurvival::timepoints(expected)))
      state <- adapter$extract_state(model, list(time = d$time,
        prediction_horizon = as.numeric(median(d$time)), time_unit = "month"))
      replay <- adapter$predict_state(unserialize(serialize(state, NULL)), query)
      expect_identical(replay$risk, observed$risk)
      expect_identical(.Random.seed, before)
      expect_identical(getOption("mc.cores"), 7L)
    }
  }
  expect_identical(methods::getMethod("predict", "bnnSurvivalEnsemble", where = asNamespace("bnnSurvival")), old_method)
  singular <- x[c(1L, 1L)]
  names(singular) <- c("a", "duplicate")
  model <- adapter$fit(singular, d$time, d$event)
  query <- singular[1:2, ]; internal <- sb_bagged_knn_internal_data(query)$data
  expected <- tryCatch(sb_with_bagged_knn_single_core(native(model$engine, internal)), error = identity)
  observed <- tryCatch(sb_predict_bagged_knn_survival_curves(model, query), error = identity)
  expect_s3_class(expected, "error")
  expect_identical(class(observed), class(expected))
  expect_identical(conditionMessage(observed), conditionMessage(expected))
})

test_that("Akritas native CDF and batched times retain native boundary behavior", {
  skip_if_not_installed("distr6"); skip_if_not_installed("survivalmodels")
  withr::local_seed(9281)
  for (p in c(1L, 2L, 20L)) {
    basis <- matrix(rnorm(60L * p), ncol = p)
    basis[2L, ] <- basis[1L, ]
    query <- rbind(basis[1:3, , drop = FALSE], rep(-1e6, p), rep(1e6, p))
    expected <- as.numeric(sb_akritas_distribution(basis)$cdf(data = query))
    expect_identical(sb_akritas_empirical_cdf(basis, query), expected)
  }
  d <- sb_example_survival(72L, 42L)
  d$time <- rep(seq_len(12L), length.out = nrow(d))
  x <- d[c("gene_a", "gene_b", "gene_c")]
  before <- .Random.seed
  for (p in c(1L, 3L)) {
    train <- x[seq_len(p)]
    model <- sb_fit_akritas_survival(train, d$time, d$event, "akritas_survival", sb_akritas_params())
    native <- utils::getFromNamespace("C_Akritas", "survivalmodels")
    distribution <- sb_akritas_distribution(unname(as.matrix(train)))
    query <- train[c(1, 2, 7), , drop = FALSE]
    q <- as.numeric(distribution$cdf(data = unname(as.matrix(query))))
    for (h in c(min(d$time) / 2, min(d$time), median(d$time), median(d$time) + 0.25, max(d$time))) {
      state <- sb_extract_akritas_state(model, list(time = d$time,
        prediction_horizon = as.numeric(h), time_unit = "month"))
      expected <- as.numeric(native(model$truth, as.numeric(h), model$grid,
        model$cdf_train, q, model$params$lambda))
      replay <- sb_predict_akritas_state(unserialize(serialize(state, NULL)), query)
      expect_identical(replay$survival_probability, expected)
      expect_identical(replay$risk, sb_predict_akritas_risk(model, query))
    }
  }
  expect_identical(.Random.seed, before)
})
