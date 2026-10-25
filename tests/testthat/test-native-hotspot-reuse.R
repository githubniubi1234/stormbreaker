test_that("Lasso selection reuses the full analysis path already fitted by CV", {
  skip_if_not_installed("glmnet")
  d <- sb_example_survival(120L, 11032L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  native <- glmnet::glmnet
  calls <- 0L
  local_mocked_bindings(glmnet = function(...) {
    calls <<- calls + 1L
    native(...)
  }, .package = "glmnet")
  actual <- sb_fit_lasso_nonzero(x, d$time, d$event, 3L, 17L,
    capture_history = TRUE)
  expect_identical(calls, 4L)
  model <- sb_glmnet_cox_learner_adapter(1, actual$lambda)$fit(x, d$time, d$event)
  coefficient <- drop(as.matrix(stats::coef(model, s = actual$lambda)))
  expect_identical(actual$scores$coefficient, unname(coefficient))
  expect_identical(actual$coefficient_path, list(lambda = model$lambda, beta = model$beta))
  expect_identical(actual$selected_features, names(coefficient)[coefficient != 0])
})

test_that("RSF selector avoids proximity without changing native selection", {
  skip_if_not_installed("randomForestSRC")
  d <- sb_example_survival(90L, 11032L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  native <- randomForestSRC::rfsrc
  requested <- list()
  local_mocked_bindings(rfsrc = function(...) {
    args <- list(...)
    requested[[length(requested) + 1L]] <<- args$proximity
    native(...)
  }, .package = "randomForestSRC")
  withr::local_seed(112L)
  before <- .Random.seed
  actual <- sb_fit_rsf_minimal_depth(x, d$time, d$event, 5L, 17L)
  after <- .Random.seed
  expect_identical(requested[[1L]], FALSE)
  original <- sb_fit_rsf_minimal_depth
  body(original) <- parse(text = sub("proximity = FALSE", "proximity = TRUE",
    paste(deparse(body(original)), collapse = "\n"), fixed = TRUE))[[1L]]
  assign(".Random.seed", before, .GlobalEnv)
  expected <- original(x, d$time, d$event, 5L, 17L)
  expect_identical(actual, expected)
  expect_identical(.Random.seed, after)
})

test_that("RSF native depth traversal can use private reference containers", {
  skip_if_not_installed("randomForestSRC")
  d <- sb_example_survival(90L, 11032L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  data <- sb_additional_survival_formula_data(x, d$time, d$event,
    random_forest_src = TRUE)
  engine <- withr::with_options(list(rf.cores = 1L, mc.cores = 1L),
    randomForestSRC::rfsrc(data$formula, data$data, ntree = 30L,
      nodesize = 5L, splitrule = "logrank", importance = TRUE, seed = 17L))
  ns <- asNamespace("randomForestSRC")
  parser <- get("rfsrcParseTree", ns)
  maximum <- randomForestSRC::max.subtree
  before <- .Random.seed
  expected <- withr::with_options(list(rf.cores = 1L, mc.cores = 1L),
    randomForestSRC::var.select(engine, conservative = "high", verbose = FALSE))
  after <- .Random.seed
  assign(".Random.seed", before, .GlobalEnv)
  actual <- withr::with_options(list(rf.cores = 1L, mc.cores = 1L),
    sb_rsf_minimal_depth_selection(engine))
  expect_identical(actual, expected)
  expect_identical(.Random.seed, after)
  expect_identical(get("rfsrcParseTree", ns), parser)
  expect_identical(randomForestSRC::max.subtree, maximum)
})
