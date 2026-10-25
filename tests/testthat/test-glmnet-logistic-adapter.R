test_that("glmnet logistic can use a converged prefix of a truncated path", {
  skip_if_not_installed("glmnet")
  x <- withr::with_seed(5106, as.data.frame(matrix(stats::rnorm(60 * 6), 60, 6)))
  y <- rep(0:1, 30L)
  native <- glmnet::glmnet
  messages <- character()
  reference <- withCallingHandlers(
    native(as.matrix(x), y, family = "binomial", maxit = 20L),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_gt(length(messages), 0L)
  expect_gte(length(reference$lambda), 3L)
  lambda <- reference$lambda[[3L]]
  local_mocked_bindings(glmnet = function(...) {
    args <- list(...)
    args$maxit <- 20L
    do.call(native, args)
  }, .package = "glmnet")
  adapter <- sb_glmnet_logistic_learner_adapter(alpha = 1, lambda = lambda)
  fitted <- adapter$fit(x, y)
  expect_identical(attr(fitted, "stormbreaker_fit_warnings"), messages)
  expect_identical(adapter$predict(fitted, x), unname(drop(stats::predict(
    reference, newx = as.matrix(x), s = lambda, type = "response"))))
  outside <- sb_glmnet_logistic_learner_adapter(
    alpha = 1, lambda = min(reference$lambda) / 10)
  expect_s3_class(tryCatch(outside$fit(x, y), error = identity),
    "stormbreaker_glmnet_lambda_outside_path")
})

test_that("the glmnet logistic adapter matches a direct path fit", {
  skip_if_not_installed("glmnet")

  d <- make_sb_binary_fixture(n = 60)
  analysis_rows <- 1:45
  assessment_rows <- 46:60
  predictors <- c("gene_a", "gene_b", "gene_c")
  x_analysis <- d[analysis_rows, predictors, drop = FALSE]
  x_assessment <- d[assessment_rows, predictors, drop = FALSE]
  manual <- glmnet::glmnet(
    x = as.matrix(x_analysis),
    y = d$outcome[analysis_rows],
    family = "binomial",
    alpha = 0.5,
    standardize = TRUE
  )
  lambda <- manual$lambda[min(20L, length(manual$lambda))]
  learner <- sb_glmnet_logistic_learner_adapter(
    alpha = 0.5,
    lambda = lambda,
    standardize = TRUE
  )
  model <- learner$fit(x_analysis, d$outcome[analysis_rows])
  observed <- learner$predict(model, x_assessment)
  expected <- unname(drop(stats::predict(
    manual,
    newx = as.matrix(x_assessment),
    s = lambda,
    type = "response"
  )))

  expect_identical(learner$id, "glmnet_logistic")
  expect_equal(observed, expected, tolerance = 1e-12)
  expect_identical(all(observed >= 0 & observed <= 1), TRUE)

  coefficients <- learner$extract_coefficients(model)
  expected_coefficients <- as.matrix(stats::coef(manual, s = lambda))
  expect_identical(
    coefficients$term,
    c("(Intercept)", predictors)
  )
  expect_equal(
    coefficients$coefficient,
    unname(expected_coefficients[, 1L]),
    tolerance = 1e-12
  )
  design <- cbind("(Intercept)" = 1, as.matrix(x_assessment))
  reconstructed <- unname(drop(stats::plogis(
    design %*% coefficients$coefficient[
      match(colnames(design), coefficients$term)
    ]
  )))
  expect_equal(observed, reconstructed, tolerance = 1e-12)
})

test_that("the glmnet logistic adapter locks schema and settings", {
  skip_if_not_installed("glmnet")

  d <- make_sb_binary_fixture(n = 60)
  manual <- glmnet::glmnet(
    x = as.matrix(d[c("gene_a", "gene_b")]),
    y = d$outcome,
    family = "binomial"
  )
  lambda <- manual$lambda[min(10L, length(manual$lambda))]
  learner <- sb_glmnet_logistic_learner_adapter(1, lambda)
  model <- learner$fit(d[c("gene_a", "gene_b")], d$outcome)

  schema_error <- capture_sb_error(
    learner$predict(model, d[c("gene_b", "gene_a")])
  )
  expect_s3_class(schema_error, "stormbreaker_learner_schema_mismatch")

  alpha_error <- capture_sb_error(
    sb_glmnet_logistic_learner_adapter(1.2, lambda)
  )
  expect_s3_class(alpha_error, "stormbreaker_bad_glmnet_alpha")

  lambda_error <- capture_sb_error(
    sb_glmnet_logistic_learner_adapter(1, 0)
  )
  expect_s3_class(lambda_error, "stormbreaker_bad_glmnet_lambda")
})

test_that("the glmnet logistic adapter exposes dependency and warnings", {
  d <- make_sb_binary_fixture(n = 40)
  local_mocked_bindings(sb_glmnet_available = function() FALSE)
  learner <- sb_glmnet_logistic_learner_adapter(1, 0.05)
  missing_error <- capture_sb_error(
    learner$fit(d[c("gene_a", "gene_b")], d$outcome)
  )
  expect_s3_class(missing_error, "stormbreaker_missing_glmnet")

  skip_if_not_installed("glmnet")
  local_mocked_bindings(sb_glmnet_available = function() TRUE)
  original_glmnet <- glmnet::glmnet
  local_mocked_bindings(
    glmnet = function(...) {
      warning("synthetic logistic glmnet warning", call. = FALSE)
      original_glmnet(...)
    },
    .package = "glmnet"
  )
  warning_error <- capture_sb_error(
    learner$fit(d[c("gene_a", "gene_b")], d$outcome)
  )
  expect_s3_class(
    warning_error,
    "stormbreaker_glmnet_logistic_fit_warning"
  )
  expect_match(
    conditionMessage(warning_error),
    "synthetic logistic glmnet warning"
  )
})

test_that("glmnet logistic rejects an ambiguous intercept feature", {
  skip_if_not_installed("glmnet")
  d <- make_sb_binary_fixture(n = 40)
  x <- d[c("gene_a", "gene_b")]
  names(x)[1L] <- "(Intercept)"
  learner <- sb_glmnet_logistic_learner_adapter(1, 0.05)

  error <- capture_sb_error(learner$fit(x, d$outcome))

  expect_s3_class(error, "stormbreaker_reserved_intercept_name")
})

test_that("glmnet logistic tolerates only roundoff at path endpoints", {
  skip_if_not_installed("glmnet")

  d <- make_sb_binary_fixture(n = 60)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  native <- glmnet::glmnet(as.matrix(x), d$outcome, family = "binomial")
  endpoints <- range(native$lambda)
  near <- endpoints * c(1 - 2 * .Machine$double.eps,
    1 + 2 * .Machine$double.eps)
  for (lambda in near) {
    learner <- sb_glmnet_logistic_learner_adapter(1, lambda)
    model <- learner$fit(x, d$outcome)
    expect_identical(model$lambda, native$lambda)
    expect_identical(attr(model, "stormbreaker_lambda"), lambda)
    expected <- as.numeric(predict(native, as.matrix(x), s = lambda,
      type = "response"))
    expect_equal(learner$predict(model, x), expected, tolerance = 0)
    state <- learner$extract_state(model)
    expect_equal(learner$predict_state(state, x)$probability, expected,
      tolerance = 1e-12)
  }
  for (lambda in endpoints * c(1 - 1e-6, 1 + 1e-6)) {
    learner <- sb_glmnet_logistic_learner_adapter(1, lambda)
    error <- capture_sb_error(learner$fit(x, d$outcome))
    expect_s3_class(error, "stormbreaker_glmnet_lambda_outside_path")
  }
  for (scale in c(1e-200, 1, 1e200)) {
    path <- c(0.01, 1) * scale
    expect_true(sb_glmnet_lambda_in_path(path[[1L]] *
      (1 - 2 * .Machine$double.eps), path))
    expect_true(sb_glmnet_lambda_in_path(path[[2L]] *
      (1 + 2 * .Machine$double.eps), path))
    expect_false(sb_glmnet_lambda_in_path(path[[1L]] * (1 - 1e-6), path))
    expect_false(sb_glmnet_lambda_in_path(path[[2L]] * (1 + 1e-6), path))
  }
})
