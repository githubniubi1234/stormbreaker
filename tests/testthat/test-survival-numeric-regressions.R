make_survival_numeric_fixture <- function() {
  withr::with_seed(20260930L, {
    n <- 600L
    x <- as.data.frame(matrix(stats::rnorm(n * 24L), nrow = n))
    names(x) <- sprintf("gene%02d", seq_len(ncol(x)))
    beta <- c(.65, -.60, .55, -.50, .45, -.40, .30, -.25, rep(0, 16L))
    eta <- drop(as.matrix(x) %*% beta)
    event_time <- 365 * (-log(stats::runif(n)) / exp(eta))^(1 / 1.3)
    censor_time <- stats::rexp(n, rate = 1 / 900)
    list(x = x, time = pmin(event_time, censor_time),
         event = as.integer(event_time <= censor_time))
  })
}

test_that("gamboost fits more than 500 rows without sparse QR warnings", {
  skip_if_not_installed("mboost")
  d <- make_survival_numeric_fixture()
  old <- getOption("mboost_useMatrix")
  ids <- c("gamboost_cox", "gamboost_weibull_aft",
    "gamboost_loglogistic_aft", "gamboost_lognormal_aft",
    "gamboost_gehan_aft", "gamboost_unoc")
  for (id in ids) {
    adapter <- sb_resolve_pipeline_runtime(sb_pipeline("none", id))$learner
    model <- adapter$fit(d$x, d$time, d$event)
    risk <- adapter$predict(model, d$x[seq_len(12L), , drop = FALSE])
    expect_identical(all(is.finite(risk)), TRUE)
    expect_gt(length(unique(risk)), 1L)
    expect_equal(mboost::mstop(model$engine), 100)
    expect_identical(getOption("mboost_useMatrix"), old)
  }
})

test_that("Gompertz initialization is stable across time units", {
  skip_if_not_installed("flexsurv")
  d <- make_survival_numeric_fixture()
  adapter <- sb_flexsurv_parametric_learner_adapter("flexsurv_gompertz")
  days <- adapter$fit(d$x, d$time, d$event)
  years <- adapter$fit(d$x, d$time / 365, d$event)
  assessment <- d$x[seq_len(12L), , drop = FALSE]
  p_days <- sb_flexsurv_parametric_direct_probability(days, assessment, 365)
  p_years <- sb_flexsurv_parametric_direct_probability(years, assessment, 1)
  expect_equal(p_days, p_years, tolerance = 1e-5)
  expect_identical(days$engine$opt$convergence, 0L)
  expect_identical(all(is.finite(days$engine$cov)), TRUE)
  native_x <- sb_flexsurv_numerical_predictors(assessment, days$numerical_fitting)
  names(native_x) <- days$engine_features
  native <- as.numeric(stats::predict(days$engine, newdata = native_x,
    type = "survival", times = 365, conf.int = FALSE)$.pred_survival)
  expect_equal(p_days, native, tolerance = 1e-12)
  expect_equal(adapter$predict(years, assessment) - adapter$predict(days, assessment),
    rep(log(365), nrow(assessment)), tolerance = 1e-5)
})

test_that("LPC explains the native fifty-feature automatic-search limitation", {
  skip_if_not_installed("lpc")
  d <- make_survival_numeric_fixture()
  endpoint <- list(task_type = "survival", time = d$time, event = d$event)
  error <- tryCatch(sb_fit_lpc_selector(d$x, endpoint, "lpc_cox_nlogn"), error = identity)
  expect_s3_class(error, "stormbreaker_invalid_lpc_selector")
  expect_s3_class(error, "stormbreaker_lpc_dimension_unsupported")
  expect_match(conditionMessage(error), "at least 50 predictors", fixed = TRUE)
})

test_that("PLS-Cox resolves late-component numerical warnings without lowering the cap", {
  skip_if_not_installed("plsRcox")
  d <- make_survival_numeric_fixture()
  panels <- list(c(1:9, 11, 13, 16, 20, 21, 24), c(1:9, 16))
  for (panel in panels) {
    x <- d$x[panel]
    adapter <- sb_plsrcox_cv_learner_adapter(inner_folds = 5L, seed = 1L,
      maximum_components = 10L)
    model <- adapter$fit(x, d$time, d$event)
    expect_length(model$audit$component_search_warnings, 0L)
    expect_identical(model$audit$maximum_components, 10L)
    expect_equal(model$audit$cox_control$eps, 1e-11)
    expect_identical(model$audit$cox_control$toler.chol,
      survival::coxph.control()$toler.chol)
    expected <- as.numeric(stats::predict(model$engine,
      newdata = x[1:12, , drop = FALSE], type = "lp"))
    actual <- adapter$predict(model, x[1:12, , drop = FALSE])
    expect_equal(actual - mean(actual), expected - mean(expected), tolerance = 1e-8)
    expect_gt(length(unique(actual)), 1L)
  }
})

test_that("PLS-Cox still rejects genuine selected-component search warnings", {
  skip_if_not_installed("plsRcox")
  d <- make_survival_numeric_fixture()
  testthat::local_mocked_bindings(cv.plsRcox = function(...) {
    warning("Loglik converged before variable  1 ; coefficient may be infinite. ")
    list(lambda.min5 = 1L)
  }, .package = "plsRcox")
  error <- tryCatch(sb_plsrcox_cv_learner_adapter()$fit(d$x, d$time, d$event),
    error = identity)
  expect_s3_class(error, "stormbreaker_plsrcox_cv_warning")
  expect_s3_class(error$parent, "warning")
})
