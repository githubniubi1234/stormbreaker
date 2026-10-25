make_cox_control_fixture <- function() {
  withr::with_seed(12413L, {
    n <- 160L
    x <- as.data.frame(matrix(stats::rnorm(n * 4L), nrow = n))
    names(x) <- paste0("gene", seq_len(ncol(x)))
    event_time <- stats::rexp(n, exp(1.1 * x$gene1 - 0.7 * x$gene2))
    censor_time <- stats::rexp(n, 0.2)
    list(x = x, time = pmin(event_time, censor_time),
      event = as.integer(event_time <= censor_time))
  })
}

test_that("StepCox tightens likelihood precision without changing warning tolerance", {
  d <- make_cox_control_fixture()
  fit <- sb_fit_stepcox(d$x, d$time, d$event, "both")
  native <- survival::coxph.control()
  expect_equal(fit$cox_control$eps, 1e-11)
  expect_identical(fit$cox_control$iter.max, 50L)
  expect_identical(fit$cox_control$toler.inf, native$toler.inf)
  expect_identical(fit$cox_control$toler.chol, native$toler.chol)
  expect_gt(length(fit$selected_features), 0L)
})

test_that("PLS-Cox preserves the native infinite-coefficient warning threshold", {
  skip_if_not_installed("plsRcox")
  d <- make_cox_control_fixture()
  fit <- sb_plsrcox_cv_learner_adapter(inner_folds = 3L, seed = 1L,
    maximum_components = 3L)$fit(d$x, d$time, d$event)
  native <- survival::coxph.control()
  expect_equal(fit$audit$cox_control$eps, 1e-11)
  expect_identical(fit$audit$cox_control$iter.max, 50L)
  expect_identical(fit$audit$cox_control$toler.inf, native$toler.inf)
  expect_identical(fit$audit$cox_control$toler.chol, native$toler.chol)
})

test_that("StepCox still rejects complete survival separation", {
  x <- data.frame(separated = rep(c(1, 0), each = 20L))
  error <- tryCatch(sb_fit_stepcox(x, seq_len(40L), rep(1L, 40L), "both"),
    error = identity)
  expect_s3_class(error, "stormbreaker_stepcox_warning")
  expect_match(conditionMessage(error), "coefficient may be infinite", fixed = TRUE)
})
