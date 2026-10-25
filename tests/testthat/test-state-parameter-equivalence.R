test_that("flexsurv states accept equivalent parameter representations", {
  skip_if_not_installed("flexsurv")
  withr::local_seed(741L)
  x <- data.frame(x1 = rnorm(60L), x2 = rnorm(60L))
  time <- rgamma(60L, shape = 2, rate = exp(0.2 * x$x1) / 10)
  event <- rep(c(1L, 1L, 0L), 20L)
  for (kind in c("parametric", "spline")) {
    id <- if (kind == "parametric") "flexsurv_gamma" else "flexsurvspline_hazard_k1"
    stem <- if (kind == "parametric") "flexsurv_parametric" else "flexsurvspline"
    params <- get(paste0("sb_", stem, "_fixed_params"))(id)
    model <- get(paste0("sb_fit_", stem, "_fixed"))(x, time, event, id, params)
    state <- get(paste0("sb_extract_", stem, "_state"))(model)
    predict <- get(paste0("sb_predict_", stem, "_state"))
    expected <- predict(state, x)
    changed <- state
    changed$params <- rev(state$params)
    expect_identical(predict(changed, x), expected)
    if (kind == "spline") {
      changed$params$k <- as.numeric(changed$params$k)
      expect_identical(predict(changed, x), expected)
      changed$params$k <- 2
    } else {
      changed$params$distribution <- "gompertz"
    }
    error <- tryCatch(predict(changed, x), error = identity)
    expect_s3_class(error, "stormbreaker_bad_learner_state")
  }
})

test_that("mboost states compare iteration values rather than storage types", {
  skip_if_not_installed("mboost")
  withr::local_seed(742L)
  x <- data.frame(x1 = rnorm(40L), x2 = rnorm(40L))
  time <- rgamma(40L, shape = 2, rate = exp(0.2 * x$x1) / 10)
  event <- rep(c(1L, 1L, 1L, 0L), 10L)
  for (kind in c("aft", "ranking")) {
    id <- if (kind == "aft") "glmboost_lognormal_aft" else "glmboost_gehan_aft"
    stem <- paste0("mboost_", kind)
    params <- get(paste0("sb_", stem, "_fixed_params"))(id)
    params$mstop <- 5L
    model <- get(paste0("sb_fit_", stem, "_fixed"))(x, time, event, id, params)
    state <- get(paste0("sb_extract_", stem, "_state"))(model)
    predict <- get(paste0("sb_predict_", stem, "_state"))
    expected <- predict(state, x)
    changed <- state
    changed$params$mstop <- 5
    expect_identical(predict(changed, x), expected)
    for (wrong in list(6, 5.01, TRUE, NA_real_)) {
      changed$params$mstop <- wrong
      error <- tryCatch(predict(changed, x), error = identity)
      expect_s3_class(error, "stormbreaker_bad_learner_state")
    }
  }
})
