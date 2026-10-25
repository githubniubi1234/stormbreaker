test_that("Bujar views preserve original coefficients and residual event steps", {
  skip_if_not_installed("bujar")
  d <- sb_example_survival(160L, seed = 8612L)
  features <- c("gene_c", "gene_a", "gene_b")
  x <- as.matrix(d[features])
  colnames(x) <- paste0("x", seq_along(features))
  native <- withr::with_seed(1L, bujar::bujar(y = log(d$time), cens = d$event, x = x,
    degree = 1, learner = "linear.regression", center = TRUE, mimpu = NULL,
    iter.bj = 20L, max.cycle = 10L, nu = .1, mstop = 50L, twin = FALSE,
    mstop2 = 100L, tuning = FALSE, cv = FALSE, vimpint = FALSE, rel.inf = FALSE,
    tol = 1e-4, n.cores = 1L, rng = 1L, trace = FALSE))
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
    prediction_horizon = stats::median(d$time))
  fit <- sb_fit_final(task, sb_pipeline(learner = "bujar_linear_boost"))
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  co <- sb_bujar_diagnostics(fit)
  expect_identical(co$features, features)
  expect_equal(co$table$coefficient, unname(native$coef.bj[-1]), tolerance = 0)
  expect_equal(co$intercept, unname(native$coef.bj[1]), tolerance = 0)
  lp <- as.numeric(stats::predict(native, newx = x))
  expect_equal(as.numeric(co$intercept + x %*% co$table$coefficient), lp, tolerance = 0)
  residual <- log(d$time) - lp
  km <- survival::survfit(survival::Surv(residual, d$event) ~ 1)
  keep <- km$n.event > 0
  curve <- sb_bujar_diagnostics(fit, "residual_survival")
  expect_equal(curve$table$residual, unname(km$time[keep]), tolerance = 0)
  expect_equal(curve$table$survival, unname(km$surv[keep]), tolerance = 0)
  expect_gt(sum(curve$table$residual < 0), 0)
  bundle <- sb_native_diagnostics(fit, complete = TRUE)
  expect_identical(bundle$status$status, c("available", "available"))
  for (v in bundle$diagnostics) {
    p <- sb_plot_algorithm_diagnostic(v)
    b <- ggplot2::ggplot_build(p)
    expect_s3_class(b, "ggplot_built")
    expect_null(p$labels$caption)
    expect_null(p$labels$subtitle)
    if (v$type == "bujar_residual_survival") {
      expect_equal(b$data[[1]]$x, v$table$residual, tolerance = 0)
      expect_equal(b$data[[1]]$y, v$table$survival, tolerance = 0)
    } else expect_equal(b$data[[2]]$x, v$table$coefficient, tolerance = 0)
  }
  expect_identical(serialize(fit, NULL, version = 2L), before)
  no_horizon <- sb_fit_final(sb_surv_task(d, "patient_id", "time", "event", features, "month"),
    sb_pipeline(learner = "bujar_linear_boost"))
  expect_identical(no_horizon$status, "success")
  expect_identical(sb_native_diagnostics(no_horizon, complete = TRUE)$status$status,
    c("available", "unavailable"))
})
