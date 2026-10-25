make_coxboost_diagnostic_task <- function() {
  set.seed(8901L)
  x <- as.data.frame(matrix(stats::rnorm(90L * 12L), nrow = 90L))
  names(x) <- sprintf("gene_%02d", seq_len(12L))
  lp <- x[[1]] - x[[2]] + .6 * x[[3]]
  et <- stats::rexp(90L, rate = exp(lp) / 15)
  ct <- stats::rexp(90L, rate = 1 / 25)
  d <- data.frame(id = seq_len(90L), time = pmax(pmin(et, ct), .05),
    event = as.integer(et <= ct), x)
  sb_surv_task(d, "id", "time", "event", rev(names(x)), "month")
}

test_that("both CoxBoost variants expose their exact captured coefficient paths", {
  skip_if_not_installed("CoxBoost")
  for (id in c("coxboost", "coxboost_cv")) {
    pipeline <- if (id == "coxboost_cv") sb_pipeline(learner = id,
      learner_tuning_folds = 3L, learner_seed = 8902L) else sb_pipeline(learner = id)
    fit <- sb_fit_final(make_coxboost_diagnostic_task(), pipeline, capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    original <- serialize(fit, NULL)
    model <- fit$diagnostic_artifacts$learner
    engine <- model$engine
    q <- sb_boosting_diagnostics(fit, "coefficient_path")
    native <- as.matrix(stats::coef(engine, at.step = seq.int(0L, engine$stepno), scaled = FALSE))
    expect_equal(unname(q$native_table), unname(native))
    expected <- sweep(native, 2L, engine$sdx, "/")
    expect_equal(q$table$coefficient, as.vector(expected))
    endpoint <- q$table[q$table$step == engine$stepno, ]
    expect_identical(endpoint$feature, fit$learner_state$features)
    expect_equal(endpoint$coefficient, fit$learner_state$coefficient)
    centered <- sweep(as.matrix(fit$diagnostic_artifacts$predictors), 2L, engine$meanx, "-")
    expect_equal(as.numeric(centered %*% endpoint$coefficient),
      as.numeric(engine$linear.predictor[nrow(engine$linear.predictor), ]), tolerance = 1e-10)
    count <- sb_boosting_diagnostics(fit, "nonzero_path")
    expect_equal(count$table$nonzero, rowSums(native != 0))
    expect_equal(count$table$step, seq.int(0L, engine$stepno))
    expect_identical(serialize(fit, NULL), original)
    if (requireNamespace("ggplot2", quietly = TRUE)) for (z in list(q, count)) {
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    }
    fit$diagnostic_artifacts <- NULL
    err <- tryCatch(sb_boosting_diagnostics(fit, "coefficient_path"), error = identity)
    expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
    cv <- tryCatch(sb_boosting_diagnostics(fit, "cv_curve"), error = identity)
    if (id == "coxboost_cv") {
      expect_identical(cv$table, model$audit$tuning_curve)
      expect_equal(cv$stopping_iteration, model$audit$optimal_step)
      expect_equal(cv$stopping_iteration, cv$table$step[which.max(cv$table$mean_log_likelihood)])
      if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(
        ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(cv)), "ggplot_built")
    } else expect_s3_class(cv, "stormbreaker_unavailable_diagnostic")
  }
})

test_that("single-feature CoxBoost paths retain a matrix and all zero-step rows", {
  skip_if_not_installed("CoxBoost")
  d <- make_sb_fixture(n = 100L, seed = 8612L)
  fit <- sb_fit_final(sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month"),
    sb_pipeline(learner = "coxboost", coxboost_stepno = 1L), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  q <- sb_boosting_diagnostics(fit, "coefficient_path")
  expect_equal(dim(q$native_table), c(2L, 1L))
  expect_equal(q$table$step, c(0L, 1L))
  expect_equal(q$table$coefficient[1], 0)
  expect_equal(q$table$coefficient[2], fit$learner_state$coefficient)
  broken <- fit
  broken$learner_state$coefficient <- broken$learner_state$coefficient + 1
  err <- tryCatch(sb_boosting_diagnostics(broken, "coefficient_path"), error = identity)
  expect_s3_class(err, "stormbreaker_invalid_diagnostic")
})
