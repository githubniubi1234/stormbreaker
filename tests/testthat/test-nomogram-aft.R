test_that("AFT nomograms reconstruct native scores and distribution survival", {
  skip_if_not_installed("survival")
  d <- sb_example_survival(160L, seed = 8612L)
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month", prediction_horizon = 7)
  for (id in sb_survreg_learner_ids()) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id),
      capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    n <- sb_nomogram(fit, d)
    v <- sb_nomogram_values(n, d)
    engine <- fit$diagnostic_artifacts$learner$engine
    eta <- drop(cbind(1, as.matrix(d[features])) %*% stats::coef(engine))
    expect_equal(v$score, -eta, tolerance = 1e-12)
    expect_equal(n$intercept, -unname(stats::coef(engine)[1]), tolerance = 0)
    expect_equal(n$features$coefficient, -unname(stats::coef(engine)[-1]), tolerance = 0)
    native_distribution <- survival::survreg.distributions[[engine$dist]]
    if (!is.null(native_distribution$dist)) native_distribution <-
      survival::survreg.distributions[[native_distribution$dist]]
    expected <- native_distribution$density((log(7) - eta) / engine$scale)[, 2]
    expect_equal(v$probability, unname(expected), tolerance = 1e-12)
    expect_identical(v$probability_type, "survival_at_saved_horizon")
    expect_null(n$baseline)
    expect_identical(n$score_scale, "negative_log_time_location")
    probabilities <- c(.001, .1, .5, .9, .999)
    expect_equal(sb_nomogram_aft_probability(n$aft, sb_nomogram_aft_score(n$aft, probabilities)),
      probabilities, tolerance = 1e-12)
    p <- sb_plot_nomogram(n)
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    expect_identical(tail(attr(p, "stormbreaker_nomogram_axes")$name, 1), "Survival at 7 month")
    compact <- fit
    compact$diagnostic_artifacts <- NULL
    expect_identical(sb_nomogram(compact, d), n)
    no_horizon <- sb_fit_final(sb_surv_task(d, "patient_id", "time", "event", features, "month"),
      sb_pipeline(learner = id))
    plain <- sb_nomogram(no_horizon, d)
    expect_null(sb_nomogram_values(plain, d)$probability)
    expect_equal(sb_nomogram_values(plain, d)$score, -eta, tolerance = 1e-12)
    expect_identical(tail(attr(sb_plot_nomogram(plain), "stormbreaker_nomogram_axes")$name, 1),
      "Negative log-time location")
    constant <- d
    for (feature in features) constant[[feature]] <- d[[feature]][1]
    flat <- sb_nomogram(fit, constant)
    expect_equal(flat$total_points, 0)
    expect_equal(sb_nomogram_values(flat, constant)$probability, rep(v$probability[1], nrow(d)), tolerance = 1e-12)
    expect_s3_class(ggplot2::ggplot_build(sb_plot_nomogram(flat)), "ggplot_built")
    outside <- d[1, , drop = FALSE]
    outside$gene_a <- max(d$gene_a) + 10
    vv <- sb_nomogram_values(n, outside)
    expect_identical(vv$outside_display_range, TRUE)
    expect_equal(vv$score, -unname(drop(cbind(1, as.matrix(outside[features])) %*% stats::coef(engine))), tolerance = 1e-12)
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
