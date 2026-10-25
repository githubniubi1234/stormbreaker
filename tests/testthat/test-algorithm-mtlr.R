test_that("MTLR diagnostics preserve native time weights and exact input-unit potentials", {
  skip_if_not_installed("MTLR")
  d <- sb_example_survival(160L, seed = 8612L)
  features <- c("gene_c", "gene_a", "gene_b")
  d$gene_c <- d$gene_c * 2 + 4
  fit <- sb_fit_final(sb_surv_task(d, "patient_id", "time", "event", features, "month"),
    sb_pipeline(learner = "mtlr_survival"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  e <- fit$learner_state$engine
  view <- sb_mtlr_diagnostics(fit)
  expect_equal(view$native_weights, e$weight_matrix, tolerance = 0)
  expect_equal(view$table$weight, as.numeric(e$weight_matrix[, -1]), tolerance = 0)
  expect_equal(view$time_points, e$time_points, tolerance = 0)
  expect_equal(view$table$time, rep(e$time_points, length(features)), tolerance = 0)
  expect_identical(view$normalization$feature, features)
  expect_equal(view$normalization$center, unname(e$scale$center), tolerance = 0)
  expect_equal(view$normalization$sd, unname(e$scale$sd), tolerance = 0)
  expect_equal(view$requested_intervals, ceiling(sqrt(nrow(d))))
  expect_equal(view$effective_time_points, nrow(e$weight_matrix))
  expect_lt(view$effective_time_points, view$requested_intervals)
  x <- as.matrix(d[features])
  normalized <- scale(x, center = e$scale$center, scale = e$scale$sd)
  expect_equal(unname(cbind(1, x) %*% t(view$original_unit_weights)),
    unname(cbind(1, normalized) %*% t(e$weight_matrix)), tolerance = 1e-12)
  newdata <- d[features]
  names(newdata) <- fit$learner_state$internal_features
  native_predict <- get("predict.mtlr", asNamespace("MTLR"))
  raw_engine <- e
  raw_engine$weight_matrix <- view$original_unit_weights
  raw_engine$scale <- NULL
  expect_equal(native_predict(raw_engine, newdata, type = "survivalcurve", add_zero = FALSE),
    native_predict(e, newdata, type = "survivalcurve", add_zero = FALSE), tolerance = 1e-12)
  biases <- sb_mtlr_diagnostics(fit, "biases")
  expect_equal(biases$table$bias, as.numeric(e$weight_matrix[, 1]), tolerance = 0)
  compact <- fit
  compact$diagnostic_artifacts <- NULL
  expect_identical(sb_mtlr_diagnostics(compact), view)
  # Synthetic zero-SD state exercises the native ignored-column prediction branch.
  zero <- compact
  zero$learner_state$engine$scale$sd[1] <- 0
  zero_view <- sb_mtlr_diagnostics(zero)
  expect_equal(as.numeric(zero_view$original_unit_weights[, 2]), rep(0, length(e$time_points)))
  zero_raw <- zero$learner_state$engine
  zero_raw$weight_matrix <- zero_view$original_unit_weights
  zero_raw$scale <- NULL
  expect_equal(native_predict(zero_raw, newdata, type = "survivalcurve", add_zero = FALSE),
    native_predict(zero$learner_state$engine, newdata, type = "survivalcurve", add_zero = FALSE), tolerance = 1e-12)
  bundle <- sb_native_diagnostics(compact, complete = TRUE)
  expect_identical(bundle$status$status, rep("available", 2L))
  for (q in bundle$diagnostics) {
    p <- sb_plot_algorithm_diagnostic(q)
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    expect_null(p$labels$subtitle)
    expect_null(p$labels$caption)
  }
  expect_identical(serialize(fit, NULL, version = 2L), before)
})
