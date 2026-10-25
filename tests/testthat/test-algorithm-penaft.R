test_that("penAFT diagnostics preserve native CV and input-unit path predictions", {
  skip_if_not_installed("penAFT")
  d <- sb_example_survival(100L, seed = 8612L)
  d$gene_c <- d$gene_c * 2 + 4
  features <- c("gene_c", "gene_a", "gene_b")
  x <- as.matrix(d[features])
  fold <- sb_make_stratified_fold_id(d$event, 5L, 1L, stratum_label = "event")
  invisible(capture.output(native <- withr::with_seed(1L, penAFT::penAFT.cv(
    X = x, logY = log(d$time), delta = d$event, nlambda = 50L,
    lambda.ratio.min = .1, penalty = "EN", alpha = .5, tol.abs = 1e-8,
    tol.rel = 2.5e-4, standardize = TRUE,
    cv.index = lapply(1:5, function(i) which(fold == i)), admm.max.iter = 10000L, quiet = TRUE))))
  fit <- sb_fit_final(sb_surv_task(d, "patient_id", "time", "event", features, "month"),
    sb_pipeline(learner = "penaft_gehan_elastic_net"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  cv <- sb_penaft_diagnostics(fit, "cv_curve")
  path <- sb_penaft_diagnostics(fit, "coefficient_path")
  co <- sb_penaft_diagnostics(fit)
  f <- native$full.fit
  expect_equal(cv$table$lambda, f$lambda, tolerance = 0)
  expect_equal(cv$table$score, native$cv.err.linPred, tolerance = 0)
  expect_equal(which(cv$table$selected), which.min(native$cv.err.linPred))
  expect_equal(path$native_path$full_fit$beta, f$beta, tolerance = 0)
  raw <- sweep(as.matrix(f$beta), 1, f$X.sd, "/")
  expect_equal(path$table$coefficient, as.numeric(t(raw)), tolerance = 0)
  expect_equal(path$table$feature, rep(features, each = length(f$lambda)))
  expect_equal(path$intercept_path, -as.numeric(crossprod(f$X.mean, raw)), tolerance = 0)
  native_lp <- scale(x, center = f$X.mean, scale = f$X.sd) %*% as.matrix(f$beta)
  replay <- sweep(x %*% raw, 2, path$intercept_path, "+")
  expect_equal(unname(replay), unname(native_lp), tolerance = 1e-12)
  j <- which.min(native$cv.err.linPred)
  expect_equal(co$table$coefficient, as.numeric(raw[, j]), tolerance = 0)
  expect_equal(co$intercept, path$intercept_path[j], tolerance = 0)
  expect_null(fit$learner_state$diagnostic_path)
  compact <- fit
  compact$diagnostic_artifacts <- NULL
  expect_identical(sb_penaft_diagnostics(compact), co)
  expect_identical(sb_native_diagnostics(compact, complete = TRUE)$status$status,
    c("available", "missing_capture", "missing_capture"))
  bundle <- sb_native_diagnostics(fit, complete = TRUE)
  expect_identical(bundle$status$status, rep("available", 3L))
  for (v in bundle$diagnostics) {
    p <- sb_plot_algorithm_diagnostic(v)
    b <- ggplot2::ggplot_build(p)
    expect_s3_class(b, "ggplot_built")
    expect_null(p$labels$caption)
    expect_null(p$labels$subtitle)
    if (v$type == "penaft_coefficients") expect_equal(b$data[[2]]$x, v$table$coefficient, tolerance = 0)
    if (v$type == "penaft_cv_curve") expect_equal(b$data[[2]]$x, log10(f$lambda[j]), tolerance = 0)
    if (v$type == "penaft_coefficient_path") expect_equal(sort(b$data[[1]]$y), sort(v$table$coefficient), tolerance = 0)
  }
  expect_identical(serialize(fit, NULL, version = 2L), before)
})
