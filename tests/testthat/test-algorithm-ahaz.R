test_that("ahaz diagnostics preserve native input-unit coefficients and complete BIC paths", {
  skip_if_not_installed("ahaz")
  d <- sb_example_survival(160L, seed = 8612L)
  d$gene_c <- d$gene_c * 2 + 4
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
  for (id in sb_ahaz_survival_learner_ids()) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    x <- as.matrix(d[features])
    colnames(x) <- paste0("x", seq_along(features))
    native <- if (id == "ahaz_lin_ying") {
      ahaz::ahaz(survival::Surv(d$time, d$event), x, robust = FALSE)
    } else ahaz::tune.ahazpen(survival::Surv(d$time, d$event), x,
      standardize = TRUE, penalty = if (id == "ahaz_lasso_bic")
        ahaz::lasso.control(alpha = 1) else ahaz::sscad.control(a = 3.7, nsteps = 1L),
      tune = ahaz::bic.control(), dfmax = 3L, nlambda = 100L)
    co <- sb_ahaz_diagnostics(fit)
    expect_identical(co$table$feature, features)
    expect_equal(co$table$coefficient, as.numeric(stats::coef(native)), tolerance = 1e-12)
    expect_equal(as.numeric(x %*% co$table$coefficient),
      as.numeric(x %*% stats::coef(native)), tolerance = 1e-12)
    expect_null(fit$learner_state$diagnostic_path)
    compact <- fit
    compact$diagnostic_artifacts <- NULL
    expect_identical(sb_ahaz_diagnostics(compact), co)
    if (id != "ahaz_lin_ying") {
      bic <- sb_ahaz_diagnostics(fit, "bic_path")
      path <- sb_ahaz_diagnostics(fit, "coefficient_path")
      expect_equal(bic$table$lambda, native$lambda, tolerance = 0)
      expect_equal(bic$table$criterion, native$tunem, tolerance = 0)
      expect_equal(bic$table$df, native$df, tolerance = 0)
      expect_equal(bic$table$lambda[bic$table$selected], native$lambda.min, tolerance = 0)
      expect_equal(path$table$coefficient, as.numeric(t(as.matrix(native$fit$beta))), tolerance = 0)
      expect_equal(path$table$feature, rep(features, each = length(native$lambda)))
      expect_equal(path$table$lambda, rep(native$lambda, length(features)), tolerance = 0)
      expect_equal(nrow(path$table), 3L * length(native$lambda))
      expect_identical(sb_native_diagnostics(compact, complete = TRUE)$status$status,
        c("available", "missing_capture", "missing_capture"))
    }
    bundle <- sb_native_diagnostics(fit, complete = TRUE)
    expect_identical(bundle$status$status, rep("available", if (id == "ahaz_lin_ying") 1L else 3L))
    for (view in bundle$diagnostics) {
      p <- sb_plot_algorithm_diagnostic(view)
      built <- ggplot2::ggplot_build(p)
      expect_s3_class(built, "ggplot_built")
      expect_null(p$labels$caption)
      expect_null(p$labels$subtitle)
      if (view$type == "ahaz_coefficients") expect_equal(built$data[[2]]$x, view$table$coefficient, tolerance = 0)
      if (view$type == "ahaz_bic_path") expect_equal(built$data[[2]]$x,
        log10(view$table$lambda[view$table$selected]), tolerance = 0)
      if (view$type == "ahaz_coefficient_path") expect_equal(sort(built$data[[1]]$y), sort(view$table$coefficient), tolerance = 0)
    }
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
