test_that("BART posterior views use original draws and original score scales", {
  skip_if_not_installed("BART")
  for (id in c("bart_binary", "bart_survival", "bart_aft_survival")) {
    d <- if (id == "bart_binary") sb_example_binary(100L, seed = 108L) else sb_example_survival(100L, seed = 8612L)
    features <- c("gene_c", "gene_a", "gene_b")
    task <- if (id == "bart_binary") sb_binary_task(d, "patient_id", "outcome", features) else
      sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    query <- d[c(1L, 7L), ]
    query$gene_c <- query$gene_c + .17
    z <- sb_bart_posterior(fit, query)
    e <- if (is.null(fit$learner_state$engine_raw)) fit$learner_state$engine else unserialize(fit$learner_state$engine_raw)
    matrix_x <- as.matrix(query[features])
    if (id == "bart_aft_survival") {
      invisible(capture.output(native_draws <- BART::pwbart(matrix_x[, e$kept_columns, drop = FALSE],
        e$treedraws, mu = e$offset, mc.cores = 1L, dodraws = TRUE)))
    } else {
      if (id == "bart_survival") matrix_x <- cbind(t = rep(e$times, 2),
        matrix_x[rep(1:2, each = length(e$times)), , drop = FALSE])
      if (length(e$rm.const)) matrix_x <- matrix_x[, e$rm.const, drop = FALSE]
      fn <- get(if(id == "bart_binary") "predict.pbart" else "predict.survbart", asNamespace("BART"))
      invisible(capture.output(original <- fn(e, matrix_x, mc.cores = 1L, openmp = FALSE)))
      native_draws <- if(id == "bart_binary") original$prob.test else original$surv.test
    }
    expect_equal(z$draws, native_draws, tolerance = 1e-12)
    expect_equal(z$table$mean, unname(colMeans(z$draws)))
    expect_equal(z$table$lower, unname(apply(z$draws, 2L, quantile, .025)))
    expect_equal(z$table$upper, unname(apply(z$draws, 2L, quantile, .975)))
    if (id == "bart_binary") {
      expect_equal(z$table$mean, sb_predict_final(fit, query, "patient_id")$predictions$probability, tolerance = 1e-12)
    } else if (id == "bart_aft_survival") {
      expect_equal(-z$table$mean, sb_predict_final(fit, query, "patient_id")$predictions$risk, tolerance = 1e-12)
    } else {
      model <- sb_restore_additional_survival_model(fit$learner_state, id)
      native <- sb_bart_survival_curves(model, query[features])
      expect_equal(z$table$mean, as.numeric(t(native$survival)), tolerance = 1e-12)
      expect_equal(z$table$time, rep(native$time, 2))
    }
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    expect_identical(sb_native_diagnostics(fit)$status$status, "unavailable")
    expect_identical(sb_native_diagnostics(fit, newdata = query)$status$status, "available")
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})
