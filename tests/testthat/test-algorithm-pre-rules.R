test_that("PRE diagnostics preserve native terms, predictions and coefficient CV", {
  skip_if_not_installed("pre")
  skip_if_not_installed("mboost")
  for (kind in c("binary", "survival")) {
    d <- make_sb_binary_fixture(n = 150L, seed = 108L)
    d$time <- exp(2 + .8 * d$gene_a - .4 * d$gene_b)
    d$event <- rep(c(1L, 1L, 0L), length.out = nrow(d))
    features <- c("gene_a", "gene_b", "gene_c")
    task <- if (kind == "binary") sb_binary_task(d, "patient_id", "outcome", features) else
      sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = paste0("pre_", kind, "_rule_ensemble")))
    expect_identical(fit$status, "success")
    expect_null(fit$diagnostic_artifacts)
    original <- fit
    q <- sb_rule_ensemble(fit)
    engine <- unserialize(fit$learner_state$engine_raw)
    native <- stats::coef(engine, penalty.par.val = fit$learner_state$params$penalty_rule)
    rownames(native) <- NULL
    expect_identical(q$table[names(native)], native)
    expect_equal(q$table$selected, native$coefficient != 0)
    expect_equal(q$table$term_type[q$table$rule %in% engine$rules$rule],
      rep("rule", sum(q$table$rule %in% engine$rules$rule)))
    # These expressions come only from this synthetic native fit, not user text.
    reconstructed <- Reduce(`+`, lapply(seq_len(nrow(q$table)), function(i) {
      value <- if (q$table$term_type[i] == "linear") {
        pmin(pmax(d[[q$table$rule[i]]], q$table$winsor_lower[i]), q$table$winsor_upper[i])
      } else eval(parse(text = q$table$description[i]), envir = d)
      q$table$coefficient[i] * value
    }))
    native_link <- as.numeric(stats::predict(engine, newdata = d[features], type = "link"))
    expect_equal(as.numeric(reconstructed), native_link, tolerance = 1e-10)
    cv <- sb_rule_ensemble(fit, "rule_cv_curve")
    expect_equal(cv$table$lambda, engine$glmnet.fit$lambda)
    expect_equal(cv$table$mean_value, unname(engine$glmnet.fit$cvm))
    expect_equal(cv$table$standard_error, unname(engine$glmnet.fit$cvsd))
    expect_equal(cv$table$nonzero, unname(engine$glmnet.fit$nzero))
    expect_equal(cv$table$lambda[cv$table$selected], engine$glmnet.fit$lambda.1se)
    expect_identical(fit, original)
    expect_equal(length(grepRaw("B001", serialize(q, NULL), fixed = TRUE)), 0)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      for (z in list(q, cv)) {
        p <- sb_plot_algorithm_diagnostic(z)
        expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
      }
      p <- sb_plot_algorithm_diagnostic(q)
      expect_equal(nrow(p$data), sum(q$table$selected & q$table$term_type != "intercept"))
      empty <- q
      empty$table$coefficient[] <- 0
      empty$table$selected[] <- FALSE
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(empty)), "ggplot_built")
    }
    bad <- fit
    engine$glmnet.fit$cvsd[1] <- -1
    bad$learner_state$engine_raw <- serialize(engine, NULL)
    err <- tryCatch(sb_rule_ensemble(bad, "rule_cv_curve"), error = identity)
    expect_s3_class(err, "stormbreaker_invalid_diagnostic")
    err <- tryCatch(sb_rule_ensemble(fit, "surrogate_tree"), error = identity)
    expect_s3_class(err, "stormbreaker_unsupported_diagnostic")
  }
})
