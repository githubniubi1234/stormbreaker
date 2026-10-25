test_that("Naive Bayes curves preserve fitted densities and posterior decomposition", {
  skip_if_not_installed("klaR")
  d <- sb_example_binary(n = 100L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  for (id in c("gaussian_naive_bayes", "kernel_naive_bayes_binary")) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    q <- sb_naive_bayes_diagnostics(fit)
    expect_null(fit$diagnostic_artifacts)
    expect_identical(q$features, features)
    p <- sb_plot_algorithm_diagnostic(q)
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    expect_null(p$labels$caption)
    expect_null(p$labels$subtitle)
    expect_identical(p$facet$params$free, list(x = TRUE, y = TRUE))
    if (id == "gaussian_naive_bayes") {
      m <- fit$learner_state$payload
      for (feature in features) for (class in c("0", "1")) {
        tab <- q$table[q$table$feature == feature & q$table$class == class, ]
        expect_equal(tab$value, dnorm(tab$x, m$means[class, feature], sqrt(m$variances[class, feature])))
      }
      ratio <- sb_naive_bayes_diagnostics(fit, "log_likelihood_ratio")
      expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(ratio)), "ggplot_built")
      expect_identical(sb_plot_algorithm_diagnostic(ratio)$facet$params$free, list(x = TRUE, y = FALSE))
      # Pick one stored grid point per feature to reconstruct the exact model
      # posterior, including class priors and the original variance smoothing.
      row <- setNames(as.data.frame(as.list(vapply(features, function(f) ratio$table$x[ratio$table$feature == f][40L], numeric(1)))), features)
      terms <- vapply(features, function(f) ratio$table$value[ratio$table$feature == f][40L], numeric(1))
      original <- sb_predict_naive_bayes_learner_state(fit$learner_state, row)
      expect_equal(plogis(ratio$prior_log_odds + sum(terms)), as.numeric(original), tolerance = 1e-12)
    } else {
      e <- unserialize(fit$learner_state$engine_raw)
      for (feature in features) for (class in c("0", "1")) {
        tab <- q$table[q$table$feature == feature & q$table$class == class, ]
        expect_identical(tab$x, e$tables[[feature]][[class]]$x)
        expect_identical(tab$value, e$tables[[feature]][[class]]$y)
      }
      expect_error(sb_naive_bayes_diagnostics(fit, "log_likelihood_ratio"), class = "stormbreaker_unavailable_diagnostic")
    }
    single <- sb_naive_bayes_diagnostics(fit, features = features[2L])
    expect_identical(unique(single$table$feature), features[2L])
    expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(single)), "ggplot_built")
    expect_error(sb_naive_bayes_diagnostics(fit, features = "missing"), class = "stormbreaker_invalid_diagnostic")
    expect_true(identical(serialize(fit, NULL, version = 2L), before))
  }
})
