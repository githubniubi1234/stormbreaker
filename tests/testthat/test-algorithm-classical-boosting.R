test_that("classical boosting views retain original ensembles and exact native counts", {
  skip_if_not_installed("fastAdaboost")
  skip_if_not_installed("caTools")
  d <- make_sb_binary_fixture(n = 100L, seed = 108L)
  d$unused <- 0
  features <- c("gene_c", "unused", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  for (id in sb_classical_binary_boosting_learner_ids()) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    original <- serialize(fit, NULL)
    engine <- unserialize(fit$learner_state$model_raw)
    usage <- sb_boosting_diagnostics(fit, "split_usage")
    expect_identical(usage$native_table, fit$learner_state$feature_importance)
    selected <- if (id == "adaboost_binary") unlist(lapply(engine$trees,
      function(tree) tree$frame$var[tree$frame$var != "<leaf>"])) else features[engine$Stump[, 1]]
    expect_equal(usage$table$count, as.numeric(tabulate(match(selected, features), nbins = length(features))))
    expect_equal(usage$table$count[usage$table$feature == "unused"], 0)
    if (id == "adaboost_binary") {
      weights <- sb_boosting_diagnostics(fit, "tree_weights")
      expect_identical(weights$native_table, engine$weights)
      expect_equal(weights$table$weight, as.numeric(engine$weights))
      sizes <- sb_boosting_diagnostics(fit, "tree_size")
      expect_equal(sizes$table$terminal_nodes, unname(vapply(engine$trees,
        function(tree) as.integer((nrow(tree$frame) + 1L) / 2L), integer(1))))
      expect_equal(nrow(weights$table), length(engine$trees))
      plots <- list(usage, weights, sizes)
    } else {
      stumps <- sb_boosting_diagnostics(fit, "stump_sequence")
      expect_identical(stumps$native_table, engine$Stump)
      expect_equal(stumps$table$feature, features[engine$Stump[, 1]])
      expect_equal(stumps$table$threshold, as.numeric(engine$Stump[, 2]))
      # Reconstruct class-1 probabilities independently from the returned rules.
      values <- d[features]
      margin <- numeric(nrow(values))
      for (i in seq_len(nrow(stumps$table))) {
        rule <- stumps$table[i, ]
        margin <- margin + ifelse(values[[rule$feature]] <= rule$threshold,
          rule$contribution_at_or_below, rule$contribution_above)
      }
      expected <- stats::predict(engine, as.matrix(values), type = "raw")[, "1"]
      expect_equal(stats::plogis(margin), unname(expected))
      plots <- list(usage, stumps)
    }
    if (requireNamespace("ggplot2", quietly = TRUE)) for (q in plots) expect_s3_class(
      ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
    err <- tryCatch(sb_boosting_diagnostics(fit, "training_loss"), error = identity)
    expect_s3_class(err, "stormbreaker_unavailable_diagnostic")
    expect_identical(serialize(fit, NULL), original)
  }
})

test_that("AdaBoost early convergence does not produce invented extra trees", {
  skip_if_not_installed("fastAdaboost")
  y <- rep(c(0L, 1L), each = 10L)
  d <- data.frame(id = seq_along(y), outcome = y, separating = ifelse(y == 1L, 1, -1), noise = seq_along(y))
  fit <- sb_fit_final(sb_binary_task(d, "id", "outcome", c("separating", "noise")),
    sb_pipeline(learner = "adaboost_binary", adaboost_iterations = 5L))
  expect_identical(fit$status, "success")
  q <- sb_boosting_diagnostics(fit, "tree_weights")
  expect_equal(q$fitted_iterations, 1L)
  expect_equal(q$requested_iterations, 5L)
  expect_equal(nrow(q$table), 1L)
  expect_equal(q$table$weight, as.numeric(unserialize(fit$learner_state$model_raw)$weights))
  if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(
    ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
})
