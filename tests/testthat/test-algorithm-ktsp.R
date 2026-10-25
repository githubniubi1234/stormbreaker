test_that("k-TSP diagnostics preserve ordered pair rules and tie votes", {
  skip_if_not_installed("switchBox")
  d <- withr::with_seed(1701L, {
    x <- matrix(rnorm(72 * 8), nrow = 72)
    colnames(x) <- paste0("gene_", seq_len(8))
    y <- as.integer(x[, 1] - .7 * x[, 2] + .3 * x[, 3] + rnorm(72, sd = .6) > 0)
    data.frame(patient_id = paste0("NM", seq_len(72)), outcome = y, x)
  })
  features <- paste0("gene_", c(3, 1, 2, 4:8))
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features),
    sb_pipeline(learner = "ktsp_binary"), capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  before <- serialize(fit, NULL, version = 2L)
  e <- fit$learner_state$engine
  view <- sb_ktsp_diagnostics(fit)
  expect_identical(view$table$first, e$pairs[, 1])
  expect_identical(view$table$second, e$pairs[, 2])
  expect_equal(view$table$native_score, e$pair_scores, tolerance = 0)
  expect_equal(view$table$pair, seq_len(e$selected_pairs))
  expect_identical(view$unused_features, setdiff(features, as.vector(e$pairs)))
  expect_identical(view$vote_aggregation, "unweighted_mean")
  expect_equal(unname(view$vote_rule), c(1, .5, 0), tolerance = 0)
  new <- d[1:3, , drop = FALSE]
  for (feature in features) new[[feature]] <- 0
  for (first in e$pairs[, 1]) new[[first]] <- c(1, 0, -1)
  expect_equal(sb_predict_final(fit, new, "patient_id")$predictions$probability, c(1, .5, 0), tolerance = 0)
  compact <- fit
  compact$diagnostic_artifacts <- NULL
  expect_identical(sb_ktsp_diagnostics(compact), view)
  expect_identical(sb_native_diagnostics(compact)$status$status, "available")
  p <- sb_plot_algorithm_diagnostic(view)
  expect_identical(levels(p$data$rule), rev(paste0("Pair ", view$table$pair, ": ", view$table$class_1_rule)))
  expect_equal(ggplot2::ggplot_build(p)$data[[1]]$x, view$table$native_score, tolerance = 0)
  expect_null(p$labels$subtitle)
  expect_null(p$labels$caption)
  expect_identical(serialize(fit, NULL, version = 2L), before)
})

test_that("a multi-pair k-TSP fit keeps every rule and averages rather than weights votes", {
  skip_if_not_installed("switchBox")
  d <- withr::with_seed(1001L, {
    y <- sample(rep(c(0L, 1L), 120L))
    x <- matrix(rnorm(240 * 12), nrow = 240)
    for (j in seq(1, 12, 2)) x[, j] <- x[, j] + .9 * y
    colnames(x) <- paste0("gene_", seq_len(12))
    data.frame(patient_id = paste0("K", seq_len(240)), outcome = y, x)
  })
  features <- paste0("gene_", seq_len(12))
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features), sb_pipeline(learner = "ktsp_binary"))
  view <- sb_ktsp_diagnostics(fit)
  expect_gt(nrow(view$table), 1L)
  expect_equal(view$table$native_score, fit$learner_state$engine$pair_scores, tolerance = 0)
  expect_equal(nrow(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(view))$data[[1]]), nrow(view$table))
  new <- d[1:3, , drop = FALSE]
  for (feature in features) new[[feature]] <- 0
  for (first in view$table$first) new[[first]] <- c(-1, 0, 1)
  new[[view$table$first[1]]][1] <- 1
  expect_equal(sb_predict_final(fit, new, "patient_id")$predictions$probability,
    c(1 / nrow(view$table), .5, 1), tolerance = 1e-12)
})
