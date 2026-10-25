test_that("randomGLM diagnostics preserve original bags and probability averaging", {
  skip_if_not_installed("randomGLM")
  d <- sb_example_binary(n = 160L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  fit <- sb_fit_final(sb_binary_task(d, "patient_id", "outcome", features),
    sb_pipeline(learner = "random_glm_binary"))
  expect_identical(fit$status, "success")
  expect_null(fit$diagnostic_artifacts)
  before <- serialize(fit, NULL, version = 2L)
  state <- fit$learner_state
  counts <- sb_randomglm_diagnostics(fit)
  expected_counts <- vapply(seq_along(features), function(j)
    sum(vapply(state$bag_features, function(index) j %in% index, logical(1))), numeric(1))
  expect_equal(counts$table$bags, expected_counts)
  expect_equal(counts$table$fraction, expected_counts / length(state$intercepts))
  expect_identical(counts$table$feature, features)
  expect_identical(counts$intercepts, state$intercepts)
  coefficients <- sb_randomglm_diagnostics(fit, "bag_coefficients")
  expect_equal(nrow(coefficients$table), length(state$intercepts) * length(features))
  beta <- matrix(coefficients$table$coefficient, nrow = length(state$intercepts))
  for (i in seq_along(state$intercepts)) {
    expect_equal(unname(beta[i, state$bag_features[[i]]]), state$bag_coefficients[[i]])
  }
  score <- sweep(as.matrix(d[features]) %*% t(beta), 2L, state$intercepts, "+")
  expect_equal(rowMeans(plogis(score)), sb_predict_final(fit, d, "patient_id")$predictions$probability,
    tolerance = 1e-12)
  sizes <- sb_randomglm_diagnostics(fit, "bag_size")
  expect_equal(sizes$table$selected_features, unname(lengths(state$bag_features)))
  expect_identical(sizes$table$bag, seq_along(state$bag_features))
  expect_error(sb_nomogram(fit, d), class = "stormbreaker_unsupported_coefficients")
  expect_error(sb_randomglm_diagnostics(fit, "cv_curve"), class = "stormbreaker_unsupported_diagnostic")
  bundle <- sb_native_diagnostics(fit, complete = TRUE)
  expect_identical(bundle$status$status, rep("available", 3L))
  if (requireNamespace("ggplot2", quietly = TRUE)) for (view in bundle$diagnostics) {
    plot <- sb_plot_algorithm_diagnostic(view)
    expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
    expect_null(plot$labels$subtitle)
    expect_null(plot$labels$caption)
  }
  expect_identical(serialize(fit, NULL, version = 2L), before)

  synthetic <- fit
  synthetic$learner_state$intercepts <- c(0, 0, 0)
  synthetic$learner_state$bag_features <- list(integer(), 2L, 1L)
  synthetic$learner_state$bag_coefficients <- list(numeric(), 0, 1)
  synthetic$learner_state$params$n_bags <- 3L
  expect_equal(sb_randomglm_diagnostics(synthetic)$table$bags, c(1, 1, 0))
  expect_equal(sb_randomglm_diagnostics(synthetic, "bag_size")$table$selected_features, c(0, 1, 1))
  tab <- sb_randomglm_diagnostics(synthetic, "bag_coefficients")$table
  expect_true(tab$included[tab$bag == 2L & tab$feature == features[2L]])
  expect_equal(tab$coefficient[tab$bag == 2L & tab$feature == features[2L]], 0)
  expect_false(any(tab$included[tab$bag == 1L]))
})
