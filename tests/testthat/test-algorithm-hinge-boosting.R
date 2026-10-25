test_that("captured hinge paths reproduce native updates without changing portable state", {
  skip_if_not_installed("bst")
  d <- make_sb_binary_fixture(n = 100L, seed = 108L)
  features <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  pipeline <- sb_pipeline(learner = "bst_hinge_binary")
  compact <- sb_fit_final(task, pipeline)
  fit <- sb_fit_final(task, pipeline, capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  expect_identical(compact$status, "success")
  expect_identical(fit$learner_state, compact$learner_state)
  expect_null(fit$learner_state$engine$diagnostic_path)
  err <- tryCatch(sb_boosting_diagnostics(compact, "training_loss"), error = identity)
  expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
  original <- serialize(fit, NULL)
  params <- pipeline$learner_params
  native <- bst::bst(as.matrix(d[features]), ifelse(d$outcome == 1, 1, -1),
    cost = params$cost, family = params$family, learner = params$base_learner,
    ctrl = bst::bst_control(mstop = params$iterations, nu = params$learning_rate,
      center = params$center, intercept = params$fit_intercept, trace = FALSE))
  loss <- sb_boosting_diagnostics(fit, "training_loss")
  usage <- sb_boosting_diagnostics(fit, "component_usage")
  path <- sb_boosting_diagnostics(fit, "coefficient_path")
  expect_equal(loss$table$loss, native$risk)
  expect_equal(loss$table$iteration, seq_len(params$iterations))
  expect_equal(usage$table$updates, tabulate(native$ensemble, nbins = length(features)))
  expect_equal(path$selections$native_feature_index, native$ensemble)
  expect_equal(path$selections$native_update, native$coef)
  expect_equal(path$table$feature[path$table$iteration == 0], c("(Intercept)", features))
  reconstructed_loss <- vapply(seq_len(params$iterations), function(i) {
    beta <- path$table$coefficient[path$table$iteration == i]
    margin_i <- beta[1] + as.numeric(as.matrix(d[features]) %*% beta[-1])
    mean(bst::loss(ifelse(d$outcome == 1, 1, -1), margin_i, cost = params$cost, family = "hinge"))
  }, numeric(1))
  expect_equal(reconstructed_loss, native$risk, tolerance = 1e-12)
  endpoint <- path$table[path$table$iteration == params$iterations, ]$coefficient
  margin <- endpoint[1] + as.numeric(as.matrix(d[features]) %*% endpoint[-1])
  expect_equal(margin, as.numeric(native$yhat), tolerance = 1e-12)
  expect_equal(stats::plogis(2 * margin), sb_predict_omics_specialist_binary_state(
    fit$learner_state, d[features])$probability, tolerance = 1e-12)
  if (requireNamespace("ggplot2", quietly = TRUE)) for (q in list(loss, usage, path)) expect_s3_class(
    ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
  expect_identical(serialize(fit, NULL), original)
})
