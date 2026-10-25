test_that("PLS-Cox CV plots the original higher-is-better criterion without capture", {
  skip_if_not_installed("plsRcox")
  task <- make_component_survival_task()
  fit <- sb_fit_final(task, sb_pipeline(learner = "plsrcox_cv", learner_tuning_folds = 3L,
    learner_seed = 8903L))
  expect_identical(fit$status, "success")
  original <- serialize(fit, NULL, version = 2L)
  q <- sb_pls_diagnostics(fit, "cv_curve")
  a <- fit$learner_state$audit
  expect_null(fit$diagnostic_artifacts)
  expect_identical(q$table$iauc, a$cv_error)
  expect_identical(q$table$native_spread, a$cv_standard_error)
  expect_identical(q$table$component, seq.int(0L, a$maximum_components))
  expect_equal(q$table$component[q$table$selected], a$selected_components)
  expect_equal(q$selected_components, which.max(q$table$iauc) - 1L)
  expect_identical(q$folds, 3L)
  p <- sb_plot_algorithm_diagnostic(q)
  expect_match(p$labels$y, "higher is better", fixed = TRUE)
  expect_match(p$labels$caption, "not conventional SE", fixed = TRUE)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  expect_true(identical(serialize(fit, NULL, version = 2L), original))
  # Unavailable components must not be imputed or cause the remaining curve to
  # be reselected; synthetic corruption and missingness exercise this boundary.
  missing <- fit
  missing$learner_state$audit$cv_error[2L] <- NA_real_
  missing$learner_state$audit$cv_standard_error[2L] <- NA_real_
  q <- sb_pls_diagnostics(missing, "cv_curve")
  expect_true(is.na(q$table$iauc[2L]))
  expect_identical(q$selected_components, a$selected_components)
  expect_warning(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), NA)
  bad <- fit
  bad$learner_state$audit$cv_standard_error[1L] <- -1
  expect_error(sb_pls_diagnostics(bad, "cv_curve"), class = "stormbreaker_invalid_diagnostic")
  bad <- fit
  bad$learner_state$audit$cv_error <- 0.5
  expect_error(sb_pls_diagnostics(bad, "cv_curve"), class = "stormbreaker_invalid_diagnostic")
  single <- fit
  single$learner_state$audit$component_selection <- "single_predictor_has_one_feasible_component"
  expect_error(sb_pls_diagnostics(single, "cv_curve"), class = "stormbreaker_unavailable_diagnostic")
  fixed <- sb_fit_final(task, sb_pipeline(learner = "plsrcox"))
  expect_identical(fixed$status, "success")
  expect_error(sb_pls_diagnostics(fixed, "cv_curve"), class = "stormbreaker_unavailable_diagnostic")
})
