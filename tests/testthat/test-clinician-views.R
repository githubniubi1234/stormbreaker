test_that("clinical assembly preserves saved DCA and native data without fitting", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_binary(n = 100L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"),
    decision_thresholds = c(.1, .3, .6))
  pipeline <- sb_pipeline(learner = "gaussian_naive_bayes")
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 41L)
  expect_identical(result$status, "success")
  fit <- sb_fit_final(task, pipeline)
  native <- sb_naive_bayes_diagnostics(fit)
  before <- serialize(result, NULL, version = 2L)
  views <- sb_clinician_views(result, c(.2, .5, .8), .5, native = list(Density = native), detail = "full")
  expect_identical(views$plots[["Decision curve"]]$data, sb_plot_dca(result)$data)
  expect_identical(views$plots[["Model: Density"]]$data, sb_plot_algorithm_diagnostic(native)$data)
  for (plot in views$plots) {
    expect_null(plot$labels$subtitle)
    expect_null(plot$labels$caption)
  }
  expect_equal(sort(unique(views$plots[["Threshold performance"]]$data$threshold)), c(.2, .5, .8))
  expect_identical(serialize(result, NULL, version = 2L), before)
  expect_equal(nrow(views$status), 8L)
  expect_false(any(views$status$status == "error"))
  expect_true(all(c("Clinical decisions", "Model diagnostics") %in% views$groups))
  expect_true(anyNA(views$plots[["Threshold performance"]]$data$value))
  file <- tempfile(fileext = ".html")
  sb_render_clinician_report(views, file)
  html <- paste(readLines(file), collapse = "\n")
  expect_identical(grepl("View availability", html, fixed = TRUE), FALSE)
  expect_match(html, "Model diagnostics", fixed = TRUE)
  expect_true(grepl("original Gaussian means", html, fixed = TRUE))
  expect_error(sb_clinician_views(result, .2, .5), class = "stormbreaker_bad_report_input")
  fit_before <- serialize(fit, NULL, version = 2L)
  local_mocked_bindings(sb_fit_final = function(...) stop("Report must not fit models"))
  requested <- sb_clinician_views(result, native = sb_native_diagnostics(fit))
  model_plots <- requested$plots[requested$groups == "Model diagnostics"]
  expect_length(model_plots, 1L)
  expect_identical(model_plots[[1L]]$data, sb_plot_algorithm_diagnostic(native)$data)
  detailed <- sb_clinician_views(result, native = sb_native_diagnostics(fit, complete = TRUE))
  expect_equal(sum(detailed$groups == "Model diagnostics"), 2L)
  expect_identical(serialize(fit, NULL, version = 2L), fit_before)
  expect_identical(serialize(result, NULL, version = 2L), before)
  result$predictions$probability[1L] <- NA_real_
  partial <- sb_clinician_views(result)
  expect_identical(partial$status$status[partial$status$view == "ROC"], "unavailable")
  expect_match(partial$status$reason[partial$status$view == "ROC"], "invalid_predictions", fixed = TRUE)
})

test_that("clinical figures move explanations to captions without changing their evidence", {
  skip_if_not_installed("ggplot2")
  result <- make_sb_curve_view_fixture("survival")$cv
  original <- sb_plot_calibration(result) + ggplot2::labs(
    subtitle = "Saved cohort <A>", caption = "Intervals use saved estimates.")
  local_mocked_bindings(sb_plot_calibration = function(...) original)
  views <- sb_clinician_views(result)
  plot <- views$plots[["Calibration"]]
  expect_identical(plot$data, original$data)
  expect_identical(plot$layers, original$layers)
  expect_identical(plot$labels$title, original$labels$title)
  expect_identical(views$captions[match("Calibration", names(views$plots))],
    "Saved cohort <A>\nIntervals use saved estimates.")
  expect_null(plot$labels$subtitle)
  expect_null(plot$labels$caption)
  file <- tempfile(fileext = ".html")
  sb_render_clinician_report(views, file)
  html <- paste(readLines(file), collapse = "\n")
  expect_match(html, "Saved cohort &lt;A&gt;", fixed = TRUE)
  expect_match(html, "Intervals use saved estimates.", fixed = TRUE)
})

test_that("explicit query diagnostics can be added without automatic extraction", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_binary(n = 60L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "knn_binary"))
  result <- make_sb_curve_view_fixture("binary")$cv
  missing <- sb_clinician_views(result, native = sb_native_diagnostics(fit))
  expect_identical(missing$status$status[missing$status$group == "Model diagnostics"], "unavailable")
  query <- d[c(2L, 9L), ]
  supplied <- sb_clinician_views(result, native = sb_native_diagnostics(fit, newdata = query))
  plot <- supplied$plots[supplied$groups == "Model diagnostics"][[1L]]
  expect_identical(plot$data, sb_plot_algorithm_diagnostic(sb_knn_neighborhood(fit, query))$data)
})

test_that("default and expanded reports never extract or plot technical diagnostics", {
  skip_if_not_installed("ggplot2")
  result <- make_sb_curve_view_fixture("binary")$cv
  unexpected <- function(...) stop("Technical diagnostics were not requested")
  local_mocked_bindings(sb_native_diagnostics = unexpected,
    sb_selector_diagnostics = unexpected, sb_plot_algorithm_diagnostic = unexpected,
    sb_plot_selector_diagnostic = unexpected)
  views <- sb_clinician_views(result)
  expect_equal(sum(views$groups %in% c("Model diagnostics", "Feature selection")), 0L)
  for (complete in c(FALSE, TRUE)) {
    file <- tempfile(fileext = ".html")
    sb_render_clinician_report(views, file, complete = complete)
    html <- paste(readLines(file), collapse = "\n")
    expect_identical(grepl("Model diagnostics|Feature selection", html), FALSE)
  }
})

test_that("survival assembly does not reinterpret censoring as binary outcomes", {
  skip_if_not_installed("ggplot2")
  result <- make_sb_curve_view_fixture("survival")$cv
  views <- sb_clinician_views(result)
  expect_identical(views$status$status[views$status$view == "ROC"], "unavailable")
  expect_null(views$plots[["Probability distribution"]])
  expect_identical(views$plots[["Decision curve"]]$data, sb_plot_dca(result)$data)
  expect_identical(views$plots[["Calibration"]]$data, sb_plot_calibration(result)$data)
})

test_that("missing clinical curves remain explicit and are never synthesized", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_binary(n = 80L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  result <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3L, seed = 41L)
  views <- sb_clinician_views(result)
  expect_identical(views$status$status[views$status$view == "Decision curve"], "unavailable")
  expect_null(views$plots[["Decision curve"]])
  expect_true(length(views$diagnostics[["Decision curve"]]) > 0L)
})
