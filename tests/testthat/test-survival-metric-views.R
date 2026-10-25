test_that("horizon metric tables preserve native CV and external values", {
  fixture <- make_sb_curve_view_fixture("survival", comparison = TRUE)
  before <- fixture
  for (name in c("cv", "external")) {
    x <- fixture[[name]]
    results <- if (name == "cv") x$results else x$evaluations
    for (metric in c("time_auc", "uno_c", "ipcw_brier")) {
      tab <- sb_result_table(x, metric = metric)
      expect_identical(tab$estimate, unname(vapply(results, function(z) z$metrics[[metric]]$estimate, numeric(1))))
      expect_identical(tab$status, unname(vapply(results, function(z) z$metrics[[metric]]$status, character(1))))
      expect_identical(tab$scope, unname(vapply(results, function(z) z$metrics[[metric]]$scope, character(1))))
      expect_identical(tab$aggregation, unname(vapply(results, function(z) z$metrics[[metric]]$aggregation, character(1))))
      expect_identical(tab$prediction_horizon, rep(6, 2))
      expect_identical(tab$time_unit, rep("month", 2))
      expect_identical(tab$metric, rep(metric, 2))
    }
  }
  expect_identical(fixture, before)
})

test_that("horizon-free requests and unavailable saved results stay distinct", {
  d <- sb_example_survival(60, seed = 1)
  task <- sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month")
  x <- sb_run_cv(task, folds = 3, seed = 1)
  for (metric in c("time_auc", "uno_c", "ipcw_brier")) {
    tab <- sb_result_table(x, metric = metric)
    expect_identical(tab$status, "not_requested")
    expect_identical(tab$estimate, NA_real_)
    expect_identical(tab$prediction_horizon, NA_real_)
    expect_identical(tab$time_unit, "month")
  }
  x$task$prediction_horizon <- 6
  expect_identical(sb_result_table(x, metric = "time_auc")$status, "not_saved")
  x$status <- "failed"
  expect_identical(sb_result_table(x, metric = "time_auc")$run_status, "failed")
})

test_that("horizon views never fit or recompute metrics", {
  x <- make_sb_curve_view_fixture("survival")
  before <- x
  local_mocked_bindings(
    sb_run_cv = function(...) stop("must not run CV"),
    sb_fit_final = function(...) stop("must not fit"),
    sb_evaluate_external = function(...) stop("must not evaluate"),
    sb_survival_time_auc_metric = function(...) stop("must not recompute"),
    sb_survival_uno_c_metric = function(...) stop("must not recompute"),
    sb_survival_ipcw_brier_metric = function(...) stop("must not recompute")
  )
  for (metric in c("time_auc", "uno_c", "ipcw_brier")) {
    expect_identical(sb_result_table(x$cv, metric = metric)$estimate, x$cv$metrics[[metric]]$estimate)
    expect_identical(sb_result_table(x$external, metric = metric)$estimate, x$external$metrics[[metric]]$estimate)
  }
  expect_identical(x, before)
})

test_that("horizon plots retain values above one and distinct targets", {
  skip_if_not_installed("ggplot2")
  a <- make_sb_curve_view_fixture("survival")$cv
  a$metrics$ipcw_brier$estimate <- 1.25
  b <- a
  b$pipeline$id <- "second_target"
  b$task$prediction_horizon <- 12
  b$metrics$ipcw_brier$prediction_horizon <- 12
  b$metrics$ipcw_brier$estimate <- 0.5
  x <- structure(list(results = list(a, b)), class = "stormbreaker_cv_comparison")
  before <- x
  p <- sb_plot_performance(x, "ipcw_brier")
  built <- ggplot2::ggplot_build(p)
  expect_equal(sort(built$data[[1]]$y), c(0.5, 1.25))
  expect_identical(p$labels$y, "IPCW Brier score")
  expect_identical(p$data$prediction_horizon, c(6, 12))
  labels <- p$scales$get_scales("x")$limits
  expect_identical(any(grepl("6 month", labels, fixed = TRUE)), TRUE)
  expect_identical(any(grepl("12 month", labels, fixed = TRUE)), TRUE)
  expect_identical(x, before)
})

test_that("horizon plots omit unavailable methods while retaining diagnostics", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("survival", comparison = TRUE)$cv
  x$results[[1]]$metrics$time_auc$status <- "incomplete"
  x$results[[1]]$metrics$time_auc$estimate <- NA_real_
  p <- sb_plot_performance(x, "time_auc")
  expect_identical(nrow(p$data), 2L)
  expect_identical(length(p$scales$get_scales("x")$limits), 1L)
  audit <- sb_plot_performance(x, "time_auc", show_unavailable = TRUE)
  expect_identical(length(audit$scales$get_scales("x")$limits), 2L)
  x$results[[2]]$metrics$time_auc$estimate <- Inf
  empty <- ggplot2::ggplot_build(sb_plot_performance(x, "time_auc"))
  expect_match(empty$data[[1]]$label, "No usable results", fixed = TRUE)
})

test_that("report includes saved horizon metrics and leaves source failures intact", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("survival")
  sections <- list(CV = x$cv, External = x$external)
  before <- serialize(sections, NULL)
  failure_before <- sb_failure_table(sections)
  output <- tempfile(fileext = ".html")
  sb_render_report(sections, output)
  html <- paste(readLines(output, warn = FALSE), collapse = "\n")
  for (metric in c("time_auc", "uno_c", "ipcw_brier")) {
    expect_identical(grepl(paste0("<summary>", metric, "</summary>"), html, fixed = TRUE), TRUE)
  }
  expect_match(html, "prediction_horizon", fixed = TRUE)
  expect_match(html, "time_unit", fixed = TRUE)
  expect_identical(serialize(sections, NULL), before)
  expect_identical(sb_failure_table(sections), failure_before)
  expect_identical(file.exists(sub(".html$", "_failures.csv", output)), TRUE)
})

test_that("external missing horizon values and original errors are preserved", {
  fixture <- make_sb_curve_view_fixture("survival")
  external <- sb_example_survival(60, seed = 3)
  external$patient_id <- paste0("NEW_", external$patient_id)
  fit <- sb_fit_final(fixture$task, sb_pipeline())
  prediction <- sb_predict_final(fit, external[c("patient_id", "gene_a", "gene_b")], "patient_id")
  endpoint <- sb_surv_endpoint(external, "patient_id", "time", "event", "month")
  evaluated <- sb_evaluate_external(prediction, endpoint)
  expect_identical(is.null(evaluated$endpoint$prediction_horizon), TRUE)
  for (metric in c("time_auc", "uno_c", "ipcw_brier")) {
    expect_identical(sb_result_table(evaluated, metric = metric)$prediction_horizon, 6)
  }
  prediction$predictions$patient_id[1] <- "missing_external_patient"
  failed <- sb_evaluate_external(prediction, endpoint)
  before <- failed
  expect_identical(failed$status, "failed")
  tab <- sb_result_table(failed, metric = "time_auc")
  expect_identical(tab$status, "not_saved")
  expect_identical(tab$estimate, NA_real_)
  expect_identical(tab$error_message, failed$error_message)
  expect_identical(tab$error_class, failed$error_class)
  expect_identical(failed, before)
})

test_that("horizon metrics remain survival-only and distinguish units", {
  skip_if_not_installed("ggplot2")
  binary <- make_sb_curve_view_fixture("binary")
  for (x in binary[c("cv", "external")]) {
    error <- tryCatch(sb_result_table(x, metric = "time_auc"), error = identity)
    expect_s3_class(error, "stormbreaker_result_metric_not_saved")
  }
  a <- make_sb_curve_view_fixture("survival")$external
  b <- a
  b$metrics$time_auc$time_unit <- "day"
  x <- structure(list(evaluations = list(a, b)), class = "stormbreaker_external_comparison")
  tab <- sb_result_table(x, metric = "time_auc")
  expect_identical(tab$time_unit, c("month", "day"))
  p <- sb_plot_performance(x, "time_auc")
  expect_length(p$scales$get_scales("x")$limits, 2)
  expect_identical(nrow(ggplot2::ggplot_build(p)$data[[1]]), 2L)
})

test_that("horizon-free HTML does not invent time-specific sections", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_survival(60, seed = 1)
  task <- sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month")
  cv <- sb_run_cv(task, folds = 3, seed = 1)
  fit <- sb_fit_final(task, sb_pipeline())
  ext <- d
  ext$patient_id <- paste0("EXT_", ext$patient_id)
  pred <- sb_predict_final(fit, ext[c("patient_id", "gene_a")], "patient_id")
  endpoint <- sb_surv_endpoint(ext, "patient_id", "time", "event", "month")
  evaluated <- sb_evaluate_external(pred, endpoint)
  expect_identical(sb_result_table(evaluated, metric = "time_auc")$status, "not_requested")
  file <- tempfile(fileext = ".html")
  sb_render_report(list(CV = cv, External = evaluated), file)
  html <- paste(readLines(file, warn = FALSE), collapse = "\n")
  for (metric in c("time_auc", "uno_c", "ipcw_brier")) {
    expect_identical(grepl(paste0("<summary>", metric, "</summary>"), html, fixed = TRUE), FALSE)
  }
})
