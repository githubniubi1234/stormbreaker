test_that("repeated displays and exports preserve separate identities and failures", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture()$cv
  failed <- x
  failed$status <- "failed"
  failed$runs$status[1] <- "failed"
  failed$runs$error_message[1] <- "original, \"reason\"\nnext line"
  repeated <- structure(list(repetitions = list(x, failed)),
                        class = "stormbreaker_repeated_cv_comparison")
  before <- serialize(repeated, NULL)
  p <- sb_plot_performance(repeated, "roc_auc", show_unavailable = TRUE)
  expect_length(levels(p$data$display_label), 2)
  expect_equal(nrow(p$layers[[1]]$data), 1)
  for (fn in list(sb_plot_calibration, sb_plot_dca)) {
    p <- fn(repeated, show_unavailable = TRUE)
    status <- attr(p, "stormbreaker_curve_status")
    expect_equal(status$repeat_id, 1:2)
    expect_length(unique(status$panel), 2)
    expect_silent(ggplot2::ggplot_build(p))
  }
  failures <- sb_failure_table(repeated)
  expect_true(all(failures$repeat_id == 2L))
  expect_true(any(failures$error_message == failed$runs$error_message[1], na.rm = TRUE))
  mixed <- sb_failure_table(list(single = x, repeated = repeated))
  expect_equal(nrow(mixed), nrow(failures))
  out <- withr::local_tempdir()
  html_path <- file.path(out, "report.html")
  sb_render_report(list(repeated = repeated, ordinary = x), html_path)
  html <- paste(readLines(html_path), collapse = "\n")
  expect_true(grepl("Development repeated CV", html, fixed = TRUE))
  csv <- read.csv(file.path(out, "report_failures.csv"))
  expect_equal(csv$repeat_id, failures$repeat_id)
  expect_equal(csv$error_message, failures$error_message)
  expect_error(sb_render_report(repeated, html_path), class = "stormbreaker_report_exists")
  expect_identical(serialize(repeated, NULL), before)
})

test_that("IECV reports and curves label development validation without pooling", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("survival")$external
  iecv <- structure(list(folds = list(
    list(heldout_cohort_id = "Site <A>", comparison = x),
    list(heldout_cohort_id = "Site B", comparison = x))),
    class = "stormbreaker_internal_external_cv")
  for (fn in list(sb_plot_calibration, sb_plot_dca)) {
    p <- fn(iecv, show_unavailable = TRUE)
    status <- attr(p, "stormbreaker_curve_status")
    expect_equal(status$heldout_cohort_id, c("Site <A>", "Site B"))
    expect_true(all(status$validation == "internal_external_cv"))
    expect_true(all(grepl("not final independent", status$panel, fixed = TRUE)))
    expect_silent(ggplot2::ggplot_build(p))
  }
  path <- file.path(withr::local_tempdir(), "iecv.html")
  sb_render_report(iecv, path)
  html <- paste(readLines(path), collapse = "\n")
  expect_true(grepl("Development IECV (not final independent external validation)", html, fixed = TRUE))
  expect_true(grepl("Site &lt;A&gt;", html, fixed = TRUE))
  expect_true(grepl("heldout_cohort_id", html, fixed = TRUE))
  expect_equal(nrow(sb_failure_table(iecv)), 0)
})

test_that("composite strip labels wrap long identifiers without losing characters", {
  labels <- c(paste(rep("long_model_identifier", 8), collapse = "__"),
              "Held-out cohort A\nDevelopment IECV")
  wrapped <- sb_view_wrap_label(labels)
  expect_true(all(nchar(unlist(strsplit(wrapped, "\n", fixed = TRUE))) <= 50L))
  expect_equal(gsub("\n", "", wrapped, fixed = TRUE), gsub("\n", "", labels, fixed = TRUE))
})
