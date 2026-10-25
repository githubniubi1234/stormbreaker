test_that("combined panels retain ordered evidence and separate explanations", {
  skip_if_not_installed("patchwork")
  result <- make_sb_curve_view_fixture("binary")$cv
  x <- sb_clinician_views(result)
  before <- serialize(x, NULL)
  selected <- c("Calibration", "ROC")
  combined <- sb_plot_clinician_panel(x, selected)
  expect_s3_class(combined, "patchwork")
  for (i in seq_along(selected)) {
    original <- x$plots[[selected[i]]]
    expect_identical(combined[[i]]$data, original$data)
    expect_identical(combined[[i]]$layers, original$layers)
    expect_identical(combined[[i]]$labels$title, selected[i])
    expect_null(combined[[i]]$labels$subtitle)
    expect_null(combined[[i]]$labels$caption)
  }
  expect_identical(attr(combined, "stormbreaker_captions")$view, selected)
  expect_identical(serialize(x, NULL), before)
  expect_s3_class(patchwork::patchworkGrob(combined), "gtable")
  expect_error(sb_plot_clinician_panel(x, "Missing"), class = "stormbreaker_bad_report_input")
  expect_error(sb_plot_clinician_panel(x, c("ROC", "ROC")), class = "stormbreaker_bad_report_input")
})
