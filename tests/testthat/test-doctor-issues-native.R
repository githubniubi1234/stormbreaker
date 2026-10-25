test_that("native unavailable tables produce one issue per view without stopping runs", {
  local_mocked_bindings(sb_failure_table = function(x) sb_failure_empty())
  native <- data.frame(view = c("path", "importance"),
    status = c("unavailable", "unavailable"),
    reason = c("Original path was not captured", "No native importance"))
  x <- list(comparison = structure(list(), class = "stormbreaker_cv_comparison"),
    pipelines = list(), details = list(model = list(fit = list(status = "success"),
      figures = list(unavailable = list(native = native, extra = "No specialized plot"),
        native = list(status = native)))))
  tab <- sb_run_issues(x)
  expect_equal(nrow(tab), 3L)
  expect_setequal(tab$stage, c("plot path", "plot importance", "plot extra"))
  expect_setequal(tab$reason, c(native$reason, "No specialized plot"))
})
