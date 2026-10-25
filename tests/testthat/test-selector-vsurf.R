test_that("VSURF saved evidence plots integrate into a separate clinician group", {
  skip_if_not_installed("VSURF")
  fixture <- withr::with_seed(9211L, {
    y <- sample(rep(c(0L, 1L), length.out = 72L))
    x <- matrix(rnorm(72L * 30L), nrow = 72L)
    x[, 1] <- x[, 1] + 1.4 * y
    x[, 2] <- x[, 2] - y
    x[, 3] <- x[, 3] + .7 * y
    x <- as.data.frame(x)
    names(x) <- paste0("gene_", seq_len(30L))
    list(x = x, y = y)
  })
  state <- sb_fit_vsurf_binary_selector(fixture$x, fixture$y)
  before <- serialize(state, NULL, version = 2L)
  importance <- sb_selector_diagnostics(state, "vsurf_importance")
  stages <- sb_selector_diagnostics(state, "vsurf_stages")
  expect_identical(importance$table, state$scores)
  expect_identical(stages$table, state$selection_history)
  expect_identical(importance$selected_features, state$selected_features)
  expect_equal(nrow(importance$table), 30L)
  ip <- sb_plot_selector_diagnostic(importance)
  sp <- sb_plot_selector_diagnostic(stages)
  expect_equal(ggplot2::ggplot_build(ip)$data[[2]]$x,
    state$scores$importance - state$scores$importance_sd)
  expect_equal(ggplot2::ggplot_build(ip)$data[[2]]$xend,
    state$scores$importance + state$scores$importance_sd)
  expect_identical(levels(ip$data$feature), rev(state$scores$feature[order(state$scores$rank)]))
  expect_equal(ggplot2::ggplot_build(sp)$data[[1]]$y, state$selection_history$n_features)
  result <- make_sb_curve_view_fixture("binary")$cv
  views <- sb_clinician_views(result, selectors = list(Importance = importance, Stages = stages))
  expect_identical(views$groups[grepl("^Selector:", names(views$plots))], rep("Feature selection", 2))
  expect_identical(views$status$status[grepl("^Selector:", views$status$view)], rep("available", 2))
  expect_null(views$plots[["Selector: Importance"]]$labels$caption)
  expect_match(views$captions[match("Selector: Importance", names(views$plots))], "between-forest SD", fixed = TRUE)
  file <- tempfile(fileext = ".html")
  sb_render_clinician_report(views, file)
  html <- paste(readLines(file), collapse = "\n")
  expect_match(html, "Feature selection", fixed = TRUE)
  expect_match(html, "between-forest SD", fixed = TRUE)
  expect_identical(serialize(state, NULL, version = 2L), before)
})
