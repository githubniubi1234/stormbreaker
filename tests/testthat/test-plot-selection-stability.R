make_sb_selection_plot_fixture <- function(failed = FALSE) {
  selections <- list(
    list(fold = 1L, status = "success", selected_features = c("a", "b")),
    list(fold = 2L, status = "success", selected_features = c("b", "c")),
    list(fold = 3L, status = if (failed) "failed" else "success", selected_features = c("a", "b")))
  structure(list(status = if (failed) "failed" else "success", pipeline = list(id = "selected"),
    selections = selections, selection_stability = sb_selection_stability("test", selections, c("a", "b", "c", "d"))),
    class = "stormbreaker_binary_cv_result")
}

test_that("selection frequency copies saved denominators and labels incomplete selections", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_selection_plot_fixture(TRUE)
  before <- x
  p <- sb_plot_selection_frequency(x)
  expect_equal(p$data$selection_frequency, c(0.5, 1, 0.5, 0))
  expect_identical(unique(p$data$n_successful_folds), 2L)
  expect_match(as.character(p$data$panel[1]), "2/3", fixed = TRUE)
  expect_match(as.character(p$data$panel[1]), "incomplete", fixed = TRUE)
  expect_identical(x, before)
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 4L)
})

test_that("Jaccard heatmap uses saved fold pairs without invented diagonal", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_selection_plot_fixture()
  p <- sb_plot_selection_jaccard(x)
  expect_equal(p$data$jaccard, x$selection_stability$pairwise$jaccard)
  expect_identical(p$data$fold_i, x$selection_stability$pairwise$fold_i)
  expect_equal(nrow(p$data), 3L)
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 3L)
})

test_that("UpSet intersections come from exact saved successful fold sets", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_selection_plot_fixture()
  p <- sb_plot_selection_upset(x)
  intersections <- attr(p, "stormbreaker_intersections")
  expect_equal(sum(intersections$count), 3L)
  expect_setequal(intersections$membership, c("1,3", "1,2,3", "2"))
  expect_setequal(unlist(intersections$features, use.names = FALSE), c("a", "b", "c"))
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  p <- sb_plot_selection_upset(make_sb_selection_plot_fixture(TRUE))
  expect_equal(attr(p, "stormbreaker_selection_status")$n_successful_folds, 2L)
  expect_match(p$labels$caption, "successful folds", fixed = TRUE)
})

test_that("non-applicable stability does not invent a selected set", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_selection_plot_fixture()
  x$selection_stability$status <- "not_applicable"
  for (fun in list(sb_plot_selection_frequency, sb_plot_selection_jaccard, sb_plot_selection_upset)) {
    p <- fun(x)
    expect_match(p$labels$subtitle, "No usable", fixed = TRUE)
    expect_identical(attr(p, "stormbreaker_selection_status")$status, "not_applicable")
  }
})
