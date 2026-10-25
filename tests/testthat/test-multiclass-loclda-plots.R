test_that("local LDA retains radius ties and excluded class neighbors", {
  skip_if_not_installed("klaR")
  d <- iris
  d$id <- seq_len(nrow(d))
  f <- sb_fit_multiclass(sb_multiclass_task(d, "id", "Species", names(iris)[1:4]), "loclda", 3, 1)
  expect_error(sb_plots(f), "Supply new_data")
  out <- sb_plots(f, new_data = d[1:7, ])
  expect_equal(length(out$pages), 2L)
  tab <- out$diagnostics$neighborhoods
  for (q in 1:7) {
    distances <- sqrt(rowSums(sweep(f$engine$learn, 2, as.numeric(d[q, 1:4]), "-")^2))
    radius <- sort(distances)[f$engine$k]
    selected <- which(distances <= radius)
    actual <- tab[tab$query == q, ]
    expect_equal(actual$reference_row, selected)
    expect_equal(actual$weight, vapply(distances[selected] / radius, f$engine$weight.func, 1), ignore_attr = TRUE)
  }
  e <- list(learn = matrix(c(0, 1, -1, 2, 5), ncol = 1), k = 3,
    grouping = factor(c("a", "a", "b", "c", "c")), lev = c("a", "b", "c"), weight.func = function(x) 1 - x)
  ties <- sb_loclda_neighbor_table(e, matrix(0, 1, 1))
  expect_equal(ties$reference_row, 1:3)
  expect_equal(ties$used, c(TRUE, TRUE, FALSE))
  expect_equal(ties$effective_weight, c(1, 0, 0))
  for (p in out$pages) expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})
