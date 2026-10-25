test_that("multiclass neighbor weights recover native predictions and paginate all queries", {
  skip_if_not_installed("kknn")
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  for (method in c("knn", "kernel_weighted_knn")) {
    f <- sb_fit_multiclass(task, method, 3, 1)
    expect_error(sb_plots(f), "Supply new_data")
    q <- d[1:9, ]
    out <- sb_plots(f, new_data = q)
    tab <- out$diagnostics$neighborhoods
    predicted <- sb_predict_multiclass(f, q)
    votes <- matrix(0, nrow(q), length(f$classes))
    for (i in seq_len(nrow(tab))) {
      cl <- match(tab$class[i], f$classes)
      votes[tab$query[i], cl] <- votes[tab$query[i], cl] + tab$weight[i]
    }
    expect_equal(unname(predicted), votes, tolerance = 1e-10)
    expect_equal(rowSums(votes), rep(1, nrow(q)))
    expect_equal(sum(vapply(out$pages, function(p) nrow(p$data), 1L)), nrow(tab))
    expect_equal(length(out$pages), 2L)
    for (p in out$pages) expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  }
})
