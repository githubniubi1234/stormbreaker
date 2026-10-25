multiclass_selection_task <- function(k = 3L) {
  withr::with_seed(1717, {
    label <- rep(seq_len(k), each = 30L)
    d <- data.frame(id = seq_along(label), group = paste0("Type ", label))
    for (j in seq_len(k)) d[[paste0("gene ", j)]] <- 3 * (label == j) + rnorm(length(label))
    d$noise <- rnorm(length(label))
    sb_multiclass_task(d, "id", "group", names(d)[-(1:2)],
      classes = rev(paste0("Type ", seq_len(k))))
  })
}

test_that("multiclass forest selectors return native selections and saved plots", {
  for (p in c("Boruta", "varSelRF", "VSURF")) skip_if_not_installed(p)
  task <- multiclass_selection_task()
  methods <- c("boruta", "varselrf", "vsurf")
  analysis <- sb_analyze_methods(task, methods)
  expect_identical(analysis$status$status, rep("success", 3))
  for (method in methods) {
    result <- analysis$results[[method]]
    expect_setequal(sb_selected_features(analysis, method), result$fit$selected_features)
    expect_identical(result$fit$classes, task$classes)
    expect_identical(result$features$selected,
      task$predictors %in% result$fit$selected_features)
    expect_identical(result$plot_status$status, rep("available", nrow(result$plot_status)))
    for (p in result$plots) {
      expect_null(p$labels$caption)
      expect_null(p$labels$subtitle)
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    }
  }
  restored <- unserialize(serialize(analysis, NULL))
  expect_identical(sb_selected_features(restored, "boruta"), sb_selected_features(analysis, "boruta"))
  expect_s3_class(tryCatch(sb_compare_multiclass(task, "boruta", folds = 3), error = identity),
    "stormbreaker_bad_method")
})

test_that("ABESS uses native multinomial class order and GIC selection", {
  skip_if_not_installed("abess")
  for (k in c(3L, 4L)) {
    task <- multiclass_selection_task(k)
    a <- sb_analyze_methods(task, "abess")
    expect_identical(a$status$status, "success")
    fit <- a$results$abess$fit
    rows <- c(task$n, 1, 31)
    p <- sb_predict_multiclass(fit, task$data[rows, ])
    native <- predict(fit$engine, as.matrix(task$data[rows, task$predictors]),
      support.size = fit$engine$best.size, type = "response")[[1]]
    expect_equal(unname(p), unname(native))
    expect_identical(colnames(p), task$classes)
    expect_equal(unname(rowSums(p)), rep(1, length(rows)))
    expect_equal(unname(sb_predict_multiclass(fit, task$data[rows[1], , drop = FALSE])),
      unname(p[1, , drop = FALSE]))
    allp <- sb_predict_multiclass(fit, task$data)
    expect_gt(mean(task$classes[max.col(allp)] == task$data$group), .85)
    beta <- as.matrix(coef(fit$engine, support.size = fit$engine$best.size)[[1]])[-1, , drop = FALSE]
    expect_setequal(sb_selected_features(a, "abess"), task$predictors[rowSums(beta != 0) > 0])
    reference <- sb_method_features(a, "abess", class = tail(task$classes, 1))
    expect_equal(reference$coefficient, rep(0, length(task$predictors)))
    expect_identical(a$results$abess$plot_status$status, c("available", "available"))
    for (plot in a$results$abess$plots) expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
    restored <- unserialize(serialize(fit, NULL))
    expect_equal(sb_predict_multiclass(restored, task$data), allp)
  }
  cv <- sb_compare_multiclass(task, "abess", folds = 3)
  expect_identical(cv$comparison$status, "success")
  expect_equal(unname(rowSums(cv$results$abess$probability)), rep(1, task$n))
})
