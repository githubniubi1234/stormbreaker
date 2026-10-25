test_that("native tree variants preserve multiclass probabilities", {
  for (p in c("ranger", "C50", "partykit", "ipred", "RRF", "tree", "evtree")) skip_if_not_installed(p)
  methods <- c("ranger", "ranger_extratrees", "c50", "c50_rules", "ctree",
    "cforest", "bagging", "rrf", "tree", "evtree")
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(d)[1:4],
    classes = c("virginica", "setosa", "versicolor"))
  a <- sb_analyze_methods(task, methods)
  expect_identical(a$status$status, rep("success", length(methods)))
  for (method in methods) {
    fit <- a$results[[method]]$fit
    rows <- c(120, 2, 61)
    p <- sb_predict_multiclass(fit, d[rows, ])
    expect_identical(colnames(p), task$classes)
    expect_equal(unname(rowSums(p)), rep(1, length(rows)), tolerance = 1e-6)
    expect_equal(unname(sb_predict_multiclass(fit, d[rows[1], , drop = FALSE])),
      unname(p[1, , drop = FALSE]))
    expect_equal(sb_predict_multiclass(unserialize(serialize(fit, NULL)), d[rows, ]), p)
    allp <- sb_predict_multiclass(fit, d)
    expect_gt(mean(task$classes[max.col(allp)] == d$Species), .8)
    expect_s3_class(tryCatch(sb_selected_features(a, method), error = identity),
      "stormbreaker_selection_not_defined")
  }
  cv <- sb_compare_multiclass(task, methods, folds = 3)
  expect_identical(cv$comparison$status, rep("success", length(methods)))
})
