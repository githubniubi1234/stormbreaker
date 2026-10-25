test_that("specialist engines return aligned portable multiclass predictions", {
  for (p in c("mclust", "HDclassif", "hda", "RLT", "aorsf", "wsrf")) skip_if_not_installed(p)
  methods <- c("edda", "mclustda", "hdda", "hda", "rlt", "aorsf", "wsrf")
  for (k in c(3L, 4L)) {
    withr::local_seed(934)
    label <- rep(seq_len(k), each = 45)
    d <- data.frame(id = seq_along(label), outcome = paste0("Type ", label))
    for (j in seq_len(6)) d[[paste0("gene ", j)]] <- 2 * (label == j) + rnorm(length(label))
    task <- sb_multiclass_task(d, "id", "outcome", names(d)[-(1:2)],
      classes = rev(paste0("Type ", seq_len(k))))
    a <- sb_analyze_methods(task, methods)
    expect_identical(a$status$status, rep("success", length(methods)))
    for (method in methods) {
      fit <- a$results[[method]]$fit
      p <- sb_predict_multiclass(fit, d)
      expect_identical(colnames(p), task$classes)
      expect_equal(unname(rowSums(p)), rep(1, nrow(d)), tolerance = 1e-6)
      if (method == "hda") {
        x <- as.matrix(d[task$predictors])
        colnames(x) <- fit$feature_codes
        native <- predict(fit$engine, newdata = x, task = "c")
        expect_equal(unname(p), unname(native$posteriors[, fit$codes, drop = FALSE]))
        expect_identical(fit$codes[max.col(p)], as.character(native$prediction))
      } else expect_gt(mean(task$classes[max.col(p)] == d$outcome), .65)
      rows <- c(nrow(d), 1, 46)
      expect_equal(unname(sb_predict_multiclass(fit, d[rows, ])), unname(p[rows, , drop = FALSE]))
      expect_equal(unname(sb_predict_multiclass(fit, d[1, , drop = FALSE])), unname(p[1, , drop = FALSE]))
      expect_equal(unname(sb_predict_multiclass(unserialize(serialize(fit, NULL)), d[rows, ])),
        unname(p[rows, , drop = FALSE]))
    }
  }
  cv <- sb_compare_multiclass(task, methods, folds = 3)
  expect_identical(cv$comparison$status, rep("success", length(methods)))
})
