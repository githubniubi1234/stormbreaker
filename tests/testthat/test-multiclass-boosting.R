test_that("boosting objectives produce portable native multiclass probabilities", {
  for (p in c("xgboost", "lightgbm", "catboost")) skip_if_not_installed(p)
  methods <- c("xgboost", "xgboost_dart", "xgboost_linear", "lightgbm",
    "lightgbm_dart", "lightgbm_goss", "lightgbm_rf", "catboost")
  for (k in c(3L, 4L)) {
    withr::local_seed(971)
    code <- rep(seq_len(k), each = 40)
    d <- data.frame(id = seq_along(code), type = paste0("Type ", code))
    for (j in seq_len(5)) d[[paste0("gene ", j)]] <- 3 * (code == j) + rnorm(length(code))
    task <- sb_multiclass_task(d, "id", "type", names(d)[-(1:2)],
      classes = rev(paste0("Type ", seq_len(k))))
    a <- sb_analyze_methods(task, methods)
    expect_identical(a$status$status, rep("success", length(methods)))
    for (method in methods) {
      fit <- a$results[[method]]$fit
      p <- sb_predict_multiclass(fit, d)
      expect_identical(colnames(p), task$classes)
      expect_equal(unname(rowSums(p)), rep(1, nrow(d)), tolerance = 1e-6)
      expect_gt(mean(task$classes[max.col(p)] == d$type), .8)
      rows <- c(nrow(d), 1, 41)
      expect_equal(unname(sb_predict_multiclass(fit, d[rows, ])), unname(p[rows, , drop = FALSE]))
      expect_equal(unname(sb_predict_multiclass(fit, d[1, , drop = FALSE])), unname(p[1, , drop = FALSE]))
      expect_equal(unname(sb_predict_multiclass(unserialize(serialize(fit, NULL)), d[rows, ])),
        unname(p[rows, , drop = FALSE]))
      x <- as.matrix(d[rows, task$predictors])
      colnames(x) <- fit$feature_codes
      native <- if (fit$package == "xgboost") {
        predict(xgboost::xgb.load.raw(fit$engine), x, reshape = TRUE)
      } else if (fit$package == "lightgbm") {
        predict(lightgbm::lgb.load(model_str = fit$engine), x)
      } else catboost::catboost.predict(unserialize(fit$engine), catboost::catboost.load_pool(x),
        prediction_type = "Probability", thread_count = 1L)
      expect_equal(unname(p[rows, , drop = FALSE]), unname(native), ignore_attr = TRUE)
      for (plot in a$results[[method]]$plots) expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
    }
  }
  comparison <- sb_compare_multiclass(task, methods, folds = 3)
  expect_identical(comparison$comparison$status, rep("success", length(methods)))
})
