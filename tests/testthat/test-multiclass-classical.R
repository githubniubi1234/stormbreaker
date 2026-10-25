classical_multiclass_task <- function(k = 3L) {
  withr::with_seed(816, {
    n <- 40L * k
    code <- rep(seq_len(k), each = 40L)
    d <- data.frame(id = seq_len(n), a = rnorm(n) + code,
      b = rnorm(n) - code / 2, c = rnorm(n), label = paste0("Type ", code))
    sb_multiclass_task(d, "id", "label", c("a", "b", "c"),
      classes = rev(paste0("Type ", seq_len(k))))
  })
}

test_that("four classical adapters retain native multiclass probabilities and evidence", {
  for (package in c("MASS", "e1071", "rpart", "ggplot2")) skip_if_not_installed(package)
  methods <- c("lda", "qda", "gaussian_naive_bayes", "decision_tree")
  for (k in c(3L, 4L)) {
    task <- classical_multiclass_task(k)
    analysis <- sb_analyze_methods(task, methods)
    expect_identical(analysis$status$status, rep("success", 4))
    for (method in methods) {
      result <- analysis$results[[method]]
      fit <- result$fit
      rows <- c(task$n, 1, 41)
      newdata <- task$data[rows, task$predictors]
      p <- sb_predict_multiclass(fit, newdata)
      native <- if (method %in% c("lda", "qda")) predict(fit$engine, newdata)$posterior else
        if (method == "gaussian_naive_bayes") predict(fit$engine, newdata, type = "raw") else {
          names(newdata) <- fit$feature_codes
          predict(fit$engine, newdata, type = "prob")
        }
      expect_equal(unname(p), unname(native[, fit$codes, drop = FALSE]))
      expect_identical(colnames(p), task$classes)
      expect_equal(unname(rowSums(p)), rep(1, length(rows)))
      expect_equal(unname(sb_predict_multiclass(fit, task$data[rows[1], , drop = FALSE])),
        unname(p[1, , drop = FALSE]))
      expect_s3_class(tryCatch(sb_selected_features(analysis, method), error = identity),
        "stormbreaker_selection_not_defined")
      if (method != "decision_tree") {
        for (label in task$classes) {
          tab <- sb_method_features(analysis, method, class = label)
          original <- task$data[task$data$label == label, task$predictors]
          expect_equal(unname(tab$mean), unname(colMeans(original)))
          if (method == "gaussian_naive_bayes") expect_equal(unname(tab$sd),
            unname(vapply(original, stats::sd, numeric(1))))
        }
      } else {
        tab <- sb_method_features(analysis, method)
        expect_identical(tab$used_for_split, fit$feature_codes %in% fit$engine$frame$var)
        expect_identical(fit$parameters$xval, 0L)
      }
      for (plot in result$plots) {
        expect_null(plot$labels$subtitle)
        expect_null(plot$labels$caption)
        expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
      }
      if (method == "decision_tree") {
        expect_match(result$notes$tree_structure, "Every original node", all = FALSE)
      }
    }
  }
  local_mocked_bindings(sb_fit_multiclass_classical = function(...) stop("Unexpected refit"))
  restored <- unserialize(serialize(analysis, NULL))
  expect_identical(sb_method_features(restored, "lda", "Type 2"),
    sb_method_features(analysis, "lda", "Type 2"))
  expect_equal(sb_predict_multiclass(restored$results$qda$fit, task$data),
    sb_predict_multiclass(analysis$results$qda$fit, task$data))
  expect_equal(sb_predict_multiclass(restored$results$decision_tree$fit, task$data),
    sb_predict_multiclass(analysis$results$decision_tree$fit, task$data))
  expect_identical(attr(restored$results$decision_tree$fit$engine$terms, ".Environment"), baseenv())
})

test_that("shared CV supports new methods and retains singular QDA failures", {
  for (package in c("MASS", "e1071", "rpart")) skip_if_not_installed(package)
  task <- classical_multiclass_task()
  methods <- c("lda", "qda", "gaussian_naive_bayes", "decision_tree")
  result <- sb_compare_multiclass(task, methods, folds = 3)
  expect_identical(result$comparison$status, rep("success", 4))
  for (method in methods) expect_equal(unname(rowSums(result$results[[method]]$probability)),
    rep(1, task$n))
  task$data$c <- task$data$a
  singular <- sb_compare_multiclass(task, c("qda", "gaussian_naive_bayes"), folds = 3)
  expect_identical(singular$comparison$status, c("failed", "success"))
  expect_null(singular$results$qda$metrics)
  expect_identical(is.na(singular$comparison$accuracy[1]), TRUE)
})
