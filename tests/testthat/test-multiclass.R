multiclass_test_task <- function(grouped = FALSE) {
  d <- datasets::iris
  d$id <- seq_len(nrow(d))
  if (grouped) {
    d$patient <- d$id
    d <- d[rep(seq_len(nrow(d)), each = 2), ]
    d$id <- seq_len(nrow(d))
  }
  sb_multiclass_task(d, "id", "Species", names(d)[1:4],
    classes = c("virginica", "setosa", "versicolor"),
    group_col = if (grouped) "patient" else NULL)
}

test_that("multiclass task preserves labels and rejects ambiguous roles", {
  task <- multiclass_test_task()
  expect_identical(task$classes, c("virginica", "setosa", "versicolor"))
  expect_identical(task$data$Species, datasets::iris$Species)
  bad <- task$data
  bad$Species[1] <- NA
  expect_s3_class(tryCatch(sb_multiclass_task(bad, "id", "Species", task$predictors),
    error = identity), "stormbreaker_bad_multiclass_outcome")
  expect_s3_class(tryCatch(sb_multiclass_task(task$data, "id", "Species", "Species"),
    error = identity), "stormbreaker_bad_predictors")
  bad <- task$data[task$data$Species != "setosa", ]
  expect_s3_class(tryCatch(sb_multiclass_task(bad, "id", "Species", task$predictors),
    error = identity), "stormbreaker_bad_multiclass_outcome")
  expect_identical(sb_method_catalog("multiclass")$method[
    sb_method_catalog("multiclass")$supported], c("boruta", "varselrf", "vsurf",
      "lasso", "elastic_net", "abess", "random_forest",
      "lda", "qda", "gaussian_naive_bayes", "decision_tree", "ridge",
      "svm_linear", "svm_radial", "svm_radial_weighted", "svm_polynomial", "svm_sigmoid",
      "kernel_naive_bayes", "rda", "nnet", "knn", "ranger", "ranger_extratrees",
      "c50", "c50_rules", "ctree", "cforest", "bagging", "rrf", "tree", "evtree",
      "fda", "mda", "pda", "slda", "sda", "loclda", "pamr", "gaussian_process",
      "kernel_weighted_knn", "liblinear_l2_logistic", "liblinear_l1_logistic",
      "xgboost", "xgboost_dart", "xgboost_linear", "lightgbm", "lightgbm_dart",
      "lightgbm_goss", "lightgbm_rf", "catboost",
      "edda", "mclustda", "hdda", "hda", "rlt", "aorsf", "wsrf",
      "logistic", "stepglm", "brglm2", "gam", "pre",
      "mrmr", "jmi", "cmim", "disr", "jmim", "lpc"))
})

test_that("native three-class fits preserve probabilities and selection per class", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("ggplot2")
  task <- multiclass_test_task()
  analysis <- sb_analyze_methods(task, c("lasso", "elastic_net", "random_forest", "mcp"), inner_folds = 3)
  expect_identical(analysis$status$status, c(rep("success", 3), "failed"))
  expect_identical(analysis$results$mcp$error_class, "stormbreaker_unsupported_method_task")
  expect_identical(analysis$task_type, "multiclass")
  rows <- c(121, 2, 61)
  for (method in c("lasso", "elastic_net", "random_forest")) {
    result <- analysis$results[[method]]
    fit <- result$fit
    p <- sb_predict_multiclass(fit, task$data[rows, ])
    expect_identical(colnames(p), task$classes)
    expect_equal(unname(rowSums(p)), rep(1, 3), tolerance = 1e-8)
    one <- sb_predict_multiclass(fit, task$data[rows[1], , drop = FALSE])
    expect_equal(unname(one), unname(p[1, , drop = FALSE]))
    if (method == "random_forest") {
      native <- predict(fit$engine, task$data[rows, task$predictors], type = "prob")
      expect_equal(unname(p), unname(as.matrix(native[, fit$codes])))
      expect_identical(sb_method_features(analysis, method)$selected, rep(NA, 4))
      expect_s3_class(tryCatch(sb_selected_features(analysis, method), error = identity),
        "stormbreaker_selection_not_defined")
    } else {
      native <- predict(fit$engine, as.matrix(task$data[rows, task$predictors]),
        s = fit$lambda, type = "response")
      expect_equal(unname(p), unname(native[, fit$codes, 1]))
      beta <- coef(fit$engine, s = fit$lambda)
      union <- character()
      for (k in seq_along(task$classes)) {
        b <- as.matrix(beta[[fit$codes[k]]])
        selected <- rownames(b)[b[, 1] != 0 & rownames(b) != "(Intercept)"]
        expect_setequal(sb_selected_features(analysis, method, class = task$classes[k]), selected)
        union <- union(union, selected)
      }
      expect_setequal(sb_selected_features(analysis, method), union)
    }
    expect_identical(result$plot_status$status, rep("available", length(result$plots)))
    class_plot <- if (method == "random_forest") result$plots$class_importance else
      result$plots$coefficient_path
    layout <- ggplot2::ggplot_build(class_plot)$layout$layout
    expect_identical(as.character(layout$class), task$classes)
    expect_identical(any(grepl("confusion", names(result$plots))), FALSE)
    for (plot in result$plots) {
      expect_null(plot$labels$subtitle)
      expect_null(plot$labels$caption)
      expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
    }
  }
  local_mocked_bindings(sb_fit_multiclass = function(...) stop("Unexpected refit"))
  saved <- unserialize(serialize(analysis, NULL))
  expect_identical(sb_selected_features(saved, "lasso", "setosa"),
    sb_selected_features(analysis, "lasso", "setosa"))
  expect_equal(sb_predict_multiclass(saved$results$lasso$fit, task$data[rows, ]),
    sb_predict_multiclass(analysis$results$lasso$fit, task$data[rows, ]))
})

test_that("metrics match hand calculations and confusion counts do not overlap", {
  truth <- c("a", "b", "c", "a", "b", "c")
  p <- matrix(c(.8,.1,.1, .1,.8,.1, .1,.1,.8,
    .1,.8,.1, .1,.8,.1, .1,.8,.1), byrow = TRUE, ncol = 3,
    dimnames = list(NULL, c("a", "b", "c")))
  m <- sb_multiclass_metrics(truth, p)
  expect_equal(m$overall$accuracy, 4 / 6)
  expect_equal(m$overall$balanced_accuracy, 2 / 3)
  expect_equal(m$overall$macro_f1, 2 / 3)
  expect_equal(m$overall$log_loss, -(4 * log(.8) + 2 * log(.1)) / 6)
  expect_equal(m$overall$brier, (4 * .06 + 2 * 1.46) / 6)
  absent <- sb_multiclass_metrics(truth[1:2], p[1:2, ])
  expect_identical(is.na(absent$overall$macro_auc), TRUE)
  expect_identical(is.na(absent$overall$balanced_accuracy), TRUE)
  invalid <- p
  invalid[1, 1] <- .9
  expect_s3_class(tryCatch(sb_multiclass_metrics(truth, invalid), error = identity),
    "stormbreaker_bad_multiclass_probabilities")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    text <- ggplot2::ggplot_build(sb_plot_multiclass_confusion(m))$data[[2]]
    expect_identical(nrow(text), 9L)
    expect_identical(anyDuplicated(text[c("x", "y")]), 0L)
    expect_equal(sum(as.numeric(text$label)), length(truth))
  }
})

test_that("four classes with engine-like labels keep their original mapping", {
  skip_if_not_installed("randomForest")
  skip_if_not_installed("glmnet")
  d <- withr::with_seed(15, {
    x <- data.frame(id = seq_len(120), x1 = rnorm(120), x2 = rnorm(120), x3 = rnorm(120))
    x$label <- rep(c("Overall", "MeanDecreaseGini", "A type", "class2"), each = 30)
    x
  })
  classes <- c("class2", "Overall", "A type", "MeanDecreaseGini")
  task <- sb_multiclass_task(d, "id", "label", c("x1", "x2", "x3"), classes)
  result <- sb_analyze_methods(task, c("random_forest", "lasso"), inner_folds = 3)
  expect_identical(result$status$status, c("success", "success"))
  for (fit in lapply(result$results, `[[`, "fit")) {
    p <- sb_predict_multiclass(fit, d[c(91, 1, 61, 31), ])
    expect_identical(colnames(p), classes)
    expect_equal(unname(rowSums(p)), rep(1, 4))
  }
})

test_that("nested folds keep repeated patients together and exclude held-out rows", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("randomForest")
  task <- multiclass_test_task(grouped = TRUE)
  v <- sb_compare_multiclass(task, c("lasso", "random_forest"), folds = 3, inner_folds = 3)
  expect_identical(v$comparison$status, c("success", "success"))
  groups <- task$data$patient
  expect_identical(all(vapply(split(v$assignment, groups), function(x) length(unique(x)) == 1, logical(1))), TRUE)
  for (r in v$results) {
    expect_identical(rownames(r$probability), as.character(task$data$id))
    expect_equal(unname(rowSums(r$probability)), rep(1, task$n))
    for (fold in r$folds) {
      expect_length(intersect(groups[fold$training_rows], groups[fold$assessment_rows]), 0)
      if (!is.null(fold$inner_fold_id)) {
        expect_length(fold$inner_fold_id, length(fold$training_rows))
        sizes <- vapply(split(fold$inner_fold_id, groups[fold$training_rows]),
          function(x) length(unique(x)), integer(1))
        expect_identical(all(sizes == 1L), TRUE)
      }
    }
  }
  failed <- sb_compare_multiclass(multiclass_test_task(), c("lasso", "random_forest"),
    folds = 3, inner_folds = 100)
  expect_identical(failed$comparison$status, c("failed", "success"))
  expect_null(failed$results$lasso$metrics)
  expect_identical(is.na(failed$comparison$accuracy[1]), TRUE)
})
