test_that("SVM-RFE short entry fits, predicts and draws its native curves", {
  skip_if_not_installed("e1071")
  set.seed(20)
  d <- data.frame(id = paste0("p", 1:120), y = rep(0:1, 60),
    a = rnorm(120), b = rnorm(120), c = rnorm(120))
  d$a <- d$a + 2 * d$y
  task <- sb_binary_task(d, "id", "y", c("a", "b", "c"))
  pipeline <- bi_pipeline("Logistic", bi_selector("SVM_RFE"))
  fit <- sb_fit_final(task, pipeline, capture_diagnostics = TRUE)
  expect_identical(fit$status, "success")
  plots <- sb_plots(fit)
  expect_equal(length(plots$unavailable), 0)
  expect_equal(all(c("selector_1_cv_accuracy", "selector_1_cv_error") %in% names(plots$plots)), TRUE)
  expect_s3_class(ggplot2::ggplot_build(plots$plots$selector_1_cv_accuracy), "ggplot_built")
  expect_identical(sb_run_cv(task, pipeline, folds = 3, seed = 1)$status, "success")
  chain <- bi_pipeline("Logistic", bi_selector(c("GBM", "SVM_RFE")))
  chained <- sb_fit_final(task, chain, capture_diagnostics = TRUE)
  expect_identical(chained$status, "success")
  stages <- chained$selection$state$stages
  evidence <- stages[[2]]$evidence
  rows <- evidence$fold_id != 1L
  upstream <- sb_fit_binary_model_selector(d[rows, c("a", "b", "c")], d$y[rows],
    "gbm_binary_importance", sb_binary_model_selector_options("gbm_binary_importance"))
  expect_setequal(evidence$fold_subsets[[1]][[length(evidence$fold_subsets[[1]])]],
    upstream$selected_features)
  expect_identical(sb_run_cv(task, chain, folds = 3, seed = 1)$status, "success")
})

test_that("SVM-RFE refits nested subsets and respects patient folds", {
  skip_if_not_installed("e1071")
  set.seed(19)
  x <- as.data.frame(matrix(rnorm(120 * 4), 120, 4))
  y <- rep(rep(c(0L, 1L), each = 30), each = 2)
  x$V1 <- x$V1 + 3 * y
  groups <- rep(seq_len(60), each = 2)
  fit <- sb_fit_svm_rfe(x, y, groups)
  e <- fit$evidence
  expect_equal(vapply(split(e$fold_id, groups), function(z) length(unique(z)), integer(1)),
    setNames(rep(1L, 60), as.character(seq_len(60))))
  expect_equal(e$curve$accuracy, colSums(e$correct) / nrow(x))
  expect_equal(e$curve$error, 1 - e$curve$accuracy)
  expect_length(fit$selected_features, which.max(e$curve$accuracy))
  for (k in 2:4) expect_equal(setdiff(e$subsets[[k - 1]], e$subsets[[k]]), character())
  train <- e$fold_id != 1
  replay <- sb_svm_rfe_path(x[train, ], y[train], x[!train, ])
  expect_identical(replay$subsets, e$fold_subsets[[1]])
})
