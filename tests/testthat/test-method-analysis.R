method_analysis_task <- function() {
  withr::with_seed(881L, {
    n <- 90L
    d <- data.frame(id = seq_len(n), a = rnorm(n), b = rnorm(n), c = rnorm(n))
    d$y <- rbinom(n, 1, plogis(2 * d$a - d$b))
    sb_binary_task(d, "id", "y", c("a", "b", "c"))
  })
}

test_that("catalog distinguishes task support and selection from ranking", {
  binary <- sb_method_catalog("binary")
  survival <- sb_method_catalog("survival")
  expect_identical(nrow(binary), 9L)
  expect_identical(binary$supported, rep(TRUE, 9))
  expect_identical(survival$method[!survival$supported], c("varselrf", "vsurf"))
  expect_identical(binary$method[binary$feature_result == "ranking"], "random_forest")
  expect_identical(sb_method_names(c("LASSO", "Elastic Net", "random-forest")),
    c("lasso", "elastic_net", "random_forest"))
})

test_that("actual LASSO and forest runs retain native results and reusable plots", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("ggplot2")
  task <- method_analysis_task()
  x <- sb_analyze_methods(task, c("LASSO", "random_forest"), inner_folds = 3L)
  expect_identical(x$status$status, rep("success", 2))
  lasso <- x$results$lasso
  expect_identical(lasso$fit$pipeline$selector, "none")
  expect_identical(lasso$fit$development$candidate_features, task$predictors)
  beta <- as.matrix(stats::coef(lasso$fit$diagnostic_artifacts$learner,
    s = lasso$fit$tuning$lambda))
  expected <- rownames(beta)[beta[, 1] != 0 & rownames(beta) != "(Intercept)"]
  expect_setequal(sb_selected_features(x, "lasso"), expected)
  rf <- sb_method_features(x, "random_forest")
  native <- unserialize(x$results$random_forest$fit$learner_state$engine_raw)
  expect_equal(rf$importance, unname(native$importance[rf$feature, "MeanDecreaseGini"]))
  expect_identical(rf$selected, rep(NA, nrow(rf)))
  expect_equal(rf$rank, rank(-rf$importance, ties.method = "min"))
  for (result in x$results) {
    expect_identical(result$plot_status$status, rep("available", length(result$plots)))
    for (plot in result$plots) {
      expect_null(plot$labels$subtitle)
      expect_null(plot$labels$caption)
      expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
    }
  }
  local_mocked_bindings(sb_fit_final = function(...) stop("unexpected refit"),
    sb_analyze_one_method = function(...) stop("unexpected analysis"))
  restored <- unserialize(serialize(x, NULL))
  expect_identical(sb_selected_features(restored, "LASSO"), sb_selected_features(x, "lasso"))
  expect_identical(sb_method_features(restored, "random_forest"), rf)
  expect_s3_class(tryCatch(sb_selected_features(restored, "random_forest"), error = identity),
    "stormbreaker_selection_not_defined")
})

test_that("method failures are retained and later methods still run", {
  called <- character()
  local_mocked_bindings(sb_analyze_one_method = function(task, spec, ...) {
    called <<- c(called, spec$method)
    if (spec$method == "boruta") stop("native failure")
    list(status = "success", method = spec$method, feature_result = spec$feature_result)
  })
  x <- sb_analyze_methods(method_analysis_task(), c("boruta", "lasso"))
  expect_identical(called, c("boruta", "lasso"))
  expect_identical(x$status$status, c("failed", "success"))
  expect_identical(x$status$reason, c("native failure", ""))
  expect_identical(x$results$boruta$error_message, "native failure")
  expect_s3_class(tryCatch(sb_selected_features(x, "boruta"), error = identity), "error")
})

test_that("ambiguous requests fail before fitting and empty selections remain empty", {
  task <- method_analysis_task()
  for (methods in list(c("LASSO", "lasso"), "svm_rfe", character(), NA_character_)) {
    expect_s3_class(tryCatch(sb_analyze_methods(task, methods), error = identity),
      "stormbreaker_bad_method")
  }
  for (alpha in c(0, 1, NA_real_)) {
    expect_s3_class(tryCatch(sb_analyze_methods(task, "elastic_net", elastic_net_alpha = alpha),
      error = identity), "stormbreaker_bad_method")
  }
  x <- structure(list(results = list(lasso = list(status = "success",
    feature_result = "selection", features = data.frame(feature = "a", selected = FALSE)))),
    class = "stormbreaker_method_analysis")
  expect_identical(sb_selected_features(x, "lasso"), character())
})
