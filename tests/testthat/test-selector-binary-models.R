binary_model_selector_task <- function(grouped = FALSE) {
  withr::local_seed(5)
  d <- as.data.frame(matrix(rnorm(180 * 8), 180, 8))
  names(d) <- paste0("gene", 1:8)
  d$id <- paste0("p", 1:180)
  d$y <- rbinom(180, 1, plogis(1.5 * d$gene1 - d$gene2))
  if (grouped) {
    d$patient <- d$id
    d <- d[rep(seq_len(nrow(d)), each = 2), ]
    d$id <- paste0("row", seq_len(nrow(d)))
  }
  sb_binary_task(d, "id", "y", paste0("gene", 1:8), group_col = if (grouped) "patient" else NULL)
}

test_that("binary GamBoost selection matches native active spline terms", {
  skip_if_not_installed("mboost")
  task <- binary_model_selector_task()
  x <- task$data[task$predictors]
  names(x) <- paste0("gene ", seq_len(ncol(x)))
  p <- sb_pipeline("gamboost_binary_nonzero", "logistic")
  selected <- sb_fit_binary_model_selector(x, task$data$y, p$selector, p$selector_params)
  work <- x; names(work) <- paste0("feature_", seq_len(ncol(x)))
  native <- sb_fit_classical_binary(work, task$data$y, "gamboost_binary",
    sb_classical_binary_params("gamboost_binary"))$engine
  beta <- stats::coef(native)
  active <- names(beta)[vapply(beta, function(b) any(b != 0), logical(1))]
  index <- match(sub("^bbs\\(([^,)]+).*", "\\1", active), names(work))
  expect_setequal(selected$selected_features, names(x)[index])
  expect_identical(sb_doctor_selector_menu("binary")[["GamBoost"]], p$selector)
  expect_null(selected$engine)
  expect_identical(selected, unserialize(serialize(selected, NULL)))
})

test_that("binary model selectors round-trip and retain genuine selected genes", {
  task <- binary_model_selector_task()
  for (id in sb_binary_model_selector_ids()) {
    pipeline <- sb_pipeline(id, "gbm_binary")
    expect_identical(sb_lint_pipeline(pipeline)$status, "admitted")
    fit <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success")
    if (id == "svm_rfe_binary") {
      evidence <- fit$selection$state$evidence
      expect_identical(fit$selection$selected_features,
        evidence$subsets[[which.max(evidence$curve$accuracy)]])
    } else expect_identical(fit$selection$selected_features,
      fit$selection$scores$feature[fit$selection$scores$selected])
    expect_identical(all(fit$selection$selected_features %in% task$predictors), TRUE)
    copied <- unserialize(serialize(fit, NULL))
    pred <- sb_predict_final(copied, task$data[c("id", task$predictors)], "id")
    expect_identical(pred$status, "success")
  }
  a <- sb_pipeline("stepglm_binary_aic", "logistic", selector_options = list(direction = "forward"))
  b <- sb_pipeline("stepglm_binary_aic", "logistic", selector_options = list(direction = "backward"))
  expect_identical(a$id == b$id, FALSE)
  expect_identical(sb_lint_pipeline(a)$status, "admitted")
  bad <- tryCatch(sb_pipeline("rf_binary_importance", "logistic", selector_options = list(thresold = 1)), error = identity)
  expect_s3_class(bad, "stormbreaker_bad_selector_options")
  empty <- sb_fit_final(task, sb_pipeline("rf_binary_importance", "logistic", selector_options = list(threshold = 1e6)))
  expect_identical(empty$status, "failed")
  expect_identical(empty$error_stage, "selection")
})

test_that("selectors execute within outer and learner-inner training folds", {
  task <- binary_model_selector_task()
  original <- sb_fit_binary_model_selector
  seen <- list()
  local_mocked_bindings(sb_fit_binary_model_selector = function(x, outcome, id, params, group_ids = NULL) {
    ids <- rownames(x)
    if (all(ids %in% rownames(task$data))) ids <- task$data$id[match(ids, rownames(task$data))]
    seen[[length(seen) + 1L]] <<- ids
    original(x, outcome, id, params, group_ids)
  })
  for (id in sb_binary_model_selector_ids()) {
    seen <- list()
    cv <- sb_run_cv(task, sb_pipeline(id, "glmnet_logistic", inner_folds = 3, inner_seed = 7,
      selector_options = if (id %in% c("lasso_binary_nonzero", "enet_binary_nonzero")) list(lambda_rule = "lambda.min") else NULL), folds = 3)
    if (id == "svm_rfe_binary") {
      # This fixture selects one feature in a nested fold; native glmnet
      # cannot fit that design. Failure must remain visible, without padding.
      expect_identical(cv$status, "failed")
      failed <- cv$runs[cv$runs$status == "failed", ]
      expect_identical(failed$error_class, "stormbreaker_glmnet_logistic_tuning_failure")
      expect_match(failed$error_message, "cannot retain 1 or less variables", fixed = TRUE)
    } else expect_identical(cv$status, "success")
    expect_gt(length(seen), 3)
    expect_identical(all(lengths(seen) < task$n), TRUE)
    expect_identical(any(lengths(seen) < floor(task$n * 2/3)), TRUE)
    for (ids in seen) expect_identical(any(vapply(cv$folds,
      function(fold) all(ids %in% fold$analysis_ids), logical(1))), TRUE)
  }
})

test_that("LASSO inner folds keep repeated patients together", {
  task <- binary_model_selector_task(TRUE)
  original <- sb_fit_binary_model_selector
  grouped_calls <- 0L
  local_mocked_bindings(sb_fit_binary_model_selector = function(x, outcome, id, params, group_ids = NULL) {
    result <- original(x, outcome, id, params, group_ids)
    if (id == "lasso_binary_nonzero") {
      expect_length(group_ids, nrow(x))
      expect_identical(all(vapply(split(result$evidence$fold_id, group_ids),
        function(z) length(unique(z)) == 1L, logical(1))), TRUE)
      grouped_calls <<- grouped_calls + 1L
    }
    result
  })
  p <- sb_pipeline("lasso_binary_nonzero", "glmnet_logistic", inner_folds = 3, inner_seed = 7,
    selector_options = list(lambda_rule = "lambda.min"))
  fit <- sb_fit_final(task, p)
  expect_identical(fit$status, "success")
  cv <- sb_run_cv(task, p, folds = 3)
  expect_identical(cv$status, "success")
  expect_gt(grouped_calls, 5)
})

test_that("new selectors refit in SVM calibration and preserve patient groups", {
  task <- binary_model_selector_task(TRUE)
  original <- sb_fit_binary_model_selector
  calls <- 0L
  local_mocked_bindings(sb_fit_binary_model_selector = function(x, outcome, id, params, group_ids = NULL) {
    result <- original(x, outcome, id, params, group_ids)
    calls <<- calls + 1L
    if (id == "lasso_binary_nonzero") {
      expect_length(group_ids, nrow(x))
      expect_identical(all(vapply(split(result$evidence$fold_id, group_ids),
        function(z) length(unique(z)) == 1L, logical(1))), TRUE)
    }
    result
  })
  for (id in sb_binary_model_selector_ids()) {
    calls <- 0L
    fit <- sb_fit_final(task, sb_pipeline(id, "svm_linear_binary",
      selector_options = if (id == "lasso_binary_nonzero") list(lambda_rule = "lambda.min") else NULL))
    expect_identical(fit$status, "success")
    expect_gt(calls, 1L)
    expect_length(fit$tuning$inner_selections, 5L)
  }
})

test_that("RF scores match native importance with original gene names", {
  task <- binary_model_selector_task()
  x <- task$data[task$predictors]
  names(x)[1:2] <- c("HLA-DRA", "gene with space")
  params <- sb_binary_model_selector_options("rf_binary_importance")
  selected <- sb_fit_binary_model_selector(x, task$data$y, "rf_binary_importance", params)
  work <- x
  names(work) <- paste0("feature_", seq_len(ncol(x)))
  withr::local_seed(params$seed)
  native <- randomForest::randomForest(x = work, y = factor(task$data$y, levels = c(0, 1)),
    ntree = params$trees, importance = TRUE, keep.forest = FALSE)
  score <- as.numeric(randomForest::importance(native, type = 1, scale = FALSE)[, 1])
  expect_equal(selected$scores$score, score)
  expect_identical(selected$selected_features, names(x)[score > params$threshold])
  for (id in c("stepglm_binary_aic", "glmboost_binary_nonzero")) {
    result <- sb_fit_binary_model_selector(x, task$data$y, id, sb_binary_model_selector_options(id))
    expect_identical(all(result$selected_features %in% names(x)), TRUE)
    expect_identical(all(result$evidence$selected_features %in% names(x)), TRUE)
  }
})
