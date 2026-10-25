test_that("matching glmnet preparation shares native tuning and refits", {
  skip_if_not_installed("glmnet")
  preparations <- 0L
  # Observe the actual native lambda searches, not only the wrapper invocation.
  native_binary <- sb_select_glmnet_logistic_lambda
  native_survival <- sb_select_glmnet_cox_lambda
  local_mocked_bindings(
    sb_select_glmnet_logistic_lambda = function(...) {
      preparations <<- preparations + 1L; native_binary(...)
    },
    sb_select_glmnet_cox_lambda = function(...) {
      preparations <<- preparations + 1L; native_survival(...)
    })
  features <- c("gene_a", "gene_b", "gene_c")
  for (type in c("binary", "survival")) {
    if (type == "binary") {
      d <- sb_example_binary(120L, 1208L)
      task <- sb_binary_task(d, "patient_id", "outcome", features)
      id <- "glmnet_logistic"
    } else {
      d <- sb_example_survival(120L, 11032L)
      task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
        prediction_horizon = 7)
      id <- "glmnet_cox"
    }
    pipelines <- list(sb_pipeline("none", id, inner_folds = 3L, inner_seed = 101L),
      sb_pipeline("fixed_features", id, selector_features = features,
        feature_set_id = "all", inner_folds = 3L, inner_seed = 101L))
    preparations <- 0L
    plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 11L, reuse_selectors = FALSE)
    expect_identical(preparations, 6L)
    preparations <- 0L
    shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 11L)
    expect_identical(preparations, 3L)
    expect_identical(shared, plain)
    expect_identical(unique(shared$comparison$status), "success")
    run <- function(enabled) sb_with_selector_reuse(
      lapply(pipelines, function(p) sb_fit_final(task, p)), enabled = enabled,
      diagnostics = TRUE)
    preparations <- 0L; plain_final <- run(FALSE)
    expect_identical(preparations, 2L)
    preparations <- 0L; shared_final <- run(TRUE)
    expect_identical(preparations, 1L)
    expect_equal(shared_final$result, plain_final$result)
    expect_identical(shared_final$reuse$final_learner_hits, 1L)
    for (i in seq_along(pipelines)) expect_identical(
      sb_predict_final(unserialize(serialize(shared_final$result[[i]], NULL)),
        d[c("patient_id", features)], "patient_id"),
      sb_predict_final(plain_final$result[[i]], d[c("patient_id", features)], "patient_id"))
  }
})

test_that("preparation identity includes upstream data, selection and grouping", {
  skip_if_not_installed("glmnet")
  d <- sb_example_binary(120L, 1208L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "binary", outcome = d$outcome)
  p <- sb_pipeline("univ_binary_auc_top_k", "glmnet_logistic",
    selector_n_features = 2L, inner_folds = 3L, inner_seed = 101L)
  runtime <- sb_resolve_pipeline_runtime(p)
  selected <- sb_fit_runtime_selector_uncached(runtime$selector, x, endpoint)
  selected_x <- runtime$selector$apply(selected, x)
  sc <- list(selector = runtime$selector, fitted_selector = selected, x = x)
  native <- sb_select_glmnet_logistic_lambda
  calls <- 0L
  local_mocked_bindings(sb_select_glmnet_logistic_lambda = function(...) {
    calls <<- calls + 1L; native(...)
  })
  result <- sb_with_selector_reuse({
    prepare <- function(pipeline = p, context = sc, groups = d$patient_id) {
      sb_prepare_runtime_learner(sb_resolve_pipeline_runtime(pipeline)$learner,
        pipeline, selected_x, endpoint, d$patient_id, context, groups)
    }
    a <- prepare(); b <- prepare()
    expect_equal(a$tuning, b$tuning)
    fixed <- sb_pipeline("fixed_features", "glmnet_logistic",
      selector_features = names(selected_x), feature_set_id = "same_outer_genes",
      inner_folds = 3L, inner_seed = 101L)
    prepare(fixed)
    changed <- sc
    omitted <- setdiff(names(x), names(selected_x))
    changed$x[[omitted]] <- -changed$x[[omitted]]
    prepare(context = changed)
    prepare(groups = rev(d$patient_id))
  }, diagnostics = TRUE)
  expect_identical(calls, 4L)
  expect_identical(result$reuse$preparation_hits, 1L)
})

test_that("preparation warnings and failures are never hidden by reuse", {
  calls <- 0L
  mode <- "warning"
  local_mocked_bindings(sb_prepare_runtime_learner_uncached = function(...) {
    calls <<- calls + 1L
    if (mode == "error") stop("native tuning failure")
    list(learner = list(), tuning = list(warnings = "native retained warning"))
  })
  p <- sb_pipeline("none", "glmnet_logistic", inner_folds = 3L, inner_seed = 1L)
  learner <- sb_resolve_pipeline_runtime(p)$learner
  result <- sb_with_selector_reuse({
    prepare <- function() sb_prepare_runtime_learner(learner, p,
      data.frame(a = 1:6, b = 6:1), list(task_type = "binary", outcome = rep(0:1, 3)),
      as.character(1:6))
    for (i in 1:2) expect_identical(prepare()$tuning$warnings, "native retained warning")
    mode <- "error"
    for (i in 1:2) expect_error(prepare(), "native tuning failure")
  }, diagnostics = TRUE)
  expect_identical(calls, 4L)
  expect_identical(result$reuse$preparation_hits, 0L)
})
