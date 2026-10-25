test_that("identical learner fits share within, never across, analysis sets", {
  skip_if_not_installed("survivalsvm")
  skip_if_not_installed("mboost")
  skip_if_not_installed("plsRglm")
  native <- sb_fit_runtime_learner
  calls <- list()
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls[[length(calls) + 1L]] <<- rownames(x)
    native(learner, x, endpoint)
  })
  features <- c("gene_a", "gene_b", "gene_c")
  binary_ids <- c("logistic", "gaussian_naive_bayes", sb_specialized_binary_linear_learner_ids())
  for (id in c(binary_ids, "coxph", sb_survivalsvm_learner_ids())) {
    if (id %in% binary_ids) {
      d <- sb_example_binary(120L, 1208L)
      task <- sb_binary_task(d, "patient_id", "outcome", features)
    } else {
      d <- sb_example_survival(120L, 11032L)
      task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
        prediction_horizon = 7)
    }
    pipelines <- list(sb_pipeline("none", id), sb_pipeline("fixed_features", id,
      selector_features = features, feature_set_id = "all"))
    calls <- list()
    plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 11L, reuse_selectors = FALSE)
    expect_length(calls, 6L)
    calls <- list()
    shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 11L)
    expect_length(calls, 3L)
    expect_equal(shared, plain)
    expect_identical(unique(shared$comparison$status), "success")
    expected <- lapply(shared$results[[1]]$folds, function(f) rownames(task$data)[f$analysis_rows])
    expect_identical(calls, expected)
    run <- function(enabled) sb_with_selector_reuse(
      lapply(pipelines, function(p) sb_fit_final(task, p, capture_diagnostics = TRUE)),
      enabled = enabled)
    calls <- list(); plain_final <- run(FALSE)
    expect_length(calls, 2L)
    calls <- list(); shared_final <- run(TRUE)
    expect_length(calls, 1L)
    expect_equal(shared_final, plain_final)
    for (i in seq_along(pipelines)) expect_identical(
      sb_predict_final(unserialize(serialize(shared_final[[i]], NULL)),
        d[c("patient_id", features)], "patient_id"),
      sb_predict_final(plain_final[[i]], d[c("patient_id", features)], "patient_id"))
  }
})

test_that("CV reuse separates full input identities and retries warning or failed fits", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls <<- calls + 1L
    if (identical(learner$params$mode, "warning")) warning("native warning")
    if (identical(learner$params$mode, "error")) stop("native failure")
    list(coefficient = c(a = 1, b = 2))
  })
  learner <- sb_logistic_learner_adapter()
  x <- data.frame(a = 1:6, b = 6:1)
  ids <- paste0("p", 1:6)
  endpoint <- list(task_type = "binary", outcome = rep(0:1, 3))
  result <- sb_with_selector_reuse({
    fit <- function(l = learner, z = x, e = endpoint, p = ids) {
      sb_fit_reusable_runtime_learner(l, z, e, p, "none")
    }
    first <- fit(); changed <- first; changed$coefficient[] <- 0
    expect_identical(fit(), first)
    fit(z = x[6:1, ]); fit(z = x[c("b", "a")]); fit(p = rev(ids))
    changed_endpoint <- endpoint; changed_endpoint$outcome <- 1L - endpoint$outcome
    fit(e = changed_endpoint)
    changed_learner <- learner; changed_learner$params <- list(mode = "different")
    fit(l = changed_learner)
    for (i in 1:2) {
      changed_learner$params$mode <- "warning"
      expect_warning(fit(l = changed_learner), "native warning")
      changed_learner$params$mode <- "error"
      expect_error(fit(l = changed_learner), "native failure")
    }
    sb_fit_final_runtime_learner(learner, x, endpoint, ids, "none")
  }, diagnostics = TRUE)
  expect_identical(calls, 11L)
  expect_identical(result$reuse$cv_learner_hits, 1L)
  expect_equal(result$reuse$cv_learner_entries, 6L)
  expect_identical(result$reuse$final_learner_misses, 1L)
  expect_null(.sb_selector_reuse_context$active)
})
