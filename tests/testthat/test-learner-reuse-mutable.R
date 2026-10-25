test_that("shared boosting fits do not alias mutable native model environments", {
  skip_if_not_installed("mboost")
  d <- sb_example_survival(80L, 11032L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "survival", time = d$time, event = d$event)
  learner <- sb_mboost_cox_learner_adapter("glmboost_cox")
  result <- sb_with_selector_reuse({
    fit <- function() sb_fit_reusable_runtime_learner(learner, x, endpoint, d$patient_id, "none")
    first <- fit()
    original <- learner$predict(first, x)
    mboost::mstop(first$engine) <- 1L
    second <- fit()
    expect_identical(as.integer(mboost::mstop(second$engine)), 100L)
    expect_identical(learner$predict(second, x), original)
    mboost::mstop(second$engine) <- 2L
    third <- fit()
    expect_identical(as.integer(mboost::mstop(third$engine)), 100L)
    expect_identical(learner$predict(third, x), original)
  }, diagnostics = TRUE)
  expect_identical(result$reuse$cv_learner_hits, 2L)
})

test_that("shared XGBoost fits own independent native boosters", {
  skip_if_not_installed("xgboost")
  d <- sb_example_binary(80L, 1208L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "binary", outcome = d$outcome)
  learner <- sb_xgboost_learner_adapter("xgboost_binary")
  for (phase in c("cv", "final")) {
    result <- sb_with_selector_reuse({
      fit <- function() sb_fit_reusable_runtime_learner(learner, x, endpoint,
        d$patient_id, "none", phase = phase)
      first <- fit()
      original <- learner$predict(first, x)
      xgboost::xgb.attr(first$booster, "consumer_note") <- "first"
      second <- fit()
      expect_null(xgboost::xgb.attr(second$booster, "consumer_note"))
      expect_identical(learner$predict(second, x), original)
      xgboost::xgb.attr(second$booster, "consumer_note") <- "second"
      third <- fit()
      expect_null(xgboost::xgb.attr(third$booster, "consumer_note"))
      expect_identical(learner$predict(third, x), original)
      expect_identical(xgboost::xgb.attr(first$booster, "consumer_note"), "first")
    }, diagnostics = TRUE)
    expect_identical(result$reuse[[paste0(phase, "_learner_hits")]], 2L)
  }
})

test_that("shared LightGBM fits isolate native iteration changes", {
  skip_if_not_installed("lightgbm")
  d <- sb_example_binary(100L, 1208L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "binary", outcome = d$outcome)
  learner <- sb_lightgbm_binary_learner_adapter()
  out <- sb_with_selector_reuse({
    fit <- function() sb_fit_reusable_runtime_learner(learner, x, endpoint,
      d$patient_id, "none")
    first <- fit()
    original <- learner$predict(first, x)
    original_text <- first$booster$save_model_to_string()
    first$booster$rollback_one_iter()
    expect_false(identical(first$booster$save_model_to_string(), original_text))
    second <- fit()
    expect_identical(second$booster$save_model_to_string(), original_text)
    expect_identical(learner$predict(second, x), original)
    # Loaded boosters need not retain rollback history. Use the already changed
    # native booster to verify stale source-precision diagnostics are not reused.
    second$booster <- first$booster
    expect_false(identical(second$booster$save_model_to_string(), original_text))
    changed <- second
    changed$reused_native_importance <- NULL
    expect_identical(sb_extract_lightgbm_tree_model(second),
      sb_extract_lightgbm_tree_model(changed))
    third <- fit()
    expect_identical(third$booster$save_model_to_string(), original_text)
    expect_identical(learner$predict(third, x), original)
  }, diagnostics = TRUE)
  expect_identical(out$reuse$cv_learner_hits, 2L)
})
