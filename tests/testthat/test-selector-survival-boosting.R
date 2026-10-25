test_that("survival boosting selectors fit, predict and expose their evidence", {
  d <- sb_example_survival(150, 82)
  t <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
  for (method in c("GBM", "glmBoost", "SuperPC", "GamBoost")) {
    p <- surv_pipeline("Cox", surv_selector(method))
    expect_equal(sb_validate_pipeline_contract(p, "survival"), p)
    f <- sb_fit_final(t, p, capture_diagnostics = TRUE)
    expect_identical(f$status, "success")
    s <- f$selection$state
    expect_identical(s$selected_features, s$scores$feature[s$scores$selected])
    expect_equal(nrow(s$scores), 3)
    if (method != "SuperPC") expect_gt(nrow(s$evidence$training_loss), 0)
    expect_identical(sb_predict_final(f, d, "patient_id")$status, "success")
    plots <- sb_plots(f)
    expect_length(plots$unavailable, 0)
    expect_s3_class(plots$plots$selector_1, "ggplot")
    if (method != "SuperPC") expect_s3_class(plots$plots$selector_1_training_loss, "ggplot")
    cv <- sb_run_cv(t, p, folds = 3, seed = 82)
    expect_identical(cv$status, "success")
  }
  chain <- sb_fit_final(t, surv_pipeline("Cox", surv_selector(c("GBM", "glmBoost"))))
  expect_identical(chain$status, "success")
  expect_length(chain$selection$state$stages, 2)
})

test_that("GamBoost selects native active spline terms, not linear glmBoost terms", {
  skip_if_not_installed("mboost")
  d <- sb_example_survival(150, 82)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("gene A", "gene-B", "gene C")
  p <- sb_pipeline("gamboost_cox_nonzero", "coxph")
  adapter <- sb_survival_boosting_selector_adapter(p$selector, p$selector_params)
  fitted <- adapter$fit(x, d$time, d$event)
  work <- x
  names(work) <- paste0("feature_", seq_len(ncol(x)))
  native <- sb_fit_mboost_cox(work, d$time, d$event,
    "gamboost_cox", p$selector_params)$engine
  beta <- stats::coef(native)
  active <- names(beta)[vapply(beta, function(b) any(b != 0), logical(1))]
  native_features <- sub("^bbs\\(([^,)]+).*", "\\1", active)
  expect_setequal(fitted$selected_features, names(x)[match(native_features, names(work))])
  expect_identical(adapter$apply(fitted, x), x[fitted$selected_features])
  expect_identical(unserialize(serialize(fitted, NULL)), fitted)
  expect_null(fitted$engine)
  expect_identical(sb_selector_is_learned(p$selector), TRUE)
})

test_that("GamBoost is refitted inside glmnet outer and inner analysis folds", {
  skip_if_not_installed("mboost")
  d <- sb_example_survival(150, 82)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month")
  original <- sb_fit_mboost_cox
  seen <- list()
  local_mocked_bindings(sb_fit_mboost_cox = function(x, time, event, learner, params) {
    seen[[length(seen) + 1L]] <<- d$patient_id[match(x[[1L]], d$gene_a)]
    original(x, time, event, learner, params)
  })
  p <- sb_pipeline("gamboost_cox_nonzero", "glmnet_cox", inner_folds = 3, inner_seed = 7)
  cv <- sb_run_cv(task, p, folds = 3, seed = 82)
  expect_identical(cv$status, "success")
  expect_gt(length(seen), 3L)
  expect_identical(all(lengths(seen) < nrow(d)), TRUE)
  expect_identical(any(lengths(seen) < 100L), TRUE)
  for (ids in seen) expect_identical(any(vapply(cv$folds,
    function(fold) all(ids %in% fold$analysis_ids), logical(1))), TRUE)
})
