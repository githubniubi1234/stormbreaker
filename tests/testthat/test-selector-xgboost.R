test_that("XGBoost selectors retain exactly native positive-gain features in both tasks", {
  skip_if_not_installed("xgboost")
  d <- sb_example_survival(150, 82)
  features <- c("gene_a", "gene_b", "gene_c")
  x <- d[features]; names(x) <- c("gene A", "gene-B", "gene C")
  d$outcome <- d$event
  original <- sb_fit_xgboost
  seen <- list()
  local_mocked_bindings(sb_fit_xgboost = function(x, endpoint, learner, params) {
    seen[[length(seen) + 1L]] <<- d$patient_id[match(x[[1]], d$gene_a)]
    original(x, endpoint, learner, params)
  })
  for (surv in c(FALSE, TRUE)) {
    id <- if (surv) "xgboost_cox_importance" else "xgboost_binary_importance"
    learner <- if (surv) "xgboost_cox" else "xgboost_binary"
    params <- if (surv) sb_survival_boosting_selector_params(id) else sb_binary_model_selector_options(id)
    a <- sb_resolve_selector_adapter(id, params)
    set.seed(11); before <- .Random.seed
    s <- if (surv) a$fit(x, d$time, d$event) else a$fit(x, d$outcome)
    expect_identical(.Random.seed, before)
    work <- x; names(work) <- paste0("feature_", 1:3)
    native <- sb_fit_xgboost(work, if (surv) list(time = d$time, event = d$event) else d$outcome,
      learner, sb_xgboost_params(learner))
    imp <- xgboost::xgb.importance(feature_names = names(work), model = native$booster)
    expect_setequal(s$selected_features, names(x)[match(imp$Feature[imp$Gain > 0], names(work))])
    expect_null(s$engine)
    expect_identical(a$apply(unserialize(serialize(s, NULL)), x), x[s$selected_features])
    task <- if (surv) sb_surv_task(d, "patient_id", "time", "event", features, "month") else
      sb_binary_task(d, "patient_id", "outcome", features)
    p <- sb_pipeline(id, if (surv) "glmnet_cox" else "glmnet_logistic", inner_folds = 3, inner_seed = 7)
    seen <- list()
    cv <- sb_run_cv(task, p, folds = 3, seed = 82)
    expect_identical(cv$status, "success")
    expect_gt(length(seen), 3L)
    expect_identical(any(lengths(seen) < 100L), TRUE)
    for (ids in seen) expect_identical(any(vapply(cv$folds,
      function(fold) all(ids %in% fold$analysis_ids), logical(1))), TRUE)
    expect_identical(cv$provenance$xgboost_version, as.character(packageVersion("xgboost")))
    fit <- sb_fit_final(task, p)
    expect_identical(fit$status, "success")
    expect_identical(fit$provenance$xgboost_version, as.character(packageVersion("xgboost")))
    pred <- sb_predict_final(unserialize(serialize(fit, NULL)), d[c("patient_id", features)], "patient_id")
    expect_identical(pred$status, "success")
  }
})
