test_that("all nine XGBoost variants preserve native scores and zero-use features", {
  skip_if_not_installed("xgboost")
  for (id in sb_xgboost_learner_ids()) {
    binary <- id %in% sb_xgboost_binary_learner_ids()
    d <- if (binary) make_sb_binary_fixture(n = 100L, seed = 108L) else make_sb_fixture(n = 100L, seed = 8612L)
    d$constant <- 0
    features <- c("gene_a", "gene_b", "gene_c", "constant")
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    original <- fit$learner_state
    q <- sb_boosting_diagnostics(fit, "importance")
    native_model <- sb_load_xgboost_state_model(fit$learner_state, id)$booster
    native <- as.data.frame(xgboost::xgb.importance(feature_names = features, model = native_model))
    expect_identical(q$native_table, native)
    expect_identical(q$table$feature, features)
    score <- if (grepl("linear", id, fixed = TRUE)) "Weight" else "Gain"
    expect_equal(q$table$value[match(native$Feature, features)], native[[score]], tolerance = 0)
    expect_equal(q$table$value[q$table$feature == "constant"], 0)
    expect_identical(fit$learner_state, original)
    err <- tryCatch(sb_boosting_diagnostics(fit, "training_loss"), error = identity)
    expect_s3_class(err, "stormbreaker_unavailable_diagnostic")
    if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(
      ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(q)), "ggplot_built")
  }
})

test_that("GBM final refits preserve the selected count and native predictions", {
  skip_if_not_installed("gbm")
  for (signal in c(FALSE, TRUE)) {
    d <- make_sb_fixture(n = 100L, seed = 8612L)
    if (signal) {
      d$time <- exp(2 - d$gene_a)
      d$event[] <- 1L
    }
    x <- d[c("gene_a", "gene_b", "gene_c")]
    params <- sb_gbm_cox_cv_params(inner_folds = 3L, maximum_trees = 30L)
    model <- sb_fit_gbm_cox_cv(x, d$time, d$event, params)
    count <- model$n_trees
    expect_equal(count, which.min(model$tuning$cv_error))
    expect_equal(count == 1L, !signal)
    formula_data <- sb_additional_survival_formula_data(x, d$time, d$event)
    native <- withr::with_seed(params$seed, gbm::gbm(
      formula = formula_data$formula, data = formula_data$data,
      distribution = params$distribution, n.trees = count,
      interaction.depth = params$interaction_depth, n.minobsinnode = params$min_obs_node,
      shrinkage = params$shrinkage, bag.fraction = params$bag_fraction,
      cv.folds = if (count == 1L) 0L else params$inner_folds,
      n.cores = 1L, keep.data = FALSE, verbose = FALSE))
    expect_equal(sb_predict_gbm_cox(model, x),
      unname(as.numeric(stats::predict(native, newdata = x, n.trees = count, type = "link"))), tolerance = 0)
  }
})

test_that("GBM paths distinguish native final fits and original CV including a one-tree optimum", {
  skip_if_not_installed("gbm")
  for (id in sb_gbm_diagnostic_ids()) {
    binary <- id == "gbm_binary"
    d <- if (binary) make_sb_binary_fixture(n = 100L, seed = 108L) else make_sb_fixture(n = 100L, seed = 8612L)
    features <- c("gene_a", "gene_b", "gene_c")
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else sb_surv_task(d, "patient_id", "time", "event", features, "month")
    pipeline <- if (id == "gbm_cox_cv") sb_pipeline(learner = id, gbm_cv_max_trees = 30L, learner_tuning_folds = 3L) else sb_pipeline(learner = id)
    fit <- sb_fit_final(task, pipeline, capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    original <- fit$learner_state
    engine <- if (binary) fit$diagnostic_artifacts$learner$engine else fit$diagnostic_artifacts$learner$booster
    q <- sb_boosting_diagnostics(fit)
    expect_equal(q$table$loss, engine$train.error, tolerance = 0)
    expect_equal(nrow(q$table), engine$n.trees)
    importance <- sb_boosting_diagnostics(fit, "importance")
    expect_equal(importance$table$value, unname(gbm::relative.influence(engine,
      n.trees = engine$n.trees, scale. = FALSE, sort. = FALSE)))
    if (id == "gbm_cox_cv") {
      cv <- sb_boosting_diagnostics(fit, "cv_curve")
      expect_equal(cv$table$loss, fit$diagnostic_artifacts$learner$tuning$cv_error)
      expect_equal(cv$selected_iteration, which.min(cv$table$loss))
      expect_equal(cv$selected_iteration, 1L)
      expect_equal(nrow(cv$table), 30L)
      expect_equal(nrow(q$table), 1L)
      expect_null(engine$cv.error)
      if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(
        ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(cv)), "ggplot_built")
    }
    expect_identical(fit$learner_state, original)
    if (requireNamespace("ggplot2", quietly = TRUE)) for (z in list(q, importance)) expect_s3_class(
      ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)), "ggplot_built")
    fit$diagnostic_artifacts <- NULL
    expect_identical(sb_boosting_diagnostics(fit, "importance"), importance)
    if (binary) {
      err <- tryCatch(sb_boosting_diagnostics(fit), error = identity)
      expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
    } else expect_identical(sb_boosting_diagnostics(fit), q)
  }
})
