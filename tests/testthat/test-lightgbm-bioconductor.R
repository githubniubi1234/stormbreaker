test_that("LightGBM logloss aliases survive Bioconductor list formatting", {
  skip_if_not_installed("lightgbm")
  skip_if_not_installed("BiocGenerics")
  requireNamespace("BiocGenerics", quietly = TRUE)
  d <- iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(iris)[1:4])
  methods <- c("lightgbm", "lightgbm_dart", "lightgbm_goss", "lightgbm_rf")
  for (method in methods) {
    fit <- sb_fit_multiclass(task, method, 3, 1)
    x <- as.matrix(d[task$predictors])
    colnames(x) <- fit$feature_codes
    params <- fit$parameters[setdiff(names(fit$parameters), c("nrounds", "metric"))]
    # The objective's native default is multi_logloss, independently of aliases.
    reference <- lightgbm::lgb.train(params,
      lightgbm::lgb.Dataset(x, label = as.integer(factor(d$Species, levels = task$classes)) - 1L),
      nrounds = fit$parameters$nrounds, verbose = -1L)
    expect_match(fit$engine, "\\[metric: multi_logloss\\]")
    expect_equal(unname(sb_predict_multiclass(fit, d)), unname(predict(reference, x)), tolerance = 0)
  }
  x <- d[1:100, 1:4]
  outcome <- rep(0:1, each = 50)
  for (id in sb_lightgbm_binary_learner_ids()) {
    params <- sb_lightgbm_binary_params(learner = id)
    fit <- sb_fit_lightgbm_binary(x, outcome, params, id)
    reference <- lightgbm::lgb.train(params[setdiff(names(params), c("nrounds", "metric"))],
      lightgbm::lgb.Dataset(as.matrix(x), label = outcome), nrounds = params$nrounds, verbose = -1L)
    expect_match(fit$booster$save_model_to_string(), "\\[metric: binary_logloss\\]")
    expect_equal(predict(fit$booster, as.matrix(x)), predict(reference, as.matrix(x)), tolerance = 0)
  }
})
