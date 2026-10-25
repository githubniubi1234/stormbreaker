sb_multiclass_boosting_catalog <- function() {
  data.frame(method = c("xgboost", "xgboost_dart", "xgboost_linear", "lightgbm",
    "lightgbm_dart", "lightgbm_goss", "lightgbm_rf", "catboost"),
    label = c("XGBoost", "XGBoost DART", "XGBoost linear", "LightGBM", "LightGBM DART",
      "LightGBM GOSS", "LightGBM RF", "CatBoost"), supported = TRUE,
    engine = c(rep("xgboost", 3), rep("lightgbm", 4), "catboost"),
    feature_result = c("ranking", "ranking", "model_features", rep("ranking", 5)),
    selection_rule = "No standalone gene-selection rule")
}

sb_fit_multiclass_boosting <- function(fit, x, y) {
  method <- fit$method
  catalog <- sb_multiclass_boosting_catalog()
  package <- catalog$engine[match(method, catalog$method)]
  if (!requireNamespace(package, quietly = TRUE)) sb_abort(
    paste("Install", package, "for this method."), "stormbreaker_missing_dependency")
  fit$package <- package
  fit$feature_codes <- paste0("feature", seq_along(fit$features))
  x <- as.matrix(x)
  colnames(x) <- fit$feature_codes
  target <- as.integer(y) - 1L
  fit$importance <- NULL
  record_warning <- function(w) fit$warnings <<- c(fit$warnings, conditionMessage(w))
  withCallingHandlers(withr::with_seed(fit$seed, {
    if (package == "xgboost") {
      id <- switch(method, xgboost = "xgboost_binary",
        xgboost_dart = "xgboost_binary_dart", xgboost_linear = "xgboost_binary_linear")
      params <- sb_xgboost_params(id)
      params$objective <- "multi:softprob"
      params$num_class <- length(fit$codes)
      params$eval_metric <- "mlogloss"
      params$seed <- fit$seed
      booster <- xgboost::xgb.train(params = params[setdiff(names(params), c("nrounds", "seed"))],
        data = xgboost::xgb.DMatrix(x, label = target), nrounds = params$nrounds, verbose = 0)
      fit$engine <- xgboost::xgb.save.raw(booster)
      if (method != "xgboost_linear") {
        importance <- xgboost::xgb.importance(feature_names = fit$feature_codes, model = booster)
        fit$importance <- stats::setNames(importance$Gain, importance$Feature)
        fit$importance_label <- "Gain share"
      }
    } else if (package == "lightgbm") {
      id <- switch(method, lightgbm = "lightgbm_binary", lightgbm_dart = "lightgbm_binary_dart",
        lightgbm_goss = "lightgbm_binary_goss", lightgbm_rf = "lightgbm_binary_rf")
      params <- sb_lightgbm_binary_params(learner = id, seed = fit$seed)
      params$objective <- "multiclass"
      params$num_class <- length(fit$codes)
      params$metric <- "multi_logloss"
      training_params <- params[setdiff(names(params), "nrounds")]
      # Native alias avoids BiocGenerics format.list truncating multi_logloss.
      training_params$metric <- "multiclass"
      booster <- lightgbm::lgb.train(params = training_params,
        data = lightgbm::lgb.Dataset(x, label = target), nrounds = params$nrounds, verbose = -1L)
      fit$engine <- booster$save_model_to_string()
      importance <- lightgbm::lgb.importance(booster, percentage = FALSE)
      fit$importance <- stats::setNames(importance$Gain, importance$Feature)
      fit$importance_label <- "Total split gain"
    } else {
      params <- sb_catboost_binary_params(seed = fit$seed)
      params$loss_function <- "MultiClass"
      params$eval_metric <- "MultiClass"
      params$classes_count <- length(fit$codes)
      params$boost_from_average <- FALSE
      booster <- catboost::catboost.train(catboost::catboost.load_pool(x, label = target),
        test_pool = NULL, params = params)
      fit$engine <- serialize(booster, NULL)
      fit$importance <- stats::setNames(catboost::catboost.get_feature_importance(booster,
        type = "PredictionValuesChange", thread_count = 1L), fit$feature_codes)
      fit$importance_label <- "Prediction value change"
    }
    fit$parameters <- params
  }), warning = record_warning)
  structure(fit, class = "stormbreaker_multiclass_fit")
}

sb_predict_multiclass_boosting <- function(fit, x) {
  if (!requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "to predict with this saved fit."), "stormbreaker_missing_dependency")
  x <- as.matrix(x)
  colnames(x) <- fit$feature_codes
  if (fit$package == "xgboost") {
    booster <- xgboost::xgb.load.raw(fit$engine)
    p <- stats::predict(booster, x, reshape = TRUE)
  } else if (fit$package == "lightgbm") {
    booster <- lightgbm::lgb.load(model_str = fit$engine)
    p <- stats::predict(booster, x)
  } else {
    p <- catboost::catboost.predict(unserialize(fit$engine), catboost::catboost.load_pool(x),
      prediction_type = "Probability", thread_count = 1L)
  }
  if (length(p) != nrow(x) * length(fit$codes)) sb_abort(
    "The native booster returned an unexpected probability shape.", "stormbreaker_bad_multiclass_probabilities")
  matrix(p, nrow = nrow(x), ncol = length(fit$codes))
}

sb_multiclass_boosting_evidence <- function(fit) {
  tab <- data.frame(feature = fit$features, selected = NA)
  if (!is.null(fit$importance)) {
    tab$importance <- 0
    index <- match(fit$feature_codes, names(fit$importance))
    tab$importance[!is.na(index)] <- unname(fit$importance[index[!is.na(index)]])
    tab$rank <- rank(-tab$importance, ties.method = "min")
  }
  list(features = tab, by_class = NULL)
}

sb_multiclass_boosting_plots <- function(fit) {
  if (is.null(fit$importance)) return(list(plots = list(), tables = list()))
  tab <- sb_multiclass_boosting_evidence(fit)$features
  tab$feature <- factor(tab$feature, levels = tab$feature[order(tab$importance)])
  catalog <- sb_multiclass_boosting_catalog()
  label <- catalog$label[match(fit$method, catalog$method)]
  feature <- importance <- NULL
  plot <- ggplot2::ggplot(tab, ggplot2::aes(importance, feature)) +
    ggplot2::geom_col(fill = "#2878A0") + sb_viz_theme() +
    ggplot2::labs(title = paste0(label, ": Importance"), x = fit$importance_label, y = NULL)
  list(plots = list(importance = plot), tables = list(importance = tab))
}
