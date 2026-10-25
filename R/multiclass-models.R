sb_fit_multiclass <- function(task, method, inner_folds, inner_seed,
                              alpha = 0.5, lambda_rule = "lambda.1se", seed = 1L) {
  catalog <- sb_method_catalog("multiclass")
  if (!method %in% catalog$method[catalog$supported] ||
      method %in% sb_multiclass_standalone_selectors()) sb_abort(
    paste(method, "is not implemented for multiclass tasks."),
    "stormbreaker_unsupported_method_task")
  x <- task$data[task$predictors]
  codes <- paste0("class", seq_along(task$classes))
  y <- factor(codes[match(as.character(task$data[[task$outcome_col]]), task$classes)], levels = codes)
  fit <- list(method = method, classes = task$classes, codes = codes,
    features = task$predictors, seed = seed, fold_id = NULL, lambda = NULL,
    warnings = character())
  record_warning <- function(w) fit$warnings <<- c(fit$warnings, conditionMessage(w))
  if (method %in% c("pre", "slope")) fit$fold_id <- sb_multiclass_folds(task, inner_folds, inner_seed)
  if (method == "abess") return(sb_fit_multiclass_abess(fit, x, y))
  if (method %in% sb_multiclass_native_catalog()$method) return(sb_fit_multiclass_native(fit, x, y))
  if (method %in% sb_multiclass_classical_catalog()$method) return(
    sb_fit_multiclass_classical(fit, x, y))
  if (method == "random_forest") {
    if (!requireNamespace("randomForest", quietly = TRUE)) sb_abort(
      "Install randomForest for this method.", "stormbreaker_missing_dependency")
    fit$engine <- withCallingHandlers(withr::with_seed(seed,
      randomForest::randomForest(x, y, ntree = 500L, importance = TRUE)),
      warning = record_warning)
    fit$parameters <- list(ntree = 500L, mtry = fit$engine$mtry, seed = seed)
  } else {
    if (!requireNamespace("glmnet", quietly = TRUE)) sb_abort(
      "Install glmnet for this method.", "stormbreaker_missing_dependency")
    if (!is.character(lambda_rule) || length(lambda_rule) != 1L || is.na(lambda_rule) ||
        !lambda_rule %in% c("lambda.1se", "lambda.min")) sb_abort(
      "Choose lambda.1se or lambda.min.", "stormbreaker_bad_method")
    inner_folds <- sb_check_positive_integer(inner_folds, "inner_folds", minimum = 3L)
    alpha <- switch(method, lasso = 1, ridge = 0, elastic_net = alpha)
    fit$fold_id <- sb_multiclass_folds(task, inner_folds, inner_seed)
    fit$engine <- withCallingHandlers(glmnet::cv.glmnet(as.matrix(x), y, family = "multinomial",
      alpha = alpha, foldid = fit$fold_id,
      type.measure = "deviance", type.multinomial = "ungrouped", standardize = TRUE),
      warning = record_warning)
    fit$lambda <- fit$engine[[lambda_rule]]
    fit$parameters <- list(alpha = alpha,
      lambda_rule = lambda_rule, lambda = fit$lambda, standardize = TRUE,
      type_multinomial = "ungrouped", inner_folds = inner_folds, inner_seed = inner_seed)
  }
  structure(fit, class = "stormbreaker_multiclass_fit")
}

#' Predict probabilities from a saved native multiclass fit
#' @param fit A saved fit at `analysis$results[[method]]$fit` from
#'   [sb_analyze_methods()] with a multiclass task.
#' @param newdata Data frame containing all fitted predictor columns.
#' @return A probability matrix with rows in input order and columns in the
#'   declared class order. No label conversion, refitting or normalization occurs.
#' @export
sb_predict_multiclass <- function(fit, newdata) {
  if (!inherits(fit, "stormbreaker_multiclass_fit")) sb_abort(
    "Supply a saved multiclass fit.", "stormbreaker_bad_multiclass_fit")
  if (!is.data.frame(newdata) || !all(fit$features %in% names(newdata))) sb_abort(
    "New data must contain all fitted predictors.", "stormbreaker_missing_predictors")
  x <- sb_validate_predictor_frame(newdata[fit$features], "newdata")
  if (fit$method == "abess") {
    if (!requireNamespace("abess", quietly = TRUE)) sb_abort(
      "Install abess to predict with this saved fit.", "stormbreaker_missing_dependency")
    probability <- stats::predict(fit$engine, as.matrix(x),
      support.size = fit$engine$best.size, type = "response")[[1]]
  } else if (fit$method %in% sb_multiclass_native_catalog()$method) {
    probability <- sb_predict_multiclass_native(fit, x)
  } else if (fit$method %in% sb_multiclass_classical_catalog()$method) {
    probability <- sb_predict_multiclass_classical(fit, x)
  } else if (fit$method == "random_forest") {
    if (!requireNamespace("randomForest", quietly = TRUE)) sb_abort(
      "Install randomForest to predict with this saved fit.", "stormbreaker_missing_dependency")
    probability <- stats::predict(fit$engine, x, type = "prob")[, fit$codes, drop = FALSE]
  } else {
    if (!requireNamespace("glmnet", quietly = TRUE)) sb_abort(
      "Install glmnet to predict with this saved fit.", "stormbreaker_missing_dependency")
    array <- stats::predict(fit$engine, as.matrix(x), s = fit$lambda, type = "response")
    probability <- matrix(array[, fit$codes, 1L], nrow = nrow(x), ncol = length(fit$codes))
  }
  colnames(probability) <- fit$classes
  rownames(probability) <- rownames(newdata)
  sb_validate_multiclass_probabilities(probability, fit$classes)
  probability
}

sb_multiclass_evidence <- function(fit) {
  if (fit$method %in% sb_multiclass_native_catalog()$method) return(sb_multiclass_native_evidence(fit))
  if (fit$method == "abess") return(sb_multiclass_abess_evidence(fit))
  if (fit$method %in% sb_multiclass_classical_catalog()$method) return(
    sb_multiclass_classical_evidence(fit))
  if (fit$method == "random_forest") {
    native <- fit$engine$importance
    features <- data.frame(feature = fit$features,
      importance = unname(native[fit$features, "MeanDecreaseGini"]), selected = NA)
    features$rank <- rank(-features$importance, ties.method = "min")
    by_class <- do.call(rbind, lapply(seq_along(fit$classes), function(k) data.frame(
      feature = fit$features, class = fit$classes[k],
      importance = unname(native[fit$features, fit$codes[k]]), selected = NA)))
    by_class$rank <- stats::ave(-by_class$importance, by_class$class,
      FUN = function(x) rank(x, ties.method = "min"))
    return(list(features = features, by_class = by_class))
  }
  beta <- stats::coef(fit$engine, s = fit$lambda)
  by_class <- do.call(rbind, lapply(seq_along(fit$classes), function(k) {
    b <- as.matrix(beta[[fit$codes[k]]])[fit$features, 1L]
    data.frame(feature = fit$features, class = fit$classes[k],
      coefficient = unname(b), selected = unname(b != 0))
  }))
  features <- data.frame(feature = fit$features, selected = vapply(fit$features,
    function(gene) any(by_class$selected[by_class$feature == gene]), logical(1)))
  if (fit$method == "ridge") {
    features$selected <- NA
    by_class$selected <- NA
  }
  list(features = features, by_class = by_class)
}

sb_multiclass_plots <- function(fit, evidence, new_data = NULL) {
  if (fit$method %in% sb_multiclass_native_catalog()$method) return(sb_multiclass_native_plots(fit, new_data))
  if (fit$method == "abess") return(sb_multiclass_abess_plots(fit, evidence))
  if (fit$method %in% sb_multiclass_classical_catalog()$method) return(
    sb_multiclass_classical_plots(fit, evidence))
  feature <- importance <- class <- iteration <- error <- lambda <- coefficient <- NULL
  mean_deviance <- standard_error <- NULL
  if (fit$method == "random_forest") {
    tab <- evidence$features
    tab$feature <- factor(tab$feature, levels = tab$feature[order(tab$importance)])
    importance_plot <- ggplot2::ggplot(tab, ggplot2::aes(importance, feature)) +
      ggplot2::geom_col(fill = "#2878A0") + sb_viz_theme() +
      ggplot2::labs(title = "Random forest: Importance", x = "Mean decrease in Gini", y = NULL)
    errors <- fit$engine$err.rate
    history <- data.frame(iteration = rep(seq_len(nrow(errors)), ncol(errors)),
      series = rep(seq_len(ncol(errors)), each = nrow(errors)), error = as.vector(errors))
    history$class <- factor(history$series, levels = seq_len(ncol(errors)),
      labels = c("Overall", paste0("Class: ", fit$classes)))
    oob <- ggplot2::ggplot(history, ggplot2::aes(iteration, error, colour = class)) +
      ggplot2::geom_line() + sb_viz_theme() +
      ggplot2::labs(title = "Random forest: OOB error", x = "Trees", y = "Error", colour = NULL)
    class_data <- evidence$by_class
    class_data$class <- factor(class_data$class, levels = fit$classes)
    class_plot <- ggplot2::ggplot(class_data,
      ggplot2::aes(importance, feature)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::facet_wrap(~class) + sb_viz_theme() +
      ggplot2::labs(title = "Random forest: Class importance",
        x = "OOB permutation accuracy decrease", y = NULL)
    return(list(plots = list(importance = importance_plot, oob_error = oob,
      class_importance = class_plot), tables = list(importance = tab,
        oob_error = history, class_importance = evidence$by_class)))
  }
  engine <- fit$engine$glmnet.fit
  path <- do.call(rbind, lapply(seq_along(fit$classes), function(k) {
    b <- as.matrix(engine$beta[[fit$codes[k]]])
    data.frame(feature = rep(rownames(b), ncol(b)),
      lambda = rep(engine$lambda, each = nrow(b)), coefficient = as.vector(b),
      class = fit$classes[k])
  }))
  path$class <- factor(path$class, levels = fit$classes)
  label <- switch(fit$method, lasso = "LASSO", ridge = "Ridge", elastic_net = "Elastic Net")
  path_plot <- ggplot2::ggplot(path, ggplot2::aes(lambda, coefficient, colour = feature)) +
    ggplot2::geom_line() + ggplot2::scale_x_log10() + ggplot2::facet_wrap(~class) +
    sb_viz_theme() + ggplot2::labs(title = paste0(label, ": Class coefficient paths"),
      x = "Lambda (log scale)", y = "Coefficient", colour = "Feature")
  cv <- data.frame(lambda = fit$engine$lambda, mean_deviance = fit$engine$cvm,
    standard_error = fit$engine$cvsd)
  curve <- ggplot2::ggplot(cv, ggplot2::aes(lambda, mean_deviance)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_deviance - standard_error,
      ymax = mean_deviance + standard_error), width = 0, colour = "grey65") +
    ggplot2::geom_line(colour = "#2878A0") + ggplot2::scale_x_log10() +
    ggplot2::geom_vline(xintercept = fit$lambda, linetype = "dashed") +
    sb_viz_theme() + ggplot2::labs(title = paste0(label, ": Cross-validation"),
      x = "Lambda (log scale)", y = "Multinomial deviance")
  list(plots = list(coefficient_path = path_plot, cv_curve = curve),
    tables = list(coefficient_path = path, cv_curve = cv))
}

sb_analyze_multiclass_method <- function(task, spec, inner_folds, inner_seed,
                                         elastic_net_alpha, lambda_rule) {
  if (spec$method == "lpc") return(sb_analyze_multiclass_lpc(task, spec))
  if (spec$method %in% sb_multiclass_information_catalog()$method) return(
    sb_analyze_multiclass_information(task, spec))
  if (spec$method %in% c("boruta", "varselrf", "vsurf")) return(
    sb_analyze_multiclass_selector(task, spec))
  fit <- sb_fit_multiclass(task, spec$method, inner_folds, inner_seed,
    elastic_net_alpha, lambda_rule)
  evidence <- sb_multiclass_evidence(fit)
  drawn <- tryCatch(sb_multiclass_plots(fit, evidence, task$data[task$predictors]), error = identity)
  failed <- inherits(drawn, "error")
  figure_notes <- list()
  if (!failed) {
    drawn$plots <- lapply(drawn$plots, sb_compact_plot)
    figure_notes <- lapply(drawn$plots, attr, which = "stormbreaker_notes")
  }
  list(status = "success", method = spec$method, label = spec$label,
    feature_result = spec$feature_result, selection_rule = spec$selection_rule,
    features = evidence$features, class_features = evidence$by_class,
    evidence = evidence, parameters = fit$parameters, fit = fit,
    plots = if (failed) list() else drawn$plots,
    diagnostics = if (failed) list() else drawn$tables,
    notes = c(figure_notes, list(interpretation = if (fit$method %in% sb_multiclass_native_catalog()$method ||
        fit$method == "ridge") paste("Full-input exploratory fit, not held-out evaluation.",
      "No gene-selection rule is defined. SVM uses native one-versus-one probability coupling.",
      "Support-vector counts are model structure, not accuracy or gene importance.",
      "Only diagnostics implemented for this engine are returned; no surrogate model is fitted.") else
      if (fit$method %in% sb_multiclass_classical_catalog()$method) paste(
      "Full-input exploration, not held-out evaluation. Class means and SD",
      "are feature distributions, not importance, confidence intervals or selected genes.",
      "Tree importance is native split importance; used_for_split separately",
      "identifies actual primary split variables. No fallback model is fitted.") else paste(
      "Full-input exploration, not held-out evaluation. Class coefficients are",
      "the native multinomial parameterization, not binary odds ratios.",
      "Overall selection is the union of nonzero coefficients across classes.",
      "Forest importance is ranking only; OOB class importance is unscaled",
      "permutation accuracy decrease, not Gini importance. Native OOB sampling",
      "is row-based, not patient-group-isolated; use grouped outer evaluation",
      "for repeated patients."))),
    plot_status = data.frame(view = if (failed) "multiclass_plots" else names(drawn$plots),
      status = if (failed) "error" else rep("available", length(drawn$plots)),
      reason = if (failed) conditionMessage(drawn) else rep("", length(drawn$plots))))
}
