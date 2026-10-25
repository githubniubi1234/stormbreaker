sb_multiclass_classical_catalog <- function() {
  data.frame(method = c("lda", "qda", "gaussian_naive_bayes", "decision_tree"),
    label = c("LDA", "QDA", "Gaussian naive Bayes", "Decision tree"),
    supported = TRUE, engine = c("MASS_lda", "MASS_qda", "e1071_naiveBayes", "rpart_class"),
    feature_result = c(rep("model_features", 3), "ranking"),
    selection_rule = "No standalone gene-selection rule")
}

sb_multiclass_classical_package <- function(method) {
  switch(method, lda = , qda = "MASS", gaussian_naive_bayes = "e1071", decision_tree = "rpart")
}

sb_fit_multiclass_classical <- function(fit, x, y) {
  package <- sb_multiclass_classical_package(fit$method)
  if (!requireNamespace(package, quietly = TRUE)) sb_abort(
    paste("Install", package, "for this method."), "stormbreaker_missing_dependency")
  record_warning <- function(w) fit$warnings <<- c(fit$warnings, conditionMessage(w))
  if (fit$method %in% c("lda", "qda")) {
    fit$engine <- withCallingHandlers(
      if (fit$method == "lda") MASS::lda(x, y, method = "moment") else
        MASS::qda(x, y, method = "moment"), warning = record_warning)
    fit$parameters <- list(method = "moment", prior = fit$engine$prior)
  } else if (fit$method == "gaussian_naive_bayes") {
    fit$engine <- withCallingHandlers(e1071::naiveBayes(x, y, laplace = 0),
      warning = record_warning)
    fit$parameters <- list(laplace = 0, prediction_eps = 0, prediction_threshold = 0.001)
  } else {
    fit$feature_codes <- paste0("feature", seq_along(fit$features))
    names(x) <- fit$feature_codes
    x$.outcome <- y
    fit$engine <- withCallingHandlers(rpart::rpart(.outcome ~ ., data = x,
      method = "class", parms = list(split = "gini"),
      control = rpart::rpart.control(cp = .01, minsplit = 20L, minbucket = 7L,
        maxdepth = 30L, maxsurrogate = 0L, usesurrogate = 0L, xval = 0L)),
      warning = record_warning)
    attr(fit$engine$terms, ".Environment") <- baseenv()
    fit$parameters <- list(cp = .01, minsplit = 20L, minbucket = 7L, maxdepth = 30L,
      split = "gini", maxsurrogate = 0L, usesurrogate = 0L, xval = 0L)
  }
  structure(fit, class = "stormbreaker_multiclass_fit")
}

sb_predict_multiclass_classical <- function(fit, x) {
  package <- sb_multiclass_classical_package(fit$method)
  if (!requireNamespace(package, quietly = TRUE)) sb_abort(
    paste("Install", package, "to predict with this saved fit."), "stormbreaker_missing_dependency")
  if (fit$method %in% c("lda", "qda")) return(
    stats::predict(fit$engine, x)$posterior[, fit$codes, drop = FALSE])
  if (fit$method == "gaussian_naive_bayes") return(
    stats::predict(fit$engine, x, type = "raw", eps = fit$parameters$prediction_eps,
      threshold = fit$parameters$prediction_threshold)[, fit$codes, drop = FALSE])
  names(x) <- fit$feature_codes
  stats::predict(fit$engine, x, type = "prob")[, fit$codes, drop = FALSE]
}

sb_multiclass_classical_evidence <- function(fit) {
  features <- data.frame(feature = fit$features, selected = NA)
  if (fit$method == "decision_tree") {
    importance <- fit$engine$variable.importance
    features$importance <- 0
    present <- match(fit$feature_codes, names(importance))
    features$importance[!is.na(present)] <- unname(importance[present[!is.na(present)]])
    features$rank <- rank(-features$importance, ties.method = "min")
    features$used_for_split <- fit$feature_codes %in% fit$engine$frame$var
    return(list(features = features, by_class = NULL))
  }
  by_class <- do.call(rbind, lapply(seq_along(fit$classes), function(k) {
    tab <- data.frame(feature = fit$features, class = fit$classes[k], selected = NA)
    if (fit$method == "gaussian_naive_bayes") {
      tab$mean <- vapply(fit$features, function(gene) fit$engine$tables[[gene]][fit$codes[k], 1], numeric(1))
      tab$sd <- vapply(fit$features, function(gene) fit$engine$tables[[gene]][fit$codes[k], 2], numeric(1))
    } else tab$mean <- unname(fit$engine$means[fit$codes[k], fit$features])
    tab
  }))
  list(features = features, by_class = by_class)
}

sb_multiclass_classical_plots <- function(fit, evidence) {
  feature <- importance <- class <- mean <- sd <- NULL
  label <- sb_multiclass_classical_catalog()$label[
    match(fit$method, sb_multiclass_classical_catalog()$method)]
  if (fit$method == "decision_tree") {
    tab <- evidence$features
    tab$feature <- factor(tab$feature, levels = tab$feature[order(tab$importance)])
    plot <- ggplot2::ggplot(tab, ggplot2::aes(importance, feature)) +
      ggplot2::geom_col(fill = "#2878A0") + sb_viz_theme() +
      ggplot2::labs(title = "Decision tree: Importance", x = "Native split importance", y = NULL)
    tree <- sb_multiclass_tree_structure(fit)
    return(list(plots = list(importance = plot, tree_structure = sb_plot_rpart_tree(tree)),
      tables = list(importance = tab, tree_structure = tree$table)))
  }
  tab <- evidence$by_class
  tab$class <- factor(tab$class, levels = fit$classes)
  plot <- ggplot2::ggplot(tab, ggplot2::aes(mean, feature)) + ggplot2::geom_point(colour = "#2878A0") +
    ggplot2::facet_wrap(~class) + sb_viz_theme() +
    ggplot2::labs(title = paste0(label, ": Class means"), x = "Mean", y = NULL)
  if (fit$method == "gaussian_naive_bayes") plot <- plot +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = mean - sd, xmax = mean + sd),
      orientation = "y", width = .2, colour = "#2878A0") +
    ggplot2::labs(title = "Gaussian naive Bayes: Class distributions", x = "Mean \u00b1 SD")
  list(plots = list(class_statistics = plot), tables = list(class_statistics = tab))
}
