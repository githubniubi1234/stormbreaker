sb_multiclass_tree_native_plots <- function(fit) {
  forest <- fit$method == "rrf"
  native <- if (forest) fit$engine$importance else C50::C5imp(fit$engine, metric = "usage", pct = TRUE)
  d <- data.frame(feature = fit$features[match(rownames(native), fit$feature_codes)],
    importance = as.numeric(native[, 1L]))
  d$feature <- factor(d$feature, levels = d$feature[order(d$importance)])
  feature <- importance <- NULL
  plot <- ggplot2::ggplot(d, ggplot2::aes(importance, feature)) +
    ggplot2::geom_col(fill = "#2878A0") + sb_viz_theme() +
    ggplot2::labs(title = if (forest) "Regularized forest importance" else "C5.0 feature usage",
      x = if (forest) "Mean decrease in Gini" else "Usage (%)", y = NULL,
      caption = if (forest) "Original impurity importance, including zero-valued features; not permutation importance." else
        "Original C5.0 usage percentages; these are not mutually exclusive and need not sum to 100%.")
  result <- list(plots = list(importance = plot), tables = list(importance = d))
  if (forest) {
    errors <- fit$engine$err.rate
    names <- colnames(errors)
    at <- match(names, fit$codes)
    names[!is.na(at)] <- fit$classes[at[!is.na(at)]]
    d <- data.frame(trees = rep(seq_len(nrow(errors)), ncol(errors)),
      series = rep(names, each = nrow(errors)), error = as.vector(errors))
    trees <- error <- series <- NULL
    result$plots$oob_error <- ggplot2::ggplot(d, ggplot2::aes(trees, error, colour = series)) +
      ggplot2::geom_line() + sb_viz_theme() + ggplot2::labs(title = "Out-of-bag error",
        x = "Trees", y = "Error rate", colour = NULL,
        caption = "Original cumulative out-of-bag classification errors, overall and by class; not external validation.")
    result$tables$oob_error <- d
  }
  result
}
