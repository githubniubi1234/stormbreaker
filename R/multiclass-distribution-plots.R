sb_multiclass_distribution_plots <- function(fit) {
  m <- fit$engine
  value <- density <- class <- feature <- NULL
  if (fit$method == "kernel_naive_bayes") {
    d <- do.call(rbind, lapply(fit$features, function(feature) {
      do.call(rbind, lapply(seq_along(fit$codes), function(k) {
        native <- m$tables[[feature]][[fit$codes[k]]]
        data.frame(feature = feature, class = fit$classes[k],
          value = native$x, density = native$y)
      }))
    }))
    p <- ggplot2::ggplot(d, ggplot2::aes(value, density, colour = class)) +
      ggplot2::geom_line() + ggplot2::facet_wrap(~feature, scales = "free", ncol = 2) +
      sb_viz_theme() + ggplot2::labs(title = "Class-conditional densities", x = "Feature value",
        y = "Density", colour = NULL,
        caption = "Saved kernel density grids used by the naive Bayes model; no new smoothing or density estimation is performed.")
    return(list(plots = list(class_densities = p), tables = list(class_densities = d)))
  }
  d <- data.frame(feature = rep(rownames(m$means), ncol(m$means)),
    class = rep(fit$classes[match(colnames(m$means), fit$codes)], each = nrow(m$means)),
    value = as.vector(m$means))
  p <- ggplot2::ggplot(d, ggplot2::aes(value, feature)) + ggplot2::geom_point(colour = "#2878A0") +
    ggplot2::facet_wrap(~class) + sb_viz_theme() +
    ggplot2::labs(title = "Class means", x = "Mean", y = NULL,
      caption = "Original fitted class means in input units, not variable importance or uncertainty intervals.")
  correlations <- do.call(rbind, lapply(seq_along(m$classes), function(k) {
    covariance <- (1 - m$regularization[2]) * m$covariances[, , k] + m$regularization[2] * m$covpooled
    covariance <- (1 - m$regularization[1]) * covariance +
      m$regularization[1] * mean(diag(covariance)) * diag(nrow(covariance))
    correlation <- stats::cov2cor(covariance)
    data.frame(feature = rep(rownames(correlation), ncol(correlation)),
      other = rep(colnames(correlation), each = nrow(correlation)),
      class = fit$classes[match(m$classes[k], fit$codes)], correlation = as.vector(correlation))
  }))
  other <- correlation <- NULL
  heat <- ggplot2::ggplot(correlations, ggplot2::aes(other, feature, fill = correlation)) +
    ggplot2::geom_tile() + ggplot2::facet_wrap(~class) +
    ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#C77335", limits = c(-1, 1)) +
    sb_viz_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(title = "Regularized class correlations", x = NULL, y = NULL, fill = "Correlation",
      caption = "Within-class correlations from the covariance matrices used for prediction, after both native RDA regularizations; not pooled empirical correlations.")
  list(plots = list(class_means = p, class_correlations = heat),
    tables = list(class_means = d, class_correlations = correlations))
}
