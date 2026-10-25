sb_multiclass_pamr_plots <- function(fit) {
  m <- fit$engine
  centers <- pamr::pamr.predict(m, newx = NULL,
    threshold = fit$parameters$threshold, type = "centroid")
  contrast <- (centers - m$centroid.overall) / m$sd
  d <- data.frame(feature = rep(fit$features[match(rownames(contrast), fit$feature_codes)], ncol(contrast)),
    class = rep(fit$classes[match(colnames(contrast), fit$codes)], each = nrow(contrast)),
    contrast = as.vector(contrast))
  d$feature <- factor(d$feature, levels = rev(fit$features))
  feature <- contrast <- NULL
  p <- ggplot2::ggplot(d, ggplot2::aes(contrast, feature)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey65") +
    ggplot2::geom_col(fill = "#2878A0") + ggplot2::facet_wrap(~class) +
    sb_viz_theme() + ggplot2::labs(title = "Shrunken centroids",
      x = "Standardized centroid difference", y = NULL,
      caption = paste0("Native PAMR centroid contrasts at fitted threshold ", fit$parameters$threshold,
        "; centered on the overall centroid and divided by the native within-class scale, as in pamr.plotcen. All features, including zero contrasts, are retained. This is not a cross-validation curve; the configured fit uses one fixed threshold."))
  list(plots = list(shrunken_centroids = p), tables = list(shrunken_centroids = d))
}
