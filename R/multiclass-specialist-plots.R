sb_multiclass_specialist_plots <- function(fit) {
  m <- fit$engine
  feature <- value <- component <- class <- NULL
  if (fit$method %in% c("edda", "mclustda")) {
    means <- do.call(rbind, lapply(seq_along(m$models), function(k) {
      a <- m$models[[k]]$parameters$mean
      data.frame(feature = rep(fit$features[match(rownames(a), fit$feature_codes)], ncol(a)),
        component = rep(seq_len(ncol(a)), each = nrow(a)),
        class = fit$classes[match(names(m$models)[k], fit$codes)], value = as.vector(a))
    }))
    p <- ggplot2::ggplot(means, ggplot2::aes(value, feature, colour = factor(component))) +
      ggplot2::geom_point() + ggplot2::facet_wrap(~class) + sb_viz_theme() +
      ggplot2::labs(title = "Mixture component means", x = "Mean", y = NULL, colour = "Component",
        caption = "Original selected Gaussian mixture component means in input units. Component numbers are local to each class, not matched biological subtypes.")
    result <- list(plots = list(component_means = p), tables = list(component_means = means))
    bic <- do.call(rbind, lapply(seq_along(m$models), function(k) {
      a <- m$models[[k]]$BIC
      if (!is.matrix(a)) return(NULL)
      data.frame(components = rep(rownames(a), ncol(a)), covariance = rep(colnames(a), each = nrow(a)),
        class = fit$classes[match(names(m$models)[k], fit$codes)], BIC = as.vector(a))
    }))
    if (!is.null(bic)) {
      components <- covariance <- BIC <- NULL
      result$plots$bic_search <- ggplot2::ggplot(bic, ggplot2::aes(components, covariance,
        fill = ifelse(is.finite(BIC), BIC, NA_real_))) +
        ggplot2::geom_tile() + ggplot2::facet_wrap(~class) + sb_viz_theme() +
        ggplot2::scale_fill_gradient(low = "#EFF4F8", high = "#2878A0", na.value = "grey85") +
        ggplot2::labs(title = "Mixture model selection", x = "Components", y = "Covariance model", fill = "BIC",
          caption = "Saved candidate BIC values, with larger values preferred by mclust. Grey cells have no finite native estimate; the full table retains those unavailable candidates. This is model selection, not validation accuracy.")
      result$tables$bic_search <- bic
    }
    return(result)
  }
  if (fit$method == "hdda") {
    values <- as.matrix(m$ev)
    d <- data.frame(class = rep(fit$classes[match(m$kname, fit$codes)], ncol(values)),
      axis = rep(seq_len(ncol(values)), each = nrow(values)), value = as.vector(values))
    d$retained <- d$axis <= rep(as.numeric(m$d), ncol(values))
    axis <- retained <- NULL
    p <- ggplot2::ggplot(d, ggplot2::aes(axis, value)) +
      ggplot2::geom_line(colour = "grey60") + ggplot2::geom_point(ggplot2::aes(shape = retained), colour = "#2878A0") +
      ggplot2::facet_wrap(~class) + sb_viz_theme() +
      ggplot2::labs(title = "Class eigenvalue spectra", x = "Axis", y = "Eigenvalue", shape = "Retained",
        caption = "All saved class eigenvalues and the native retained dimensions, on the fitted scaled-input representation. No new PCA, threshold selection or truncation is performed.")
    return(list(plots = list(eigenvalue_spectrum = p), tables = list(eigenvalue_spectrum = d)))
  }
  a <- m$hda.loadings
  d <- data.frame(feature = rep(fit$features[match(rownames(a), fit$feature_codes)], ncol(a)),
    axis = rep(seq_len(ncol(a)), each = nrow(a)), value = as.vector(a))
  d$retained <- d$axis <= m$reduced.dimension
  axis <- retained <- NULL
  p <- ggplot2::ggplot(d, ggplot2::aes(factor(axis), feature, fill = value)) +
    ggplot2::geom_tile() + ggplot2::facet_wrap(~retained, scales = "free_x", labeller = ggplot2::as_labeller(
      c(`TRUE` = "Used for classification", `FALSE` = "Remaining axes"))) +
    ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#C77335") + sb_viz_theme() +
    ggplot2::labs(title = "HDA projection loadings", x = "Axis", y = NULL, fill = "Loading",
      caption = "All original projection loadings, separating the configured classification dimensions from the remaining saved axes. Axis signs are arbitrary; these are not gene effect sizes.")
  list(plots = list(projection_loadings = p), tables = list(projection_loadings = d))
}
