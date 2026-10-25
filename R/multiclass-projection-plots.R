sb_multiclass_projection_plots <- function(fit) {
  m <- fit$engine
  result <- list(plots = list(), tables = list())
  add_matrix <- function(values, rows, columns, key, title, axis, note) {
    d <- data.frame(label = rep(rows, ncol(values)), panel = rep(columns, each = nrow(values)),
      value = as.vector(values))
    d$label <- factor(d$label, levels = rev(unique(rows)))
    value <- label <- NULL
    p <- ggplot2::ggplot(d, ggplot2::aes(value, label)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey75") +
      ggplot2::geom_point(colour = "#2878A0") + ggplot2::facet_wrap(~panel) +
      sb_viz_theme() + ggplot2::labs(title = title, x = axis, y = NULL, caption = note)
    result$plots[[key]] <<- p
    result$tables[[key]] <<- d
  }
  if (fit$method %in% c("fda", "mda", "pda")) {
    rows <- rownames(m$means)
    if (fit$method == "mda") {
      for (k in seq_along(m$assign.theta)) {
        at <- m$assign.theta[[k]]
        rows[at] <- paste(fit$classes[match(names(m$assign.theta)[k], fit$codes)],
          "subclass", seq_along(at))
      }
    } else rows <- fit$classes[match(rows, fit$codes)]
    add_matrix(m$means, rows, paste("Axis", seq_len(ncol(m$means))), "projected_centroids",
      "Discriminant centroids", "Discriminant coordinate",
      "Saved class or subclass centroids in the fitted discriminant space; these are not original gene expression means.")
    if (fit$method == "mda") {
      priors <- do.call(rbind, lapply(seq_along(m$sub.prior), function(k) {
        data.frame(class = fit$classes[match(names(m$sub.prior)[k], fit$codes)],
          subclass = names(m$sub.prior[[k]]), probability = as.numeric(m$sub.prior[[k]]))
      }))
      subclass <- probability <- NULL
      result$plots$subclass_weights <- ggplot2::ggplot(priors, ggplot2::aes(subclass, probability)) +
        ggplot2::geom_col(fill = "#2878A0") + ggplot2::facet_wrap(~class) + sb_viz_theme() +
        ggplot2::labs(title = "Subclass weights", x = NULL, y = "Within-class probability",
          caption = "Saved mixture probabilities within each outcome class, not overall class prevalence.")
      result$tables$subclass_weights <- priors
    }
  } else if (fit$method == "slda") {
    add_matrix(m$scores, fit$features, paste("Axis", seq_len(ncol(m$scores))), "projection_loadings",
      "Projection loadings", "Loading", "Original stabilized-LDA projection matrix; signs depend on axis orientation and are not causal effects.")
    add_matrix(m$mylda$means, fit$classes[match(rownames(m$mylda$means), fit$codes)],
      paste("Axis", seq_len(ncol(m$mylda$means))), "projected_centroids", "Projected class means",
      "Projected coordinate", "Fitted class means after the native stabilizing projection, not original gene expression means.")
  } else {
    add_matrix(t(m$beta), fit$features[match(colnames(m$beta), fit$feature_codes)],
      fit$classes[match(rownames(m$beta), fit$codes)], "discriminant_coefficients",
      "Discriminant coefficients", "Coefficient", "Original class-specific discriminant-score coefficients, not standalone gene importance.")
    d <- data.frame(parameter = names(m$regularization), intensity = unname(m$regularization))
    parameter <- intensity <- NULL
    result$plots$shrinkage <- ggplot2::ggplot(d, ggplot2::aes(intensity, parameter)) +
      ggplot2::geom_point(colour = "#2878A0") + sb_viz_theme() +
      ggplot2::labs(title = "Shrinkage intensities", x = "Intensity", y = NULL,
        caption = "Native shrinkage estimates for correlation, variance and class frequencies; not a tuning-search history.")
    result$tables$shrinkage <- d
  }
  result
}
