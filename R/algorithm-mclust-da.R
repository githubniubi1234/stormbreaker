#' Inspect original model-based discriminant components
#' @param fit Successful EDDA or MclustDA final fit.
#' @param type centers, correlations, or mixture_weights.
#' @return Original component means, within-component correlations, or
#'   within-class mixing weights. Components are latent model structure, not
#'   clinical subgroups. Covariances are original fitted parameters, not
#'   empirical sample correlations. BIC paths absent from compact fits are not
#'   reconstructed. No mixture refitting or component reselection occurs.
#' @export
sb_mclust_da_diagnostics <- function(fit, type = "centers") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% sb_mclust_da_learner_ids() || !is.character(type) ||
      length(type) != 1L || is.na(type) || !type %in% c("centers", "correlations", "mixture_weights")) sb_abort(
    "Choose centers, correlations or mixture_weights for EDDA/MclustDA.", "stormbreaker_unsupported_diagnostic")
  e <- s$engine
  sb_validate_mclust_da_prediction_engine(e, s$params)
  p <- length(s$features)
  tab <- do.call(rbind, lapply(names(e$models), function(cls) {
    model <- e$models[[cls]]
    par <- model$parameters
    g <- length(par$pro)
    means <- matrix(par$mean, nrow = p, ncol = g)
    if (type == "centers") return(data.frame(class = cls,
      component = rep(seq_len(g), each = p), feature = rep(s$features, g), mean = as.numeric(means)))
    if (type == "mixture_weights") return(data.frame(class = cls, component = seq_len(g),
      weight = as.numeric(par$pro), class_prior = unname(e$prop[cls]), model = model$modelName))
    covariance <- par$variance$sigma
    if (p == 1L) covariance <- array(rep(par$variance$sigmasq, length.out = g), c(1L, 1L, g))
    if (!is.array(covariance) || !identical(as.integer(dim(covariance)), as.integer(c(p, p, g))) ||
        any(!is.finite(covariance))) sb_abort("Original component covariances are invalid.", "stormbreaker_invalid_diagnostic")
    do.call(rbind, lapply(seq_len(g), function(j) {
      m <- matrix(covariance[, , j], p, p)
      if (any(diag(m) <= 0)) sb_abort("Original component variances must be positive.", "stormbreaker_invalid_diagnostic")
      correlation <- stats::cov2cor(m)
      data.frame(class = cls, component = j, feature = rep(s$features, p),
        feature2 = rep(s$features, each = p), correlation = as.numeric(correlation))
    }))
  }))
  structure(list(type = paste0("mclust_da_", type), table = tab,
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "Original saved mclust discriminant mixture parameters",
    interpretation = "Latent components are not clinical subgroups. Correlations are fitted component correlations, not raw sample correlations. Mixing weights are conditional on class. No refitting or BIC reconstruction."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_mclust_da_diagnostic <- function(x) {
  mean <- feature <- feature2 <- correlation <- component <- weight <- NULL
  if (x$type == "mclust_da_centers") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(mean, feature)) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::facet_wrap(~class + component, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "Mixture component centers", x = "Mean", y = NULL)
  } else if (x$type == "mclust_da_correlations") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(feature2, feature, fill = correlation)) +
      ggplot2::geom_tile(colour = "white", linewidth = .2) +
      ggplot2::facet_wrap(~class + component, labeller = ggplot2::label_both) +
      ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#B86236", limits = c(-1, 1),
        guide = ggplot2::guide_colourbar(display = "rectangles")) +
      ggplot2::labs(title = "Fitted component correlations", x = NULL, y = NULL, fill = "Correlation")
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(factor(component), weight)) +
      ggplot2::geom_col(fill = "#2878A0", width = .6) + ggplot2::facet_wrap(~class, labeller = ggplot2::label_both) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = "Within-class mixing weights", x = "Component", y = "Weight")
  }
  p <- p + sb_viz_theme()
  if (x$type == "mclust_da_correlations") p <- p + ggplot2::coord_fixed() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank())
  attr(p, "stormbreaker_suggested_size") <- c(width = 12, height = 6)
  p
}
