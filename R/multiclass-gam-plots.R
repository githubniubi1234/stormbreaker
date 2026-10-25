sb_multiclass_gam_plots <- function(fit) {
  m <- fit$engine
  lpi <- attr(m$formula, "lpi")
  curves <- rugs <- list()
  for (j in seq_along(m$smooth)) {
    smooth <- m$smooth[[j]]
    term <- smooth$term
    index <- smooth$first.para:smooth$last.para
    lp <- which(vapply(lpi, function(z) all(index %in% z), logical(1)))
    values <- m$model[[term]]
    grid <- stats::setNames(data.frame(seq(min(values), max(values), length.out = 100)), term)
    X <- mgcv::PredictMat(smooth, grid)
    estimate <- drop(X %*% stats::coef(m)[index])
    se <- sqrt(rowSums((X %*% m$Vp[index, index, drop = FALSE]) * X))
    feature <- fit$features[match(term, fit$feature_codes)]
    curves[[j]] <- data.frame(feature = feature, class = fit$classes[lp + 1L],
      value = grid[[1]], estimate = estimate, se = se,
      lower = estimate - stats::qnorm(.975) * se, upper = estimate + stats::qnorm(.975) * se)
    rugs[[j]] <- data.frame(feature = feature, class = fit$classes[lp + 1L], value = values)
  }
  d <- do.call(rbind, curves)
  rug <- do.call(rbind, rugs)
  value <- estimate <- lower <- upper <- NULL
  plots <- list()
  for (cl in unique(d$class)) {
    plots[[paste0("smooths_", cl)]] <- ggplot2::ggplot(d[d$class == cl, ], ggplot2::aes(value, estimate)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "#DDEBF3") +
      ggplot2::geom_hline(yintercept = 0, colour = "grey65", linetype = 2) +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::geom_rug(data = rug[rug$class == cl, ], ggplot2::aes(x = value),
        inherit.aes = FALSE, sides = "b", alpha = .25) +
      ggplot2::facet_wrap(~feature, ncol = 1, scales = "free") + sb_viz_theme() +
      ggplot2::labs(title = paste(cl, "vs", fit$classes[1], "smooth effects"),
        x = "Feature value", y = "Partial log-odds",
        caption = paste0("Native centered GAM smooth terms for ", cl, " versus reference ", fit$classes[1],
          ". Pointwise approximate 95% Bayesian intervals use the fitted coefficient covariance conditional on smoothing parameters; not simultaneous bands or predicted probabilities. Rug marks retain every training value. No refitting, extrapolation or causal interpretation."))
  }
  list(plots = plots, tables = list(smooths = d, training_rug = rug))
}
