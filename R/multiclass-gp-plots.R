sb_multiclass_gp_plots <- function(fit) {
  m <- fit$engine
  pairs <- utils::combn(seq_along(m@lev), 2L)
  classes <- fit$classes[match(m@lev, fit$codes)]
  d <- do.call(rbind, lapply(seq_along(m@alpha), function(k) {
    data.frame(pair = paste(classes[pairs[, k]], collapse = " vs "),
      basis = seq_along(m@alpha[[k]]), reference_row = m@alphaindex[[k]],
      weight = as.numeric(m@alpha[[k]]))
  }))
  basis <- weight <- NULL
  p <- ggplot2::ggplot(d, ggplot2::aes(basis, weight)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey65") +
    ggplot2::geom_point(colour = "#2878A0", size = 1.2) +
    ggplot2::facet_wrap(~pair) + sb_viz_theme() +
    ggplot2::labs(title = "Pairwise kernel weights", x = "Native basis", y = "Kernel weight",
      caption = "All original kernlab Gaussian-process one-vs-one kernel coefficients. Positive latent scores favor the first named class. These weights belong to reference observations, not genes; they are not feature importance or SHAP. No posterior uncertainty band is invented, and no observations or fitted coefficients are removed.")
  list(plots = list(kernel_weights = p), tables = list(kernel_weights = d))
}
