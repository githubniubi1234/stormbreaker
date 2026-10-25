sb_multiclass_neural_plots <- function(fit) {
  m <- fit$engine
  labels <- c("Bias", fit$features, paste("Hidden", seq_len(m$n[2L])), fit$classes)
  destination <- rep(seq_len(length(m$nconn) - 1L) - 1L, diff(m$nconn))
  tab <- data.frame(source = labels[m$conn + 1L],
    destination = labels[destination + 1L], weight = m$wts,
    layer = ifelse(destination <= sum(m$n[1:2]), "Input to hidden", "Hidden to output"))
  tab$source <- factor(tab$source, levels = unique(labels))
  tab$destination <- factor(tab$destination, levels = unique(labels))
  source <- weight <- NULL
  plot <- ggplot2::ggplot(tab, ggplot2::aes(destination, source, fill = weight)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~layer, scales = "free") +
    sb_neural_weight_scale() +
    sb_viz_theme() + ggplot2::labs(title = "Neural network weights", x = NULL,
      y = NULL, fill = "Weight", caption = paste(
        "Original fitted connections, including biases. Inputs were standardized;",
        "outputs use softmax. Connection weights are not gene importance or causal effects.",
        "Weight colors use symmetric pseudo-log scaling (sigma = 1), with original weights on the legend; no coefficient clipping."))
  list(plots = list(connection_weights = plot), tables = list(connection_weights = tab))
}
