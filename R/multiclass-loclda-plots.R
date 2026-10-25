sb_multiclass_loclda_plots <- function(fit, new_data) {
  if (is.null(new_data)) stop("Supply new_data for local-LDA neighborhoods.", call. = FALSE)
  x <- as.matrix(sb_validate_predictor_frame(new_data[fit$features], "new_data"))
  tab <- sb_loclda_neighbor_table(fit$engine, x)
  tab$class <- fit$classes[match(tab$class, fit$codes)]
  distance <- weight <- class <- used <- NULL
  p <- ggplot2::ggplot(tab, ggplot2::aes(distance, weight, colour = class, shape = used)) +
    ggplot2::geom_point(size = 1.8, alpha = .8) + ggplot2::facet_wrap(~query) +
    ggplot2::scale_shape_manual(values = c("FALSE" = 4, "TRUE" = 16),
      labels = c("FALSE" = "Excluded", "TRUE" = "Used")) + sb_viz_theme() +
    ggplot2::labs(title = "Local-LDA neighborhoods", x = "Distance", y = "Local weight",
      colour = "Class", shape = NULL,
      caption = "Original distance-radius neighborhoods, including distance ties. Classes represented by at most one neighbor are excluded by the native local-LDA rule, but remain visible. All supplied queries are retained across pages. Self-neighbors remain for training queries. Local fitting weights are not probabilities, feature importance or validation performance; reference rows and classes may be sensitive.")
  list(plots = list(neighborhoods = p), tables = list(neighborhoods = tab))
}
