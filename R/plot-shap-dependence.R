#' Relate saved SHAP contributions to observed feature values
#' @param x A successful [sb_additive_shap()] or [sb_tree_shap()] result.
#' @param feature One explained feature name.
#' @param feature_values Feature data with explicit observation IDs.
#' @param id_col ID column, matching explained observations exactly.
#' @param colour_by Optional other explained feature used only for colouring.
#' @return A ggplot object containing all observations without original IDs.
#'   Colouring by another variable does not compute SHAP interaction values.
#'   No threshold, smoothing, causal effect or new explanation is inferred.
#' @export
sb_plot_shap_dependence <- function(x, feature, feature_values, id_col, colour_by = NULL) {
  m <- sb_viz_shap_matrix(x)
  valid_name <- function(z) is.character(z) && length(z) == 1L && !is.na(z) && z %in% colnames(m)
  if (!valid_name(feature) || (!is.null(colour_by) && !valid_name(colour_by))) {
    sb_abort("Choose named explained features.", "stormbreaker_bad_visualization_argument")
  }
  needed <- unique(c(feature, colour_by))
  if (!is.data.frame(feature_values) || !is.character(id_col) || length(id_col) != 1L ||
      is.na(id_col) || !all(c(id_col, needed) %in% names(feature_values))) {
    sb_abort("Supply feature values with an explicit ID column and requested features.",
      "stormbreaker_bad_visualization_argument")
  }
  ids <- as.character(feature_values[[id_col]])
  if (anyNA(ids) || anyDuplicated(ids) || length(ids) != nrow(m) || !setequal(ids, rownames(m))) {
    sb_abort("Feature-value IDs must match explained observations exactly.",
      "stormbreaker_bad_visualization_argument")
  }
  values <- feature_values[match(rownames(m), ids), needed, drop = FALSE]
  if (!all(vapply(values, function(v) is.numeric(v) && all(is.finite(v)), logical(1)))) {
    sb_abort("Feature values must be finite numeric columns.", "stormbreaker_bad_visualization_argument")
  }
  d <- data.frame(observation = seq_len(nrow(m)), feature_value = values[[feature]],
    contribution = unname(m[, feature]))
  feature_value <- contribution <- colour_value <- NULL
  p <- ggplot2::ggplot(d, ggplot2::aes(feature_value, contribution)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey65", linetype = 2)
  if (is.null(colour_by)) {
    p <- p + ggplot2::geom_point(colour = "#2878B5", alpha = .7)
  } else {
    d$colour_value <- values[[colour_by]]
    p$data <- d
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = colour_value), alpha = .7) +
      ggplot2::scale_colour_gradient(low = "#2878B5", high = "#D64E66", name = colour_by)
  }
  p + ggplot2::geom_rug(sides = "b", alpha = .2) +
    ggplot2::labs(title = paste("SHAP:", feature),
      x = feature, y = paste0("SHAP (", sb_viz_score_label(x$score_scale), ")"),
      subtitle = paste(nrow(m), "saved observations | original feature units"),
      caption = "Predictive association, not causality. No clinical cutoff inferred. Colouring is not an interaction estimate. Individual values may be sensitive.") +
    sb_shap_theme()
}
