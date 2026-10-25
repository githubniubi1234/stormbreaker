#' Inspect original weighted-subspace forest diagnostics
#' @param fit A successful weighted-subspace forest final fit.
#' @param type importance, tree_nodes or tree_oob. Importance and OOB require
#'   original capture_diagnostics; compact OOB placeholders are never used.
#' @return Native mean-decrease-information-gain-ratio importance, total node
#'   counts per tree, or individual-tree OOB errors. Individual-tree OOB is
#'   not cumulative ensemble OOB or external performance. Artificial constant
#'   bridge features, if present, remain separately identified in the table.
#'   No permutation importance, learning curve or missing history is invented.
#' @export
sb_wsrf_diagnostics <- function(fit, type = "importance") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!identical(s$learner, "wsrf_binary")) sb_abort(
    "This diagnostic requires WSRF.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("importance", "tree_nodes", "tree_oob")) sb_abort(
    "Choose importance, tree_nodes or tree_oob.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  n <- length(e$trees)
  m <- NULL
  if (type == "tree_nodes") {
    tab <- data.frame(tree = seq_len(n), nodes = lengths(e$trees))
  } else {
    m <- fit$diagnostic_artifacts$learner$diagnostic_metrics
    if (is.null(m)) sb_abort("Original WSRF metrics were not captured; compact placeholders are not diagnostic values.",
      "stormbreaker_missing_diagnostic_artifacts")
    if (type == "tree_oob") {
      if (!is.numeric(m$tree_oob) || length(m$tree_oob) != n || any(is.infinite(m$tree_oob)) ||
          any(m$tree_oob < 0 | m$tree_oob > 1, na.rm = TRUE)) sb_abort(
        "Original WSRF tree OOB values are invalid.", "stormbreaker_invalid_diagnostic")
      tab <- data.frame(tree = seq_len(n), error = m$tree_oob, available = is.finite(m$tree_oob))
    } else {
      z <- m$importance
      if (!is.matrix(z) || !"MeanDecreaseIGR" %in% colnames(z) ||
          !identical(rownames(z), s$internal_features) || any(!is.finite(z))) sb_abort(
        "Original WSRF information-gain-ratio importance is invalid.", "stormbreaker_invalid_diagnostic")
      bridge <- seq_len(nrow(z)) > length(s$features)
      tab <- data.frame(feature = c(s$features, if (any(bridge)) "Constant bridge"),
        internal_feature = s$internal_features, artificial_bridge = bridge,
        importance = unname(z[, "MeanDecreaseIGR"]))
    }
  }
  structure(list(type = paste0("wsrf_", type), table = tab, learner = s$learner,
    pipeline_id = fit$pipeline$id, forest_oob = m$forest_oob,
    source = "original WSRF node lists and separately captured pre-compaction metrics",
    interpretation = paste("Importance is native mean decrease in information gain ratio, not Gini or permutation importance.",
      "Node counts include internal and terminal nodes. OOB points are individual-tree errors, not a cumulative learning curve.",
      "Compact prediction placeholders are never reported as OOB results. No refitting or uncertainty reconstruction.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_wsrf_diagnostic <- function(x) {
  importance <- feature <- nodes <- trees <- tree <- error <- NULL
  d <- x$table
  if (x$type == "wsrf_importance") {
    d$feature <- factor(d$feature, levels = d$feature[order(d$importance)])
    p <- ggplot2::ggplot(d, ggplot2::aes(importance, feature)) +
      ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::labs(title = "Weighted-subspace forest importance", x = "Mean decrease in information gain ratio", y = NULL)
  } else if (x$type == "wsrf_tree_nodes") {
    counts <- as.data.frame(table(d$nodes))
    names(counts) <- c("nodes", "trees")
    counts$nodes <- as.numeric(as.character(counts$nodes))
    p <- ggplot2::ggplot(counts, ggplot2::aes(nodes, trees)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::labs(title = "Weighted-subspace forest sizes", x = "Total nodes per tree", y = "Number of trees")
  } else p <- ggplot2::ggplot(d[d$available, , drop = FALSE], ggplot2::aes(tree, error)) +
    ggplot2::geom_point(colour = "#2878A0", size = 1) + ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = "Individual-tree OOB errors", x = "Tree index", y = "OOB misclassification rate")
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}
