#' Inspect a fitted single rpart tree without replacing it
#' @param fit A final rpart classification or survival tree, fitted with
#'   capture_diagnostics = TRUE. Forests and boosted ensembles are not trees.
#' @return An algorithm diagnostic with every node, parent, native split rule,
#'   node count, node prediction, exact native split matrix and complexity table.
#'   Binary node probabilities are training-node proportions, not validation
#'   estimates. Survival node values are relative rates, not survival probabilities.
#' @export
sb_tree_structure <- function(fit) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!fit$pipeline$learner %in% c("decision_tree_binary", "rpart_survival")) sb_abort(
    "This adapter requires a single rpart tree, not an ensemble or another tree implementation.",
    "stormbreaker_unsupported_diagnostic")
  if (is.null(fit$diagnostic_artifacts)) sb_abort(
    "The original tree was not captured; enable capture_diagnostics when fitting.",
    "stormbreaker_missing_diagnostic_artifacts")
  engine <- fit$diagnostic_artifacts$learner$engine
  if (!inherits(engine, "rpart")) sb_abort("The captured rpart engine is missing.",
    "stormbreaker_invalid_diagnostic")
  frame <- engine$frame
  nodes <- as.numeric(rownames(frame))
  if (!length(nodes) || any(!is.finite(nodes)) || anyDuplicated(nodes) ||
      nodes[1] != 1 || any(nodes < 1) || any(nodes != floor(nodes))) sb_abort(
    "Native tree node identifiers are invalid.", "stormbreaker_invalid_diagnostic")
  leaf <- frame$var == "<leaf>"
  if (any(!leaf & (!nodes * 2 %in% nodes | !(nodes * 2 + 1) %in% nodes))) sb_abort(
    "Native tree contains an incomplete branch.", "stormbreaker_invalid_diagnostic")
  tab <- data.frame(node = nodes, parent = ifelse(nodes == 1, NA_real_, floor(nodes / 2)),
    depth = floor(log2(nodes)), terminal = leaf, n = frame$n,
    split_feature = as.character(frame$var), incoming_rule = labels(engine, digits = 7, minlength = 0),
    native_value = frame$yval)
  if (identical(engine$method, "class")) {
    levels <- attr(engine, "ylevels")
    if (!identical(as.character(levels), c("0", "1")) || ncol(frame$yval2) != 6L) sb_abort(
      "The tree's binary class mapping is not the expected 0/1 mapping.", "stormbreaker_invalid_diagnostic")
    tab$prediction <- frame$yval2[, 5L]
    tab$class_0_n <- frame$yval2[, 2L]
    tab$class_1_n <- frame$yval2[, 3L]
    scale <- "Training-node P(class 1)"
  } else {
    tab$prediction <- frame$yval
    scale <- "Native relative event rate"
  }
  xpos <- rep(NA_real_, length(nodes))
  xpos[leaf] <- seq_len(sum(leaf))
  position <- function(node) {
    i <- match(node, nodes)
    if (is.finite(xpos[i])) return(xpos[i])
    xpos[i] <<- mean(c(position(2 * node), position(2 * node + 1)))
    xpos[i]
  }
  position(1)
  tab$x <- xpos
  tab$y <- -tab$depth
  structure(list(type = "tree_structure", table = tab,
    native_splits = engine$splits, native_categorical_splits = engine$csplit,
    complexity = engine$cptable, prediction_scale = scale,
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    source = "all nodes of original rpart fit; no pruning or refitting",
    suggested_width = max(9, 1.6 * sum(leaf)),
    suggested_height = max(5, 1.2 * (max(tab$depth) + 1))),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_rpart_tree <- function(x) {
  d <- x$table
  parent <- match(d$parent, d$node)
  edges <- data.frame(x = d$x[parent], y = d$y[parent], xend = d$x, yend = d$y)
  edges <- edges[!is.na(parent), , drop = FALSE]
  edges$y <- edges$y - .33
  edges$yend <- edges$yend + .33
  d$label <- vapply(seq_len(nrow(d)), function(i) {
    if ("display_label" %in% names(d)) return(d$display_label[i])
    rule <- paste(strwrap(d$incoming_rule[i], width = 24), collapse = "\n")
    paste0("Node ", d$node[i], " | n = ", d$n[i], "\n", rule, "\n",
      if (x$learner == "decision_tree_binary") "P(1) = " else "Relative rate = ",
      format(signif(d$prediction[i], 3), trim = TRUE))
  }, "")
  node_x <- node_y <- y <- xend <- yend <- label <- terminal <- NULL
  names(d)[names(d) == "x"] <- "node_x"
  names(d)[names(d) == "y"] <- "node_y"
  p <- ggplot2::ggplot(d, ggplot2::aes(node_x, node_y)) +
    ggplot2::geom_segment(data = edges, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE, colour = "#909090") +
    ggplot2::geom_label(ggplot2::aes(label = label, fill = terminal), size = 3,
      linewidth = .2, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c("FALSE" = "#FFFFFF", "TRUE" = "#DDEBF3")) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = .75)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = .7)) +
    ggplot2::labs(title = "Single-tree decision structure",
      subtitle = x$prediction_scale,
      caption = "Every original node shown; a node's rule describes its incoming branch. Terminal nodes are shaded.\nDisplay thresholds are rounded; exact native splits remain in the diagnostic. No pruning or refitting.") +
    sb_viz_theme() + ggplot2::theme(axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank())
  attr(p, "stormbreaker_suggested_size") <- c(width = x$suggested_width, height = x$suggested_height)
  p
}
