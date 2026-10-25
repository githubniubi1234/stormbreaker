sb_multiclass_forest_structure_plots <- function(fit) {
  m <- fit$engine
  counts <- stats::setNames(integer(length(fit$features)), fit$feature_codes)
  if (fit$method == "cforest") {
    visit <- function(node) {
      kids <- partykit::kids_node(node)
      if (!length(kids)) return(1L)
      name <- names(m$data)[partykit::varid_split(partykit::split_node(node))]
      counts[name] <<- counts[name] + 1L
      sum(vapply(kids, visit, 1L))
    }
    sizes <- vapply(m$nodes, visit, 1L)
  } else if (fit$method == "bagging") {
    sizes <- vapply(m$mtrees, function(tree) {
      vars <- as.character(tree$btree$frame$var)
      counts <<- counts + tabulate(match(vars[vars != "<leaf>"], fit$feature_codes), length(counts))
      sum(vars == "<leaf>")
    }, 1L)
  } else if (fit$method == "rlt") {
    sizes <- vapply(m$FittedTrees, function(tree) {
      internal <- tree["NodeType", ] == 1
      counts <<- counts + tabulate(as.integer(tree["SplitVar1", internal]), length(counts))
      sum(!internal)
    }, 1L)
  } else if (fit$method == "aorsf") {
    sizes <- vapply(m$forest$leaf_summary, function(x) sum(x != 0), 1L)
  } else sizes <- lengths(m$trees)
  d <- data.frame(tree = seq_along(sizes), nodes = as.integer(sizes))
  frequencies <- as.data.frame(table(d$nodes))
  names(frequencies) <- c("nodes", "trees")
  frequencies$nodes <- as.integer(as.character(frequencies$nodes))
  nodes <- trees <- NULL
  integer_breaks <- function(limits) {
    ticks <- pretty(limits, n = 5)
    ticks[ticks == floor(ticks) & ticks >= 0]
  }
  p <- ggplot2::ggplot(frequencies, ggplot2::aes(nodes, trees)) + ggplot2::geom_col(fill = "#2878A0") +
    ggplot2::scale_x_continuous(breaks = integer_breaks) +
    ggplot2::scale_y_continuous(breaks = integer_breaks) +
    sb_viz_theme() + ggplot2::labs(title = "Tree-size distribution",
      x = if (fit$method == "wsrf") "Total nodes per tree" else "Terminal nodes per tree", y = "Trees",
      caption = "Every original fitted tree is counted once; node counts describe structural complexity, not predictive performance. No pruning or surrogate tree is used.")
  result <- list(plots = list(tree_size = p), tables = list(tree_size = d))
  if (fit$method %in% c("cforest", "bagging", "rlt")) {
    d <- data.frame(feature = fit$features, splits = unname(counts))
    d$feature <- factor(d$feature, levels = d$feature[order(d$splits)])
    feature <- splits <- NULL
    result$plots$split_usage <- ggplot2::ggplot(d, ggplot2::aes(splits, feature)) +
      ggplot2::geom_col(fill = "#2878A0") + sb_viz_theme() +
      ggplot2::scale_x_continuous(breaks = integer_breaks) +
      ggplot2::labs(title = "Split usage", x = "Splits", y = NULL,
        caption = "Counts of original primary axis-aligned splits, including zero-use features; not impurity, permutation importance or causal effects.")
    result$tables$split_usage <- d
  }
  if (fit$method == "wsrf") {
    d <- data.frame(tree = seq_along(m$treeOOBErrorRates), error = as.numeric(m$treeOOBErrorRates))
    tree <- error <- NULL
    result$plots$individual_tree_oob <- ggplot2::ggplot(d, ggplot2::aes(tree, error)) +
      ggplot2::geom_point(colour = "#2878A0", size = 1) + sb_viz_theme() +
      ggplot2::labs(title = "Individual-tree OOB error", x = "Tree", y = "Error rate",
        caption = "Saved individual-tree out-of-bag errors, not a cumulative ensemble learning curve or held-out validation.")
    result$tables$individual_tree_oob <- d
  }
  result
}
