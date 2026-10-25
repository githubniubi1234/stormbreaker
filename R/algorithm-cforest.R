#' Inspect original conditional inference forest structure
#' @param fit A successful binary or survival cforest final fit.
#' @param type split_usage or tree_size.
#' @return All-feature split counts or all-tree node/leaf counts. Native
#'   partynode child structure and original model-frame column mapping are
#'   used; response columns are not mistaken for predictors. Counts describe
#'   structure, not conditional permutation importance, statistical significance
#'   or held-out performance. No new predictions, permutations or refits occur.
#' @export
sb_cforest_diagnostics <- function(fit, type = "split_usage") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% c("cforest_binary", "cforest_survival")) sb_abort(
    "This diagnostic requires cforest.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("split_usage", "tree_size")) sb_abort(
    "Choose split_usage or tree_size.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  if (!inherits(e, "cforest") || length(e$nodes) != s$params$num_trees) sb_abort(
    "Original conditional forest nodes are incomplete.", "stormbreaker_bad_learner_state")
  counts <- stats::setNames(integer(length(s$features)), s$features)
  visit <- function(node) {
    kids <- partykit::kids_node(node)
    if (!length(kids)) return(c(nodes = 1L, terminal_nodes = 1L))
    index <- partykit::varid_split(partykit::split_node(node))
    feature <- names(e$data)[index]
    if (length(feature) != 1L || is.na(feature) || !feature %in% s$features) sb_abort(
      "Original conditional forest split mapping is invalid.", "stormbreaker_invalid_diagnostic")
    counts[feature] <<- counts[feature] + 1L
    c(nodes = 1L, terminal_nodes = 0L) + Reduce(`+`, lapply(kids, visit))
  }
  sizes <- t(vapply(e$nodes, visit, c(nodes = 0L, terminal_nodes = 0L)))
  tab <- if (type == "split_usage") data.frame(feature = s$features, splits = unname(counts)) else
    data.frame(tree = seq_len(nrow(sizes)), nodes = sizes[, "nodes"], terminal_nodes = sizes[, "terminal_nodes"])
  structure(list(type = paste0("cforest_", type), table = tab, learner = s$learner,
    pipeline_id = fit$pipeline$id, source = "original saved partykit cforest nodes and model-frame mapping",
    interpretation = paste("Every original tree and feature is retained. Counts are structural use, not conditional permutation importance.",
      "No new significance tests, OOB evaluation, permutations, surrogate trees or predictions are computed.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_cforest_diagnostic <- function(x) {
  feature <- splits <- terminal_nodes <- trees <- NULL
  d <- x$table
  if (x$type == "cforest_split_usage") {
    d$feature <- factor(d$feature, levels = d$feature[order(d$splits)])
    p <- ggplot2::ggplot(d, ggplot2::aes(splits, feature)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::labs(title = "Conditional forest split usage", x = "Number of original splits", y = NULL)
  } else {
    counts <- as.data.frame(table(d$terminal_nodes))
    names(counts) <- c("terminal_nodes", "trees")
    counts$terminal_nodes <- as.numeric(as.character(counts$terminal_nodes))
    p <- ggplot2::ggplot(counts, ggplot2::aes(terminal_nodes, trees)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::scale_y_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::labs(title = "Conditional forest tree sizes", x = "Terminal nodes per tree", y = "Number of trees")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}
