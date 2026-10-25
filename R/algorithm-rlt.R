#' Inspect original reinforcement-learning forest structure
#' @param fit A successful binary or survival RLT final fit.
#' @param type split_usage or tree_size.
#' @return Original axis-aligned split counts for every feature, or node and
#'   terminal-node counts for every tree. Split counts are structural usage,
#'   not permutation importance or predictive effects. The original pipeline
#'   disabled importance computation; missing OOB or reinforcement histories
#'   are not inferred. No model is retrained or replaced by a surrogate tree.
#' @export
sb_rlt_diagnostics <- function(fit, type = "split_usage") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% sb_rlt_learner_ids()) sb_abort(
    "This diagnostic requires an RLT model.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("split_usage", "tree_size")) sb_abort(
    "Choose split_usage or tree_size.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  sb_validate_rlt_engine(e, s$learner, s$features, s$params)
  rows <- lapply(seq_along(e$FittedTrees), function(i) {
    z <- e$FittedTrees[[i]]
    if (!is.matrix(z) || !all(c("NodeType", "Node", "NextLeft", "NextRight", "NumOfComb", "SplitVar1") %in% rownames(z)) ||
        any(!z["NodeType", ] %in% c(1, 2))) sb_abort(
      "Saved native RLT tree structure is invalid.", "stormbreaker_invalid_diagnostic")
    internal <- z["NodeType", ] == 1
    v <- z["SplitVar1", internal]
    if (any(z["NumOfComb", internal] != 1) || any(!is.finite(v)) ||
        any(v < 1 | v > length(s$features) | v != floor(v))) sb_abort(
      "Expected original axis-aligned RLT split indices.", "stormbreaker_invalid_diagnostic")
    list(counts = tabulate(as.integer(v), nbins = length(s$features)),
      size = data.frame(tree = i, nodes = ncol(z), terminal_nodes = sum(!internal)))
  })
  tab <- if (type == "split_usage") data.frame(feature = s$features,
    splits = Reduce(`+`, lapply(rows, `[[`, "counts"))) else do.call(rbind, lapply(rows, `[[`, "size"))
  structure(list(type = paste0("rlt_", type), table = tab, learner = s$learner,
    pipeline_id = fit$pipeline$id, source = "original saved RLT node arrays",
    interpretation = paste("Every original tree and feature retained. Counts describe structural use, not variable importance or causal effects.",
      "Original importance computation was disabled; no OOB or reinforcement-training trajectory is reconstructed.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_rlt_diagnostic <- function(x) {
  feature <- splits <- terminal_nodes <- trees <- NULL
  d <- x$table
  if (x$type == "rlt_split_usage") {
    d$feature <- factor(d$feature, levels = d$feature[order(d$splits)])
    p <- ggplot2::ggplot(d, ggplot2::aes(splits, feature)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::labs(title = "RLT split-variable usage", x = "Number of original splits", y = NULL)
  } else {
    counts <- as.data.frame(table(d$terminal_nodes))
    names(counts) <- c("terminal_nodes", "trees")
    counts$terminal_nodes <- as.numeric(as.character(counts$terminal_nodes))
    p <- ggplot2::ggplot(counts, ggplot2::aes(terminal_nodes, trees)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::scale_y_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::labs(title = "RLT tree-size distribution", x = "Terminal nodes per tree", y = "Number of trees")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}
