#' Inspect original bagged classification trees
#' @param fit A successful ipred bagging final fit.
#' @param type split_usage, tree_size or importance.
#' @return Original primary split counts, node/leaf counts, or complete
#'   tree-by-feature native rpart importance. Absent feature importance entries
#'   in a tree are zero, including stumps. The importance plot averages original
#'   rpart scores across all bags; these may credit surrogate splits and are not
#'   permutation importance or validation performance. OOB was disabled.
#' @export
sb_bagging_diagnostics <- function(fit, type = "split_usage") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!identical(s$learner, "bagging_binary")) sb_abort(
    "This diagnostic requires ipred bagging.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("split_usage", "tree_size", "importance")) sb_abort(
    "Choose split_usage, tree_size or importance.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  if (!inherits(e, "classbagg") || length(e$mtrees) != s$params$num_trees) sb_abort(
    "Saved bagged trees are incomplete.", "stormbreaker_bad_learner_state")
  counts <- stats::setNames(integer(length(s$features)), s$features)
  rows <- lapply(seq_along(e$mtrees), function(i) {
    tree <- e$mtrees[[i]]$btree
    frame <- tree$frame
    if (!inherits(tree, "rpart") || !is.data.frame(frame) || !nrow(frame)) sb_abort(
      "An original bagged rpart tree is invalid.", "stormbreaker_invalid_diagnostic")
    vars <- frame$var[frame$var != "<leaf>"]
    if (any(!vars %in% s$features)) sb_abort("A bagged split has an unknown feature.", "stormbreaker_invalid_diagnostic")
    counts <<- counts + as.integer(table(factor(vars, levels = s$features)))
    values <- stats::setNames(numeric(length(s$features)), s$features)
    native <- tree$variable.importance
    if (length(native)) {
      if (any(!names(native) %in% s$features) || any(!is.finite(native))) sb_abort(
        "Native bagged-tree importance is invalid.", "stormbreaker_invalid_diagnostic")
      values[names(native)] <- native
    }
    list(size = data.frame(tree = i, nodes = nrow(frame), terminal_nodes = sum(frame$var == "<leaf>")),
      importance = data.frame(tree = i, feature = s$features, importance = unname(values)))
  })
  tab <- if (type == "split_usage") data.frame(feature = s$features, splits = unname(counts)) else
    do.call(rbind, lapply(rows, `[[`, if (type == "tree_size") "size" else "importance"))
  structure(list(type = paste0("bagging_", type), table = tab, learner = s$learner,
    pipeline_id = fit$pipeline$id, source = "original ipred bags and native rpart frames/importance",
    interpretation = paste("All original bags and features retained. Split usage counts primary splits only.",
      "Native rpart importance may credit surrogate splits; its across-bag mean is not permutation importance or a causal effect.",
      "OOB was disabled. No OOB history, new predictions, confidence intervals or surrogate model is reconstructed.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_bagging_diagnostic <- function(x) {
  feature <- value <- terminal_nodes <- trees <- NULL
  d <- x$table
  if (x$type != "bagging_tree_size") {
    if (x$type == "bagging_importance") {
      d <- stats::aggregate(importance ~ feature, data = d, FUN = mean)
      names(d)[2] <- "value"
    } else names(d)[2] <- "value"
    d$feature <- factor(d$feature, levels = d$feature[order(d$value)])
    p <- ggplot2::ggplot(d, ggplot2::aes(value, feature)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::labs(title = if (x$type == "bagging_importance") "Bagged-tree native importance" else "Bagged-tree split usage",
        x = if (x$type == "bagging_importance") "Mean native rpart importance" else "Number of primary splits", y = NULL)
  } else {
    counts <- as.data.frame(table(d$terminal_nodes))
    names(counts) <- c("terminal_nodes", "trees")
    counts$terminal_nodes <- as.numeric(as.character(counts$terminal_nodes))
    p <- ggplot2::ggplot(counts, ggplot2::aes(terminal_nodes, trees)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::scale_y_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::labs(title = "Bagged-tree sizes", x = "Terminal nodes per tree", y = "Number of trees")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}
