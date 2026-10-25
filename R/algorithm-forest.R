#' Inspect native forest importance, tree sizes or saved OOB errors
#' @param fit A final ranger, randomForest/RRF binary or randomForestSRC survival
#'   model. Importance and tree sizes use saved prediction state. OOB error
#'   requires capture_diagnostics = TRUE for randomForest/randomForestSRC.
#' @param type One of importance, tree_size or oob_error.
#' @return An algorithm diagnostic with the full native table and source labels.
#'   Ranger importance is its saved impurity measure under the declared split
#'   rule; randomForest importance is mean decrease in Gini. Neither is
#'   permutation importance, a causal effect or directly comparable across
#'   engines. All features, including zero/negative values, remain in the table.
#'
#'   Tree sizes count original terminal nodes, without extracting a surrogate
#'   tree. OOB values are development-forest diagnostics, not validation of
#'   preprocessing/feature selection or untouched-cohort performance. Missing
#'   OOB values remain explicit. randomForestSRC may save only its final OOB
#'   error: this produces one point, not an invented learning curve. Ranger OOB
#'   computation was disabled in these pipelines and is not reconstructed.
#'   No fitting, permutation, parameter selection or prediction occurs.
#' @export
sb_forest_diagnostics <- function(fit, type = "importance") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  learner <- fit$pipeline$learner
  ranger <- learner %in% sb_ranger_learner_ids()
  rf <- learner %in% c("random_forest_binary", "rrf_binary")
  src <- identical(learner, "rsf_survival")
  if (!ranger && !rf && !src) sb_abort("This forest engine is not yet supported.",
    "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("importance", "tree_size", "oob_error")) sb_abort(
    "Choose importance, tree_size or oob_error.", "stormbreaker_unsupported_diagnostic")
  state <- fit$learner_state
  if (type == "oob_error") {
    if (ranger) sb_abort(
      "Ranger OOB computation was disabled; no original OOB error sequence exists.",
      "stormbreaker_unavailable_diagnostic")
    if (is.null(fit$diagnostic_artifacts)) sb_abort(
      "Original OOB results require capture_diagnostics when fitting; they are not in prediction state.",
      "stormbreaker_missing_diagnostic_artifacts")
    engine <- fit$diagnostic_artifacts$learner$engine
    if (rf) {
      errors <- engine$err.rate
      if (!is.matrix(errors) || !identical(colnames(errors), c("OOB", "0", "1"))) sb_abort(
        "Native OOB error columns are not the expected overall/0/1 mapping.",
        "stormbreaker_invalid_diagnostic")
      tab <- data.frame(trees = rep(seq_len(nrow(errors)), ncol(errors)),
        series = rep(c("Overall", "Class 0", "Class 1"), each = nrow(errors)),
        error = as.vector(errors))
      measure <- "OOB misclassification rate"
    } else {
      errors <- engine$err.rate
      if (!is.numeric(errors) || !is.null(dim(errors))) sb_abort(
        "Expected a single-outcome native survival OOB error sequence.",
        "stormbreaker_invalid_diagnostic")
      tab <- data.frame(trees = seq_along(errors), series = "Overall", error = errors)
      measure <- "Native OOB error (1 - concordance)"
    }
    if (!nrow(tab) || max(tab$trees) != engine$ntree ||
        any(is.infinite(tab$error)) || any(tab$error < 0 | tab$error > 1, na.rm = TRUE)) sb_abort(
      "Native forest OOB errors are incomplete or invalid.", "stormbreaker_invalid_diagnostic")
    tab$available <- is.finite(tab$error)
    source <- "original captured forest OOB results"
  } else if (type == "importance") {
    if (src) sb_abort(
      "This randomForestSRC fit disabled importance computation; none is invented.",
      "stormbreaker_unavailable_diagnostic")
    if (ranger) {
      tab <- state$feature_importance
      measure <- paste0("Native impurity importance (", state$params$splitrule, " splits)")
    } else {
      engine <- unserialize(state$engine_raw)
      importance <- engine$importance
      if (!is.matrix(importance) || !"MeanDecreaseGini" %in% colnames(importance)) sb_abort(
        "Native mean-decrease-Gini importance is missing.", "stormbreaker_invalid_diagnostic")
      tab <- data.frame(feature = rownames(importance), importance = importance[, "MeanDecreaseGini"])
      rownames(tab) <- NULL
      measure <- "Mean decrease in Gini"
    }
    if (!is.data.frame(tab) || !all(c("feature", "importance") %in% names(tab)) ||
        !is.numeric(tab$importance) || any(!is.finite(tab$importance)) ||
        anyDuplicated(tab$feature) || !setequal(tab$feature, state$features)) sb_abort(
      "Native forest importance is missing or invalid.", "stormbreaker_invalid_diagnostic")
    source <- "original forest importance retained in prediction state"
  } else {
    if (ranger) {
      forest <- unserialize(state$model_raw)
      leaves <- vapply(forest$child.nodeIDs, function(children) {
        if (length(children) != 2L || length(children[[1]]) != length(children[[2]])) sb_abort(
          "Native ranger child-node arrays are invalid.", "stormbreaker_invalid_diagnostic")
        sum(children[[1]] == 0 & children[[2]] == 0)
      }, integer(1))
      trees <- state$num_trees
    } else {
      engine <- unserialize(state$engine_raw)
      leaves <- if (rf) colSums(engine$forest$nodestatus == -1L) else engine$leaf.count
      trees <- engine$ntree
    }
    if (!is.numeric(leaves) || length(leaves) != trees || !length(leaves) ||
        any(!is.finite(leaves)) || any(leaves < 1) || any(leaves != floor(leaves))) sb_abort(
      "Native terminal-node counts are incomplete or invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(tree = seq_along(leaves), terminal_nodes = as.integer(leaves))
    measure <- "Terminal nodes per original tree"
    source <- "original forest structure retained in prediction state"
  }
  structure(list(type = paste0("forest_", type), table = tab, measure = measure,
    learner = learner, pipeline_id = fit$pipeline$id, source = source),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_forest_diagnostic <- function(x) {
  importance <- feature <- terminal_nodes <- trees <- error <- series <- NULL
  if (x$type == "forest_importance") {
    tab <- x$table
    tab$feature <- factor(tab$feature, levels = tab$feature[order(tab$importance)])
    p <- ggplot2::ggplot(tab, ggplot2::aes(importance, feature)) +
      ggplot2::geom_col(fill = "#2878A0", width = .65) +
      ggplot2::geom_vline(xintercept = 0, colour = "#909090") +
      ggplot2::labs(title = "Forest importance", subtitle = x$learner,
        x = x$measure, y = NULL,
        caption = "Native split-based importance, not permutation importance or a causal effect.\nCan favour variables with more possible splits; correlated features can share importance.") + sb_viz_theme()
    attr(p, "stormbreaker_suggested_size") <- c(width = 10, height = max(5, 2.5 + .28 * nrow(tab)))
    return(p)
  }
  if (x$type == "forest_tree_size") {
    counts <- as.data.frame(table(x$table$terminal_nodes), stringsAsFactors = FALSE)
    names(counts) <- c("terminal_nodes", "trees")
    counts$terminal_nodes <- as.integer(counts$terminal_nodes)
    return(ggplot2::ggplot(counts, ggplot2::aes(terminal_nodes, trees)) +
      ggplot2::geom_col(fill = "#2878A0", width = .8) +
      ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::scale_y_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::labs(title = "Forest tree-size distribution", subtitle = x$learner,
        x = x$measure, y = "Number of trees",
        caption = "Every fitted tree counted once. Original terminal-node counts; no pruning or surrogate.\nStructural complexity, not an estimate of predictive performance.") + sb_viz_theme())
  }
  tab <- x$table[x$table$available, , drop = FALSE]
  p <- ggplot2::ggplot(tab, ggplot2::aes(trees, error, colour = series)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_colour_manual(values = c(Overall = "#2878A0", "Class 0" = "#909090", "Class 1" = "#D77935")) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_x_continuous(limits = c(0, max(x$table$trees)),
      breaks = function(z) unique(round(pretty(z)))) +
    ggplot2::labs(title = "Forest OOB results", subtitle = x$learner,
      x = "Number of accumulated trees", y = x$measure, colour = "Outcome",
      caption = "Only native recorded values shown; missing tree counts remain in the diagnostic table.\nOOB does not validate prior preprocessing/selection and is not external-cohort performance.") + sb_viz_theme()
  if (nrow(tab) > 1L) {
    # Keep missing tree positions in the line layer, so gaps are not interpolated.
    p <- p + ggplot2::geom_line(data = x$table, na.rm = TRUE, linewidth = .4)
  }
  if (!nrow(tab)) {
    p <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
      label = "No finite native OOB error was recorded.") + ggplot2::theme_void() +
      ggplot2::labs(title = "Forest OOB results", subtitle = x$learner)
  }
  p
}
