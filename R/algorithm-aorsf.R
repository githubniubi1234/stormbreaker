#' Inspect saved accelerated oblique forest diagnostics
#' @param fit A successful aorsf binary or survival final fit.
#' @param type tree_size or oob_score.
#' @return Original terminal-node counts, or saved OOB evaluation values on
#'   the native evaluation schedule. Binary scores are AUC-ROC and survival
#'   scores are Harrell's C-index, not error rates. A single final evaluation
#'   stays a single point. Importance was disabled and is not reconstructed.
#'   OOB evaluation does not validate prior preprocessing or feature selection.
#'   Serialized R6 state is read into a separate object, without updating the fit.
#' @export
sb_aorsf_diagnostics <- function(fit, type = "tree_size") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% c("aorsf_binary", "aorsf_survival")) sb_abort(
    "This diagnostic requires an aorsf model.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("tree_size", "oob_score")) sb_abort(
    "Choose tree_size or oob_score.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  if (type == "tree_size") {
    leaves <- vapply(e$forest$leaf_summary, function(z) sum(z != 0), integer(1))
    if (length(leaves) != e$n_tree || anyNA(leaves) || any(leaves < 1)) sb_abort(
      "Original oblique forest leaf summaries are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(tree = seq_along(leaves), terminal_nodes = leaves)
    measure <- "Terminal nodes per tree"
  } else {
    measure <- e$eval_oobag$stat_type
    expected <- if (s$learner == "aorsf_binary") "AUC-ROC" else "Harrell's C-index"
    z <- e$eval_oobag$stat_values
    steps <- unique(c(seq.int(e$oobag_eval_every, e$n_tree, by = e$oobag_eval_every), e$n_tree))
    if (!identical(measure, expected) || !is.matrix(z) || ncol(z) != 1L ||
        nrow(z) != length(steps) || any(is.infinite(z)) || any(z < 0 | z > 1, na.rm = TRUE)) sb_abort(
      "Saved native OOB score or evaluation schedule is invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(trees = steps, score = as.numeric(z), available = is.finite(as.numeric(z)))
  }
  structure(list(type = paste0("aorsf_", type), table = tab, measure = measure,
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "saved native aorsf leaf summaries and OOB evaluations",
    interpretation = paste("Counts use original native leaf summaries. OOB scores retain their original metric and evaluation schedule.",
      "AUC and concordance are scores, not error rates. A final-only evaluation remains one point, not an invented learning curve.",
      "OOB is not validation of prior preprocessing/selection or external performance; no importance, intervals or new predictions are computed.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_aorsf_diagnostic <- function(x) {
  terminal_nodes <- trees <- score <- NULL
  d <- x$table
  if (x$type == "aorsf_tree_size") {
    counts <- as.data.frame(table(d$terminal_nodes))
    names(counts) <- c("terminal_nodes", "trees")
    counts$terminal_nodes <- as.numeric(as.character(counts$terminal_nodes))
    p <- ggplot2::ggplot(counts, ggplot2::aes(terminal_nodes, trees)) + ggplot2::geom_col(fill = "#2878A0") +
      ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::scale_y_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::labs(title = "Oblique forest tree sizes", x = x$measure, y = "Number of trees")
  } else p <- ggplot2::ggplot(d[d$available, , drop = FALSE], ggplot2::aes(trees, score)) +
    ggplot2::geom_point(colour = "#2878A0", size = 3) + ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_x_continuous(limits = c(0, max(d$trees))) +
    ggplot2::labs(title = "Recorded oblique forest OOB score", x = "Number of trees", y = x$measure)
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}
