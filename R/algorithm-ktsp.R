#' Inspect the original k-top-scoring-pairs decision rules
#' @param fit A successful k-TSP binary final fit.
#' @return An algorithm diagnostic preserving all selected ordered feature pairs
#'   and their original training scores. A first-greater-than-second comparison
#'   contributes one class-one vote, equality contributes one half and a smaller
#'   first value contributes zero. The package averages pair votes; this fraction
#'   is not independently calibrated probability. Native pair scores are training
#'   selection evidence, not validation accuracy, uncertainty or vote weights.
#'   Compact saved state suffices; no pairs are retrained or reordered.
#' @export
sb_ktsp_diagnostics <- function(fit) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (fit$pipeline$learner != "ktsp_binary") sb_abort(
    "This diagnostic requires the k-TSP binary learner.", "stormbreaker_unsupported_diagnostic")
  e <- fit$learner_state$engine
  if (!is.matrix(e$pairs) || ncol(e$pairs) != 2L || !nrow(e$pairs) ||
      !all(e$pairs %in% fit$learner_state$features) ||
      length(e$pair_scores) != nrow(e$pairs) || any(!is.finite(e$pair_scores)) ||
      !identical(e$labels, c("0", "1"))) sb_abort(
    "Original k-TSP pair mapping or scores are invalid.", "stormbreaker_bad_learner_state")
  tab <- data.frame(pair = seq_len(nrow(e$pairs)), first = e$pairs[, 1], second = e$pairs[, 2],
    native_score = e$pair_scores, class_1_rule = paste(e$pairs[, 1], ">", e$pairs[, 2]))
  structure(list(type = "ktsp_pair_scores", table = tab,
    vote_rule = c(greater = 1, equal = .5, smaller = 0), vote_aggregation = "unweighted_mean",
    unused_features = setdiff(fit$learner_state$features, as.vector(e$pairs)),
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    source = "original saved switchBox ordered pairs and pair selection scores",
    interpretation = paste("Each displayed comparison contributes one class-one vote when true, one half for equality, zero otherwise.",
      "Original pair order retained. Final output is the unweighted mean of votes; native scores are not voting weights.",
      "Scores describe training pair selection, not validation accuracy or uncertainty; vote fractions are not independently calibrated probabilities.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_ktsp_diagnostic <- function(x) {
  rule <- native_score <- NULL
  d <- x$table
  labels <- paste0("Pair ", d$pair, ": ", d$class_1_rule)
  d$rule <- factor(labels, levels = rev(labels))
  p <- ggplot2::ggplot(d, ggplot2::aes(native_score, rule)) +
    ggplot2::geom_point(size = 2.5, colour = "#2878A0") +
    ggplot2::labs(title = "k-TSP pair rules", x = "Native pair selection score", y = NULL) +
    sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = max(4.5, 2 + .4 * nrow(d)))
  p
}
