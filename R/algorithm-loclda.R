#' Inspect the original local-LDA prediction neighborhood
#' @param fit Successful local-LDA final fit.
#' @param newdata Explicit query rows; no implicit training queries.
#' @return All neighbors within the native kth-distance radius, including ties,
#'   their native weights, classes and whether the native class-count rule uses
#'   them. Native posterior probabilities are also returned. Reference indices
#'   and classes may be sensitive. No tuning or separate classifier is fitted.
#' @export
sb_loclda_neighborhood <- function(fit, newdata = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (s$learner != "loclda_binary") sb_abort(
    "This diagnostic requires local LDA.", "stormbreaker_unsupported_diagnostic")
  if (is.null(newdata)) sb_abort("Supply explicit query rows for local-LDA neighborhoods.", "stormbreaker_unavailable_diagnostic")
  if (!requireNamespace("klaR", quietly = TRUE)) sb_abort("Local-LDA diagnostics require klaR.", "stormbreaker_unavailable_diagnostic")
  e <- unserialize(s$engine_raw)
  x <- as.matrix(sb_select_final_feature_frame(newdata, s$features, "Query data", "query fitted features"))
  if (!nrow(x) || !inherits(e, "loclda") || !identical(colnames(e$learn), s$features)) sb_abort(
    "Original local-LDA data or query is invalid.", "stormbreaker_invalid_diagnostic")
  tab <- sb_loclda_neighbor_table(e, x)
  probability <- sb_predict_loclda_binary(e, x)
  structure(list(type = "loclda_neighborhood", table = tab,
    predictions = data.frame(query = seq_len(nrow(x)), probability = as.numeric(probability)),
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "Original learn/grouping/k/weight.func and native predict.loclda",
    interpretation = "Distance-radius neighbors include ties. Classes with at most one neighbor are excluded by the native prediction rule. Weights are local prediction weights, not feature importance. Reference indices and labels may be sensitive. No tuning or extra model fit."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_loclda_neighbor_table <- function(e, x) {
  do.call(rbind, lapply(seq_len(nrow(x)), function(i) {
    distances <- sqrt(rowSums(sweep(e$learn, 2L, x[i, ], "-")^2))
    radius <- sort(distances)[e$k]
    if (!is.finite(radius) || radius <= 0) sb_abort(
      "Native local-LDA neighborhood radius is zero or invalid.", "stormbreaker_unavailable_diagnostic")
    selected <- which(distances <= radius)
    classes <- as.character(e$grouping[selected])
    counts <- table(factor(classes, levels = e$lev))
    weights <- vapply(distances[selected] / radius, e$weight.func, numeric(1))
    used <- counts[classes] > 1L
    data.frame(query = i, reference_row = selected, class = classes,
      distance = distances[selected], radius = radius, weight = weights,
      used = as.logical(used), effective_weight = weights * as.logical(used))
  }))
}

sb_plot_loclda_neighborhood <- function(x) {
  distance <- weight <- class <- used <- NULL
  p <- ggplot2::ggplot(x$table, ggplot2::aes(distance, weight, colour = class, shape = used)) +
    ggplot2::geom_point(size = 2, alpha = .8) +
    ggplot2::facet_wrap(~query, labeller = ggplot2::label_both) +
    ggplot2::scale_colour_manual(values = c("0" = "#2878A0", "1" = "#B86236")) +
    ggplot2::scale_shape_manual(values = c("FALSE" = 4, "TRUE" = 16), labels = c("FALSE" = "Excluded", "TRUE" = "Used")) +
    ggplot2::labs(title = "Local-LDA neighborhood", x = "Distance", y = "Local weight", colour = "Class", shape = NULL) + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}
