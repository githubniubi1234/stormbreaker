#' Plot a binary confusion matrix directly from labels and probabilities
#' @param outcome Numeric vector of observed 0/1 labels in patient row order.
#' @param probability Numeric vector of aligned class-1 probabilities in `[0,1]`.
#'   Rows are matched by position, not vector names. No rows are dropped.
#' @param kind Prediction origin: "apparent" training, "external", or "oof".
#'   This records the caller's declaration; it does not verify independence.
#' @param threshold Fixed decision threshold, default 0.5. Probability greater
#'   than or equal to threshold predicts class 1. No threshold optimization.
#' @param cohort Optional short cohort title.
#' @param class_labels Two distinct display labels for classes 0 and 1.
#' @param normalize "none" for counts, "observed" for proportions within each
#'   observed class. A missing observed class has undefined proportions.
#' @return A ggplot, with observed classes on x and predicted classes on y.
#'   The stormbreaker_prediction_context attribute records kind, threshold,
#'   sample count and TP/FP/TN/FN. Plot data contain counts, not patient records.
#' @details Use a prespecified or development-derived threshold unchanged on
#'   validation cohorts. Apparent training predictions are not cross-validation
#'   predictions. This function does not fit models or enter automatic reports.
#' @export
sb_plot_prediction_confusion <- function(outcome, probability, kind,
                                          threshold = .5, cohort = NULL,
                                          class_labels = c("0", "1"),
                                          normalize = c("none", "observed")) {
  kind <- match.arg(kind, c("apparent", "external", "oof"))
  normalize <- match.arg(normalize)
  if (!is.numeric(outcome) || !is.null(dim(outcome)) || !length(outcome) ||
      anyNA(outcome) || !all(outcome %in% c(0, 1)) ||
      !is.numeric(probability) || !is.null(dim(probability)) ||
      length(probability) != length(outcome) || any(!is.finite(probability)) ||
      any(probability < 0 | probability > 1)) stop("Supply aligned 0/1 labels and finite probabilities in [0,1].")
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) || threshold < 0 || threshold > 1) stop("threshold must be in [0,1].")
  if (!is.character(class_labels) || length(class_labels) != 2L || anyNA(class_labels) ||
      any(!nzchar(trimws(class_labels))) || anyDuplicated(class_labels)) stop("Supply two distinct class labels for 0 and 1.")
  if (!is.null(cohort) && (!is.character(cohort) || length(cohort) != 1L || is.na(cohort) || !nzchar(trimws(cohort)))) stop("cohort must be a short nonempty label or NULL.")
  positive <- probability >= threshold
  counts <- data.frame(threshold = threshold,
    tp = sum(positive & outcome == 1), fp = sum(positive & outcome == 0),
    tn = sum(!positive & outcome == 0), fn = sum(!positive & outcome == 1))
  table <- cbind(counts, view_id = "predictions", pipeline_id = "Model",
    fold = NA_integer_, panel = "Confusion matrix")
  metrics <- structure(list(table = table, diagnostics = data.frame(), saved_metrics = list()),
    class = "stormbreaker_threshold_metrics")
  p <- sb_plot_confusion(metrics, threshold, normalize = normalize)
  observed <- predicted <- value <- NULL
  p$mapping <- ggplot2::aes(observed, predicted, fill = value)
  labels <- stats::setNames(class_labels, c("0", "1"))
  p <- p + ggplot2::scale_x_discrete(limits = c("0", "1"), labels = labels, drop = FALSE) +
    ggplot2::scale_y_discrete(limits = c("1", "0"), labels = labels, drop = FALSE) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = if (is.null(cohort)) "Confusion matrix" else cohort,
      subtitle = NULL, caption = NULL, x = "Observed", y = "Predicted") +
    ggplot2::theme(strip.text = ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
      legend.position = "none")
  attr(p, "stormbreaker_prediction_context") <- list(kind = kind, cohort = cohort,
    threshold = threshold, n = length(outcome), class_labels = class_labels, counts = counts)
  attr(p, "stormbreaker_suggested_size") <- c(width = 4, height = 4)
  p
}
