sb_roc_interval <- function(y, p, reps = 0L) {
  auc <- function(y, p) {
    n1 <- sum(y == 1)
    n0 <- sum(y == 0)
    if (!n1 || !n0) return(NA_real_)
    (sum(rank(p)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  values <- if (reps > 0) replicate(reps, {
    i <- sample.int(length(y), replace = TRUE)
    auc(y[i], p[i])
  }) else numeric()
  ci <- if (sum(is.finite(values)) >= 2) stats::quantile(values, c(.025, .975), na.rm = TRUE) else c(NA_real_, NA_real_)
  data.frame(auc = auc(y, p), lower = unname(ci[1]), upper = unname(ci[2]),
    bootstrap_success = sum(is.finite(values)), bootstrap_requested = reps)
}

#' Plot one binary ROC per supplied score column
#' @param outcome Numeric 0/1 outcomes, one independent patient per row.
#' @param probability Numeric score vector or matrix with named model columns.
#'   Scores need not be probabilities; higher values must indicate class 1.
#' @param cohort Short cohort label.
#' @param kind Explicit prediction origin: apparent training, external, or OOF.
#' @param ci_reps Patient-bootstrap replicates; 0 omits intervals.
#' @param seed Random seed.
#' @return A ggplot with AUC/95 percent percentile intervals in its legend and
#'   a stormbreaker_roc_statistics attribute. Intervals condition on supplied
#'   scores, do not include model development uncertainty, and use ordinary
#'   patient resampling. Single-class resamples are counted as unsuccessful.
#' @export
sb_plot_prediction_roc <- function(outcome, probability, kind, cohort = "Cohort",
                                    ci_reps = 1000L, seed = 1L) {
  kind <- match.arg(kind, c("apparent", "external", "oof"))
  sb_require_curve_plotting()
  if (is.numeric(probability) && is.null(dim(probability))) probability <- matrix(probability, ncol = 1, dimnames = list(NULL, "Model"))
  if (!is.matrix(probability) || !is.numeric(probability) || !ncol(probability) ||
      nrow(probability) != length(outcome) || any(!is.finite(probability)) ||
      !is.numeric(outcome) || anyNA(outcome) || !all(outcome %in% c(0, 1)) ||
      length(unique(outcome)) != 2L) stop("Supply aligned finite scores and both 0/1 classes; no rows are dropped.")
  if (is.null(colnames(probability)) || anyDuplicated(colnames(probability)) || any(!nzchar(colnames(probability)))) stop("Score columns need unique model names.")
  if (length(ci_reps) != 1L || !is.finite(ci_reps) || ci_reps < 0 || ci_reps != as.integer(ci_reps) || ci_reps == 1) stop("ci_reps must be 0 or an integer of at least 2.")
  withr::local_seed(seed)
  tables <- lapply(seq_len(ncol(probability)), function(i) {
    z <- sb_binary_threshold_coordinates(outcome, probability[, i], "roc")
    z$model <- colnames(probability)[i]
    z
  })
  a <- do.call(rbind, lapply(seq_len(ncol(probability)), function(i) sb_roc_interval(outcome, probability[, i], ci_reps)))
  a$model <- colnames(probability)
  labels <- paste0(a$model, ": AUC ", sprintf("%.3f", a$auc),
    if (ci_reps > 0) paste0(" (95% CI ", sprintf("%.3f", a$lower), "\u2013", sprintf("%.3f", a$upper), ")") else "")
  x <- y <- model <- NULL
  p <- ggplot2::ggplot(do.call(rbind, tables), ggplot2::aes(x, y, colour = model)) +
    ggplot2::geom_path(linewidth = .7) + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey65") +
    ggplot2::scale_colour_discrete(breaks = a$model, labels = labels) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(title = cohort, x = "1 \u2212 Specificity", y = "Sensitivity", colour = NULL) +
    sb_viz_theme() + ggplot2::theme(legend.position = "bottom")
  attr(p, "stormbreaker_roc_statistics") <- a
  attr(p, "prediction_origin") <- kind
  p
}
