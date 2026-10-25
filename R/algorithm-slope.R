#' Inspect the original SLOPE cross-validation search
#' @param fit A successful SLOPE final fit captured with capture_diagnostics.
#' @return An algorithm diagnostic retaining the full native q/gamma/alpha
#'   summary, means, standard errors, bounds and selected optimum. Paths may
#'   stop early; no uncomputed alpha values are added. Bounds are the native
#'   mean plus/minus qt(0.975, score-row-count - 1) times native SE, not external
#'   validation uncertainty. Training observations are not retained in this
#'   added diagnostic summary. Historical compact models have only their saved
#'   optimum and cannot supply a complete CV path. Original coefficients remain
#'   available through [sb_model_coefficients()], including without capture.
#' @export
sb_slope_diagnostics <- function(fit) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (fit$pipeline$learner != "slope_binary") sb_abort(
    "This diagnostic requires SLOPE.", "stormbreaker_unsupported_diagnostic")
  path <- fit$diagnostic_artifacts$learner$diagnostic_path
  if (is.null(path)) sb_abort("The original SLOPE CV summary was not captured.",
    "stormbreaker_missing_diagnostic_artifacts")
  tab <- path$summary
  columns <- c("q", "gamma", "alpha", "measure", "mean", "se", "lo", "hi")
  if (!is.data.frame(tab) || !all(columns %in% names(tab)) || !nrow(tab) ||
      any(!is.finite(as.matrix(tab[setdiff(columns, "measure")]))) || any(tab$alpha <= 0) ||
      any(tab$se < 0) || any(tab$lo > tab$mean | tab$hi < tab$mean) ||
      !identical(unique(tab$measure), "deviance")) sb_abort(
    "The original SLOPE deviance path is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
  optimum <- path$optima[1L, , drop = FALSE]
  tab$selected <- tab$q == optimum$q & tab$gamma == optimum$gamma & tab$alpha == optimum$alpha
  if (sum(tab$selected) != 1L || !isTRUE(all.equal(tab$mean[tab$selected], fit$learner_state$tuning$mean))) sb_abort(
    "The native SLOPE optimum does not match the saved prediction model.", "stormbreaker_invalid_diagnostic")
  structure(list(type = "slope_cv_curve", table = tab, native_table = path$summary,
    selected = optimum, score_rows = path$score_rows,
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    source = "original final-fit SLOPE cross-validation summary",
    interpretation = paste("Original mean CV deviance and native t-based 95% bounds; lower is better.",
      "Orange point identifies the original chosen q/gamma/alpha, not a new selection.",
      "Native early-stopped paths retained; not external-validation uncertainty or a coefficient path.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_slope_diagnostic <- function(x) {
  alpha <- mean <- lo <- hi <- NULL
  tab <- x$table
  p <- ggplot2::ggplot(tab, ggplot2::aes(alpha, mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "#2878A0", alpha = .16) +
    ggplot2::geom_line(colour = "#2878A0") +
    ggplot2::geom_point(data = tab[tab$selected, , drop = FALSE], colour = "#D77935", size = 2.5) +
    ggplot2::facet_wrap(~q + gamma, labeller = ggplot2::label_both) +
    ggplot2::scale_x_log10() + ggplot2::labs(title = "SLOPE cross-validation", x = "Alpha (log scale)",
      y = "CV deviance") + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11,
    height = max(5, 2 + 2.5 * ceiling(nrow(unique(tab[c("q", "gamma")])) / 3)))
  p
}
