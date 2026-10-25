sb_selector_score_views <- function(id) {
  if (id == "univ_binary_auc_top_k") return(c(discrimination = "AUC (direction-adjusted)"))
  if (id == "univ_cox_top_k") return(c(wald_z = "Wald Z"))
  if (id == "rsf_minimal_depth") return(c(minimal_depth = "Minimal depth", importance = "Importance"))
  if (id == "sis_binary_nlogn") return(c(rank = "Screening rank"))
  if (id == "mrmr_binary_nlogn") return(c(mrmr_score = "mRMR score"))
  if (id %in% c("jmi_binary_nlogn", "cmim_binary_nlogn", "disr_binary_nlogn", "jmim_binary_nlogn")) return(c(information_score = "Selection score"))
  if (id %in% c("lpc_binary_nlogn", "lpc_cox_nlogn")) return(c(lpc_score = "LPC score"))
  if (id %in% c("stabsel_binary_lasso", "stabsel_cox_lasso")) return(c(selection_probability = "Selection probability"))
  character()
}

#' Plot saved screening scores or stability-selection probabilities
#' @param x Saved selector state or a final fit with one selector. For a chain,
#'   pass an individual stage, or use sb_plots to draw every stage.
#' @param metric NULL uses the selector's default diagnostic. RSF also accepts
#'   importance in addition to minimal_depth.
#' @param top_n Maximum number of genes shown, ordered by native rank when
#'   available, otherwise by score (absolute LPC score, smallest RSF depth).
#' @return A ggplot with the displayed native values in its data. The complete
#'   score table is retained in attribute selector_scores. Unevaluated scores
#'   remain missing and are never replaced by zero.
#' @details SIS supplies ranks, not invented importance values. Information
#'   criteria are sequential selection scores, not classifier performance or
#'   causal effects. Stability probabilities are the saved subsampling results;
#'   the line is the original selection cutoff, not a newly optimized threshold.
#'   Univariable binary screening displays the saved direction-adjusted AUC.
#'   Univariable Cox screening displays signed Wald Z, ordered by the original
#'   p-value ranking. Neither view refits the selector.
#' @export
sb_plot_selector_scores <- function(x, metric = NULL, top_n = 20L) {
  sb_require_curve_plotting()
  s <- if (inherits(x, "stormbreaker_final_fit")) x$selection$state else x
  if (!inherits(s, "stormbreaker_selector_fit")) stop("Supply a saved selector state.")
  views <- sb_selector_score_views(s$id)
  if (is.null(metric)) metric <- names(views)[1L]
  if (length(metric) != 1L || is.na(metric) || !metric %in% names(views)) stop("No supported score view for this selector.")
  top_n <- sb_check_positive_integer(top_n, "top_n", minimum = 1L)
  all <- s$scores
  if (!all(c("feature", "selected", metric) %in% names(all))) stop("The original score table is incomplete.")
  d <- all[is.finite(all[[metric]]), , drop = FALSE]
  if (!nrow(d)) stop("No finite native scores are available.")
  priority <- if ("rank" %in% names(d) && metric != "importance") d$rank else
    if (metric == "minimal_depth") d[[metric]] else
      if (metric == "lpc_score") -abs(d[[metric]]) else -d[[metric]]
  d <- utils::head(d[order(priority, seq_len(nrow(d)), na.last = TRUE), , drop = FALSE], top_n)
  d$value <- d[[metric]]
  d$feature <- factor(d$feature, levels = rev(d$feature))
  d$selection <- ifelse(d$selected, "Selected", "Not selected")
  value <- feature <- selection <- NULL
  p <- ggplot2::ggplot(d, ggplot2::aes(value, feature, colour = selection)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(values = c(Selected = "#2878B5", `Not selected` = "grey65")) +
    ggplot2::labs(x = unname(views[metric]), y = NULL, colour = NULL) + sb_viz_theme()
  if (metric == "selection_probability") p <- p +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::geom_vline(xintercept = s$params$cutoff, linetype = 2, colour = "grey50")
  attr(p, "selector_scores") <- all
  p
}
