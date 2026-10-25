#' Fit single-gene Cox models and draw a hazard-ratio forest plot
#' @param task A survival task with one independent patient per row.
#' @param conf_level Confidence level, default 0.95.
#' @return List with plot, complete results table and confidence level.
#'   Failed fits remain in the table and are not plotted. No significance
#'   filtering is performed. P values and BH-adjusted P values are returned.
#' @details Each gene is fitted separately using Efron ties. Hazard ratios
#'   describe a one-unit increase on the supplied expression scale, not an
#'   adjusted or causal effect. Wald intervals are pointwise, not simultaneous.
#'   Grouped tasks are not supported. This analysis does not select genes for
#'   a validated pipeline; screening must be repeated within training folds.
#' @export
sb_cox_forest <- function(task, conf_level = 0.95) {
  task <- sb_validate_task_contract(task)
  if (!is.null(task$group_col)) stop("Use one independent patient per row for this Cox forest.")
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) stop("conf_level must be in (0,1).")
  sb_require_curve_plotting()
  x <- task$data[task$predictors]
  response <- survival::Surv(task$data[[task$time_col]], task$data[[task$event_col]])
  tab <- do.call(rbind, lapply(seq_along(x), function(i) {
    sb_score_univariable_cox_feature(names(x)[i], x, response,
      rep(1, nrow(x)), survival::coxph.control(), i)
  }))
  z <- stats::qnorm((1 + conf_level) / 2)
  tab$hr <- exp(tab$coefficient)
  tab$lower <- exp(tab$coefficient - z * tab$standard_error)
  tab$upper <- exp(tab$coefficient + z * tab$standard_error)
  tab$p_adjusted <- stats::p.adjust(tab$p_value, method = "BH")
  tab$plot_available <- tab$status == "success" & is.finite(tab$hr) &
    is.finite(tab$lower) & is.finite(tab$upper) & tab$lower > 0
  shown <- tab[tab$plot_available, ]
  shown$feature <- factor(shown$feature, levels = rev(tab$feature))
  shown$direction <- ifelse(shown$hr >= 1, "HR > 1", "HR < 1")
  shown$direction[shown$hr == 1] <- "HR = 1"
  feature <- hr <- lower <- upper <- direction <- NULL
  p <- ggplot2::ggplot(shown, ggplot2::aes(hr, feature, colour = direction)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper, yend = feature), linewidth = .6) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c("HR > 1" = "#C87558", "HR < 1" = "#287E8E", "HR = 1" = "grey50")) +
    ggplot2::labs(x = "Hazard ratio", y = NULL) + sb_viz_theme() +
    ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank())
  list(plot = p, table = tab, conf_level = conf_level)
}
