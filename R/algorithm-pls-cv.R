sb_pls_cv_diagnostic <- function(fit) {
  if (fit$pipeline$learner != "plsrcox_cv") sb_abort(
    "This PLS learner has no component-count CV search.", "stormbreaker_unavailable_diagnostic")
  audit <- fit$learner_state$audit
  if (identical(audit$component_selection, "single_predictor_has_one_feasible_component")) sb_abort(
    "One predictor has one feasible component; no CV search was performed.", "stormbreaker_unavailable_diagnostic")
  value <- audit$cv_error
  spread <- audit$cv_standard_error
  maximum <- audit$maximum_components
  selected <- audit$selected_components
  if (!is.numeric(maximum) || length(maximum) != 1L || !is.finite(maximum) || maximum < 1L ||
      maximum != floor(maximum) || !is.numeric(value) || length(value) != maximum + 1L ||
      !is.numeric(spread) || length(spread) != length(value) || any(is.infinite(value)) ||
      any(is.infinite(spread)) || any(spread < 0, na.rm = TRUE) ||
      !is.numeric(selected) || length(selected) != 1L || !is.finite(selected) ||
      selected < 1L || selected > maximum || selected != floor(selected) ||
      !is.finite(value[selected + 1L])) sb_abort(
    "Saved plsRcox CV criterion or component count is invalid.", "stormbreaker_invalid_diagnostic")
  tab <- data.frame(component = seq.int(0L, maximum), iauc = value,
    native_spread = spread, selected = seq.int(0L, maximum) == selected)
  structure(list(type = "pls_cv_curve", table = tab,
    native_table = list(cv_error = value, cv_standard_error = spread),
    selected_components = as.integer(selected), folds = length(audit$folds),
    search_warnings = audit$component_search_warnings,
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    interpretation = paste0("Original component search, including zero-component baseline; not whole-pipeline validation.\n",
      "Bars: saved native cv.se5 spread (fold SD / fold count), not conventional SE or confidence intervals.\n",
      "Missing values remain missing; fold-level completion cannot be recovered from historical summaries."),
    source = "saved native plsRcox criterion 5 (iAUC_sh); no fitting"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_pls_cv_diagnostic <- function(x) {
  component <- iauc <- native_spread <- NULL
  tab <- x$table
  p <- ggplot2::ggplot(tab, ggplot2::aes(component, iauc)) +
    ggplot2::geom_line(colour = "#2878A0", na.rm = TRUE) +
    ggplot2::geom_point(colour = "#2878A0", na.rm = TRUE) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = iauc - native_spread,
      ymax = iauc + native_spread), width = .15, colour = "#2878A0", na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = x$selected_components, colour = "#D77935", linetype = 2) +
    ggplot2::scale_x_continuous(breaks = function(limits) unique(round(pretty(limits)))) +
    ggplot2::labs(title = "PLS-Cox component-count selection",
      subtitle = paste0(x$folds, " original folds; selected components: ", x$selected_components,
        if (length(x$search_warnings)) paste0("; retained search warnings: ", length(x$search_warnings)) else ""),
      x = "Number of components", y = "Integrated AUC (iAUC_sh; higher is better)",
      caption = x$interpretation) + sb_viz_theme()
  p
}
