#' Inspect original accelerated failure-time coefficients
#' @param fit A successful Weibull, exponential, log-normal or log-logistic
#'   survreg final fit.
#' @param type log_time_coefficients or time_ratios.
#' @return An algorithm diagnostic with original feature coefficients and
#'   exp(coefficient) time ratios per input unit. Intercept and distribution scale
#'   are retained separately. A ratio above one multiplies conditional survival
#'   time quantiles upward; it is not a hazard ratio or causal treatment effect.
#'   Only point estimates are provided; missing covariance is not reconstructed.
#'   Compact saved parameters suffice and no model is refitted.
#' @export
sb_aft_diagnostics <- function(fit, type = "time_ratios") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  state <- fit$learner_state
  if (!state$learner %in% sb_survreg_learner_ids()) sb_abort(
    "This diagnostic requires one of the four survreg AFT models.",
    "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("log_time_coefficients", "time_ratios")) sb_abort(
    "Choose log_time_coefficients or time_ratios.", "stormbreaker_unsupported_diagnostic")
  sb_restore_additional_survival_model(state, state$learner)
  tab <- data.frame(feature = state$features, coefficient = state$coefficient,
    time_ratio = exp(state$coefficient))
  structure(list(type = paste0("aft_", type), table = tab,
    intercept = state$intercept, distribution_scale = state$scale,
    distribution = state$params$distribution, learner = state$learner,
    pipeline_id = fit$pipeline$id,
    source = "original saved survreg log-time location coefficients",
    interpretation = paste("Point estimates per original input unit, with all fitted features retained.",
      "exp(coefficient) multiplies conditional time quantiles; it is not a hazard ratio or a causal effect.",
      "The intercept and distribution scale are retained separately, not ranked as feature effects.",
      "No uncertainty is fabricated; ratios outside floating-point range remain in the table and are omitted from the ratio plot.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_aft_diagnostic <- function(x) {
  feature <- value <- NULL
  tab <- x$table
  ratio <- x$type == "aft_time_ratios"
  tab$value <- if (ratio) tab$time_ratio else tab$coefficient
  tab$value[!is.finite(tab$value) | (ratio & tab$value <= 0)] <- NA_real_
  if (!any(is.finite(tab$value))) sb_abort("No finite values are available on this axis; use log_time_coefficients.",
    "stormbreaker_unavailable_diagnostic")
  tab$feature <- factor(tab$feature, levels = rev(tab$feature))
  p <- ggplot2::ggplot(tab, ggplot2::aes(value, feature)) +
    ggplot2::geom_vline(xintercept = if (ratio) 1 else 0, colour = "grey65", linewidth = .4) +
    ggplot2::geom_point(colour = "#2878A0", size = 2.5, na.rm = TRUE) +
    ggplot2::labs(title = if (ratio) "AFT time ratios" else "AFT log-time coefficients",
      x = if (ratio) "Time ratio per input unit" else "Log-time coefficient per input unit", y = NULL) +
    sb_viz_theme()
  if (ratio) p <- p + ggplot2::scale_x_log10()
  attr(p, "stormbreaker_suggested_size") <- c(width = 9, height = max(5, 2 + .3 * nrow(tab)))
  p
}
