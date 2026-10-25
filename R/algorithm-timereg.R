#' Inspect original additive-hazard and proportional-odds processes
#' @param fit A successful timereg Aalen or proportional-odds final fit.
#' @param type effects or baseline.
#' @return A diagnostic preserving every saved time and cumulative coefficient
#'   for Aalen, or original log-odds coefficients for proportional odds. Baselines
#'   retain the original cumulative intercept or cumulative odds, respectively.
#'   Steps are not smoothed or differentiated; negative/nonmonotone estimates
#'   are retained. No confidence bands or resampling tests are reconstructed.
#'   The zero-predictor baseline is not a representative-patient prediction.
#' @export
sb_timereg_diagnostics <- function(fit, type = "effects") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% sb_timereg_survival_learner_ids()) sb_abort(
    "This diagnostic requires an Aalen or proportional-odds timereg fit.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) || !type %in% c("effects", "baseline")) sb_abort(
    "Choose effects or baseline.", "stormbreaker_unsupported_diagnostic")
  aalen <- s$learner == "aalen_additive_survival"
  zero <- as.data.frame(as.list(stats::setNames(rep(0, length(s$features)), s$features)), check.names = FALSE)
  sb_restore_timereg_survival_model(s, zero)
  if (type == "effects" && aalen) {
    tab <- data.frame(time = rep(s$time, length(s$features)),
      feature = rep(s$features, each = length(s$time)),
      value = as.numeric(s$cumulative_coefficient))
    scale <- "Cumulative coefficient"
  } else if (type == "effects") {
    tab <- data.frame(feature = s$features, value = s$coefficient)
    scale <- "Log cumulative-odds coefficient"
  } else {
    tab <- data.frame(time = s$time, value = if (aalen) s$cumulative_intercept else s$baseline_cumulative_odds)
    scale <- if (aalen) "Cumulative intercept" else "Baseline cumulative odds"
  }
  structure(list(type = paste0("timereg_", if (aalen) "aalen_" else "odds_", type),
    table = tab, value_scale = scale, time_unit = fit$development$endpoint$time_unit,
    features = s$features, time = s$time, learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "original compact timereg coefficient and baseline processes",
    interpretation = paste(if (aalen) "Aalen effects are cumulative coefficients, not instantaneous hazard effects or hazard ratios." else
      "Proportional-odds coefficients act on cumulative event odds, not hazards.",
      "All saved times and coefficients retained. Steps are right-continuous; no smoothing, differentiation, clipping or monotonicity correction.",
      "Baseline refers to all predictors equal to zero, not an average or representative patient.",
      "Feature panels use independent vertical scales. No unavailable variance, confidence bands or tests reconstructed; conditional model parameters are not causal effects.",
      if (aalen && type == "baseline" && (any(tab$value < 0) || any(diff(tab$value) < 0)))
        "This unconstrained cumulative-intercept estimate includes negative values or decreases; do not interpret it as a valid monotone patient cumulative-risk curve." else "")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_timereg_diagnostic <- function(x) {
  time <- value <- feature <- NULL
  d <- x$table
  if (x$type == "timereg_odds_effects") {
    d$feature <- factor(d$feature, levels = rev(x$features))
    p <- ggplot2::ggplot(d, ggplot2::aes(value, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = .4) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::labs(title = "Proportional-odds coefficients", x = x$value_scale, y = NULL)
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(time, value)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = .4) +
      ggplot2::geom_step(direction = "hv", colour = "#2878A0", linewidth = .6) +
      ggplot2::labs(title = if (x$type == "timereg_aalen_effects") "Aalen cumulative effects" else x$value_scale,
        x = paste0("Time (", x$time_unit, ")"), y = x$value_scale)
    if (x$type == "timereg_aalen_effects") {
      p$data$feature <- factor(p$data$feature, levels = x$features)
      p <- p + ggplot2::facet_wrap(~feature, ncol = 1, scales = "free_y")
    }
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11,
    height = if (x$type == "timereg_aalen_effects") max(5, 2.5 * length(x$features)) else 5)
  p
}
