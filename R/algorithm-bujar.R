#' Inspect saved Buckley-James linear boosting effects
#' @param fit A successful Bujar linear final fit.
#' @param type coefficients or residual_survival. The latter requires an
#'   originally requested prediction horizon and its saved residual curve.
#' @return Original log-time coefficients, or the saved event-step residual
#'   Kaplan-Meier distribution used by the model's probability mapping.
#'   The intercept is retained separately. Coefficients are not log hazard
#'   ratios; the package risk score is negative predicted log time. The residual
#'   distribution is training-derived, not held-out calibration or a patient
#'   survival curve. No variance, censored-only steps or training iteration
#'   history is reconstructed from compact state.
#' @export
sb_bujar_diagnostics <- function(fit, type = "coefficients") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!identical(s$learner, sb_bujar_linear_learner_ids())) sb_abort(
    "This diagnostic requires the Bujar linear learner.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("coefficients", "residual_survival")) sb_abort(
    "Choose coefficients or residual_survival.", "stormbreaker_unsupported_diagnostic")
  if (type == "coefficients") {
    if (length(s$coefficient) != length(s$features) || any(!is.finite(s$coefficient)) ||
        length(s$intercept) != 1L || !is.finite(s$intercept)) sb_abort(
      "Original Bujar coefficients are invalid.", "stormbreaker_bad_learner_state")
    tab <- data.frame(feature = s$features, coefficient = s$coefficient)
  } else {
    if (is.null(s$survival_probability)) sb_abort(
      "No residual distribution was saved: the original fit did not request a prediction horizon.",
      "stormbreaker_unavailable_diagnostic")
    probability <- sb_validate_bujar_probability_state(s)
    tab <- data.frame(residual = probability$residual_time, survival = probability$residual_survival)
  }
  structure(list(type = paste0("bujar_", type), table = tab,
    features = s$features, intercept = s$intercept, learner = s$learner,
    pipeline_id = fit$pipeline$id,
    source = "original Bujar compact coefficients and saved analysis residual Kaplan-Meier mapping",
    interpretation = paste("Coefficients act on predicted log time, not log hazard ratios; risk is negative predicted log time.",
      "The intercept is retained separately. Residual event steps are training-derived, not held-out calibration or patient survival.",
      "No refitting, smoothing, invented iteration history, censored-only steps or uncertainty bands.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_bujar_diagnostic <- function(x) {
  feature <- coefficient <- residual <- survival <- NULL
  d <- x$table
  if (x$type == "bujar_coefficients") {
    d$feature <- factor(d$feature, levels = rev(x$features))
    p <- ggplot2::ggplot(d, ggplot2::aes(coefficient, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = .4) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::labs(title = "Buckley-James coefficients", x = "Log-time coefficient", y = NULL)
  } else p <- ggplot2::ggplot(d, ggplot2::aes(residual, survival)) +
    ggplot2::geom_step(direction = "hv", colour = "#2878A0") +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = "Residual survival distribution", x = "Log-time residual", y = "Residual survival")
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11,
    height = if (x$type == "bujar_coefficients") max(5, 2 + .3 * length(x$features)) else 5)
  p
}
