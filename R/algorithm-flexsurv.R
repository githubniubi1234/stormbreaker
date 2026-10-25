#' Inspect original flexible parametric survival parameters
#' @param fit A successful fixed flexsurv or flexsurvspline final fit.
#' @param type location_coefficients, or spline_baseline for spline models.
#' @return Original regression coefficients on their actual parameter scale,
#'   or the saved spline's survival curve with all input predictors set to zero.
#'   This mathematical baseline is not an average or representative patient and
#'   zero can be outside the training predictor range. The curve is confined to
#'   saved boundary knots; internal knot positions are shown without refitting.
#'   Compact parameters suffice. No hazard-ratio conversion or uncertainty is
#'   invented. Distribution parameters and full-precision knots remain separate.
#' @export
sb_flexsurv_diagnostics <- function(fit, type = "location_coefficients") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  spline <- s$learner %in% sb_flexsurvspline_learner_ids()
  if (!spline && !s$learner %in% sb_flexsurv_parametric_learner_ids()) sb_abort(
    "This diagnostic requires a fixed flexsurv or flexsurvspline model.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("location_coefficients", "spline_baseline") || (type == "spline_baseline" && !spline)) sb_abort(
    "Choose location_coefficients, or spline_baseline for a spline model.", "stormbreaker_unsupported_diagnostic")
  if (spline) sb_restore_flexsurvspline_model(s, s$learner) else sb_restore_flexsurv_parametric_model(s, s$learner)
  scale <- if (spline) switch(s$params$scale, hazard = "Log cumulative hazard",
    odds = "Log cumulative odds", normal = "Probit event probability") else
    if (s$params$distribution == "gengamma") "Log-time location (mu)" else "Log rate"
  if (type == "location_coefficients") {
    tab <- data.frame(feature = s$features, coefficient = s$coefficient)
    interpretation <- paste("Original location-regression coefficients per input unit on the stated native scale.",
      "All fitted features retained; ancillary parameters and intercept/baseline coefficients remain separate.",
      "Not a universal hazard ratio, SHAP importance, confidence interval or causal effect.")
  } else {
    times <- sort(unique(c(exp(seq(min(s$knots), max(s$knots), length.out = 201)), exp(s$knots))))
    probability <- flexsurv::psurvspline(times, gamma = s$gamma, knots = s$knots,
      scale = s$params$scale, timescale = s$params$timescale, spline = s$params$spline, lower.tail = FALSE)
    tab <- data.frame(time = times, survival = as.numeric(probability))
    interpretation <- paste("Original spline baseline with every input predictor set to zero, not a typical patient or cohort average.",
      "Zero may lie outside the training predictor range. Time is restricted to the saved boundary knots.",
      "Vertical markers are original internal knots. No uncertainty or clinical validation is implied.")
  }
  structure(list(type = paste0("flexsurv_", type), table = tab, parameter_scale = scale,
    intercept = s$intercept, ancillary = s$ancillary, gamma = s$gamma, knots = s$knots,
    params = s$params, learner = s$learner, pipeline_id = fit$pipeline$id,
    time_unit = fit$development$endpoint$time_unit,
    source = "original compact flexsurv prediction parameters", interpretation = interpretation),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_flexsurv_diagnostic <- function(x) {
  feature <- coefficient <- time <- survival <- NULL
  d <- x$table
  if (x$type == "flexsurv_location_coefficients") {
    d$feature <- factor(d$feature, levels = rev(d$feature))
    p <- ggplot2::ggplot(d, ggplot2::aes(coefficient, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey65", linewidth = .4) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::labs(title = "Location coefficients", x = paste0(x$parameter_scale, " coefficient"), y = NULL)
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(time, survival)) +
      ggplot2::geom_vline(xintercept = exp(x$knots[-c(1, length(x$knots))]), colour = "grey65", linetype = 2) +
      ggplot2::geom_line(colour = "#2878A0", linewidth = .7) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = "Spline baseline",
        x = if (length(x$time_unit)) paste0("Time (", x$time_unit, ")") else "Time (model units)", y = "Survival probability")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 10, height = 5)
  p
}
