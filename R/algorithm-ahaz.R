#' Inspect original constant-coefficient additive-hazard models
#' @param fit A successful Lin-Ying, LASSO or SSCAD ahaz final fit.
#' @param type coefficients, bic_path or coefficient_path. Paths require
#'   capture_diagnostics = TRUE during the original penalized final fit.
#' @return Original input-unit additive-hazard coefficients or native penalized
#'   coefficient/BIC paths. All path points, zero coefficients and original
#'   selected lambda are retained. A separately saved Cox-Breslow calibration
#'   is not the raw additive-hazard equation and is not plotted here. No missing
#'   path, uncertainty, baseline or likelihood is reconstructed.
#' @export
sb_ahaz_diagnostics <- function(fit, type = "coefficients") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% sb_ahaz_survival_learner_ids()) sb_abort(
    "This diagnostic requires an ahaz learner.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("coefficients", "bic_path", "coefficient_path")) sb_abort(
    "Choose coefficients, bic_path or coefficient_path.", "stormbreaker_unsupported_diagnostic")
  if (!is.numeric(s$coefficient) || length(s$coefficient) != length(s$features) || any(!is.finite(s$coefficient))) sb_abort(
    "Saved additive-hazard coefficients are invalid.", "stormbreaker_bad_learner_state")
  path <- NULL
  if (type == "coefficients") {
    tab <- data.frame(feature = s$features, coefficient = s$coefficient)
  } else {
    if (s$learner == "ahaz_lin_ying") sb_abort(
      "Unpenalized Lin-Ying has no penalized BIC path.", "stormbreaker_unsupported_diagnostic")
    path <- fit$diagnostic_artifacts$learner$diagnostic_path
    if (is.null(path)) sb_abort("Original ahaz paths were not captured.", "stormbreaker_missing_diagnostic_artifacts")
    n <- length(path$lambda)
    if (!n || !identical(path$criterion_type, "BIC") || any(!is.finite(path$lambda)) || any(path$lambda <= 0) ||
        length(path$criterion) != n || length(path$df) != n ||
        any(!is.finite(path$criterion)) || any(!is.finite(path$df)) ||
        !identical(dim(path$beta), c(as.integer(length(s$features)), as.integer(n)))) sb_abort(
      "Original additive-hazard path is invalid.", "stormbreaker_invalid_diagnostic")
    selected <- which(path$lambda == path$selected_lambda)
    if (length(selected) != 1L || !isTRUE(all.equal(as.numeric(path$beta[, selected]), s$coefficient, tolerance = 1e-12))) sb_abort(
      "Original selected path does not match the saved model.", "stormbreaker_invalid_diagnostic")
    tab <- if (type == "bic_path") data.frame(lambda = path$lambda, criterion = path$criterion,
      df = path$df, selected = seq_len(n) == selected) else data.frame(
        feature = rep(s$features, each = n), lambda = rep(path$lambda, length(s$features)),
        coefficient = as.numeric(t(as.matrix(path$beta))))
  }
  structure(list(type = paste0("ahaz_", type), table = tab, native_path = path,
    features = s$features, selected_lambda = s$audit$selected_lambda,
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "original ahaz additive-hazard coefficients and saved native BIC path",
    interpretation = paste("Original input-unit additive-hazard coefficients, not log hazard ratios or causal effects.",
      "Native BIC is the original selection criterion, not held-out performance; marked lambda is the original choice.",
      "Full path and zero coefficients retained, without refitting or uncertainty reconstruction.",
      "Any separate saved Cox-Breslow probability calibration is not represented by these raw coefficients.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_ahaz_diagnostic <- function(x) {
  feature <- coefficient <- lambda <- criterion <- NULL
  d <- x$table
  if (x$type == "ahaz_coefficients") {
    d$feature <- factor(d$feature, levels = rev(x$features))
    p <- ggplot2::ggplot(d, ggplot2::aes(coefficient, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = .4) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::labs(title = "Additive-hazard coefficients", x = "Coefficient per input unit", y = NULL)
  } else if (x$type == "ahaz_bic_path") {
    p <- ggplot2::ggplot(d, ggplot2::aes(lambda, criterion)) +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::geom_point(data = d[d$selected, , drop = FALSE], colour = "#C96942", size = 2.5) +
      ggplot2::scale_x_log10() + ggplot2::labs(title = "Additive-hazard BIC path", x = "Lambda", y = "Native BIC criterion")
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(lambda, coefficient, colour = feature)) +
      ggplot2::geom_line() + ggplot2::geom_vline(xintercept = x$selected_lambda, colour = "grey60", linetype = 2) +
      ggplot2::scale_x_log10() + ggplot2::labs(title = "Additive-hazard coefficient path", x = "Lambda", y = "Coefficient", colour = NULL)
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11,
    height = if (x$type == "ahaz_coefficients") max(5, 2 + .3 * length(x$features)) else 5)
  p
}
