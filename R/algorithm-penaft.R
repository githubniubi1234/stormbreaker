#' Inspect original Gehan elastic-net AFT coefficients and tuning
#' @param fit A successful penAFT final fit.
#' @param type coefficients, cv_curve or coefficient_path. Paths require
#'   capture_diagnostics = TRUE during the original final fit.
#' @return Original input-unit log-time coefficients and centered intercept,
#'   or the complete original path and linear-predictor CV selection score.
#'   Standardized native path coefficients are converted using saved means/SDs;
#'   the native path remains available separately. The CV score is not AUC or
#'   held-out performance. No standard errors, surrogate fit or missing path
#'   is reconstructed. Package risk is negative centered log-time prediction.
#' @export
sb_penaft_diagnostics <- function(fit, type = "coefficients") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!identical(s$learner, "penaft_gehan_elastic_net")) sb_abort(
    "This diagnostic requires penAFT.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("coefficients", "cv_curve", "coefficient_path")) sb_abort(
    "Choose coefficients, cv_curve or coefficient_path.", "stormbreaker_unsupported_diagnostic")
  path <- NULL
  intercept_path <- NULL
  if (type == "coefficients") tab <- data.frame(feature = s$features, coefficient = s$coefficient) else {
    path <- fit$diagnostic_artifacts$learner$diagnostic_path
    if (is.null(path)) sb_abort("Original penAFT paths were not captured.", "stormbreaker_missing_diagnostic_artifacts")
    f <- path$full_fit
    beta <- as.matrix(f$beta)
    n <- length(f$lambda)
    j <- path$selected_index
    if (!n || !identical(dim(beta), c(as.integer(length(s$features)), as.integer(n))) ||
        length(path$cv_error) != n || any(!is.finite(path$cv_error)) ||
        any(!is.finite(beta)) || any(!is.finite(f$lambda)) || any(f$lambda <= 0) ||
        length(j) != 1L || !is.finite(j) || j < 1 || j > n) sb_abort(
      "Original penAFT path is invalid.", "stormbreaker_invalid_diagnostic")
    intercept_path <- rep(0, n)
    if (isTRUE(f$standardize)) {
      if (length(f$X.sd) != nrow(beta) || length(f$X.mean) != nrow(beta) ||
          any(!is.finite(f$X.sd)) || any(f$X.sd <= 0) || any(!is.finite(f$X.mean))) sb_abort(
        "Original penAFT normalization is invalid.", "stormbreaker_invalid_diagnostic")
      beta <- sweep(beta, 1, f$X.sd, "/")
      intercept_path <- -as.numeric(crossprod(f$X.mean, beta))
    }
    if (!isTRUE(all.equal(as.numeric(beta[, j]), s$coefficient, tolerance = 1e-12)) ||
        !isTRUE(all.equal(intercept_path[j], s$intercept, tolerance = 1e-12)) ||
        f$lambda[j] != s$selected_lambda || path$cv_error[j] != s$selected_cv_error) sb_abort(
      "Original selected penAFT path does not match the saved prediction state.", "stormbreaker_invalid_diagnostic")
    tab <- if (type == "cv_curve") data.frame(lambda = f$lambda, score = path$cv_error,
      selected = seq_len(n) == j) else data.frame(feature = rep(s$features, each = n),
        lambda = rep(f$lambda, length(s$features)), coefficient = as.numeric(t(beta)))
  }
  structure(list(type = paste0("penaft_", type), table = tab, native_path = path,
    intercept = s$intercept, intercept_path = intercept_path, features = s$features,
    selected_lambda = s$selected_lambda, learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "original penAFT Gehan elastic-net model and captured linear-predictor CV path",
    interpretation = paste("Original input-unit centered log-time coefficients, not log hazard ratios; risk has the opposite sign.",
      "Centering intercepts are retained separately. Native standardized coefficients and all path points are preserved.",
      "The marked lambda is the original selection; CV score is not AUC or untouched-test performance. No invented uncertainty.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_penaft_diagnostic <- function(x) {
  feature <- coefficient <- lambda <- score <- NULL
  d <- x$table
  if (x$type == "penaft_coefficients") {
    d$feature <- factor(d$feature, levels = rev(x$features))
    p <- ggplot2::ggplot(d, ggplot2::aes(coefficient, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = .4) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::labs(title = "Gehan elastic-net coefficients", x = "Log-time coefficient", y = NULL)
  } else if (x$type == "penaft_cv_curve") p <- ggplot2::ggplot(d, ggplot2::aes(lambda, score)) +
    ggplot2::geom_line(colour = "#2878A0") +
    ggplot2::geom_point(data = d[d$selected, , drop = FALSE], colour = "#C96942", size = 2.5) +
    ggplot2::scale_x_log10() + ggplot2::labs(title = "Gehan elastic-net selection", x = "Lambda", y = "CV linear-predictor score") else
      p <- ggplot2::ggplot(d, ggplot2::aes(lambda, coefficient, colour = feature)) +
        ggplot2::geom_line() + ggplot2::geom_vline(xintercept = x$selected_lambda, linetype = 2, colour = "grey60") +
        ggplot2::scale_x_log10() + ggplot2::labs(title = "Gehan elastic-net coefficient path", x = "Lambda", y = "Log-time coefficient", colour = NULL)
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11,
    height = if (x$type == "penaft_coefficients") max(5, 2 + .3 * length(x$features)) else 5)
  p
}
