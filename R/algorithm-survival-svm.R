#' Inspect native survival-SVM structure without refitting
#' @param fit A final survival-SVM fit with capture_diagnostics = TRUE.
#' @param type score_distribution, basis_weights or coefficients (linear kernel).
#' @return An algorithm diagnostic. All twelve formulation/kernel combinations
#'   retain their native score scale, before the package's Cox calibration.
#'   Regression scores follow native predicted-time units; ranking/hybrid scores
#'   are not probabilities or log hazards. Basis weights refer to kernel rows,
#'   not features or patient importance. Linear coefficients include a separately
#'   returned intercept and reproduce native scores in original input units.
#'   Score and basis-row tables omit identifiers but remain sensitive training
#'   data. No fitting, tuning, calibration or inverse probability calculation.
#' @export
sb_survival_svm_diagnostics <- function(fit, type = "score_distribution") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!fit$pipeline$learner %in% sb_survivalsvm_learner_ids()) sb_abort(
    "This diagnostic requires a survival-SVM learner.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("score_distribution", "basis_weights", "coefficients")) sb_abort(
    "Choose score_distribution, basis_weights or coefficients.", "stormbreaker_unsupported_diagnostic")
  a <- fit$diagnostic_artifacts
  model <- a$learner
  if (is.null(model)) sb_abort("Original survival-SVM artifacts were not captured; enable capture_diagnostics at final fitting.",
    "stormbreaker_missing_diagnostic_artifacts")
  features <- model$features
  formulation <- model$params$type
  kernel <- model$params$kernel
  native <- sb_survivalsvm_native_components(model$engine, formulation, length(features))
  intercept <- NULL
  if (type == "basis_weights") {
    tab <- data.frame(basis_row = seq_along(native$row_weight), value = native$row_weight)
    label <- if (formulation == "regression") "Native support-basis weights" else "Native combined kernel-row weights"
    interpretation <- "Signed weights of the original kernel expansion; not feature or patient importance.\nRanking/hybrid basis rows are not an e1071-style sparse support-vector set. Zero weights are retained."
  } else if (type == "coefficients") {
    if (kernel != "lin_kernel") sb_abort("This nonlinear kernel has no single exact input-feature coefficient vector.",
      "stormbreaker_unavailable_diagnostic")
    w <- as.numeric(crossprod(native$basis, native$row_weight))
    beta <- w / model$audit$scale
    intercept <- native$intercept - sum(w * model$audit$center / model$audit$scale)
    tab <- data.frame(feature = features, value = as.numeric(beta))
    label <- "Native survival-SVM original-unit weights"
    interpretation <- "Original native score equation before Cox calibration; intercept is returned separately.\nNot calibrated log-hazard coefficients, hazard ratios, uncertainty or feature importance."
  } else {
    if (!identical(names(a$predictors), features)) sb_abort("Captured predictor order does not match survival-SVM features.",
      "stormbreaker_invalid_diagnostic")
    if (!requireNamespace("survivalsvm", quietly = TRUE)) sb_abort("Package survivalsvm is required for native score evaluation.",
      "stormbreaker_missing_dependency")
    x <- as.data.frame(sweep(sweep(as.matrix(a$predictors), 2L, model$audit$center, "-"),
      2L, model$audit$scale, "/"))
    score <- as.numeric(stats::predict(model$engine, newdata = x)$predicted)
    if (length(score) != nrow(x)) sb_abort("Native survival-SVM score length is invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(observation = seq_along(score), value = score)
    label <- if (formulation == "regression") "Native predicted-time score" else "Native ranking / hybrid score"
    interpretation <- "Original development-score distribution before Cox calibration; not held-out discrimination.\nScores are not survival probabilities or calibrated log hazards. Censoring is not treated as a binary class."
  }
  if (any(!is.finite(tab$value)) || (!is.null(intercept) && !is.finite(intercept))) sb_abort(
    "Original survival-SVM diagnostic values are nonfinite.", "stormbreaker_invalid_diagnostic")
  structure(list(type = paste0("survival_svm_", type), table = tab, intercept = unname(intercept),
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    formulation = formulation, kernel = kernel, label = label,
    calibration_coefficient = model$audit$calibration_coefficient,
    sensitive_training_data = type != "coefficients", interpretation = interpretation,
    source = "original captured survival-SVM; no fitting"), class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_survival_svm_diagnostic <- function(x) {
  feature <- value <- NULL
  tab <- x$table
  if (x$type == "survival_svm_coefficients") {
    tab$feature <- factor(tab$feature, levels = rev(tab$feature))
    p <- ggplot2::ggplot(tab, ggplot2::aes(value, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
      ggplot2::geom_point(colour = "#2878A0") + ggplot2::labs(y = NULL)
  } else {
    p <- ggplot2::ggplot(tab, ggplot2::aes(value)) +
      ggplot2::stat_ecdf(geom = "step", pad = FALSE, colour = "#2878A0") +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(y = "Cumulative fraction")
    if (length(unique(tab$value)) == 1L) p <- p +
      ggplot2::geom_point(data = data.frame(value = tab$value[1L]), ggplot2::aes(y = 1), colour = "#2878A0")
  }
  p + ggplot2::labs(title = x$label, subtitle = x$learner, x = x$label,
    caption = x$interpretation) + sb_viz_theme()
}
