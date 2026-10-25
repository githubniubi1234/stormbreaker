#' Inspect original binary SVM and DWD model structure
#' @param fit A successful final binary e1071 SVM, LIBLINEAR SVC, LS-SVM or DWD fit.
#' @param type support_counts, decision_distribution or coefficients (linear only).
#' @return An algorithm diagnostic. Counts and linear weights use portable state.
#'   Decision distributions require capture_diagnostics = TRUE and evaluate the
#'   original model on captured development predictors without fitting. Scores
#'   are oriented toward encoded class 1, not calibrated probabilities or
#'   geometric distances. Zero is the native classification boundary, not
#'   necessarily the calibrated probability 0.5 boundary. Row-level scores remain
#'   sensitive even though patient identifiers are omitted. Linear coefficients
#'   and intercept are converted back to original predictor units; they are not
#'   odds ratios or feature importance. No tuning or support-vector search is run.
#' @export
sb_svm_diagnostics <- function(fit, type = "support_counts") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (fit$pipeline$learner == "lssvm_radial_binary") return(sb_lssvm_diagnostic(fit, type))
  if (fit$pipeline$learner %in% sb_linear_margin_diagnostic_ids()) return(sb_linear_margin_diagnostic(fit, type))
  ids <- c("svm_linear_binary", "svm_radial_binary",
    "svm_radial_weighted_binary", "svm_polynomial_binary", "svm_sigmoid_binary")
  id <- fit$pipeline$learner
  if (!id %in% ids) sb_abort("This diagnostic requires a supported binary LIBSVM learner.",
    "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("support_counts", "decision_distribution", "coefficients")) sb_abort(
    "Choose support_counts, decision_distribution or coefficients.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(fit$learner_state$engine_raw)
  features <- fit$learner_state$features
  if (!inherits(e, "svm") || e$nclasses != 2L || !setequal(e$levels, c("0", "1")) ||
      !setequal(e$labels, 1:2) || length(e$nSV) != 2L || any(!is.finite(e$nSV)) ||
      any(e$nSV < 0) || any(e$nSV != floor(e$nSV)) || sum(e$nSV) != nrow(e$SV) ||
      !identical(colnames(e$SV), features)) sb_abort(
    "Original LIBSVM structure or class mapping is invalid.", "stormbreaker_invalid_diagnostic")
  classes <- e$levels[e$labels]
  orientation <- if (classes[1L] == "1") 1 else -1
  sensitive <- FALSE
  intercept <- NULL
  if (type == "support_counts") {
    tab <- data.frame(class = classes, support_vectors = as.integer(e$nSV))
    native <- list(nSV = e$nSV, labels = e$labels, levels = e$levels)
    interpretation <- "Original support-vector counts by encoded class; not class prevalence or feature importance."
  } else if (type == "coefficients") {
    if (id != "svm_linear_binary" || e$kernel != 0) sb_abort(
      "Nonlinear kernels have no single exact input-feature coefficient vector.", "stormbreaker_unavailable_diagnostic")
    w <- as.numeric(crossprod(e$coefs, e$SV))
    center <- rep(0, length(features))
    scale <- rep(1, length(features))
    if (any(e$scaled)) {
      center[e$scaled] <- e$x.scale[["scaled:center"]]
      scale[e$scaled] <- e$x.scale[["scaled:scale"]]
    }
    if (length(w) != length(features) || any(!is.finite(w)) || any(!is.finite(scale)) ||
        any(scale <= 0) || any(!is.finite(center)) || length(e$rho) != 1L || !is.finite(e$rho)) sb_abort(
      "Original linear LIBSVM scaling or coefficients are invalid.", "stormbreaker_invalid_diagnostic")
    beta <- orientation * w / scale
    intercept <- orientation * (-e$rho - sum(w * center / scale))
    tab <- data.frame(feature = features, coefficient = beta)
    native <- list(coefs = e$coefs, rho = e$rho, scaled_weights = w, scale = scale, center = center)
    interpretation <- "Exact decision-score weights in original predictor units, oriented toward class 1.\nIntercept is returned separately; not log-odds, odds ratios, uncertainty or feature importance."
  } else {
    a <- fit$diagnostic_artifacts
    if (is.null(a$predictors)) sb_abort("Development predictors were not captured; enable capture_diagnostics at final fitting.",
      "stormbreaker_missing_diagnostic_artifacts")
    if (!identical(names(a$predictors), features)) sb_abort(
      "Captured predictor order does not match the original SVM.", "stormbreaker_invalid_diagnostic")
    if (!requireNamespace("e1071", quietly = TRUE)) sb_abort("Package e1071 is required to evaluate native decision values.",
      "stormbreaker_missing_dependency")
    pred <- stats::predict(e, newdata = data.matrix(a$predictors), decision.values = TRUE, probability = FALSE)
    raw <- attr(pred, "decision.values")
    outcome <- a$endpoint$outcome
    if (!is.matrix(raw) || ncol(raw) != 1L || nrow(raw) != nrow(a$predictors) || any(!is.finite(raw)) ||
        length(outcome) != nrow(raw) || anyNA(outcome) || !all(outcome %in% c(0, 1))) sb_abort(
      "Native decision values or captured outcomes are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(observation = seq_len(nrow(raw)), class = as.character(outcome),
      decision = orientation * as.numeric(raw))
    native <- unname(raw)
    sensitive <- TRUE
    interpretation <- "Original development-sample decision scores, oriented toward class 1; not held-out performance.\nZero: native classification boundary, not necessarily calibrated probability 0.5. Scores are not probabilities."
  }
  structure(list(type = paste0("svm_", type), table = tab, native_table = native,
    intercept = intercept, orientation = orientation, learner = id, pipeline_id = fit$pipeline$id,
    sensitive_training_data = sensitive, interpretation = interpretation,
    source = "original binary LIBSVM model; no fitting"), class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_svm_diagnostic <- function(x) {
  class <- support_vectors <- feature <- coefficient <- decision <- NULL
  tab <- x$table
  if (x$type == "svm_support_counts") {
    tab$class <- factor(tab$class, levels = c("0", "1"))
    p <- ggplot2::ggplot(tab, ggplot2::aes(class, support_vectors)) +
      ggplot2::geom_col(fill = "#2878A0", width = .5) +
      ggplot2::scale_y_continuous(breaks = function(limits) unique(round(pretty(limits)))) +
      ggplot2::labs(title = "Support-vector counts", x = "Encoded class", y = "Number of support vectors")
  } else if (x$type == "svm_coefficients") {
    tab$feature <- factor(tab$feature, levels = rev(tab$feature))
    p <- ggplot2::ggplot(tab, ggplot2::aes(coefficient, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
      ggplot2::geom_point(colour = "#2878A0") +
      ggplot2::labs(title = "Original-unit decision weights", x = "Decision-score coefficient", y = NULL)
  } else {
    p <- ggplot2::ggplot(tab, ggplot2::aes(decision, colour = class)) +
      ggplot2::stat_ecdf(geom = "step", pad = FALSE) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey50", linetype = 2) +
      ggplot2::scale_colour_manual(values = c("0" = "#2878A0", "1" = "#D77935")) +
      ggplot2::labs(title = "Development decision-score distributions", x = "Decision score toward class 1",
        y = "Cumulative fraction within class", colour = "Observed class")
  }
  p + ggplot2::labs(subtitle = x$learner, caption = x$interpretation) + sb_viz_theme()
}
