sb_lssvm_diagnostic <- function(fit, type) {
  if (!is.character(type) || length(type) != 1L || is.na(type)) sb_abort(
    "Choose decision_distribution for LS-SVM.", "stormbreaker_unsupported_diagnostic")
  if (type != "decision_distribution") sb_abort(
    "Radial LS-SVM retains a dense kernel basis, not a sparse support-vector set or input-feature coefficients.",
    "stormbreaker_unavailable_diagnostic")
  a <- fit$diagnostic_artifacts
  if (is.null(a$predictors)) sb_abort("Development predictors were not captured; enable capture_diagnostics at final fitting.",
    "stormbreaker_missing_diagnostic_artifacts")
  e <- unserialize(fit$learner_state$engine_raw)
  features <- fit$learner_state$features
  if (!identical(names(a$predictors), features) || !identical(e$class_labels, c("0", "1")) ||
      !identical(names(e$center), features) || !identical(names(e$scale), features) ||
      any(!is.finite(e$scale)) || any(e$scale <= 0) || !is.matrix(e$alpha) || ncol(e$alpha) != 2L ||
      !is.matrix(e$scaled_training) || nrow(e$scaled_training) != nrow(e$alpha) ||
      ncol(e$scaled_training) != length(features) || length(e$intercept) != 2L ||
      length(e$sigma) != 1L || !is.finite(e$sigma) || e$sigma <= 0) sb_abort(
    "Original LS-SVM kernel basis or class mapping is invalid.", "stormbreaker_invalid_diagnostic")
  z <- sweep(sweep(data.matrix(a$predictors), 2L, e$center, "-"), 2L, e$scale, "/")
  # Same saved kernel expansion, not a refit or inverse probability transform.
  distance <- outer(rowSums(z^2), rowSums(e$scaled_training^2), "+") - 2 * tcrossprod(z, e$scaled_training)
  scores <- sweep(exp(-e$sigma * pmax(distance, 0)) %*% e$alpha, 2L, e$intercept, "+")
  outcome <- a$endpoint$outcome
  if (any(!is.finite(scores)) || length(outcome) != nrow(scores) || anyNA(outcome) ||
      !all(outcome %in% c(0, 1))) sb_abort("Original LS-SVM scores or captured outcomes are invalid.",
    "stormbreaker_invalid_diagnostic")
  tab <- data.frame(observation = seq_len(nrow(scores)), class = as.character(outcome),
    decision = as.numeric(scores[, 2L] - scores[, 1L]))
  structure(list(type = "svm_decision_distribution", table = tab, native_table = unname(scores),
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id, sensitive_training_data = TRUE,
    interpretation = "Original development class-score contrast (class 1 minus class 0); not held-out performance.\nZero is the class boundary. The package uses plogis(contrast), not a fitted calibration; scores are not probabilities.",
    source = "original saved LS-SVM kernel expansion; no fitting"), class = "stormbreaker_algorithm_diagnostic")
}
