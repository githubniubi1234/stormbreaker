sb_linear_margin_diagnostic_ids <- function() {
  c("liblinear_l2_l2loss_svc_binary", "liblinear_l2_l1loss_svc_binary",
    "liblinear_l1_l2loss_svc_binary", "dwd_binary")
}

sb_linear_margin_diagnostic <- function(fit, type) {
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("coefficients", "decision_distribution", "support_counts")) sb_abort(
    "Choose coefficients or decision_distribution for this linear margin model.", "stormbreaker_unsupported_diagnostic")
  if (type == "support_counts") sb_abort(
    "This engine does not retain an e1071-style support-vector set; counts are unavailable.",
    "stormbreaker_unavailable_diagnostic")
  state <- fit$learner_state
  features <- state$features
  id <- fit$pipeline$learner
  if (id == "dwd_binary") {
    e <- state$engine
    beta <- e$coefficient
    intercept <- e$intercept
    native <- list(coefficient = beta, intercept = intercept)
    orientation <- 1
  } else {
    e <- unserialize(state$engine_raw)
    m <- e$model
    if (!inherits(m, "LiblineaR") || !m$Type %in% c(1L, 3L, 5L) ||
        !setequal(m$ClassNames, c(0L, 1L)) || !is.matrix(m$W) ||
        !identical(dim(m$W), c(1L, length(features) + 1L)) ||
        !identical(names(e$center), features) || !identical(names(e$scale), features) ||
        any(!is.finite(e$scale)) || any(e$scale <= 0) || any(!is.finite(e$center)) ||
        length(m$Bias) != 1L || !is.finite(m$Bias) || m$Bias < 0) sb_abort(
      "Original LIBLINEAR weights, bias or scaling are invalid.", "stormbreaker_invalid_diagnostic")
    orientation <- if (m$ClassNames[1L] == 1L) 1 else -1
    w <- as.numeric(m$W[1L, seq_along(features)])
    beta <- stats::setNames(orientation * w / e$scale, features)
    intercept <- orientation * (m$W[1L, length(features) + 1L] * m$Bias - sum(w * e$center / e$scale))
    native <- list(W = m$W, Bias = m$Bias, ClassNames = m$ClassNames,
      center = e$center, scale = e$scale)
  }
  if (!is.numeric(beta) || !identical(names(beta), features) || any(!is.finite(beta)) ||
      length(intercept) != 1L || !is.finite(intercept)) sb_abort(
    "Original linear decision function is invalid.", "stormbreaker_invalid_diagnostic")
  if (type == "coefficients") {
    tab <- data.frame(feature = features, coefficient = as.numeric(beta))
    interpretation <- "Exact original-unit decision weights toward class 1; intercept returned separately.\nNot log-odds, odds ratios, confidence intervals or feature importance; no refit."
  } else {
    a <- fit$diagnostic_artifacts
    if (is.null(a$predictors)) sb_abort("Development predictors were not captured; enable capture_diagnostics at final fitting.",
      "stormbreaker_missing_diagnostic_artifacts")
    if (!identical(names(a$predictors), features) || length(a$endpoint$outcome) != nrow(a$predictors) ||
        anyNA(a$endpoint$outcome) || !all(a$endpoint$outcome %in% c(0, 1))) sb_abort(
      "Captured predictor order or outcomes are invalid.", "stormbreaker_invalid_diagnostic")
    decision <- as.numeric(data.matrix(a$predictors) %*% beta + intercept)
    if (any(!is.finite(decision))) sb_abort("Original linear decision scores are nonfinite.",
      "stormbreaker_invalid_diagnostic")
    tab <- data.frame(observation = seq_along(decision), class = as.character(a$endpoint$outcome), decision = decision)
    native <- NULL
    mapping <- if (id == "dwd_binary") "plogis(2 * score)" else "plogis(score)"
    interpretation <- paste0("Original development-sample decision scores toward class 1; not held-out performance.\n",
      "Zero is the class boundary. The package maps these scores with ", mapping, "; this is not calibration.")
  }
  structure(list(type = paste0("svm_", type), table = tab, native_table = native,
    intercept = unname(intercept), orientation = orientation, learner = id, pipeline_id = fit$pipeline$id,
    sensitive_training_data = type == "decision_distribution", interpretation = interpretation,
    source = "original saved linear margin model; no fitting"), class = "stormbreaker_algorithm_diagnostic")
}
