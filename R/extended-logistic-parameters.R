sb_extended_logistic_learner_ids <- function() c(sb_specialized_binary_linear_learner_ids(),
  "brglm2_logistic", "bayesglm_logistic_cauchy",
  "liblinear_l2_logistic_binary", "liblinear_l1_logistic_binary",
  "robust_logistic_mqle")

sb_extended_logistic_parameters <- function(state) {
  if (identical(state$learner, "brglm2_logistic") && sb_is_constant_sum_state(state)) {
    return(sb_constant_sum_linear_parameters(state))
  }
  if (state$learner %in% c("liblinear_l2_logistic_binary",
                          "liblinear_l1_logistic_binary")) {
    e <- unserialize(state$engine_raw)
    m <- e$model
    features <- state$features
    expected_type <- if (state$learner == "liblinear_l2_logistic_binary") 0L else 6L
    if (!inherits(m, "LiblineaR") || !identical(as.integer(m$Type), expected_type) ||
        !identical(sort(as.integer(m$ClassNames)), c(0L, 1L)) ||
        !is.matrix(m$W) || !identical(dim(m$W), c(1L, length(features) + 1L)) ||
        any(!is.finite(m$W)) || !identical(names(e$center), features) ||
        !identical(names(e$scale), features) || any(!is.finite(e$center)) ||
        any(!is.finite(e$scale)) || any(e$scale <= 0) ||
        length(m$Bias) != 1L || !is.finite(m$Bias) || m$Bias < 0) sb_abort(
      "Original LIBLINEAR logistic weights, class order or scaling are invalid.",
      "stormbreaker_bad_learner_state")
    orientation <- if (m$ClassNames[1L] == 1L) 1 else -1
    w <- as.numeric(m$W[1L, seq_along(features)])
    return(list(coefficient = stats::setNames(orientation * w / e$scale, features),
      intercept = unname(orientation * (m$W[1L, length(features) + 1L] * m$Bias -
        sum(w * e$center / e$scale)))))
  }
  if (state$learner %in% sb_specialized_binary_linear_learner_ids()) {
    model <- sb_restore_specialized_binary_linear_model(state, state$learner)
    return(list(coefficient = model$coefficient,
      intercept = model$intercept - sum(model$center * model$coefficient)))
  }
  if (!state$learner %in% c("brglm2_logistic", "bayesglm_logistic_cauchy",
                          "robust_logistic_mqle")) sb_abort(
    "This learner has no extended logistic parameter adapter.", "stormbreaker_unsupported_coefficients")
  p <- if (state$learner == "robust_logistic_mqle") unserialize(state$engine_raw) else state$engine
  if (!is.numeric(p$coefficient) || !identical(names(p$coefficient), state$features) ||
      any(!is.finite(p$coefficient)) || !is.numeric(p$intercept) || length(p$intercept) != 1L ||
      !is.finite(p$intercept)) sb_abort("Saved logistic parameters are invalid.", "stormbreaker_bad_learner_state")
  list(coefficient = p$coefficient, intercept = p$intercept)
}

sb_linear_coefficients_diagnostic <- function(fit) {
  parameters <- sb_model_coefficients(fit)
  structure(list(type = "linear_coefficients", table = parameters$table,
    parameters = parameters, learner = parameters$learner, pipeline_id = parameters$pipeline_id,
    source = "original saved linear predictor in input units",
    interpretation = paste("Signed coefficients on", parameters$score_scale,
      "; original intercept retained separately. Not SHAP, causal effects or confidence intervals.")),
    class = "stormbreaker_algorithm_diagnostic")
}
