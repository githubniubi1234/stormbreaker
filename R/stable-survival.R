# Preserve a finite scale already defined by the fitted ordinary Cox model.
# No clipping, external-data refitting, or conversion of arbitrary risk scores.
sb_final_cox_log_cumulative_hazard <- function(fit, values) {
  if (is.null(values) || !identical(fit$learner_state$learner, "coxph") ||
      !any(values$survival_probability %in% c(0, 1))) return(NULL)
  h0 <- fit$learner_state$baseline_hazard$cumulative_baseline_hazard
  if (!is.numeric(h0) || length(h0) != 1L || !is.finite(h0) || h0 <= 0) return(NULL)
  stable <- log(h0) + values$risk
  if (!all(is.finite(stable))) return(NULL)
  # Emit only when this additional representation agrees with the saved
  # probability. Do not make previously valid predictions fail here.
  if (any(abs(exp(-exp(stable)) - values$survival_probability) >
      64 * .Machine$double.eps)) return(NULL)
  sb_validate_stable_survival_scale(stable, values$survival_probability)
}

sb_validate_stable_survival_scale <- function(value, probability) {
  if (!is.numeric(value) || length(value) != length(probability) ||
      anyNA(value) || any(!is.finite(value)) || !is.numeric(probability) ||
      anyNA(probability) || any(!is.finite(probability)) ||
      any(probability < 0 | probability > 1)) {
    sb_abort("The saved log cumulative hazard must be finite and patient-aligned.",
      "stormbreaker_bad_external_prediction")
  }
  reconstructed <- exp(-exp(value))
  if (any(abs(reconstructed - probability) > 64 * .Machine$double.eps)) {
    sb_abort("Saved log cumulative hazard and survival probability disagree.",
      "stormbreaker_bad_external_prediction")
  }
  unname(as.numeric(value))
}

sb_external_stable_scale <- function(rows, probability) {
  if (!"log_cumulative_hazard" %in% names(rows)) return(NULL)
  value <- sb_external_column(rows, "log_cumulative_hazard")
  sb_validate_stable_survival_scale(value, probability)
}
