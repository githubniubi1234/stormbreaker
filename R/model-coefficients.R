#' Inspect the parameters of a saved linear model
#'
#' @param fit Successful portable linear final fit.
#' @param reference_data Optional original-training feature data. Only fitted
#'   numeric features are read. At least two rows are needed for sample SD.
#' @return A `stormbreaker_coefficients` list with a feature table, intercept
#'   (NULL for Cox), score scale and reference sample count. Zero coefficients
#'   remain in the table. Scaled coefficients are beta times sample SD, not a
#'   replacement model and not a SHAP importance or causal effect.
#'   Stepwise GLM, GLMBoost and plsRglm retain their fitted slopes and convert
#'   the original centered intercept to input units. Bias-reduced and Cauchy-prior
#'   logistic fits retain their original estimates, not a new ordinary GLM.
#'   LDA retains the exact posterior log-odds equation from its original class
#'   priors, means and discriminant scaling; no logistic model is fitted.
#'   SLOPE uses its original selected penalized logistic coefficients in input units.
#'   LIBLINEAR L1/L2 logistic models undo their original standardization and
#'   orient native weights toward class 1, including the original bias term.
#'   Their coefficients are penalized log-odds, without ordinary GLM intervals.
#'   Robust Mqle logistic models retain their original robust estimates;
#'   neither ordinary logistic refitting nor robust weight reconstruction occurs.
#'   Survreg AFT models use the original ranking score, negative log-time
#'   location: both slopes and intercept are the negatives of native AFT values.
#' @export
sb_model_coefficients <- function(fit, reference_data = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  state <- fit$learner_state
  supported <- c(
    "logistic", "glmnet_logistic", "lda_binary", "slope_binary", sb_abess_binary_learner_ids(),
    sb_ncvreg_binary_learner_ids(), sb_linear_survival_learner_ids(), sb_extended_logistic_learner_ids(),
    sb_survreg_learner_ids()
  )
  if (!state$learner %in% supported) {
    sb_abort(
      "This learner has no supported linear coefficient view.",
      "stormbreaker_unsupported_coefficients"
    )
  }
  reference_sd <- rep(NA_real_, length(state$features))
  reference_n <- 0L
  if (!is.null(reference_data)) {
    reference <- sb_select_final_feature_frame(
      reference_data, state$features,
      "Reference data", "reference fitted features"
    )
    if (nrow(reference) < 2L) {
      sb_abort(
        "Reference data need at least two rows.",
        "stormbreaker_bad_reference_data"
      )
    }
    reference_sd <- unname(vapply(reference, stats::sd, numeric(1)))
    reference_n <- nrow(reference)
  }
  survival <- state$learner %in% sb_linear_survival_learner_ids()
  aft <- state$learner %in% sb_survreg_learner_ids()
  if (aft) sb_restore_additional_survival_model(state, state$learner)
  parameters <- if (state$learner == "lda_binary") sb_discriminant_parameters(state) else
    if (state$learner %in% sb_extended_logistic_learner_ids()) sb_extended_logistic_parameters(state) else
    if (aft) list(coefficient = -state$coefficient, intercept = -state$intercept) else
    list(coefficient = state$coefficient, intercept = if (survival) NULL else state$intercept)
  structure(list(
    table = data.frame(
      feature = state$features, coefficient = parameters$coefficient,
      reference_sd = reference_sd,
      scaled_coefficient = parameters$coefficient * reference_sd
    ),
    intercept = parameters$intercept,
    score_scale = if (aft) "negative_log_time_location" else
      if (survival) "zero_reference_log_hazard_ratio" else "class_1_log_odds",
    reference_n = reference_n, learner = state$learner,
    pipeline_id = fit$pipeline$id
  ), class = "stormbreaker_coefficients")
}

#' Plot signed coefficients, optionally per reference standard deviation
#'
#' @param x A result from [sb_model_coefficients()].
#' @param scale Raw input-unit coefficients or coefficients multiplied by the
#'   supplied reference sample SD. Predictive models are never modified.
#' @param top_n Maximum number of terms shown, or NULL for all.
#' @param include_zero Whether to display zero coefficients (default FALSE).
#' @return A ggplot object. The table retains all coefficients regardless of
#'   display choices. Direction describes conditional association, not causation.
#' @export
sb_plot_coefficients <- function(x, scale = c("raw", "reference_sd"),
                                 top_n = NULL, include_zero = FALSE) {
  sb_require_curve_plotting()
  if (!inherits(x, "stormbreaker_coefficients")) {
    sb_abort(
      "Supply an sb_model_coefficients result.", "stormbreaker_bad_coefficients"
    )
  }
  scale <- match.arg(scale)
  if (!is.logical(include_zero) || length(include_zero) != 1L || is.na(include_zero)) {
    sb_abort("include_zero must be TRUE or FALSE.", "stormbreaker_bad_coefficients")
  }
  if (scale == "reference_sd" && x$reference_n < 2L) {
    sb_abort(
      "Supply original-training reference data to use SD scaling.",
      "stormbreaker_bad_reference_data"
    )
  }
  d <- x$table
  d$value <- if (scale == "raw") d$coefficient else d$scaled_coefficient
  if (!include_zero) d <- d[d$coefficient != 0, , drop = FALSE]
  d <- d[order(-abs(d$value), seq_len(nrow(d))), , drop = FALSE]
  if (!is.null(top_n)) d <- utils::head(d, sb_viz_positive_integer(top_n, "top_n"))
  d$feature_axis <- factor(d$feature, levels = rev(d$feature))
  d$direction <- ifelse(d$value >= 0, "Positive", "Negative")
  feature_axis <- value <- direction <- NULL
  intercept <- if (is.null(x$intercept)) {
    "Survival score: no intercept; baseline defines absolute survival."
  } else {
    paste0("Intercept = ", format(signif(x$intercept, 5)), " (", x$score_scale, ").")
  }
  score_label <- sb_viz_score_label(x$score_scale)
  ggplot2::ggplot(d, ggplot2::aes(x = value, y = feature_axis, fill = direction)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#9AA6AF", linewidth = 0.4) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::scale_fill_manual(values = c(Negative = "#24547A", Positive = "#B86236")) +
    ggplot2::labs(
      title = "Coefficients", y = NULL, fill = NULL,
      subtitle = paste(nrow(d), "displayed /", nrow(x$table), "fitted terms |", x$learner),
      x = paste0("Coefficient (", score_label,
        if (scale == "raw") "/unit)" else "/SD)"),
      caption = sb_viz_wrap(paste(
        intercept,
        "Conditional associations, not causal effects or SHAP importance.",
        if (scale == "reference_sd") {
          paste("Reference:", x$reference_n, "caller-supplied training rows.")
        } else {
          "Input units are those of the supplied model features."
        }
      ))
    ) +
    sb_viz_theme()
}
