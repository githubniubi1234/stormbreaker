#' Recover original GAM smooth effects and conditional uncertainty
#' @param fit A final mgcv GAM fit with capture_diagnostics = TRUE.
#' @param features Predictor names; NULL displays every fitted predictor.
#' @param n Grid points spanning each predictor's observed training range.
#' @param level Pointwise normal-approximation interval level.
#' @return An algorithm diagnostic with original centered smooth contributions,
#'   conditional standard errors, and observed predictor values without IDs.
#'   Values can still be sensitive. Other predictors are held at training
#'   medians; this adapter supports only the package's additive univariate GAMs.
#'   These are link-scale associations, not probabilities or causal effects.
#' @export
sb_gam_effects <- function(fit, features = NULL, n = 100L, level = .95) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!fit$pipeline$learner %in% sb_mgcv_gam_learner_ids()) sb_abort(
    "Smooth-effect diagnostics require an original mgcv GAM learner.",
    "stormbreaker_unsupported_diagnostic")
  artifact <- fit$diagnostic_artifacts
  if (is.null(artifact)) sb_abort(
    "Original GAM and reference predictors were not captured; enable capture_diagnostics at final fitting.",
    "stormbreaker_missing_diagnostic_artifacts")
  model <- artifact$learner
  if (!inherits(model, "stormbreaker_mgcv_gam")) sb_abort(
    "Captured GAM model is invalid.", "stormbreaker_invalid_diagnostic")
  n <- sb_viz_positive_integer(n, "n")
  if (n < 2 || !is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1) sb_abort(
    "Use n >= 2 and an interval level strictly between zero and one.",
    "stormbreaker_bad_visualization_argument")
  if (is.null(features)) features <- model$features
  if (!is.character(features) || !length(features) || anyNA(features) ||
      anyDuplicated(features) || !all(features %in% model$features)) sb_abort(
    "Specify distinct fitted predictor names.", "stormbreaker_bad_visualization_argument")
  ref <- sb_select_final_feature_frame(artifact$predictors, model$features,
    "Captured development data", "GAM reference")
  engine <- sb_validate_mgcv_engine(model$engine, model$learner, model$engine_features)
  medians <- vapply(ref, stats::median, numeric(1))
  tab <- support <- list()
  z <- stats::qnorm((1 + level) / 2)
  for (feature in features) {
    i <- match(feature, model$features)
    grid <- as.data.frame(as.list(medians), check.names = FALSE)[rep(1L, n), , drop = FALSE]
    grid[[feature]] <- seq(min(ref[[feature]]), max(ref[[feature]]), length.out = n)
    values <- grid[[feature]]
    names(grid) <- model$engine_features
    term <- engine$smooth[[i]]$label
    predicted <- mgcv::predict.gam(engine, newdata = grid, type = "terms", se.fit = TRUE,
      terms = term, unconditional = FALSE)
    effect <- as.numeric(predicted$fit[, term])
    se <- as.numeric(predicted$se.fit[, term])
    if (length(effect) != n || length(se) != n || any(!is.finite(effect)) ||
        any(!is.finite(se)) || any(se < 0)) sb_abort(
      "The original GAM could not return finite smooth contributions and uncertainty.",
      "stormbreaker_invalid_diagnostic")
    tab[[feature]] <- data.frame(feature = feature, value = values, effect = effect,
      standard_error = se, lower = effect - z * se, upper = effect + z * se)
    support[[feature]] <- data.frame(feature = feature, value = ref[[feature]])
  }
  structure(list(type = "gam_effects", table = do.call(rbind, tab),
    support = do.call(rbind, support), level = level, learner = model$learner,
    pipeline_id = fit$pipeline$id, score_scale = if (fit$task_type == "binary")
      "Centered smooth contribution to log odds" else "Centered smooth contribution to log hazard",
    source = "original mgcv predict(type='terms'); conditional on fitted smoothing parameters"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_gam_effects <- function(x) {
  value <- effect <- lower <- upper <- feature <- NULL
  ggplot2::ggplot(x$table, ggplot2::aes(value, effect)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "#B9D7E5", alpha = .6) +
    ggplot2::geom_hline(yintercept = 0, colour = "#909090", linetype = 2) +
    ggplot2::geom_line(colour = "#2878A0") +
    ggplot2::geom_rug(data = x$support, ggplot2::aes(x = value), inherit.aes = FALSE,
      sides = "b", alpha = .25) +
    ggplot2::facet_wrap(~feature, scales = "free_x", ncol = 2) +
    ggplot2::labs(title = "GAM smooth effects", x = "Predictor value (original units)",
      y = x$score_scale,
      caption = paste0(round(x$level * 100, 1), "% pointwise normal-approximation bands, conditional on fitted smoothing parameters.\nRugs: all captured training values. Centered terms exclude the intercept; not probabilities or causal effects.")) +
    sb_viz_theme()
}
