sb_coxboost_diagnostic <- function(fit, type) {
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("coefficient_path", "nonzero_path", "cv_curve")) {
    sb_abort("Choose coefficient_path, nonzero_path or cv_curve for CoxBoost.",
      "stormbreaker_unavailable_diagnostic")
  }
  id <- fit$pipeline$learner
  native <- NULL
  if (type == "cv_curve") {
    if (id != "coxboost_cv") sb_abort("Fixed CoxBoost has no CV search.",
      "stormbreaker_unavailable_diagnostic")
    audit <- fit$learner_state$audit
    tab <- audit$tuning_curve
    stop <- audit$optimal_step
    if (!is.data.frame(tab) || !all(c("step", "mean_log_likelihood", "standard_error") %in% names(tab)) ||
        !nrow(tab) || !all(vapply(tab, function(z) is.numeric(z) && all(is.finite(z)), logical(1))) ||
        any(tab$standard_error < 0) || !identical(as.integer(tab$step), seq_len(nrow(tab)) - 1L) ||
        length(stop) != 1L || !stop %in% tab$step) sb_abort(
      "Original CoxBoost CV curve is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
    native <- tab
    measure <- "Cross-validated partial log-likelihood"
    interpretation <- paste0("Original ", audit$tuning_type,
      " CV search; ribbon is +/- one standard error, not a confidence interval.\n",
      "Penalty was chosen before this step-count search; not validation of the complete pipeline.")
  } else {
    model <- fit$diagnostic_artifacts$learner
    if (is.null(model)) sb_abort(
      "Original CoxBoost path was not captured; enable capture_diagnostics at final fitting.",
      "stormbreaker_missing_diagnostic_artifacts")
    engine <- model$engine
    if (!inherits(engine, "CoxBoost") || !identical(model$features, engine$xnames)) sb_abort(
      "Captured CoxBoost engine or feature ordering is invalid.", "stormbreaker_invalid_diagnostic")
    stop <- engine$stepno
    native <- as.matrix(engine$coefficients)
    scale <- engine$sdx
    if (length(stop) != 1L || !is.finite(stop) || stop < 1 || stop != floor(stop) ||
        !identical(dim(native), c(as.integer(stop + 1L), length(model$features))) ||
        any(!is.finite(native)) || length(scale) != ncol(native) || any(!is.finite(scale)) || any(scale <= 0)) sb_abort(
      "Original CoxBoost coefficient path is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
    # Use the same conversion as the original package fitting adapter, not a
    # newly fitted model or a different coefficient convention from coef().
    beta <- sweep(native, 2L, scale, "/")
    if (!isTRUE(all.equal(as.numeric(beta[nrow(beta), ]),
      as.numeric(fit$learner_state$coefficient), tolerance = 1e-12))) sb_abort(
      "Path endpoint disagrees with the saved prediction coefficients.", "stormbreaker_invalid_diagnostic")
    tab <- if (type == "coefficient_path") {
      data.frame(step = rep(seq.int(0L, stop), times = ncol(beta)),
        feature = rep(model$features, each = nrow(beta)), coefficient = as.vector(beta))
    } else {
      data.frame(step = seq.int(0L, stop), nonzero = rowSums(beta != 0))
    }
    measure <- if (type == "coefficient_path") "Original-input CoxBoost coefficient path" else "CoxBoost nonzero predictor count"
    interpretation <- "Original final-fit path including step zero; no additional training or change in stopping.\nCoefficients follow the saved prediction scale, not hazard ratios or uncertainty intervals."
  }
  structure(list(type = paste0("coxboost_", type), table = tab, native_table = native,
    stopping_iteration = as.integer(stop), learner = id, pipeline_id = fit$pipeline$id,
    measure = measure, interpretation = interpretation,
    source = "original CoxBoost artifacts; no refitting"), class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_coxboost_diagnostic <- function(x) {
  step <- feature <- coefficient <- nonzero <- mean_log_likelihood <- standard_error <- NULL
  tab <- x$table
  if (x$type == "coxboost_cv_curve") {
    p <- ggplot2::ggplot(tab, ggplot2::aes(step, mean_log_likelihood)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = mean_log_likelihood - standard_error,
        ymax = mean_log_likelihood + standard_error), fill = "#2878A0", alpha = .15) +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::labs(y = "Partial log-likelihood (higher is better)")
  } else if (x$type == "coxboost_coefficient_path") {
    p <- ggplot2::ggplot(tab, ggplot2::aes(step, coefficient, colour = feature)) +
      ggplot2::geom_line() + ggplot2::labs(y = "Coefficient on saved prediction scale", colour = "Feature")
  } else {
    p <- ggplot2::ggplot(tab, ggplot2::aes(step, nonzero)) +
      ggplot2::geom_step(colour = "#2878A0") +
      ggplot2::scale_y_continuous(breaks = function(limits) unique(round(pretty(limits))), limits = c(0, NA)) +
      ggplot2::labs(y = "Number of nonzero predictors")
  }
  p + ggplot2::geom_vline(xintercept = x$stopping_iteration, colour = "#D77935", linetype = 2) +
    ggplot2::labs(title = x$measure, subtitle = paste0(x$learner, "; original stopping step: ", x$stopping_iteration),
      x = "Boosting step", caption = x$interpretation) + sb_viz_theme()
}
