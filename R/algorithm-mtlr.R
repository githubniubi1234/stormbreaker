#' Inspect the original multi-task logistic survival weights
#' @param fit A successful MTLR final fit.
#' @param type feature_weights or biases.
#' @return An algorithm diagnostic preserving every effective native time point
#'   and original weight. Feature weights act on the native normalized inputs;
#'   biases are displayed separately. These coupled interval potentials are not
#'   independent logistic models, hazard ratios or feature-attribution scores.
#'   Exact original-input-unit weights are retained separately for inspection,
#'   with centering absorbed into each bias and native zero-SD behavior retained.
#'   Neither an optimization trace nor uncertainty is reconstructed.
#' @export
sb_mtlr_diagnostics <- function(fit, type = "feature_weights") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (fit$pipeline$learner != "mtlr_survival") sb_abort(
    "This diagnostic requires an MTLR final fit.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("feature_weights", "biases")) sb_abort(
    "Choose feature_weights or biases.", "stormbreaker_unsupported_diagnostic")
  s <- fit$learner_state
  sb_restore_mtlr_survival_model(s)
  e <- s$engine
  w <- e$weight_matrix
  times <- e$time_points
  if (!is.matrix(w) || nrow(w) != length(times) ||
      !identical(colnames(w), c("Bias", s$internal_features)) || any(!is.finite(w)) ||
      !length(times) || any(!is.finite(times)) || any(diff(times) <= 0)) sb_abort(
    "Original MTLR weight/time mapping is invalid.", "stormbreaker_bad_learner_state")
  center <- rep(0, length(s$features))
  sd <- rep(1, length(s$features))
  if (!is.null(e$scale)) {
    center <- unname(e$scale$center)
    sd <- unname(e$scale$sd)
    if (length(center) != length(s$features) || length(sd) != length(s$features) ||
        any(!is.finite(center)) || any(!is.finite(sd)) || any(sd < 0)) sb_abort(
      "Original MTLR normalization is invalid.", "stormbreaker_bad_learner_state")
  }
  slopes <- sweep(w[, -1, drop = FALSE], 2, ifelse(sd == 0, 1, sd), "/")
  slopes[, sd == 0] <- 0
  raw <- cbind(w[, 1] - drop(slopes %*% center), slopes)
  dimnames(raw) <- dimnames(w)
  tab <- if (type == "feature_weights") data.frame(
    time_index = rep(seq_along(times), length(s$features)),
    time = rep(times, length(s$features)), feature = rep(s$features, each = length(times)),
    weight = as.numeric(w[, -1, drop = FALSE])) else data.frame(
      time_index = seq_along(times), time = times, bias = as.numeric(w[, 1]))
  structure(list(type = paste0("mtlr_", type), table = tab, native_weights = w,
    original_unit_weights = raw, time_points = times,
    normalization = data.frame(feature = s$features, center = center, sd = sd),
    requested_intervals = s$requested_intervals, effective_time_points = s$effective_time_points,
    time_unit = fit$development$endpoint$time_unit, learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "original compact MTLR effective time grid, weights and normalization",
    interpretation = paste("Native normalized-input feature weights and biases are separate; all effective time points retained.",
      "Columns represent ordered native time points, not equal-length time intervals.",
      "MTLR couples these potentials to produce a survival distribution; individual weights are not hazard ratios, independent logits or SHAP.",
      "Exact input-unit equivalents are retained separately; no refitting, optimization history or uncertainty is added.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_mtlr_diagnostic <- function(x) {
  time_index <- feature <- weight <- time <- bias <- NULL
  d <- x$table
  if (x$type == "mtlr_feature_weights") {
    d$feature <- factor(d$feature, levels = rev(x$normalization$feature))
    ticks <- unique(round(seq(1, length(x$time_points), length.out = min(8, length(x$time_points)))))
    bound <- max(abs(d$weight))
    if (bound == 0) bound <- 1
    p <- ggplot2::ggplot(d, ggplot2::aes(time_index, feature, fill = weight)) +
      ggplot2::geom_tile(colour = "white", linewidth = .2) +
      ggplot2::scale_x_continuous(breaks = ticks, labels = format(signif(x$time_points[ticks], 3), trim = TRUE), expand = c(0, 0)) +
      ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#C96942", limits = c(-bound, bound),
        guide = ggplot2::guide_colourbar(display = "rectangles")) +
      ggplot2::labs(title = "MTLR feature weights", x = paste0("Ordered time points (", x$time_unit, ")"), y = NULL, fill = "Weight")
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(time, bias)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = .4) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::labs(title = "MTLR biases", x = paste0("Time (", x$time_unit, ")"), y = "Native bias")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11,
    height = if (x$type == "mtlr_feature_weights") max(5, 2 + .3 * nrow(x$normalization)) else 5)
  p
}
