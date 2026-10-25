#' Inspect original Cancerclass centroids and probability mapping
#' @param fit Successful Cancerclass final fit.
#' @param type centroids or score_mapping.
#' @param score_range Finite increasing native-score plotting range, default
#'   c(-1, 1). This is a display range, not an observed training-score range.
#' @return Original centroids or the saved inner-OOF-fitted logistic mapping.
#'   The mapping is not a calibration-reliability curve or a validation result.
#'   Original OOF scores are not present in compact fits and are not recreated.
#'   No new calibration fit or surrogate model is introduced.
#' @export
sb_cancerclass_diagnostics <- function(fit, type = "centroids", score_range = c(-1, 1)) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (s$learner != "cancerclass_binary" || !is.character(type) ||
      length(type) != 1L || is.na(type) || !type %in% c("centroids", "score_mapping")) sb_abort(
    "Choose centroids or score_mapping for Cancerclass.", "stormbreaker_unsupported_diagnostic")
  if (!is.matrix(s$predictor) || !setequal(rownames(s$predictor), s$features) ||
      !all(c("0", "1") %in% colnames(s$predictor))) sb_abort(
    "Original Cancerclass centroids are invalid.", "stormbreaker_invalid_diagnostic")
  if (type == "centroids") {
    m <- s$predictor[s$features, c("0", "1"), drop = FALSE]
    tab <- data.frame(feature = rep(s$features, 2L), class = rep(c("0", "1"), each = length(s$features)), mean = as.numeric(m))
  } else {
    if (!is.numeric(score_range) || length(score_range) != 2L ||
        any(!is.finite(score_range)) || score_range[2L] <= score_range[1L]) sb_abort(
      "score_range must contain two finite increasing values.", "stormbreaker_invalid_diagnostic")
    score <- seq(score_range[1L], score_range[2L], length.out = 201L)
    b <- s$calibration_coefficients
    tab <- data.frame(score = score, probability = stats::plogis(b[[1L]] + b[[2L]] * score))
  }
  structure(list(type = paste0("cancerclass_", type), table = tab,
    calibration_coefficients = s$calibration_coefficients,
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "Original saved Cancerclass centroids and existing OOF-fitted score mapping",
    interpretation = "The score-mapping curve is the existing probability formula, not a calibration-reliability assessment. Its domain is an explicit display range, not observed support. No OOF reconstruction, refit or surrogate."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_cancerclass_diagnostic <- function(x) {
  mean <- feature <- class <- score <- probability <- NULL
  if (x$type == "cancerclass_centroids") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(mean, feature, colour = class)) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::scale_colour_manual(values = c("0" = "#2878A0", "1" = "#B86236")) +
      ggplot2::labs(title = "Cancerclass centroids", x = "Centroid value", y = NULL, colour = "Class")
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(score, probability)) +
      ggplot2::geom_line(colour = "#2878A0", linewidth = .7) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = "Cancerclass probability mapping", x = "Native z score", y = "Class-1 probability")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 10, height = 4)
  p
}
