#' Construct an exact points representation of an existing linear model
#' @param fit Successful final model supported by [sb_model_coefficients()].
#' @param reference_data Original-scale training predictors defining display
#'   ranges. No outcomes are used and no model is refitted.
#' @return A `stormbreaker_nomogram` with the unchanged model coefficients,
#'   display ranges, point conversion and optional saved survival baseline.
#'   Models without an absolute survival baseline retain a risk-score axis.
#'   AFT models use their original distribution and scale at the saved horizon,
#'   not a Cox baseline. Without a saved horizon they retain their original
#'   negative-log-time-location score axis.
#' @examples
#' d <- bi_example(40, seed = 42)
#' genes <- c("gene_a", "gene_b")
#' fit <- fit_final(bi_task(d, "patient_id", "outcome", genes),
#'   bi_pipeline("Logistic"))
#' points <- nomogram(fit, d[genes])
#' points$features
#' @export
sb_nomogram <- function(fit, reference_data) {
  parameters <- sb_model_coefficients(fit)
  aft <- NULL
  if (identical(parameters$score_scale, "negative_log_time_location") &&
      !is.null(fit$learner_state$survival_probability)) {
    aft <- sb_validate_additional_survival_probability_state(fit$learner_state, fit$learner_state$learner)
    aft$distribution <- fit$learner_state$params$distribution
    aft$scale <- fit$learner_state$scale
  }
  tab <- parameters$table
  features <- tab$feature
  values <- sb_select_final_feature_frame(reference_data, features,
    "Reference data", "nomogram feature ranges")
  if (!nrow(values)) sb_abort("Reference data must contain observations.",
    "stormbreaker_bad_reference_data")
  tab$minimum <- vapply(values, min, numeric(1))
  tab$maximum <- vapply(values, max, numeric(1))
  tab$minimum_contribution <- pmin(tab$minimum * tab$coefficient, tab$maximum * tab$coefficient)
  tab$span <- abs(tab$coefficient) * (tab$maximum - tab$minimum)
  largest <- max(tab$span)
  unit <- if (largest > 0) largest / 100 else 1
  intercept <- if (is.null(parameters$intercept)) 0 else parameters$intercept
  tab$maximum_points <- tab$span / unit
  structure(list(features = tab, intercept = intercept,
    score_origin = intercept + sum(tab$minimum_contribution), score_per_point = unit,
    total_points = sum(tab$maximum_points), score_scale = parameters$score_scale,
    binary = identical(parameters$score_scale, "class_1_log_odds"),
    baseline = fit$learner_state$baseline_hazard, aft = aft, pipeline_id = parameters$pipeline_id,
    reference_n = nrow(values)), class = "stormbreaker_nomogram")
}

#' Recover the original prediction from nomogram points
#' @param x A result from [sb_nomogram()].
#' @param new_data Original-scale predictors. No outcomes are read.
#' @return List with per-feature points, total points, original score, and
#'   probability when identified. Outside-display-range values are flagged,
#'   never clipped; such points cannot be read from the displayed axes.
#' @examples
#' d <- bi_example(40, seed = 42)
#' genes <- c("gene_a", "gene_b")
#' fit <- fit_final(bi_task(d, "patient_id", "outcome", genes),
#'   bi_pipeline("Logistic"))
#' points <- nomogram(fit, d[genes])
#' nom_values(points, d[1:3, genes])
#' @export
sb_nomogram_values <- function(x, new_data) {
  if (!inherits(x, "stormbreaker_nomogram")) sb_abort(
    "Supply an sb_nomogram result.", "stormbreaker_bad_visualization_argument")
  tab <- x$features
  values <- as.matrix(sb_select_final_feature_frame(new_data, tab$feature,
    "New data", "nomogram predictors"))
  contribution <- sweep(values, 2L, tab$coefficient, "*")
  points <- sweep(contribution, 2L, tab$minimum_contribution, "-") / x$score_per_point
  total <- unname(rowSums(points))
  score <- x$score_origin + total * x$score_per_point
  outside <- unname(rowSums(sweep(values, 2L, tab$minimum, "<") |
    sweep(values, 2L, tab$maximum, ">")) > 0)
  probability <- NULL
  if (x$binary) probability <- stats::plogis(score) else if (!is.null(x$aft)) {
    probability <- sb_nomogram_aft_probability(x$aft, score)
  } else if (!is.null(x$baseline)) {
    hazard <- x$baseline$cumulative_baseline_hazard
    if (!is.numeric(hazard) || length(hazard) != 1L || !is.finite(hazard) || hazard < 0) {
      sb_abort("Saved survival baseline is invalid.", "stormbreaker_invalid_baseline_hazard")
    }
    probability <- if (hazard == 0) rep(1, length(score)) else exp(-exp(log(hazard) + score))
  }
  list(points = points, total_points = total, score = score,
    probability = probability, probability_type = if (x$binary) "class_1" else
      if (!is.null(probability)) "survival_at_saved_horizon" else NULL,
    outside_display_range = outside)
}

#' Draw a nomogram without fitting a surrogate model
#' @param x A result from [sb_nomogram()].
#' @param ticks Approximate ticks per feature axis.
#' @param show_zero Show axes for exactly zero coefficients. By default these
#'   contribute no points and are omitted from the display only; their names
#'   remain in the plot metadata and the original nomogram table.
#' @param reference_data Optional original-scale reference observations for
#'   density overlays. They must lie within the nomogram's existing ranges.
#' @param patient Optional one-row feature frame for exact patient point and
#'   probability markers. Neither option fits or changes the model.
#' @return A ggplot with separate feature, points, total-points and outcome
#'   axes. Total points use their own labelled scale, as in a conventional
#'   nomogram. Short feature axes stagger endpoint labels to avoid collisions.
#'   Instructions are retained in the stormbreaker_notes attribute, not drawn
#'   as a subtitle or figure footer.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' d <- bi_example(40, seed = 42)
#' genes <- c("gene_a", "gene_b")
#' fit <- fit_final(bi_task(d, "patient_id", "outcome", genes),
#'   bi_pipeline("Logistic"))
#' points <- nomogram(fit, d[genes])
#' plot_nomogram(points, patient = d[1, genes, drop = FALSE])
#' @export
sb_plot_nomogram <- function(x, ticks = 5L, show_zero = FALSE,
                             reference_data = NULL, patient = NULL) {
  sb_require_curve_plotting()
  if (!inherits(x, "stormbreaker_nomogram")) sb_abort(
    "Supply an sb_nomogram result.", "stormbreaker_bad_visualization_argument")
  ticks <- sb_viz_positive_integer(ticks, "ticks")
  if (!is.logical(show_zero) || length(show_zero) != 1L || is.na(show_zero)) sb_abort(
    "show_zero must be TRUE or FALSE.", "stormbreaker_bad_visualization_argument")
  tab <- x$features
  omitted <- if (show_zero) character() else tab$feature[tab$coefficient == 0]
  if (!show_zero) tab <- tab[tab$coefficient != 0, , drop = FALSE]
  k <- nrow(tab)
  axis <- data.frame(name = c("Points", tab$feature, "Total points", "Model score"),
    y = rev(seq_len(k + 3)), start = 0,
    end = c(100, tab$maximum_points, 100, 100))
  if (identical(x$score_scale, "negative_log_time_location")) axis$name[k + 3] <- "Negative log-time location"
  marks <- list(data.frame(x = seq(0, 100, 20), y = k + 3,
    label = as.character(seq(0, 100, 20)), offset = -.25))
  for (j in seq_len(k)) {
    count <- max(2L, min(ticks, floor(tab$maximum_points[j] / 15) + 1L))
    v <- unique(seq(tab$minimum[j], tab$maximum[j], length.out = count))
    pos <- (v * tab$coefficient[j] - tab$minimum_contribution[j]) / x$score_per_point
    if (tab$coefficient[j] == 0) {
      v <- paste(format(signif(tab$minimum[j], 3)), "to", format(signif(tab$maximum[j], 3)), "(0 points)")
      pos <- 0
    } else v <- format(signif(v, 3), trim = TRUE)
    offset <- if (length(pos) == 2L && abs(diff(pos)) < 12) c(.25, -.3) else -.25
    marks[[length(marks) + 1L]] <- data.frame(x = pos, y = k + 3 - j,
      label = v, offset = offset)
  }
  total <- unique(seq(0, x$total_points, length.out = max(2L, ticks)))
  pos <- if (x$total_points > 0) total / x$total_points * 100 else rep(0, length(total))
  marks[[length(marks) + 1L]] <- data.frame(x = pos, y = 2,
    label = format(signif(total, 3), trim = TRUE), offset = -.25)
  score <- x$score_origin + total * x$score_per_point
  label <- format(signif(score, 3), trim = TRUE)
  if (x$binary || !is.null(x$baseline) || !is.null(x$aft)) {
    probability <- c(.001, .01, .05, .1, .2, .4, .5, .6, .8, .9, .95, .99, .999)
    desired_score <- if (x$binary) stats::qlogis(probability) else if (!is.null(x$aft))
      sb_nomogram_aft_score(x$aft, probability) else
      log(-log(probability)) - log(x$baseline$cumulative_baseline_hazard)
    pts <- (desired_score - x$score_origin) / x$score_per_point
    keep <- is.finite(pts) & pts >= 0 & pts <= x$total_points
    if (any(keep) && x$total_points > 0) {
      pos <- pts[keep] / x$total_points * 100
      label <- format(probability[keep], trim = TRUE)
      # Keep readable tick spacing, prioritizing the middle probability.
      priority <- order(abs(probability[keep] - .5))
      chosen <- integer()
      for (i in priority) {
        if (!length(chosen) || all(abs(pos[i] - pos[chosen]) >= 8)) chosen <- c(chosen, i)
      }
      chosen <- chosen[order(pos[chosen])]
      pos <- pos[chosen]
      label <- label[chosen]
    } else {
      pos <- 0
      p <- if (x$binary) stats::plogis(x$score_origin) else if (!is.null(x$aft))
        sb_nomogram_aft_probability(x$aft, x$score_origin) else
        exp(-x$baseline$cumulative_baseline_hazard * exp(x$score_origin))
      label <- format(signif(p, 4), trim = TRUE)
    }
    axis$name[k + 3] <- if (x$binary) "Probability (class 1)" else if (!is.null(x$aft))
      paste("Survival at", x$aft$prediction_horizon, x$aft$time_unit) else
      paste("Survival at", x$baseline$prediction_horizon, x$baseline$time_unit)
  }
  marks[[length(marks) + 1L]] <- data.frame(x = pos, y = 1, label = label, offset = -.25)
  data <- do.call(rbind, marks)
  data$label_y <- data$y + data$offset
  start <- end <- y <- label <- label_y <- NULL
  p <- ggplot2::ggplot(data, ggplot2::aes(x, y)) +
    ggplot2::geom_segment(data = axis, ggplot2::aes(x = start, xend = end, y = y, yend = y), inherit.aes = FALSE) +
    ggplot2::geom_segment(ggplot2::aes(xend = x, yend = y - .10)) +
    ggplot2::geom_text(ggplot2::aes(y = label_y, label = label), size = 3) +
    ggplot2::scale_y_continuous(breaks = axis$y, labels = axis$name, limits = c(.4, k + 3.5)) +
    ggplot2::scale_x_continuous(limits = c(-3, 103)) +
    ggplot2::labs(title = "Nomogram", x = NULL, y = NULL) +
    sb_viz_theme() + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank())
  attr(p, "stormbreaker_nomogram_axes") <- axis
  attr(p, "stormbreaker_nomogram_omitted_zero_terms") <- omitted
  attr(p, "stormbreaker_notes") <- c(
    "Sum feature points, then read total-points and outcome axes. Short axes have staggered endpoint labels.",
    "Original score and probability mapping; no refitting. Display ranges come from reference data.",
    if (length(omitted)) paste(length(omitted),
      "exactly zero-coefficient terms omitted from axes; retained in model/table.")
  )
  if (!is.null(reference_data) || !is.null(patient)) p <- sb_nomogram_overlays(p, x, reference_data, patient)
  p
}
