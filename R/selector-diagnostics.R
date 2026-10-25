#' Extract saved feature-selector evidence without fitting a classifier
#' @param x A final fit or saved selector state.
#' @param type One of boruta_summary, boruta_hits, boruta_history, boruta_boxplot, or
#'   selection_history (varSelRF's saved OOB sequence), vsurf_importance or
#'   vsurf_stages. VSURF ranges are native between-forest SD, not confidence
#'   intervals. Stage counts do not reconstruct unsaved stage membership.
#' @return A selector diagnostic with the original evidence table and audit.
#'   Boruta history requires final fitting with capture_diagnostics = TRUE.
#'   Boruta's negative-infinity sentinel after rejection is preserved in
#'   raw_importance and represented as missing in the plotted importance.
#'   Missing importance iterations remain missing. Boruta importance and hits
#'   are selection evidence, not classifier accuracy or causal importance.
#' @export
sb_selector_diagnostics <- function(x, type = "boruta_summary") {
  if (identical(type, "boruta_boxplot")) {
    result <- sb_selector_diagnostics(x, "boruta_history")
    result$type <- "boruta_boxplot"
    return(result)
  }
  state <- if (inherits(x, "stormbreaker_final_fit")) x$selection$state else x
  if (!inherits(state, "stormbreaker_selector_fit")) sb_abort(
    "Supply a final fit with saved selector state, or the selector state itself.",
    "stormbreaker_bad_visualization_argument")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("boruta_summary", "boruta_hits", "boruta_history", "selection_history", "vsurf_importance", "vsurf_stages")) sb_abort(
    "Unsupported selector diagnostic type.", "stormbreaker_unsupported_diagnostic")
  if (type %in% c("boruta_summary", "boruta_hits", "boruta_history")) {
    if (!state$id %in% c(sb_boruta_selector_ids(), "boruta_multiclass_roughfix")) sb_abort(
      "Boruta diagnostics require an actual Boruta selector.", "stormbreaker_unsupported_diagnostic")
    tab <- state$scores
    required <- c("feature", "median_importance", "minimum_importance", "maximum_importance",
      "normalized_hits", "decision", "selected")
    if (!is.data.frame(tab) || !all(required %in% names(tab)) ||
        !identical(tab$feature, state$source_features)) sb_abort(
      "Saved Boruta feature statistics are incomplete.", "stormbreaker_invalid_diagnostic")
    if (type == "boruta_history") {
      history <- state$diagnostic_history$importance
      if (is.null(history)) sb_abort(
        "Boruta history was not captured; enable capture_diagnostics in the original final fit.",
        "stormbreaker_missing_diagnostic_artifacts")
      if (!is.matrix(history) || !is.numeric(history) || !nrow(history) ||
          is.null(colnames(history)) || anyDuplicated(colnames(history)) ||
          !all(state$source_features %in% colnames(history)) || any(history == Inf, na.rm = TRUE)) sb_abort(
        "Saved Boruta importance history is invalid.", "stormbreaker_invalid_diagnostic")
      feature <- rep(colnames(history), each = nrow(history))
      raw <- as.vector(history)
      stopped <- !is.na(raw) & raw == -Inf
      final_decision <- state$scores$decision[match(feature, state$scores$feature)]
      if (any(stopped & (is.na(final_decision) | final_decision != "Rejected"))) sb_abort(
        "Boruta rejection sentinels do not match final rejected features.", "stormbreaker_invalid_diagnostic")
      tab <- data.frame(iteration = rep(seq_len(nrow(history)), ncol(history)),
        feature = feature, raw_importance = raw,
        importance = ifelse(stopped, NA_real_, raw), not_evaluated_after_rejection = stopped,
        kind = ifelse(feature %in% state$source_features, "Predictor", "Shadow reference"),
        final_decision = final_decision)
    }
  } else if (type %in% c("vsurf_importance", "vsurf_stages")) {
    if (!state$id %in% c("vsurf_binary_prediction", "vsurf_multiclass_prediction")) sb_abort(
      "VSURF diagnostics require the original VSURF selector.", "stormbreaker_unsupported_diagnostic")
    tab <- if (type == "vsurf_importance") state$scores else state$selection_history
    required <- if (type == "vsurf_importance") c("feature", "importance", "importance_sd", "rank", "selected") else c("stage", "n_features")
    if (!is.data.frame(tab) || !nrow(tab) || !all(required %in% names(tab))) sb_abort(
      "Saved VSURF evidence is incomplete.", "stormbreaker_invalid_diagnostic")
    if (type == "vsurf_importance" && (!identical(tab$feature, state$source_features) ||
        any(!is.finite(tab$importance)) || any(!is.finite(tab$importance_sd)) || any(tab$importance_sd < 0))) sb_abort(
      "Saved VSURF importance is invalid.", "stormbreaker_invalid_diagnostic")
    if (type == "vsurf_stages" && (!identical(tab$stage, c("threshold", "interpretation", "prediction")) ||
        any(!is.finite(tab$n_features)) || any(tab$n_features < 0))) sb_abort(
      "Saved VSURF stages are invalid.", "stormbreaker_invalid_diagnostic")
  } else {
    if (!state$id %in% c(sb_varselrf_selector_ids(), "varselrf_multiclass_oob_1se")) sb_abort(
      "This selection-history adapter requires varSelRF.", "stormbreaker_unsupported_diagnostic")
    tab <- state$selection_history
    required <- c("number_of_variables", "oob_error", "oob_standard_error")
    if (!is.data.frame(tab) || !nrow(tab) || !all(required %in% names(tab)) ||
        !all(vapply(tab[required], function(z) is.numeric(z) && all(is.finite(z)), logical(1)))) sb_abort(
      "Saved varSelRF selection history is incomplete.", "stormbreaker_invalid_diagnostic")
  }
  structure(list(type = type, selector = state$id, table = tab, audit = state$audit,
    selected_features = state$selected_features,
    decision_before_rough_fix = state$diagnostic_history$decision_before_rough_fix,
    decision_after_rough_fix = state$diagnostic_history$decision_after_rough_fix,
    source = "original selector fit; no fitting or predictive performance calculation"),
    class = "stormbreaker_selector_diagnostic")
}

#' Plot original feature-selection evidence
#' @param x An [sb_selector_diagnostics()] result.
#' @param ncol Columns for Boruta history facets; all predictors and shadow
#'   references are shown, without top-N selection.
#' @return A ggplot. Boruta ranges are iteration minima/maxima, not confidence
#'   intervals. varSelRF error bars are the native OOB standard deviations.
#' @export
sb_plot_selector_diagnostic <- function(x, ncol = 3L) {
  sb_require_curve_plotting()
  if (!inherits(x, "stormbreaker_selector_diagnostic")) sb_abort(
    "Supply a selector diagnostic.", "stormbreaker_bad_visualization_argument")
  ncol <- sb_viz_positive_integer(ncol, "ncol")
  d <- x$table
  feature <- median_importance <- minimum_importance <- maximum_importance <- NULL
  normalized_hits <- decision <- iteration <- importance <- kind <- NULL
  number_of_variables <- oob_error <- oob_standard_error <- importance_sd <- selected <- stage <- n_features <- NULL
  if (x$type %in% c("boruta_summary", "boruta_hits")) {
    d$feature <- factor(d$feature, levels = d$feature[order(d$median_importance)])
    p <- ggplot2::ggplot(d, ggplot2::aes(y = feature, colour = decision))
    if (x$type == "boruta_summary") {
      p <- p + ggplot2::geom_segment(ggplot2::aes(x = minimum_importance,
        xend = maximum_importance, yend = feature)) +
        ggplot2::geom_point(ggplot2::aes(x = median_importance), size = 2) +
        ggplot2::labs(x = "Native importance", title = "Boruta importance distribution summary",
          caption = "Point: median; line: observed iteration minimum to maximum, not a confidence interval.\nAll original features; decisions after any declared tentative rough fix.")
    } else {
      p <- p + ggplot2::geom_point(ggplot2::aes(x = normalized_hits), size = 2) +
        ggplot2::scale_x_continuous(limits = c(0, 1)) +
        ggplot2::labs(x = "Native normalized hit rate", title = "Boruta feature hits",
          caption = "Selection evidence from original Boruta attStats; not prediction accuracy.\nAll original features; decisions after any declared tentative rough fix.")
    }
    p <- p + ggplot2::scale_colour_manual(values = c(Confirmed = "#2878A0",
      Tentative = "#D77935", Rejected = "#909090"), drop = FALSE) +
      ggplot2::labs(y = NULL, colour = "Final decision")
  } else if (x$type == "boruta_boxplot") {
    d <- d[is.finite(d$importance), , drop = FALSE]
    if (!nrow(d)) stop("No observed Boruta importance history is available.")
    d$decision <- ifelse(d$kind == "Shadow reference", "Shadow", d$final_decision)
    medians <- tapply(d$importance, d$feature, stats::median)
    d$feature <- factor(d$feature, levels = names(sort(medians)))
    p <- ggplot2::ggplot(d, ggplot2::aes(importance, feature, fill = decision)) +
      ggplot2::geom_boxplot(width = .65, outlier.size = .7) +
      ggplot2::scale_fill_manual(values = c(Confirmed = "#2CA25F", Tentative = "#F0C75E",
        Rejected = "#D95F59", Shadow = "#4C78A8")) +
      ggplot2::labs(x = "Importance", y = NULL, fill = NULL, title = "Boruta")
  } else if (x$type == "boruta_history") {
    p <- ggplot2::ggplot(d, ggplot2::aes(iteration, importance, colour = kind)) +
      ggplot2::geom_line(na.rm = TRUE) +
      ggplot2::facet_wrap(~feature, ncol = ncol) +
      ggplot2::scale_colour_manual(values = c(Predictor = "#2878A0", "Shadow reference" = "#909090")) +
      ggplot2::labs(x = "Original iteration", y = "Native importance",
        title = "Boruta importance history", colour = NULL,
        caption = "Original predictor and shadow histories; iterations after rejection are not plotted as scores.\nMissing entries remain missing. Selection process evidence, not classifier performance.")
  } else if (x$type == "selection_history") {
    p <- ggplot2::ggplot(d, ggplot2::aes(number_of_variables, oob_error)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = oob_error - oob_standard_error,
        ymax = oob_error + oob_standard_error), width = 0, colour = "#909090") +
      ggplot2::geom_line(colour = "#2878A0") + ggplot2::geom_point(colour = "#2878A0") +
      ggplot2::geom_vline(xintercept = length(x$selected_features), linetype = 2) +
      ggplot2::labs(x = "Number of variables", y = "OOB error",
        title = "varSelRF elimination history",
        caption = "Bars: native OOB standard deviation, not confidence intervals. Dashed: selected subset size.\nSelection-stage OOB error, not untouched-cohort performance.")
  } else if (x$type == "vsurf_importance") {
    d$feature <- factor(d$feature, levels = rev(d$feature[order(d$rank)]))
    p <- ggplot2::ggplot(d, ggplot2::aes(importance, feature, colour = selected)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey75", linewidth = .3) +
      ggplot2::geom_segment(ggplot2::aes(x = importance - importance_sd,
        xend = importance + importance_sd, yend = feature), linewidth = .4) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_colour_manual(values = c("FALSE" = "#909090", "TRUE" = "#2878A0"),
        breaks = c(FALSE, TRUE), labels = c("Not retained", "Retained")) +
      ggplot2::labs(title = "VSURF importance", x = "Native importance", y = NULL, colour = NULL,
        caption = "Points: original forest-averaged importance; bars: plus/minus one native between-forest SD, not confidence intervals.\nAll input features shown; colour indicates the final prediction-stage subset, not stage-specific membership or clinical importance.")
  } else if (x$type == "vsurf_stages") {
    d$stage <- factor(d$stage, levels = c("threshold", "interpretation", "prediction"),
      labels = c("Threshold", "Interpretation", "Prediction"))
    p <- ggplot2::ggplot(d, ggplot2::aes(stage, n_features)) +
      ggplot2::geom_col(fill = "#2878A0", width = .6) +
      ggplot2::scale_y_continuous(breaks = function(limits) unique(round(pretty(limits)))) +
      ggplot2::labs(title = "VSURF selection stages", x = NULL, y = "Retained features",
        caption = "Original threshold, interpretation and prediction-stage subset sizes.\nCounts alone do not identify each intermediate subset; no missing membership or OOB sequence is reconstructed.")
  } else sb_abort("Unsupported selector diagnostic.", "stormbreaker_unsupported_diagnostic")
  p + sb_viz_theme()
}
