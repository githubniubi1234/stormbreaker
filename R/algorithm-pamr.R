#' Inspect original and shrunken PAMR centroids
#' @param fit A successful binary or survival PAMR final fit.
#' @return All feature/class centroids before and after the model's original
#'   fixed threshold. Raw input-unit centroids and deviations from the overall
#'   centroid divided by the saved PAMR SD are retained. This is the native
#'   PAMR scaling, not a new feature-wise z-score. Survival classes are native
#'   censoring-aware early/late groups, not external patient subgroups or
#'   horizon-specific survival probabilities. No threshold tuning is performed.
#' @export
sb_pamr_diagnostics <- function(fit) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% c("pamr_binary", "pamr_survival")) sb_abort(
    "This diagnostic requires PAMR.", "stormbreaker_unsupported_diagnostic")
  sb_require_pamr_survival_dependency()
  e <- if (s$learner == "pamr_survival") s$engine else unserialize(s$engine_raw)
  p <- length(s$features)
  original <- e$centroids
  if (!is.matrix(original) || nrow(original) != p || ncol(original) != 2L ||
      !identical(rownames(original), s$features) ||
      !identical(colnames(original), if (s$learner == "pamr_survival") c("1", "2") else c("0", "1")) ||
      any(!is.finite(original)) || length(e$centroid.overall) != p ||
      length(e$sd) != p || any(!is.finite(e$sd)) || any(e$sd <= 0)) sb_abort(
    "Saved PAMR centroids or scaling are invalid.", "stormbreaker_bad_learner_state")
  dummy <- matrix(0, nrow = p, ncol = 1L)
  shrunken <- pamr::pamr.predict(e, newx = dummy, threshold = s$params$threshold, type = "centroid")
  if (!identical(dim(shrunken), dim(original)) || any(!is.finite(shrunken))) sb_abort(
    "Native shrunken centroids are invalid.", "stormbreaker_invalid_diagnostic")
  classes <- colnames(original)
  labels <- if (s$learner == "pamr_survival") c("Early", "Late") else classes
  make <- function(z, stage) data.frame(feature = rep(s$features, 2L),
    class = rep(classes, each = p), column = rep(paste(stage, labels), each = p),
    stage = stage, original_value = as.numeric(z),
    value = as.numeric((z - e$centroid.overall) / e$sd))
  tab <- rbind(make(original, "Before"), make(shrunken, "After"))
  structure(list(type = "pamr_centroids", table = tab, features = s$features,
    columns = unique(tab$column), threshold = s$params$threshold,
    overall = e$centroid.overall, sd = e$sd, learner = s$learner,
    pipeline_id = fit$pipeline$id, source = "saved native PAMR centroids and pamr.predict centroid transformation",
    interpretation = paste("All features retained before and after the original fixed threshold; no retuning.",
      "Colour is deviation from overall centroid divided by saved PAMR SD; original-unit values remain in the table.",
      "Survival columns denote native early/late survival groups, not patient subgroup analysis or survival probabilities.",
      "No training history, confidence intervals or predictive performance is inferred.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_pamr_diagnostic <- function(x) {
  column <- feature <- value <- NULL
  d <- x$table
  d$column <- factor(d$column, levels = x$columns)
  d$feature <- factor(d$feature, levels = rev(x$features))
  p <- ggplot2::ggplot(d, ggplot2::aes(column, feature, fill = value)) +
    ggplot2::geom_tile(colour = "#D9DFE3", linewidth = .2) +
    ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#C96942", midpoint = 0,
      guide = ggplot2::guide_colourbar(display = "rectangles")) +
    ggplot2::labs(title = "PAMR centroid shrinkage", x = NULL, y = NULL, fill = "PAMR-scaled\ndeviation") + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = max(5, 2 + .3 * length(x$features)))
  p
}
