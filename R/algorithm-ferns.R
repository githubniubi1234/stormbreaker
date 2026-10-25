#' Inspect original random-fern splits
#' @param fit Successful random-fern final fit.
#' @param type split_usage or thresholds.
#' @return Original split records and counts, without refitting. Split usage
#'   describes randomly generated structure, not predictive importance.
#'   Threshold plots currently require real-valued predictors.
#' @export
sb_ferns_diagnostics <- function(fit, type = "split_usage") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (s$learner != "random_ferns_binary" || length(type) != 1L ||
      !is.character(type) || is.na(type) || !type %in% c("split_usage", "thresholds")) {
    sb_abort("Choose split_usage or thresholds for random ferns.", "stormbreaker_unsupported_diagnostic")
  }
  e <- unserialize(s$engine_raw)
  depth <- unname(e$parameters["depth"])
  nferns <- unname(e$parameters["ferns"])
  idx <- e$model$splitAttIdxs
  features <- names(e$isStruct$predictorLevels)
  if (!inherits(e, "rFerns") || !identical(features, s$features) ||
      length(idx) != depth * nferns || anyNA(idx) ||
      any(idx < 0L | idx >= length(features))) {
    sb_abort("Original random-fern splits are invalid.", "stormbreaker_invalid_diagnostic")
  }
  splits <- data.frame(fern = rep(seq_len(nferns), each = depth),
    split = rep(seq_len(depth), nferns), feature = features[idx + 1L])
  if (type == "thresholds") {
    if (any(lengths(e$isStruct$predictorLevels) > 0L) ||
        any(e$isStruct$integerPredictors)) {
      sb_abort("Threshold plots require original real-valued predictors; integer and categorical encodings are not converted.",
        "stormbreaker_unavailable_diagnostic")
    }
    if (length(e$model$threReal) != nrow(splits) || any(!is.finite(e$model$threReal))) {
      sb_abort("Original real-valued thresholds are invalid.", "stormbreaker_invalid_diagnostic")
    }
    splits$threshold <- e$model$threReal
    tab <- splits
  } else {
    tab <- data.frame(feature = features, splits = tabulate(idx + 1L, nbins = length(features)))
  }
  structure(list(type = paste0("ferns_", type), table = tab, splits = splits,
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "Original saved rFerns split indices and thresholds",
    interpretation = "Random split usage is not predictive importance. Thresholds describe fitted structure, not clinical cutoffs. No OOB performance is reconstructed."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_ferns_diagnostic <- function(x) {
  feature <- splits <- threshold <- NULL
  if (x$type == "ferns_split_usage") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(splits, feature)) +
      ggplot2::geom_col(fill = "#2878A0", width = .65) +
      ggplot2::labs(title = "Random-fern split usage", x = "Number of splits", y = NULL)
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(threshold, feature)) +
      ggplot2::geom_boxplot(fill = "#A9CEDF", width = .55, outlier.size = .8) +
      ggplot2::labs(title = "Random-fern split thresholds", x = "Original threshold", y = NULL)
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 10, height = max(4, length(unique(x$table$feature)) * .3 + 2))
  p
}
