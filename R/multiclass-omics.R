sb_multiclass_omics_catalog <- function() {
  data.frame(method = "lpc", label = "LPC", supported = TRUE, engine = "lpc",
    feature_result = "selection",
    selection_rule = "Top min(p, max(2, floor(n/log(n)))) absolute native LPC scores")
}

sb_analyze_multiclass_lpc <- function(task, spec) {
  if (!requireNamespace("lpc", quietly = TRUE)) sb_abort(
    "Install lpc for this selector.", "stormbreaker_missing_dependency")
  x <- task$data[task$predictors]
  # Native automatic threshold selection averages ranks 1:50 unconditionally.
  if (ncol(x) < 50L) sb_abort(
    "Native LPC automatic threshold selection requires at least 50 candidate genes.",
    "stormbreaker_lpc_too_few_predictors")
  y <- as.integer(factor(task$data[[task$outcome_col]], levels = task$classes))
  native <- NULL
  invisible(utils::capture.output(native <- withr::with_seed(1L,
    lpc::LPC(t(as.matrix(x)), y, type = "multiclass", soft.thresh = NULL))))
  scores <- as.numeric(native$lpcscores)
  if (length(scores) != ncol(x) || any(!is.finite(scores)) ||
      length(native$soft.thresh) != 1L || !is.finite(native$soft.thresh)) sb_abort(
    "LPC did not return finite gene scores and a selected threshold.",
    "stormbreaker_invalid_lpc_selector")
  size <- sb_high_dimensional_nlogn_size(nrow(x), ncol(x))
  ranking <- order(-abs(scores), seq_along(scores))
  selected <- utils::head(ranking, size)
  tab <- data.frame(feature = names(x), lpc_score = scores,
    base_score = as.numeric(native$tscores), rank = match(seq_along(x), ranking),
    selected = seq_along(x) %in% selected)
  params <- list(native_type = "multiclass", seed = 1L, screen_size = size,
    soft_threshold = native$soft.thresh, internal_resampling = "native_row_based_threshold_search")
  fit <- structure(list(id = "lpc_multiclass_nlogn", source_features = names(x),
    selected_features = names(x)[selected], scores = tab, params = params,
    classes = task$classes), class = "stormbreaker_selector_fit")
  feature <- lpc_score <- NULL
  shown <- tab[selected, ]
  shown$feature <- factor(shown$feature, levels = rev(shown$feature))
  plot <- ggplot2::ggplot(shown, ggplot2::aes(lpc_score, feature)) +
    ggplot2::geom_point(colour = "#2878A0") + ggplot2::theme_minimal() +
    ggplot2::labs(title = "LPC", x = "LPC score", y = NULL)
  list(status = "success", method = spec$method, label = spec$label,
    feature_result = "selection", selection_rule = spec$selection_rule,
    features = tab, class_features = NULL, evidence = tab, parameters = params,
    fit = fit, diagnostics = list(selection_scores = tab), plots = list(selection_scores = plot),
    notes = list(interpretation = paste("Full-input exploratory feature screening.",
      "Native threshold search is row-based, not patient-group-isolated.",
      "LPC scores do not give patient risk probabilities.")),
    plot_status = data.frame(view = "selection_scores", status = "available", reason = ""))
}
