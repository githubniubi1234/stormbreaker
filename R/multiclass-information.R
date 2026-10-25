sb_multiclass_information_catalog <- function() {
  data.frame(method = c("mrmr", "jmi", "cmim", "disr", "jmim"),
    label = c("mRMR", "JMI", "CMIM", "DISR", "JMIM"), supported = TRUE,
    engine = "praznik", feature_result = "selection",
    selection_rule = "Native sequential information ranking; at most min(p, max(2, floor(n/log(n)))) genes")
}

sb_multiclass_standalone_selectors <- function() {
  c("boruta", "varselrf", "vsurf", sb_multiclass_information_catalog()$method,
    sb_multiclass_omics_catalog()$method)
}

sb_analyze_multiclass_information <- function(task, spec) {
  if (!requireNamespace("praznik", quietly = TRUE)) sb_abort(
    "Install praznik for this selector.", "stormbreaker_missing_dependency")
  x <- task$data[task$predictors]
  x[] <- lapply(x, as.numeric)
  y <- factor(task$data[[task$outcome_col]], levels = task$classes)
  size <- sb_high_dimensional_nlogn_size(nrow(x), ncol(x))
  args <- list(X = x, Y = y, k = size, threads = 1L)
  if (spec$method == "mrmr") args$positive <- FALSE
  native <- do.call(getExportedValue("praznik", toupper(spec$method)), args)
  selected <- as.integer(native$selection)
  if (length(selected) > size || anyNA(selected) || anyDuplicated(selected) ||
      any(!selected %in% seq_along(x)) || length(native$score) != length(selected) ||
      any(!is.finite(native$score))) sb_abort(
    "The information selector returned an invalid feature subset.",
    "stormbreaker_invalid_information_selector")
  tab <- data.frame(feature = names(x), score = NA_real_,
    rank = match(seq_along(x), selected), selected = seq_along(x) %in% selected)
  tab$score[selected] <- as.numeric(native$score)
  params <- list(screen_size = size, screen_size_rule = spec$selection_rule,
    threads = 1L, native_method = toupper(spec$method), selected_size = length(selected),
    discretization = "native_equal_width_bins_for_numeric_predictors")
  if (spec$method == "mrmr") params$positive <- FALSE
  fit <- structure(list(id = paste0(spec$method, "_multiclass_nlogn"),
    source_features = names(x), selected_features = names(x)[selected],
    scores = tab, params = params, classes = task$classes), class = "stormbreaker_selector_fit")
  feature <- score <- NULL
  shown <- tab[selected, ]
  shown$feature <- factor(shown$feature, levels = rev(shown$feature))
  plot <- ggplot2::ggplot(shown, ggplot2::aes(score, feature)) +
    ggplot2::geom_point(colour = "#2878A0") + ggplot2::theme_minimal() +
    ggplot2::labs(title = spec$label, x = "Selection score", y = NULL)
  list(status = "success", method = spec$method, label = spec$label,
    feature_result = "selection", selection_rule = spec$selection_rule,
    features = tab, class_features = NULL, evidence = tab, parameters = params,
    fit = fit, diagnostics = list(selection_scores = tab), plots = list(selection_scores = plot),
    notes = list(interpretation = paste("Full-input feature selection, not held-out evaluation.",
      "Scores are native sequential selection criteria, not patient risks or marginal gene effects.")),
    plot_status = data.frame(view = "selection_scores", status = "available", reason = ""))
}
