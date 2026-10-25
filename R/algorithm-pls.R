#' Inspect native PLS component structure
#' @param fit A final plsRglm, PLS-DA, sparse PLS-DA or plsRcox model fitted with
#'   capture_diagnostics = TRUE for components. CV curves use saved tuning and
#'   do not require capture. Historical compact binary fits lack components.
#' @param type loadings, score_distribution or cv_curve (CV plsRcox only).
#' @return An algorithm diagnostic with all original components and features.
#'   mixOmics returns its native X component weights; plsRglm/plsRcox return X
#'   reconstruction loadings (pp). These definitions are not interchangeable.
#'   Scores describe the development samples, not held-out discrimination.
#'   Score tables omit patient identifiers but remain sensitive row-level data.
#'   No new decomposition, rotation, sign alignment or fitting is performed.
#'   CV criterion 5 is integrated AUC (iAUC_sh), higher is better, despite the
#'   historical field name cv_error. Native cv.se5 is retained as native spread,
#'   not called a conventional standard error or confidence interval. The native
#'   implementation uses fold standard deviation divided by the fold count.
#' @export
sb_pls_diagnostics <- function(fit, type = "loadings") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  id <- fit$pipeline$learner
  if (!id %in% c("plsrglm_binary", "plsda_binary", "splsda_binary",
      "plsrcox", "plsrcox_cv")) sb_abort("This diagnostic requires a supported native PLS learner.",
    "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("loadings", "score_distribution", "cv_curve")) sb_abort(
    "Choose loadings, score_distribution or cv_curve.", "stormbreaker_unsupported_diagnostic")
  if (type == "cv_curve") return(sb_pls_cv_diagnostic(fit))
  model <- fit$diagnostic_artifacts$learner
  if (is.null(model)) sb_abort("Original PLS components were not captured; enable capture_diagnostics at final fitting.",
    "stormbreaker_missing_diagnostic_artifacts")
  native <- if (id %in% c("plsrcox", "plsrcox_cv")) list(
    loadings = model$engine$pp, scores = model$engine$tt,
    loading_definition = "Native plsRcox X reconstruction loadings") else model$diagnostic_components
  if (is.null(native)) sb_abort("This historical model did not retain native PLS components.",
    "stormbreaker_missing_diagnostic_artifacts")
  mat <- if (type == "loadings") native$loadings else native$scores
  expected <- if (type == "loadings") length(model$features) else nrow(fit$diagnostic_artifacts$predictors)
  if (!is.matrix(mat) || !is.numeric(mat) || nrow(mat) != expected || ncol(mat) < 1L ||
      any(!is.finite(mat))) sb_abort("Native PLS component matrix is invalid.", "stormbreaker_invalid_diagnostic")
  if (type == "loadings" && !identical(rownames(mat), model$features)) sb_abort(
    "Native PLS loading order does not match the fitted features.", "stormbreaker_invalid_diagnostic")
  tab <- if (type == "loadings") data.frame(feature = rep(model$features, times = ncol(mat)),
    component = rep(seq_len(ncol(mat)), each = nrow(mat)), loading = as.vector(mat)) else
    data.frame(observation = rep(seq_len(nrow(mat)), times = ncol(mat)),
      component = rep(seq_len(ncol(mat)), each = nrow(mat)), score = as.vector(mat))
  structure(list(type = paste0("pls_", type), table = tab,
    native_table = if (type == "loadings") mat else unname(mat),
    learner = id, pipeline_id = fit$pipeline$id, features = model$features,
    loading_definition = native$loading_definition,
    sensitive_training_data = type == "score_distribution",
    interpretation = "Original development-fit components; native signs and component order retained.\nNot feature importance, hazard ratios, probability effects or evidence of held-out class separation.",
    source = "original captured PLS components; no fitting or decomposition"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_pls_diagnostic <- function(x) {
  if (x$type == "pls_cv_curve") return(sb_plot_pls_cv_diagnostic(x))
  feature <- component <- loading <- score <- NULL
  tab <- x$table
  tab$component <- factor(tab$component, levels = sort(unique(tab$component)))
  if (x$type == "pls_loadings") {
    tab$feature <- factor(tab$feature, levels = rev(x$features))
    extent <- max(abs(tab$loading))
    if (extent == 0) extent <- 1
    p <- ggplot2::ggplot(tab, ggplot2::aes(component, feature, fill = loading)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#D77935",
        midpoint = 0, limits = c(-extent, extent),
        guide = ggplot2::guide_colourbar(display = "rectangles")) +
      ggplot2::labs(title = x$loading_definition, x = "Original component", y = NULL, fill = "Loading")
  } else {
    p <- ggplot2::ggplot(tab, ggplot2::aes(component, score)) +
      ggplot2::geom_boxplot(fill = "#B5D3E1", outlier.shape = NA, width = .45) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = .12, height = 0, seed = 1L),
        colour = "#2878A0", size = .8, alpha = .5) +
      ggplot2::labs(title = "Development-sample component scores",
        x = "Original component", y = "Native component score")
  }
  caption <- x$interpretation
  if (x$type == "pls_score_distribution") caption <- paste(caption,
    "Box: median and interquartile range; whiskers: observations within 1.5 IQR.\nAll observations shown; boxes and whiskers are not confidence intervals.", sep = "\n")
  p + ggplot2::labs(subtitle = x$learner, caption = caption) + sb_viz_theme()
}
