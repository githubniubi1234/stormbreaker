#' Inspect original supervised principal-component diagnostics
#' @param fit A final SuperPC model fitted with capture_diagnostics = TRUE.
#' @param type feature_scores, loadings or cv_curve (CV variant only).
#' @return An algorithm diagnostic retaining native feature scores, selected
#'   singular-vector loadings, or the full component-by-threshold CV table.
#'   Feature scores are training associations, not model coefficients. Loadings
#'   retain the native SVD sign and are not hazard ratios or variable importance.
#'   CV scores are geometric means of held-out likelihood-ratio statistics;
#'   bounds use plus/minus one standard error on the log scale, not confidence
#'   intervals. The original pipeline selects thresholds using one component,
#'   even though native CV evaluates more components. No refitting is performed.
#'   Compare thresholds within a component count; likelihood-ratio statistics
#'   with different degrees of freedom are not a model ranking across counts.
#' @export
sb_superpc_diagnostics <- function(fit, type = "feature_scores") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  id <- fit$pipeline$learner
  if (!id %in% c("superpc", "superpc_cv")) sb_abort(
    "This diagnostic requires an original SuperPC learner.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("feature_scores", "loadings", "cv_curve")) sb_abort(
    "Choose feature_scores, loadings or cv_curve.", "stormbreaker_unsupported_diagnostic")
  model <- fit$diagnostic_artifacts$learner
  if (is.null(model)) sb_abort("Original SuperPC objects were not captured; enable capture_diagnostics on final fitting.",
    "stormbreaker_missing_diagnostic_artifacts")
  engine <- model$engine
  features <- model$features
  selected <- engine$prediction$which.features
  threshold <- if (id == "superpc_cv") model$audit$selected_threshold else model$params$threshold
  if (!is.logical(selected) || length(selected) != length(features) || anyNA(selected)) sb_abort(
    "Original SuperPC feature selection is invalid.", "stormbreaker_invalid_diagnostic")
  folds <- NULL
  if (type == "feature_scores") {
    native <- engine$trained$feature.scores
    if (!is.numeric(native) || length(native) != length(features) || any(!is.finite(native)) ||
        !identical(names(native), features) || !identical(unname(abs(native) >= threshold), unname(selected))) sb_abort(
      "Native SuperPC scores disagree with the original selection.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(feature = features, score = as.numeric(native), selected = unname(selected))
    measure <- "Original SuperPC feature-screening scores"
    interpretation <- "Signed native survival-screening associations, not fitted coefficients or causal effects.\nDashed lines: original absolute-score threshold. All candidate features retained."
  } else if (type == "loadings") {
    native <- engine$prediction$u
    mat <- if (is.null(dim(native))) matrix(native, ncol = 1L) else native
    if (!is.matrix(mat) || nrow(mat) != sum(selected) ||
        ncol(mat) != engine$prediction$n.components || any(!is.finite(mat))) sb_abort(
      "Original SuperPC singular vectors are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(feature = rep(features[selected], times = ncol(mat)),
      component = rep(seq_len(ncol(mat)), each = nrow(mat)), loading = as.vector(mat))
    measure <- "Original SuperPC singular-vector loadings"
    interpretation <- "Selected features only; native SVD signs retained. Not coefficients, hazard ratios or importance.\nExcluded input features remain in the feature-screening table, not fabricated loading rows."
  } else {
    if (id != "superpc_cv") sb_abort("Fixed SuperPC has no CV threshold search.",
      "stormbreaker_unavailable_diagnostic")
    cv <- engine$cross_validation
    native <- list(thresholds = cv$thresholds, scores = cv$scor, lower = cv$scor.lower, upper = cv$scor.upper)
    if (!is.numeric(cv$thresholds) || !length(cv$thresholds) || any(!is.finite(cv$thresholds)) ||
        !is.matrix(cv$scor) || !nrow(cv$scor) || ncol(cv$scor) != length(cv$thresholds) ||
        !identical(dim(cv$scor), dim(cv$scor.lower)) || !identical(dim(cv$scor), dim(cv$scor.upper))) sb_abort(
      "Original SuperPC CV dimensions are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(threshold = rep(cv$thresholds, each = nrow(cv$scor)),
      component = rep(seq_len(nrow(cv$scor)), times = ncol(cv$scor)),
      score = as.vector(cv$scor), lower = as.vector(cv$scor.lower), upper = as.vector(cv$scor.upper))
    tab$selected <- tab$component == model$params$prediction_components &
      rep(seq_along(cv$thresholds) == model$audit$selected_threshold_index, each = nrow(cv$scor))
    if (sum(tab$selected) != 1L || !is.finite(tab$score[tab$selected]) ||
        !identical(as.numeric(tab$threshold[tab$selected]), as.numeric(threshold)) ||
        any(vapply(tab[c("score", "lower", "upper")], function(z) any(is.infinite(z)), logical(1)))) sb_abort(
      "Original SuperPC CV scores or selected threshold are invalid.", "stormbreaker_invalid_diagnostic")
    folds <- c(requested = model$audit$requested_folds, effective = cv$n.fold)
    measure <- "Original SuperPC threshold cross-validation"
    interpretation <- paste0("Geometric mean held-out likelihood-ratio statistic; bands are +/- one log-scale SE, not CI.\n",
      "Compare thresholds within a component count, not across counts with different degrees of freedom.\n",
      "Original one-component threshold selection retained; folds requested/effective: ", folds[1], "/", folds[2], ".")
  }
  structure(list(type = paste0("superpc_", type), table = tab, native_table = native,
    learner = id, pipeline_id = fit$pipeline$id, threshold = threshold, folds = folds,
    measure = measure, interpretation = interpretation, source = "original captured SuperPC fit; no refitting"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_superpc_diagnostic <- function(x) {
  feature <- score <- selected <- loading <- component <- threshold <- lower <- upper <- NULL
  tab <- x$table
  if (x$type == "superpc_feature_scores") {
    tab$feature <- factor(tab$feature, levels = tab$feature[order(tab$score)])
    p <- ggplot2::ggplot(tab, ggplot2::aes(score, feature, fill = selected)) +
      ggplot2::geom_col(width = .65) +
      ggplot2::geom_vline(xintercept = c(-x$threshold, x$threshold), linetype = 2, colour = "#D77935") +
      ggplot2::scale_fill_manual(values = c("FALSE" = "#B8B8B8", "TRUE" = "#2878A0")) +
      ggplot2::labs(x = "Native screening score", y = NULL, fill = "Selected")
  } else if (x$type == "superpc_loadings") {
    tab$component <- factor(tab$component)
    p <- ggplot2::ggplot(tab, ggplot2::aes(loading, feature, fill = component)) +
      ggplot2::geom_col(position = "dodge", width = .65) +
      ggplot2::labs(x = "Native singular-vector loading", y = NULL, fill = "Component")
  } else {
    selected_component <- as.character(tab$component[tab$selected])
    tab$component <- factor(tab$component)
    tab$panel <- ifelse(as.character(tab$component) == selected_component,
      paste("Component", tab$component, "(used)"), paste("Component", tab$component))
    p <- ggplot2::ggplot(tab, ggplot2::aes(threshold, score)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "#2878A0", alpha = .12, colour = NA, na.rm = TRUE) +
      ggplot2::geom_line(colour = "#2878A0", na.rm = TRUE) +
      ggplot2::geom_point(data = tab[tab$selected, ], shape = 21, colour = "black", size = 3) +
      ggplot2::geom_vline(data = tab[tab$selected, ], ggplot2::aes(xintercept = threshold), colour = "#D77935", linetype = 2) +
      ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
      ggplot2::labs(x = "Feature-screening threshold", y = "CV likelihood-ratio statistic")
  }
  p + ggplot2::labs(title = x$measure, subtitle = x$learner, caption = x$interpretation) + sb_viz_theme()
}
