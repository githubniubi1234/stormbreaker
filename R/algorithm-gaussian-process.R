#' Inspect native Gaussian-process classifier kernel terms
#' @param fit Successful Gaussian-process binary final fit.
#' @param type basis_weights or query_contributions.
#' @param newdata Explicit queries for query_contributions.
#' @return Original basis weights, or per-query kernel terms and native
#'   probabilities. Terms are oriented toward class 1 and sum to its native
#'   latent score; they are reference-row contributions, not feature SHAP.
#'   No posterior interval is generated: this classifier does not save the
#'   regression posterior covariance required by the native variance method.
#'   Reference indices can be sensitive. No refitting or calibration occurs.
#' @export
sb_gaussian_process_diagnostics <- function(fit, type = "basis_weights", newdata = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (s$learner != "gaussian_process_binary" || !is.character(type) ||
      length(type) != 1L || is.na(type) || !type %in% c("basis_weights", "query_contributions")) {
    sb_abort("Choose basis_weights or query_contributions for a Gaussian-process classifier.",
      "stormbreaker_unsupported_diagnostic")
  }
  if (!requireNamespace("kernlab", quietly = TRUE)) sb_abort(
    "Gaussian-process diagnostics require kernlab.", "stormbreaker_unavailable_diagnostic")
  e <- unserialize(s$engine_raw)
  if (!methods::is(e, "gausspr") || e@type != "classification" ||
      !identical(e@lev, c("0", "1")) || length(e@alpha) != 1L ||
      !identical(colnames(e@xmatrix), s$features)) {
    sb_abort("Original binary Gaussian-process model is invalid.", "stormbreaker_invalid_diagnostic")
  }
  indices <- e@alphaindex[[1L]]
  weights <- -as.numeric(e@alpha[[1L]])
  if (length(indices) != length(weights) || any(!is.finite(weights))) {
    sb_abort("Original Gaussian-process basis weights are invalid.", "stormbreaker_invalid_diagnostic")
  }
  summary <- NULL
  if (type == "basis_weights") {
    tab <- data.frame(basis = seq_along(weights), reference_row = indices, weight = weights)
  } else {
    if (is.null(newdata)) sb_abort("Supply explicit queries for kernel contributions.", "stormbreaker_unavailable_diagnostic")
    x <- as.matrix(sb_select_final_feature_frame(newdata, s$features, "Query data", "query fitted features"))
    if (!nrow(x)) sb_abort("Queries must contain at least one row.", "stormbreaker_invalid_diagnostic")
    native <- kernlab::predict(e, x, type = "probabilities")
    if (is.list(e@scaling)) {
      selected <- e@scaling$scaled
      x[, selected] <- scale(x[, selected, drop = FALSE],
        center = e@scaling$x.scale[["scaled:center"]], scale = e@scaling$x.scale[["scaled:scale"]])
    }
    k <- as.matrix(kernlab::kernelMatrix(e@kernelf, x, e@xmatrix[indices, , drop = FALSE]))
    terms <- sweep(k, 2L, weights, "*")
    tab <- data.frame(query = rep(seq_len(nrow(x)), each = length(weights)),
      basis = rep(seq_along(weights), nrow(x)), reference_row = rep(indices, nrow(x)),
      similarity = as.numeric(t(k)), contribution = as.numeric(t(terms)))
    summary <- data.frame(query = seq_len(nrow(x)), latent_score = rowSums(terms),
      probability = as.numeric(native[, "1"]))
  }
  structure(list(type = paste0("gaussian_process_", type), table = tab, summary = summary,
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "Original kernlab gausspr basis coefficients and saved kernel",
    interpretation = "Class-1-oriented reference-row kernel terms, not feature effects or SHAP. No uncertainty band, refit or calibration. Reference indices may be sensitive."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_gaussian_process_diagnostic <- function(x) {
  basis <- weight <- contribution <- NULL
  if (x$type == "gaussian_process_basis_weights") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(basis, weight)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey65", linewidth = .3) +
      ggplot2::geom_point(colour = "#2878A0", size = 1.3) +
      ggplot2::labs(title = "Gaussian-process basis weights", x = "Native basis", y = "Class-1 weight")
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(basis, contribution)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey65", linewidth = .3) +
      ggplot2::geom_point(colour = "#2878A0", size = 1.3) +
      ggplot2::facet_wrap(~query, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "Kernel contributions", x = "Native basis", y = "Class-1 latent contribution")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}
