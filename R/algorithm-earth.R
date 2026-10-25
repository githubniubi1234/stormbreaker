#' Inspect the original MARS logistic basis and pruning criterion
#' @param fit A successful earth binary final fit.
#' @param type basis_coefficients or pruning_gcv.
#' @return An algorithm diagnostic retaining original selected logistic basis
#'   coefficients, or the native basis-pruning GCV sequence. Logistic parameters
#'   come from glm.coefficients, not the separate least-squares coefficients.
#'   Native basis labels are retained for display; full-precision cuts and
#'   directions remain available separately and must not be reconstructed from
#'   rounded labels. The GCV sequence needs original captured artifacts and is
#'   the basis-selection criterion, not held-out binomial deviance, calibration
#'   or the final logistic model's cross-validation performance. Term counts
#'   include the intercept. No refitting, surrogate model or uncertainty is added.
#' @export
sb_earth_diagnostics <- function(fit, type = "basis_coefficients") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (fit$pipeline$learner != "earth_binary") sb_abort(
    "This diagnostic requires the earth binary learner.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("basis_coefficients", "pruning_gcv")) sb_abort(
    "Choose basis_coefficients or pruning_gcv.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(fit$learner_state$engine_raw)
  if (type == "basis_coefficients") {
    beta <- e$glm.coefficients
    if (!is.matrix(beta) || ncol(beta) != 1L || nrow(beta) != length(e$selected.terms) ||
        is.null(rownames(beta)) || any(!is.finite(beta))) sb_abort(
      "Original MARS logistic basis coefficients are invalid.", "stormbreaker_bad_learner_state")
    tab <- data.frame(term = rownames(beta), coefficient = as.numeric(beta),
      basis_index = e$selected.terms, intercept = e$selected.terms == 1L)
    interpretation <- paste("Original selected logistic basis coefficients, including the intercept; not the least-squares coefficients.",
      "Basis labels may round knot locations; exact cuts/directions are retained separately.",
      "Not raw-feature effects, SHAP or uncertainty.")
  } else {
    native <- fit$diagnostic_artifacts$learner$engine
    if (is.null(native$gcv.per.subset)) sb_abort("Original MARS pruning diagnostics were not captured.",
      "stormbreaker_missing_diagnostic_artifacts")
    gcv <- native$gcv.per.subset
    if (!is.numeric(gcv) || !length(gcv) || !any(is.finite(gcv))) sb_abort(
      "The native pruning sequence has no finite criterion values.", "stormbreaker_unavailable_diagnostic")
    tab <- data.frame(terms = seq_along(gcv), gcv = as.numeric(gcv),
      selected = seq_along(gcv) == length(e$selected.terms))
    interpretation <- paste("Native basis-pruning GCV; smaller is better. Term count includes the intercept.",
      "Orange point marks the original selected basis size, not a new optimum.",
      "This is not final logistic-model held-out deviance or a confidence interval; nonfinite criteria remain in the table and leave plot gaps.")
  }
  structure(list(type = paste0("earth_", type), table = tab,
    directions = e$dirs, cuts = e$cuts, selected_terms = e$selected.terms,
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    source = "original saved earth logistic basis and native pruning criterion", interpretation = interpretation),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_additive_shap_earth_terms <- function(state, x) {
  e <- unserialize(state$engine_raw)
  directions <- e$dirs[e$selected.terms, , drop = FALSE]
  if (any(rowSums(directions != 0) > 1L)) sb_abort_unsupported_additive_shap(state$learner)
  basis <- get("model.matrix.earth", asNamespace("earth"))(e, x)
  beta <- as.numeric(e$glm.coefficients)
  contributions <- matrix(0, nrow(x), ncol(x), dimnames = list(NULL, state$features))
  constant <- 0
  for (j in seq_len(ncol(basis))) {
    feature <- which(directions[j, ] != 0)
    if (length(feature)) contributions[, feature] <- contributions[, feature] + basis[, j] * beta[j] else
      constant <- constant + beta[j]
  }
  list(constant = constant, terms = contributions, score_scale = "class_1_log_odds")
}

sb_plot_earth_diagnostic <- function(x) {
  term <- coefficient <- terms <- gcv <- NULL
  tab <- x$table
  if (x$type == "earth_basis_coefficients") {
    tab$term <- factor(tab$term, levels = rev(tab$term))
    p <- ggplot2::ggplot(tab, ggplot2::aes(coefficient, term)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = .4) +
      ggplot2::geom_col(fill = "#2878A0", width = .65) +
      ggplot2::labs(title = "MARS logistic basis coefficients", x = "Logistic coefficient per basis unit", y = NULL)
  } else {
    tab$gcv[!is.finite(tab$gcv)] <- NA_real_
    p <- ggplot2::ggplot(tab, ggplot2::aes(terms, gcv)) +
      ggplot2::geom_line(colour = "#2878A0", na.rm = TRUE) +
      ggplot2::geom_point(colour = "#2878A0", na.rm = TRUE) +
      ggplot2::geom_point(data = tab[tab$selected, , drop = FALSE], colour = "#D77935", size = 2.5, na.rm = TRUE) +
      ggplot2::labs(title = "MARS basis pruning", x = "Number of basis terms", y = "GCV")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11,
    height = if (x$type == "earth_basis_coefficients") max(5, 2 + .3 * nrow(tab)) else 5)
  p
}
