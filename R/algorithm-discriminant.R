sb_discriminant_parameters <- function(state) {
  id <- state$learner
  if (!id %in% c("lda_binary", "qda_binary", "rda_binary")) sb_abort(
    "This learner has no supported discriminant parameter adapter.", "stormbreaker_unsupported_diagnostic")
  if (sb_is_constant_sum_state(state)) return(sb_constant_sum_linear_parameters(state))
  if (!is.raw(state$engine_raw) || !length(state$engine_raw)) sb_abort(
    "Saved discriminant parameters are missing.", "stormbreaker_bad_learner_state")
  e <- unserialize(state$engine_raw)
  expected <- switch(id, lda_binary = "lda", qda_binary = "qda", rda_binary = "rda")
  features <- state$features
  p <- length(features)
  classes <- if (id == "rda_binary") e$classes else e$lev
  means <- if (id == "rda_binary") t(e$means) else e$means
  if (!inherits(e, expected) || !identical(classes, c("0", "1")) ||
      !is.matrix(means) || !identical(dim(means), c(2L, as.integer(p))) ||
      !identical(colnames(means), features) || any(!is.finite(means)) ||
      !identical(names(e$prior), classes) || any(!is.finite(e$prior)) || any(e$prior <= 0)) sb_abort(
    "Original discriminant class means or priors are invalid.", "stormbreaker_bad_learner_state")
  prior <- e$prior
  covariance <- NULL
  if (id == "lda_binary") {
    if (!is.matrix(e$scaling) || nrow(e$scaling) != p || !ncol(e$scaling) ||
        !identical(rownames(e$scaling), features) || any(!is.finite(e$scaling))) sb_abort(
      "Original LDA scaling is invalid.", "stormbreaker_bad_learner_state")
    center <- colSums(prior * means)
    projected <- sweep(means, 2L, center) %*% e$scaling
    beta <- drop(e$scaling %*% (projected[2L, ] - projected[1L, ]))
    intercept <- log(prior[2L] / prior[1L]) - sum(center * beta) -
      0.5 * (sum(projected[2L, ]^2) - sum(projected[1L, ]^2))
    quadratic <- matrix(0, p, p)
  } else {
    if (id == "qda_binary") {
      precision <- lapply(seq_len(2L), function(i) tcrossprod(matrix(e$scaling[, , i], p, p)))
      covariance <- lapply(precision, solve)
      logdet <- e$ldet
    } else {
      gamma <- e$regularization[1L]
      lambda <- e$regularization[2L]
      if (length(c(gamma, lambda)) != 2L || any(!is.finite(c(gamma, lambda))) ||
          any(c(gamma, lambda) < 0 | c(gamma, lambda) > 1)) sb_abort(
        "Original RDA regularization is invalid.", "stormbreaker_bad_learner_state")
      covariance <- lapply(seq_len(2L), function(i) {
        pooled <- (1 - lambda) * matrix(e$covariances[, , i], p, p) + lambda * e$covpooled
        (1 - gamma) * pooled + gamma * mean(diag(pooled)) * diag(p)
      })
      precision <- lapply(covariance, solve)
      logdet <- vapply(covariance, function(z) as.numeric(determinant(z, logarithm = TRUE)$modulus), 0)
    }
    quadratic <- 0.5 * (precision[[1L]] - precision[[2L]])
    beta <- drop(precision[[2L]] %*% means[2L, ] - precision[[1L]] %*% means[1L, ])
    intercept <- log(prior[2L] / prior[1L]) - 0.5 * (logdet[2L] - logdet[1L] +
      sum(means[2L, ] * (precision[[2L]] %*% means[2L, ])) -
      sum(means[1L, ] * (precision[[1L]] %*% means[1L, ])))
    covariance <- lapply(covariance, function(z) {
      dimnames(z) <- list(features, features)
      z
    })
    names(covariance) <- classes
  }
  if (any(!is.finite(c(beta, intercept, quadratic)))) sb_abort(
    "Original discriminant equation is not finite.", "stormbreaker_bad_learner_state")
  names(beta) <- features
  dimnames(quadratic) <- list(features, features)
  list(coefficient = beta, intercept = unname(intercept), quadratic = quadratic,
    means = means, priors = prior, covariance = covariance)
}

#' Inspect original LDA, QDA and regularized discriminant parameters
#' @param fit A successful LDA, QDA or fixed-regularization klaR RDA final fit.
#' @param type class_means or class_correlation (QDA/RDA only).
#' @return An algorithm diagnostic from compact saved model parameters, without
#'   capture, refitting, covariance estimation or surrogate models. Class means
#'   use original feature units; they are not confidence intervals. Correlations
#'   describe the fitted within-class covariance used for prediction. QDA uses
#'   its native precision factors; RDA applies both native regularizations.
#'   LDA's rank-reduced discriminant scaling cannot identify the full covariance,
#'   so this view is explicitly unavailable for LDA. Its original linear log-odds
#'   equation is available through [sb_model_coefficients()] and [sb_nomogram()].
#'   The result retains the exact equation as `parameters`: intercept plus
#'   x times coefficient plus x' quadratic x. Off-diagonal matrix entries occur
#'   twice in that quadratic form. QDA/RDA are not treated as additive models.
#'   Source and interpretation belong in a separate figure legend, not panels.
#' @export
sb_discriminant_diagnostics <- function(fit, type = "class_means") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("class_means", "class_correlation")) sb_abort(
    "Choose class_means or class_correlation.", "stormbreaker_unsupported_diagnostic")
  parameters <- sb_discriminant_parameters(fit$learner_state)
  features <- fit$learner_state$features
  if (type == "class_means") {
    tab <- data.frame(class = rep(c("0", "1"), length(features)),
      feature = rep(features, each = 2L), value = as.vector(parameters$means))
    interpretation <- "Original fitted class means in input units; not uncertainty intervals or variable importance."
  } else {
    if (is.null(parameters$covariance)) sb_abort(
      "Rank-reduced LDA scaling does not retain a recoverable full within-class covariance.",
      "stormbreaker_unavailable_diagnostic")
    tab <- do.call(rbind, lapply(c("0", "1"), function(class) {
      correlation <- stats::cov2cor(parameters$covariance[[class]])
      data.frame(class = class, feature = rep(features, length(features)),
        feature2 = rep(features, each = length(features)), value = as.vector(correlation))
    }))
    interpretation <- paste("Correlations implied by the original class-specific prediction covariance.",
      "RDA includes its original pooling and spherical shrinkage; not raw-data correlation, importance or SHAP interaction.")
  }
  structure(list(type = paste0("discriminant_", type), table = tab, features = features,
    parameters = parameters, learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    source = "original saved discriminant model parameters", interpretation = interpretation),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_discriminant_diagnostic <- function(x) {
  feature <- feature2 <- class <- value <- NULL
  tab <- x$table
  tab$feature <- factor(tab$feature, levels = x$features)
  if (x$type == "discriminant_class_means") {
    p <- ggplot2::ggplot(tab, ggplot2::aes(value, class, colour = class)) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::facet_wrap(~feature, scales = "free_x", ncol = min(3L, length(x$features))) +
      ggplot2::scale_colour_manual(values = c("0" = "#2878A0", "1" = "#D77935"), guide = "none") +
      ggplot2::labs(title = "Fitted class means", x = "Mean (original units)", y = "Encoded class") + sb_viz_theme()
    size <- c(width = 11, height = max(4, 2 + 2 * ceiling(length(x$features) / 3)))
  } else {
    tab$feature2 <- factor(tab$feature2, levels = rev(x$features))
    p <- ggplot2::ggplot(tab, ggplot2::aes(feature, feature2, fill = value)) +
      ggplot2::geom_tile() + ggplot2::facet_wrap(~class, nrow = 1L, labeller = ggplot2::label_both) +
      ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#D77935", limits = c(-1, 1),
        guide = ggplot2::guide_colourbar(display = "rectangles")) +
      ggplot2::coord_fixed() + ggplot2::labs(title = "Fitted within-class correlations", x = NULL, y = NULL,
        fill = "Correlation") + sb_viz_theme() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    size <- c(width = max(11, 4 + .6 * length(x$features)), height = max(5, 2 + .3 * length(x$features)))
  }
  attr(p, "stormbreaker_suggested_size") <- size
  p
}
