#' Inspect saved shrinkage and subspace discriminant models
#' @param fit Successful SDA or SLDA final fit.
#' @param type coefficients or regularization for SDA; loadings or
#'   class_means for SLDA. NULL chooses the first applicable view.
#' @return Native shrinkage parameters, class-1 log-score contrast, subspace
#'   loadings, or class means in that subspace. SDA's native posterior method
#'   rounds probabilities; coefficients describe its pre-rounding log odds.
#'   SLDA axes are fitted subspace coordinates, not original features.
#'   No covariance reconstruction, feature selection or refitting occurs.
#' @export
sb_shrinkage_discriminant_diagnostics <- function(fit, type = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  is_sda <- identical(s$learner, "sda_binary")
  if (!s$learner %in% c("sda_binary", "slda_binary")) sb_abort(
    "This diagnostic requires SDA or SLDA.", "stormbreaker_unsupported_diagnostic")
  choices <- if (is_sda) c("coefficients", "regularization") else c("loadings", "class_means")
  if (is.null(type)) type <- choices[1L]
  if (!is.character(type) || length(type) != 1L || is.na(type) || !type %in% choices) sb_abort(
    paste("Choose", paste(choices, collapse = " or ")), "stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  intercept <- NULL
  if (is_sda) {
    if (!inherits(e, "sda") || !identical(colnames(e$beta), s$features) ||
        !identical(rownames(e$beta), c("0", "1")) || any(!is.finite(e$beta))) sb_abort(
      "Original SDA coefficients are invalid.", "stormbreaker_invalid_diagnostic")
    if (type == "coefficients") {
      tab <- data.frame(feature = s$features, coefficient = as.numeric(e$beta["1", ] - e$beta["0", ]))
      intercept <- unname(e$alpha["1"] - e$alpha["0"])
    } else {
      tab <- data.frame(parameter = names(e$regularization), value = as.numeric(e$regularization))
    }
  } else {
    if (!inherits(e, "slda") || !is.matrix(e$scores) || nrow(e$scores) != length(s$features) ||
        any(!is.finite(e$scores)) || !identical(e$mylda$lev, c("0", "1"))) sb_abort(
      "Original SLDA subspace is invalid.", "stormbreaker_invalid_diagnostic")
    if (type == "loadings") {
      tab <- data.frame(feature = rep(s$features, ncol(e$scores)),
        axis = rep(seq_len(ncol(e$scores)), each = nrow(e$scores)), loading = as.numeric(e$scores))
    } else {
      m <- e$mylda$means
      tab <- data.frame(class = rep(rownames(m), ncol(m)),
        axis = rep(seq_len(ncol(m)), each = nrow(m)), mean = as.numeric(m))
    }
  }
  structure(list(type = paste0("shrinkage_discriminant_", type), table = tab,
    intercept = intercept, learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "Original SDA alpha/beta/shrinkage or SLDA scores/mylda",
    interpretation = "SDA coefficients describe pre-rounding class-1 log odds; native probabilities remain unchanged. SLDA loading signs are basis conventions, not importance or clinical effects. No refitting."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_shrinkage_discriminant_diagnostic <- function(x) {
  feature <- coefficient <- parameter <- value <- axis <- loading <- class <- mean <- NULL
  if (x$type == "shrinkage_discriminant_coefficients") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(coefficient, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey65", linewidth = .3) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::labs(title = "SDA class-1 coefficients", x = "Log-odds coefficient", y = NULL)
  } else if (x$type == "shrinkage_discriminant_regularization") {
    labels <- c(lambda = "Correlation", lambda.var = "Variance", lambda.freqs = "Class frequency")
    p <- ggplot2::ggplot(x$table, ggplot2::aes(value, parameter)) +
      ggplot2::geom_col(fill = "#2878A0", width = .6) +
      ggplot2::scale_y_discrete(labels = labels) + ggplot2::scale_x_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = "SDA shrinkage", x = "Shrinkage intensity", y = NULL)
  } else if (x$type == "shrinkage_discriminant_loadings") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(loading, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey65", linewidth = .3) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::facet_wrap(~axis, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "SLDA subspace loadings", x = "Loading", y = NULL)
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(mean, class)) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::facet_wrap(~axis, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "SLDA class centers", x = "Subspace coordinate", y = "Class")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 10, height = 4)
  p
}
