sb_multiclass_logistic_plots <- function(fit) {
  coefficients <- stats::coef(fit$engine)
  terms <- colnames(coefficients)
  at <- match(terms, fit$feature_codes)
  terms[!is.na(at)] <- fit$features[at[!is.na(at)]]
  classes <- fit$classes[match(rownames(coefficients), fit$codes)]
  d <- data.frame(feature = rep(terms, each = nrow(coefficients)),
    contrast = rep(paste(classes, "vs", fit$classes[1L]), ncol(coefficients)),
    coefficient = as.vector(coefficients))
  d$feature <- factor(d$feature, levels = rev(terms))
  coefficient <- feature <- NULL
  plot <- ggplot2::ggplot(d, ggplot2::aes(coefficient, feature)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
    ggplot2::geom_point(colour = "#2878A0") +
    ggplot2::facet_wrap(~contrast) + sb_viz_theme() +
    ggplot2::labs(title = "Multinomial coefficients", x = "Log-odds coefficient", y = NULL,
      caption = "Original fitted coefficients, including the intercept, relative to the reference class. No refitting or coefficient thresholding is performed.")
  result <- list(plots = list(coefficients = plot), tables = list(coefficients = d))
  if (fit$method == "stepglm") {
    path <- fit$engine$anova
    path$step <- seq_len(nrow(path)) - 1L
    step <- AIC <- NULL
    result$plots$selection_path <- ggplot2::ggplot(path, ggplot2::aes(step, AIC)) +
      ggplot2::geom_line(colour = "#2878A0") + ggplot2::geom_point(colour = "#2878A0") +
      sb_viz_theme() + ggplot2::labs(title = "Stepwise selection", x = "Step", y = "AIC",
        caption = "Actual accepted stepAIC steps, including the initial model; not a cross-validation performance curve.")
    result$tables$selection_path <- path
  }
  result
}
