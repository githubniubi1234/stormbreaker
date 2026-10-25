sb_glmnet_tuning_diagnostic <- function(fit, type) {
  if (!fit$pipeline$learner %in% c("glmnet_logistic", "glmnet_cox")) sb_abort(
    "This diagnostic requires a saved glmnet inner-validation curve.",
    "stormbreaker_unsupported_diagnostic")
  tuning <- fit$tuning
  tab <- tuning$curve
  value_col <- if (identical(tuning$measure, "deviance")) "mean_deviance" else "mean_measure"
  required <- c("lambda", value_col, "standard_error", "nonzero")
  valid <- identical(tuning$status, "success") && is.data.frame(tab) && nrow(tab) > 0 &&
    all(required %in% names(tab)) && all(vapply(tab[required], function(z) {
      is.numeric(z) && all(is.finite(z))
    }, logical(1))) && all(tab$lambda > 0) && !anyDuplicated(tab$lambda) &&
    all(tab$standard_error >= 0) && all(tab$nonzero >= 0) &&
    all(tab$nonzero == floor(tab$nonzero)) &&
    is.numeric(tuning$lambda) && length(tuning$lambda) == 1L &&
    is.finite(tuning$lambda) && tuning$lambda %in% tab$lambda
  if (!valid) sb_abort("Saved glmnet tuning curve is missing or invalid.",
    "stormbreaker_invalid_diagnostic")
  tab$mean_value <- tab[[value_col]]
  tab$selected <- tab$lambda == tuning$lambda
  labels <- c(deviance = "Deviance", auc = "ROC AUC", C = "Concordance",
    class = "Classification error", mse = "Mean squared error", mae = "Mean absolute error")
  measure <- unname(labels[tuning$measure])
  if (length(measure) != 1L || is.na(measure)) measure <- tuning$measure
  structure(list(type = type, table = tab, learner = fit$pipeline$learner,
    pipeline_id = fit$pipeline$id, measure = measure,
    selected_lambda = tuning$lambda, selection_rule = tuning$selection_rule,
    alpha = tuning$alpha, source = "saved final-development inner cross-validation",
    uncertainty = "one saved standard error; not a confidence interval"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_glmnet_tuning_diagnostic <- function(x) {
  lambda <- mean_value <- standard_error <- nonzero <- NULL
  p <- ggplot2::ggplot(x$table, ggplot2::aes(lambda))
  if (x$type == "cv_curve") {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_value - standard_error,
      ymax = mean_value + standard_error), width = 0, colour = "#A6A6A6") +
      ggplot2::geom_line(ggplot2::aes(y = mean_value), colour = "#2878A0") +
      ggplot2::geom_point(ggplot2::aes(y = mean_value), size = 1, colour = "#2878A0") +
      ggplot2::labs(title = "CV error", y = x$measure,
        caption = "Bars: saved mean +/- one standard error, not a confidence interval.\nDevelopment inner validation, not untouched-cohort performance; no refitting.")
  } else {
    p <- p + ggplot2::geom_step(ggplot2::aes(y = nonzero), colour = "#2878A0") +
      ggplot2::geom_point(ggplot2::aes(y = nonzero), size = 1, colour = "#2878A0") +
      ggplot2::scale_y_continuous(limits = c(0, max(1, x$table$nonzero)),
        breaks = function(limits) unique(round(pretty(limits, n = 5)))) +
      ggplot2::labs(title = "Model size",
        y = "Nonzero predictors",
        caption = "Counts retained from the original glmnet tuning object.\nThis is model complexity, not validation performance; no refitting.")
  }
  p + ggplot2::geom_vline(xintercept = x$selected_lambda, linetype = 2, colour = "#D77935") +
    ggplot2::scale_x_log10() + ggplot2::labs(x = "Lambda (log scale)",
      subtitle = paste0("Original rule: ", x$selection_rule, "; alpha = ", x$alpha,
        "; selected lambda = ", format(signif(x$selected_lambda, 4)))) + sb_viz_theme()
}
