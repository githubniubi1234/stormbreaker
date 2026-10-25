sb_multiclass_pre_plots <- function(fit) {
  m <- fit$engine
  native <- stats::coef(m, penalty.par.val = "lambda.1se")
  tab <- do.call(rbind, lapply(seq_along(fit$codes), function(k) {
    data.frame(feature = native$rule, description = native$description,
      class = fit$classes[k], coefficient = native[[paste0("coefficient.", fit$codes[k])]])
  }))
  for (i in order(nchar(fit$feature_codes), decreasing = TRUE))
    tab$description <- gsub(fit$feature_codes[i], fit$features[i], tab$description, fixed = TRUE)
  tab$term_type <- ifelse(tab$feature == "(Intercept)", "intercept",
    ifelse(tab$feature %in% m$rules$rule, "rule", "linear"))
  cv <- m$glmnet.fit
  tuning <- data.frame(lambda = cv$lambda, mean_value = cv$cvm, standard_error = cv$cvsd)
  plots <- list(coefficient_cv = sb_plot_pre_diagnostic(list(type = "rule_cv_curve",
    table = tuning, selected_lambda = cv$lambda.1se, selection_rule = "lambda.1se", measure = cv$name)))
  coefficient <- feature <- term_type <- NULL
  for (cl in fit$classes) {
    d <- tab[tab$class == cl & tab$coefficient != 0 & tab$term_type != "intercept", ]
    if (!nrow(d)) next
    d$feature <- factor(d$feature, levels = d$feature[order(d$coefficient)])
    plots[[paste0("nonzero_terms_", cl)]] <- ggplot2::ggplot(d, ggplot2::aes(coefficient, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey65") +
      ggplot2::geom_col(ggplot2::aes(fill = term_type)) +
      ggplot2::scale_fill_manual(values = c(rule = "#2878A0", linear = "#D77935"), name = "Term") +
      sb_viz_theme() + ggplot2::labs(title = paste(cl, "nonzero terms"),
        x = "Class-score coefficient", y = "Native term ID",
        caption = "Nonzero non-intercept terms at the saved lambda.1se. Complete class coefficients, intercepts, zero terms and rule definitions remain in the diagnostic table. Native multinomial scores use softmax, not one-reference-class log odds. Rule and linear terms have different units; coefficient magnitude is not feature importance. No refitting.")
  }
  list(plots = plots, tables = list(coefficients = tab, coefficient_cv = tuning))
}
