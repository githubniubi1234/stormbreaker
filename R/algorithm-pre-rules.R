#' Inspect an original PRE rule ensemble and its coefficient tuning
#' @param fit A final binary or survival PRE model. Its saved prediction state
#'   already contains the native engine; extra diagnostic capture is not needed.
#' @param type Either rule_coefficients or rule_cv_curve.
#' @return An algorithm diagnostic. For coefficients, table contains every
#'   original term, its exact native description, coefficient and term type,
#'   including zero coefficients and the intercept where present. Linear terms
#'   retain native winsorization bounds and coefficient rescaling. Their native
#'   interval description denotes clipping, not an interval-indicator rule.
#'   Coefficient size is
#'   not comparable feature importance. The plot excludes zero coefficients and
#'   the intercept, but the complete table remains unchanged.
#'
#'   For tuning, table contains the saved coefficient-only cross-validation
#'   sequence, with its original standard errors and selected lambda. Rules were
#'   induced on the complete development data before this internal validation;
#'   this is not validation of the complete rule-learning pipeline. No model is
#'   refitted. Rule descriptions are returned as text, not evaluated.
#' @export
sb_rule_ensemble <- function(fit, type = "rule_coefficients") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!fit$pipeline$learner %in% sb_pre_rule_ensemble_learner_ids()) {
    sb_abort("This diagnostic requires an original PRE rule ensemble.",
      "stormbreaker_unsupported_diagnostic")
  }
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("rule_coefficients", "rule_cv_curve")) {
    sb_abort("Choose rule_coefficients or rule_cv_curve.",
      "stormbreaker_unsupported_diagnostic")
  }
  sb_require_rule_aft_dependency(fit$pipeline$learner)
  state <- fit$learner_state
  engine <- unserialize(state$engine_raw)
  if (!inherits(engine, "pre") || !inherits(engine$glmnet.fit, "cv.glmnet")) {
    sb_abort("The saved PRE engine is incomplete.", "stormbreaker_invalid_diagnostic")
  }
  rule <- state$params$penalty_rule
  if (!identical(rule, "lambda.1se")) {
    sb_abort("The saved PRE selection rule is unsupported.", "stormbreaker_invalid_diagnostic")
  }
  cv <- engine$glmnet.fit
  selected <- cv[[rule]]
  if (!is.numeric(selected) || length(selected) != 1L || !is.finite(selected) ||
      selected <= 0 || !selected %in% cv$lambda) {
    sb_abort("The original selected PRE lambda is invalid.", "stormbreaker_invalid_diagnostic")
  }
  if (type == "rule_coefficients") {
    tab <- stats::coef(engine, penalty.par.val = rule)
    if (!is.data.frame(tab) || !all(c("rule", "coefficient", "description") %in% names(tab)) ||
        anyNA(tab) || !is.numeric(tab$coefficient) || any(!is.finite(tab$coefficient)) ||
        anyDuplicated(tab$rule)) {
      sb_abort("Original PRE term definitions or coefficients are invalid.",
        "stormbreaker_invalid_diagnostic")
    }
    tab$term_type <- ifelse(tab$rule == "(Intercept)", "intercept",
      ifelse(tab$rule %in% engine$rules$rule, "rule", "linear"))
    bound_index <- match(tab$rule, engine$wins_points$varname)
    tab$winsor_lower <- engine$wins_points$lb[bound_index]
    tab$winsor_upper <- engine$wins_points$ub[bound_index]
    tab$selected <- tab$coefficient != 0
    rownames(tab) <- NULL
  } else {
    values <- list(lambda = cv$lambda, mean_value = cv$cvm,
      standard_error = cv$cvsd, nonzero = cv$nzero)
    if (!length(cv$lambda) || !all(vapply(values, function(z) {
      is.numeric(z) && length(z) == length(cv$lambda) && all(is.finite(z))
    }, logical(1))) || any(cv$lambda <= 0) || anyDuplicated(cv$lambda) ||
        any(cv$cvsd < 0) || any(cv$nzero < 0) || any(cv$nzero != floor(cv$nzero))) {
      sb_abort("The saved PRE coefficient-validation path is invalid.",
        "stormbreaker_invalid_diagnostic")
    }
    tab <- as.data.frame(values)
    tab$selected <- tab$lambda == selected
  }
  structure(list(type = type, table = tab, learner = fit$pipeline$learner,
    pipeline_id = fit$pipeline$id, selected_lambda = selected, selection_rule = rule,
    measure = cv$name, scale = if (engine$family == "binomial") "Log odds" else "Log hazard",
    source = "original saved PRE engine; no refitting",
    validation_scope = "coefficient-only CV after rule induction on complete development data"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_pre_diagnostic <- function(x) {
  if (x$type == "rule_cv_curve") {
    lambda <- mean_value <- standard_error <- NULL
    return(ggplot2::ggplot(x$table, ggplot2::aes(lambda, mean_value)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_value - standard_error,
        ymax = mean_value + standard_error), width = 0, colour = "#A6A6A6") +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::geom_point(size = 1, colour = "#2878A0") +
      ggplot2::geom_vline(xintercept = x$selected_lambda, linetype = 2, colour = "#D77935") +
      ggplot2::scale_x_log10() +
      ggplot2::labs(title = "PRE coefficient-only validation",
        subtitle = paste("Original selection:", x$selection_rule),
        x = "Lambda (log scale)", y = x$measure,
        caption = "Rules learned before internal coefficient CV; not validation of complete rule learning.\nBars: mean +/- one standard error, not a confidence interval. Dashed line: selected lambda.") +
      sb_viz_theme())
  }
  tab <- x$table[x$table$selected & x$table$term_type != "intercept", , drop = FALSE]
  tab$rule <- factor(tab$rule, levels = rev(tab$rule[order(abs(tab$coefficient), decreasing = TRUE)]))
  coefficient <- rule <- term_type <- NULL
  intercept <- x$table$coefficient[x$table$term_type == "intercept"]
  subtitle <- if (length(intercept)) paste("Original intercept:", format(signif(intercept, 4))) else
    "Cox ensemble: no intercept"
  p <- ggplot2::ggplot(tab, ggplot2::aes(coefficient, rule)) +
    ggplot2::labs(title = "PRE nonzero terms", subtitle = subtitle,
      x = paste(x$scale, "coefficient (native term units)"), y = "Native term ID",
      caption = "Complete term definitions and zero coefficients are retained in the diagnostic table.\nRule coefficients and linear slopes have different units; this is not feature importance.") +
    sb_viz_theme()
  if (nrow(tab)) {
    p <- p + ggplot2::geom_vline(xintercept = 0, colour = "#9A9A9A") +
      ggplot2::geom_col(ggplot2::aes(fill = term_type), width = .65) +
      ggplot2::scale_fill_manual(values = c(rule = "#2878A0", linear = "#D77935"), name = "Term")
  } else {
    p <- p + ggplot2::annotate("text", x = 0, y = 1,
      label = "No nonzero rule or linear terms at the original selected lambda.", size = 3.2) +
      ggplot2::scale_x_continuous(limits = c(-1, 1)) +
      ggplot2::scale_y_continuous(limits = c(0, 2)) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
  }
  attr(p, "stormbreaker_suggested_size") <- c(width = 10, height = max(5, 2.5 + .28 * nrow(tab)))
  p
}
