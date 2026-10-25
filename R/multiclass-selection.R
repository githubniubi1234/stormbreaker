sb_multiclass_selector_inputs <- function(x, y) {
  x <- sb_validate_predictor_frame(x, "x")
  if (!is.factor(y) || nlevels(y) < 3L || length(y) != nrow(x) ||
      anyNA(y) || any(table(y) == 0L)) sb_abort(
    "Supply a factor with at least three observed classes.", "stormbreaker_bad_multiclass_outcome")
  list(x = x, endpoint = list(outcome = y))
}

sb_analyze_multiclass_selector <- function(task, spec) {
  x <- task$data[task$predictors]
  y <- factor(task$data[[task$outcome_col]], levels = task$classes)
  method <- spec$method
  fit <- switch(method,
    boruta = sb_fit_boruta_selector(x, list(outcome = y),
      "boruta_multiclass_roughfix", capture_history = TRUE),
    varselrf = sb_fit_varselrf_selector(x, y, "varselrf_multiclass_oob_1se"),
    vsurf = sb_fit_vsurf_classification_selector(x, y, "vsurf_multiclass_prediction"))
  fit$classes <- task$classes
  types <- switch(method, boruta = c("boruta_summary", "boruta_hits", "boruta_history"),
    varselrf = "selection_history", vsurf = c("vsurf_importance", "vsurf_stages"))
  diagnostics <- plots <- notes <- list()
  status <- lapply(types, function(type) {
    error <- tryCatch({
      diagnostic <- sb_selector_diagnostics(fit, type)
      plot <- sb_plot_selector_diagnostic(diagnostic)
      notes[[type]] <<- list(labels = plot$labels, source = diagnostic$source)
      diagnostics[[type]] <<- diagnostic
      plots[[type]] <<- sb_method_plot(plot, type, spec$label)
      NULL
    }, error = identity)
    data.frame(view = type, status = if (is.null(error)) "available" else "error",
      reason = if (is.null(error)) "" else conditionMessage(error))
  })
  notes$interpretation <- paste("Full-input multiclass feature selection, not held-out evaluation.",
    "These are selectors, not probability-producing classifiers.",
    "Native forest OOB resampling is row-based, not patient-group-isolated.")
  list(status = "success", method = method, label = spec$label,
    feature_result = "selection", selection_rule = spec$selection_rule,
    features = fit$scores, class_features = NULL, evidence = fit$scores,
    parameters = fit$params, fit = fit, diagnostics = diagnostics, plots = plots,
    notes = notes, plot_status = do.call(rbind, status))
}

sb_fit_multiclass_abess <- function(fit, x, y) {
  if (!requireNamespace("abess", quietly = TRUE)) sb_abort(
    "Install abess for this method.", "stormbreaker_missing_dependency")
  support <- 0:min(10L, ncol(x))
  args <- list(x = as.matrix(x), y = y, family = "multinomial",
      tune.path = "sequence", tune.type = "gic", support.size = support,
      normalize = 2L, fit.intercept = TRUE, lambda = 0, c.max = 2L,
      splicing.type = 2L, max.splicing.iter = 20L, newton = "approx",
      max.newton.iter = 10L, early.stop = FALSE, warm.start = TRUE,
      important.search = ncol(x), num.threads = 1L, seed = fit$seed)
  run <- function(args) {
    warnings <- character()
    engine <- withCallingHandlers(withr::with_seed(args$seed,
      do.call(abess::abess, args)), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    list(engine = engine, warnings = warnings)
  }
  environment(run) <- baseenv()
  if (.Platform$OS.type == "windows") {
    # Native multinomial ABESS can terminate a mixed-engine Windows session.
    # A fresh process retains the exact algorithm while containing native faults.
    worker <- parallel::makePSOCKcluster(1L)
    on.exit(parallel::stopCluster(worker), add = TRUE)
    result <- parallel::clusterCall(worker, function(run, args, paths) {
      .libPaths(paths)
      run(args)
    }, run, args, .libPaths())[[1L]]
  } else result <- run(args)
  fit$engine <- result$engine
  fit$warnings <- c(fit$warnings, result$warnings)
  fit$engine$call <- NULL
  fit$parameters <- list(family = "multinomial", tune_type = "gic", support_size = support,
    selected_support_size = fit$engine$best.size, normalize = 2L, lambda = 0,
    newton = "approx", max_newton_iter = 10L,
    reference_class = utils::tail(fit$classes, 1L), seed = fit$seed)
  structure(fit, class = "stormbreaker_multiclass_fit")
}

sb_multiclass_abess_evidence <- function(fit) {
  beta <- as.matrix(stats::coef(fit$engine, support.size = fit$engine$best.size)[[1]])
  # ABESS drops the last response column as reference. Its 0.4.11 coefficient
  # column labels instead drop the first; use the declared factor-level order.
  beta <- cbind(beta[fit$features, , drop = FALSE], 0)
  by_class <- do.call(rbind, lapply(seq_along(fit$classes), function(k)
    data.frame(feature = fit$features, class = fit$classes[k],
      coefficient = unname(beta[, k]), selected = unname(beta[, k] != 0),
      reference = k == length(fit$classes))))
  features <- data.frame(feature = fit$features, selected = rowSums(beta != 0) > 0,
    row.names = NULL)
  list(features = features, by_class = by_class)
}

sb_multiclass_abess_plots <- function(fit, evidence) {
  support_size <- gic <- coefficient <- feature <- NULL
  criterion <- data.frame(support_size = fit$engine$support.size, gic = fit$engine$tune.value)
  curve <- ggplot2::ggplot(criterion, ggplot2::aes(support_size, gic)) +
    ggplot2::geom_line(colour = "#2878A0") + ggplot2::geom_point(colour = "#2878A0") +
    ggplot2::geom_vline(xintercept = fit$engine$best.size, linetype = "dashed") +
    sb_viz_theme() + ggplot2::labs(title = "ABESS: Subset selection",
      x = "Features", y = "GIC")
  tab <- evidence$by_class
  tab$class <- factor(tab$class, levels = fit$classes)
  coefficients <- ggplot2::ggplot(tab, ggplot2::aes(coefficient, feature)) +
    ggplot2::geom_col(fill = "#2878A0") + ggplot2::facet_wrap(~class) +
    sb_viz_theme() + ggplot2::labs(title = "ABESS: Class coefficients",
      x = "Coefficient", y = NULL)
  list(plots = list(support_criterion = curve, class_coefficients = coefficients),
    tables = list(support_criterion = criterion, class_coefficients = tab))
}
