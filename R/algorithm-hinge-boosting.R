sb_hinge_boosting_diagnostic <- function(fit, type) {
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("training_loss", "component_usage", "coefficient_path")) sb_abort(
    "Hinge boosting supports training_loss, component_usage and coefficient_path.",
    "stormbreaker_unavailable_diagnostic")
  model <- fit$diagnostic_artifacts$learner
  path <- model$diagnostic_path
  if (is.null(path)) sb_abort(
    "Original hinge-boosting iterations were not captured; enable capture_diagnostics on a new final fit. Historical compact fits contain only the final risk.",
    "stormbreaker_missing_diagnostic_artifacts")
  n <- model$params$iterations
  features <- model$features
  if (!identical(path$family, "hinge") || !identical(path$learner, "ls") ||
      length(path$risk) != n || length(path$selected) != n || length(path$updates) != n ||
      any(!is.finite(path$risk)) || any(!is.finite(path$updates)) ||
      any(!is.finite(path$selected)) || any(path$selected != floor(path$selected)) ||
      any(path$selected < 1L) || any(path$selected > length(features)) ||
      length(path$meanx) != length(features) || any(!is.finite(path$meanx)) ||
      length(path$offset) != 1L || !is.finite(path$offset) ||
      length(path$learning_rate) != 1L || !is.finite(path$learning_rate) ||
      !is.null(path$intercept_updates)) sb_abort(
    "Captured native hinge-boosting path is invalid or uses unsupported intercept updates.",
    "stormbreaker_invalid_diagnostic")
  selections <- data.frame(iteration = seq_len(n), feature = features[path$selected],
    native_feature_index = path$selected, native_update = path$updates)
  if (type == "training_loss") {
    tab <- data.frame(iteration = seq_len(n), loss = as.numeric(path$risk))
    measure <- "Native mean hinge loss after each update"
  } else if (type == "component_usage") {
    tab <- data.frame(feature = features, updates = tabulate(path$selected, nbins = length(features)))
    measure <- "Native hinge-boosting component updates"
  } else {
    beta <- matrix(0, nrow = n + 1L, ncol = length(features))
    for (i in seq_len(n)) {
      beta[i + 1L, ] <- beta[i, ]
      j <- path$selected[i]
      beta[i + 1L, j] <- beta[i + 1L, j] + path$learning_rate * path$updates[i]
    }
    intercept <- as.numeric(path$offset - beta %*% path$meanx)
    if (!isTRUE(all.equal(as.numeric(beta[n + 1L, ]), as.numeric(model$engine$coefficient), tolerance = 1e-12)) ||
        !isTRUE(all.equal(intercept[n + 1L], model$engine$intercept, tolerance = 1e-12))) sb_abort(
      "Hinge path endpoint disagrees with the original margin equation.", "stormbreaker_invalid_diagnostic")
    beta <- cbind(intercept, beta)
    tab <- data.frame(iteration = rep(seq.int(0L, n), times = ncol(beta)),
      feature = rep(c("(Intercept)", features), each = nrow(beta)), coefficient = as.vector(beta))
    measure <- "Original hinge-margin coefficient path"
  }
  structure(list(type = paste0("hinge_", type), table = tab, selections = selections,
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id, measure = measure,
    source = "captured original hinge fit; no additional training",
    interpretation = "Native hinge-margin scale, not probability or log-odds. Original stopping retained.\nTraining diagnostics, not external validation, uncertainty intervals or causal effects."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_hinge_boosting_diagnostic <- function(x) {
  iteration <- loss <- feature <- coefficient <- NULL
  if (x$type == "hinge_component_usage") {
    z <- x
    z$type <- "boosting_component_usage"
    return(sb_plot_boosting_diagnostic(z) + ggplot2::labs(title = x$measure))
  }
  p <- if (x$type == "hinge_training_loss") {
    ggplot2::ggplot(x$table, ggplot2::aes(iteration, loss)) +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::labs(y = "Mean native hinge loss", x = "Completed boosting update (starts at 1)")
  } else {
    ggplot2::ggplot(x$table, ggplot2::aes(iteration, coefficient, colour = feature)) +
      ggplot2::geom_line() + ggplot2::labs(y = "Coefficient in original margin equation",
        x = "Boosting update (zero = initial equation)", colour = "Term")
  }
  p + ggplot2::labs(title = x$measure, subtitle = x$learner, caption = x$interpretation) + sb_viz_theme()
}
