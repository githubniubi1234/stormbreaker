sb_captured_selector_adapter <- function(id, params) {
  if (id %in% sb_survival_weight_selector_ids()) return(sb_survival_weight_selector_adapter(id, params, TRUE))
  if (id == "ordered_chain") return(sb_chain_selector_adapter(params, TRUE))
  if (id %in% sb_boruta_selector_ids()) return(sb_forest_selector_adapter(id, capture_history = TRUE, options = params))
  if (id %in% sb_binary_model_selector_ids()) return(sb_binary_model_selector_adapter(id, params, TRUE))
  if (id %in% sb_abess_selector_ids()) return(sb_abess_selector_adapter(id, TRUE))
  if (id == "coxboost_nonzero") return(sb_coxboost_nonzero_selector_adapter(params$inner_folds, params$seed, TRUE))
  if (id == "coxboost_then_stepcox") {
    adapter <- sb_sequential_survival_selector_adapter(
      sb_coxboost_nonzero_selector_adapter(params$inner_folds, params$seed, TRUE),
      sb_stepcox_selector_adapter(params$direction), id)
    adapter$params <- params
    return(adapter)
  }
  if (id %in% c("lasso_nonzero", "enet_cox_nonzero")) return(
    sb_lasso_nonzero_selector_adapter(params$inner_folds, params$seed, id, TRUE))
  if (id == "lasso_then_stepcox") {
    adapter <- sb_sequential_survival_selector_adapter(
      sb_lasso_nonzero_selector_adapter(params$inner_folds, params$seed, capture_history = TRUE),
      sb_stepcox_selector_adapter(params$direction), id)
    adapter$params <- params
    return(adapter)
  }
  sb_resolve_selector_adapter(id, params)
}

sb_selector_path_plot <- function(path, selected_lambda) {
  if (is.null(path) || is.null(path$beta) || !length(path$lambda)) {
    stop("Coefficient path was not retained. Refit once with capture_diagnostics=TRUE.")
  }
  b <- as.matrix(path$beta)
  d <- data.frame(lambda = rep(path$lambda, each = nrow(b)),
    feature = rep(rownames(b), ncol(b)), coefficient = as.vector(b))
  k <- unique(round(seq(1, ncol(b), length.out = min(7L, ncol(b)))))
  lambda <- feature <- coefficient <- NULL
  ggplot2::ggplot(d, ggplot2::aes(lambda, coefficient, colour = feature)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = selected_lambda, linetype = 2, colour = "grey50") +
    ggplot2::scale_x_log10(sec.axis = ggplot2::dup_axis(breaks = path$lambda[k],
      labels = colSums(b != 0)[k], name = "Nonzero coefficients")) +
    ggplot2::labs(x = "Lambda (log scale)", y = "Coefficient", colour = "Gene") + sb_viz_theme()
}

sb_selector_oob_plot <- function(d) {
  if (is.null(d)) stop("OOB history was not retained. Refit once with capture_diagnostics=TRUE.")
  classes <- setdiff(names(d), "trees")
  long <- data.frame(trees = rep(d$trees, length(classes)),
    series = rep(classes, each = nrow(d)), error = unlist(d[classes], use.names = FALSE))
  trees <- series <- error <- NULL
  colours <- stats::setNames(grDevices::hcl.colors(length(classes), "Dark 3"), classes)
  known <- c("0" = "#2878A0", "1" = "#B86236", "oob" = "#686868", "OOB" = "#686868")
  common <- intersect(names(colours), names(known))
  colours[common] <- known[common]
  ggplot2::ggplot(long, ggplot2::aes(trees, error, colour = series)) +
    ggplot2::geom_line(linewidth = .5) + ggplot2::scale_colour_manual(values = colours) +
    ggplot2::labs(x = "Trees", y = "OOB error", colour = NULL) + sb_viz_theme()
}
