#' Inspect original flexible, penalized and mixture discriminant geometry
#' @param fit Successful FDA, PDA or MDA final fit.
#' @param type centers, query_variates, or mixture_weights (MDA only).
#' @param newdata Explicit query rows for query_variates.
#' @return Original class/subclass centers, within-class mixture weights,
#'   or native query discriminant coordinates. Axes are discriminant variates,
#'   not original features. Mixture subclasses are fitted latent components,
#'   not clinical subgroups. No refitting, new projection or covariance estimate.
#' @export
sb_fda_diagnostics <- function(fit, type = "centers", newdata = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% c("fda_binary", "pda_binary", "mda_binary") ||
      !is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("centers", "query_variates", "mixture_weights")) {
    sb_abort("Choose centers, query_variates or mixture_weights for FDA/PDA/MDA.",
      "stormbreaker_unsupported_diagnostic")
  }
  if (!requireNamespace("mda", quietly = TRUE)) sb_abort(
    "These diagnostics require mda.", "stormbreaker_unavailable_diagnostic")
  e <- unserialize(s$engine_raw)
  mixture <- s$learner == "mda_binary"
  if (!inherits(e, if (mixture) "mda" else "fda") || !is.matrix(e$means) ||
      any(!is.finite(e$means))) sb_abort(
    "Original discriminant centers are invalid.", "stormbreaker_invalid_diagnostic")
  mapping <- if (mixture) do.call(rbind, lapply(names(e$assign.theta), function(cls) {
    data.frame(center = e$assign.theta[[cls]], class = cls,
      component = seq_along(e$assign.theta[[cls]]))
  })) else data.frame(center = seq_len(nrow(e$means)), class = rownames(e$means), component = 1L)
  mapping <- mapping[order(mapping$center), ]
  if (!identical(as.integer(mapping$center), seq_len(nrow(e$means)))) sb_abort(
    "Original subclass assignments are invalid.", "stormbreaker_invalid_diagnostic")
  if (type == "centers") {
    tab <- do.call(rbind, lapply(seq_len(ncol(e$means)), function(j) {
      data.frame(mapping, axis = j, coordinate = e$means[, j])
    }))
  } else if (type == "mixture_weights") {
    if (!mixture) sb_abort("Only MDA has fitted mixture subclasses.", "stormbreaker_unavailable_diagnostic")
    tab <- mapping
    tab$weight <- vapply(seq_len(nrow(tab)), function(i) {
      unname(e$sub.prior[[tab$class[i]]][tab$component[i]])
    }, numeric(1))
  } else {
    if (is.null(newdata)) sb_abort("Supply explicit query rows for discriminant variates.", "stormbreaker_unavailable_diagnostic")
    x <- sb_select_final_feature_frame(newdata, s$features, "Query data", "query fitted features")
    if (!nrow(x)) sb_abort("Queries must contain at least one row.", "stormbreaker_invalid_diagnostic")
    m <- sb_predict_fda_native(e, x, type = "variates", dimension = e$dimension)
    if (!is.matrix(m) || nrow(m) != nrow(x) || any(!is.finite(m))) sb_abort(
      "Native query variates are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(query = rep(seq_len(nrow(m)), ncol(m)),
      axis = rep(seq_len(ncol(m)), each = nrow(m)), coordinate = as.numeric(m))
  }
  structure(list(type = paste0("fda_", type), table = tab,
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "Original mda-package centers, subclass priors and native variates",
    interpretation = "Discriminant axes are not input features. MDA components are latent model structure, not clinical subgroups. Axis signs are basis conventions. No refitting or new projection."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_fda_diagnostic <- function(x) {
  coordinate <- class <- component <- weight <- query <- NULL
  if (x$type == "fda_centers") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(coordinate, factor(component), colour = class)) +
      ggplot2::geom_point(size = 2.5) + ggplot2::facet_wrap(~axis, labeller = ggplot2::label_both) +
      ggplot2::scale_colour_manual(values = c("0" = "#2878A0", "1" = "#B86236")) +
      ggplot2::labs(title = "Discriminant centers", x = "Discriminant coordinate", y = "Component", colour = "Class")
  } else if (x$type == "fda_mixture_weights") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(factor(component), weight)) +
      ggplot2::geom_col(fill = "#2878A0", width = .6) +
      ggplot2::facet_wrap(~class, labeller = ggplot2::label_both) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = "Within-class mixture weights", x = "Component", y = "Weight")
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(coordinate, factor(query))) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::facet_wrap(~axis, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "Discriminant query coordinates", x = "Discriminant coordinate", y = "Query")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}
