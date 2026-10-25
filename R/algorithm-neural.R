sb_neural_parameters <- function(state) {
  if (!state$learner %in% c("nnet_binary", "neuralnet_multilayer_binary")) sb_abort(
    "This diagnostic requires a supported neural network.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(state$engine_raw)
  m <- e$model
  if (state$learner == "nnet_binary") {
    sizes <- m$n
    weights <- list()
    offset <- c(0L, cumsum(sizes))
    for (layer in 1:2) {
      w <- matrix(0, sizes[layer] + 1L, sizes[layer + 1L])
      for (j in seq_len(ncol(w))) {
        destination <- offset[layer + 1L] + j
        at <- seq.int(m$nconn[destination + 1L] + 1L, m$nconn[destination + 2L])
        sources <- m$conn[at]
        row <- ifelse(sources == 0, 1L, sources - offset[layer] + 1L)
        if (length(at) != nrow(w) || anyDuplicated(row) || any(row < 1L | row > nrow(w))) sb_abort(
          "Saved nnet connections do not match the supported feed-forward topology.", "stormbreaker_invalid_diagnostic")
        w[row, j] <- m$wts[at]
      }
      weights[[layer]] <- w
    }
  } else {
    weights <- m$weights[[1L]]
    sizes <- c(nrow(weights[[1L]]) - 1L, vapply(weights, ncol, integer(1)))
  }
  if (sizes[1L] != length(state$features) || utils::tail(sizes, 1L) != 1L ||
      any(!is.finite(unlist(weights))) || any(vapply(seq_along(weights), function(i)
        nrow(weights[[i]]) != sizes[i] + 1L, logical(1)))) sb_abort(
    "Saved neural weight dimensions or values are invalid.", "stormbreaker_bad_learner_state")
  list(weights = weights, sizes = sizes, center = e$center, scale = e$scale,
    activation = if (state$learner == "nnet_binary") "nnet_saturated_logistic" else "logistic",
    saturation = if (state$learner == "nnet_binary") c(lower = -15, upper = 15) else NULL)
}

#' Inspect original neural-network connections
#' @param fit A successful nnet or multilayer neuralnet final fit.
#' @param type topology or weights.
#' @return An algorithm diagnostic from compact original connection weights.
#'   Bias connections and all hidden units are retained. Topology edges indicate
#'   connectivity, not importance or causality. Weight matrices use the model's
#'   original standardized inputs; weights are not log odds, SHAP or input
#'   importance, and are not comparable as feature effects across layers.
#'   Bias nodes are repeated visually per receiving layer for readability and
#'   represent the same constant value one. No new fit or optimization-history
#'   reconstruction is performed. The saved `parameters` permit native forward
#'   calculation with logistic activations at each hidden and output layer.
#'   nnet's native activation is exactly zero below -15 and one above 15;
#'   its saturation thresholds are retained rather than replaced by plogis.
#' @export
sb_neural_diagnostics <- function(fit, type = "topology") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("topology", "weights")) sb_abort(
    "Choose topology or weights.", "stormbreaker_unsupported_diagnostic")
  p <- sb_neural_parameters(fit$learner_state)
  layers <- length(p$weights)
  labels <- c(list(fit$learner_state$features), lapply(seq_len(layers), function(i)
    if (i == layers) "Class 1" else paste0("H", i, ".", seq_len(p$sizes[i + 1L]))))
  nodes <- do.call(rbind, lapply(seq_along(labels), function(i) {
    n <- length(labels[[i]])
    data.frame(id = paste0("L", i - 1L, "_", seq_len(n)), label = labels[[i]], layer = i - 1L,
      position = rev(seq_len(n)) - (n + 1) / 2, bias = FALSE)
  }))
  bias <- data.frame(id = paste0("B", seq_len(layers)), label = "Bias", layer = seq_len(layers) - 1L,
    position = max(p$sizes) / 2 + 1, bias = TRUE)
  nodes <- rbind(nodes, bias)
  edges <- do.call(rbind, lapply(seq_len(layers), function(i) {
    w <- p$weights[[i]]
    data.frame(stage = paste0("Layer ", i - 1L, " -> ", i),
      source = rep(c(paste0("B", i), paste0("L", i - 1L, "_", seq_len(nrow(w) - 1L))), ncol(w)),
      target = rep(paste0("L", i, "_", seq_len(ncol(w))), each = nrow(w)),
      source_label = rep(c("Bias", labels[[i]]), ncol(w)),
      target_label = rep(labels[[i + 1L]], each = nrow(w)), weight = as.vector(w))
  }))
  structure(list(type = paste0("neural_", type), table = edges, nodes = nodes, parameters = p,
    learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    source = "original saved neural-network weights, connections and input scaling",
    interpretation = paste("Original feed-forward logistic network, retaining every connection and bias.",
      "Input coordinates use the fitted training centering/scaling; weights are not feature importance, log odds or SHAP.",
      "Repeated bias nodes represent constant one. No optimization history or new fit.",
      if (p$activation == "nnet_saturated_logistic") "Native nnet activations saturate to 0 below -15 and 1 above 15." else
        "Original neuralnet logistic activation retained.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_neural_diagnostic <- function(x) {
  layer <- position <- label <- xend <- yend <- source_label <- target_label <- weight <- NULL
  start <- end <- NULL
  nodes <- x$nodes
  tab <- x$table
  if (x$type == "neural_topology") {
    from <- match(tab$source, nodes$id)
    to <- match(tab$target, nodes$id)
    tab$layer <- nodes$layer[from]
    tab$position <- nodes$position[from]
    tab$xend <- nodes$layer[to]
    tab$yend <- nodes$position[to]
    tab$start <- tab$layer + .12
    tab$end <- tab$xend - .12
    p <- ggplot2::ggplot(tab, ggplot2::aes(layer, position)) +
      ggplot2::geom_segment(ggplot2::aes(xend = start, yend = position), colour = "#9AA6AF", alpha = .55, linewidth = .35) +
      ggplot2::geom_segment(ggplot2::aes(x = start, xend = end, yend = yend), colour = "#9AA6AF", alpha = .55, linewidth = .35) +
      ggplot2::geom_segment(ggplot2::aes(x = end, y = yend, xend = xend, yend = yend), colour = "#9AA6AF", alpha = .55, linewidth = .35) +
      ggplot2::geom_point(data = nodes, shape = 21, fill = "white", colour = "#2878A0", size = 4) +
      ggplot2::geom_label(data = nodes, ggplot2::aes(label = label), nudge_y = .3,
        size = 3, linewidth = 0, fill = "white", label.padding = grid::unit(.08, "lines")) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = .35)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = .6)) +
      ggplot2::labs(title = "Neural-network structure") + sb_viz_theme() +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
  } else {
    tab$source_label <- factor(tab$source_label, levels = rev(unique(tab$source_label)))
    tab$target_label <- factor(tab$target_label, levels = unique(tab$target_label))
    p <- ggplot2::ggplot(tab, ggplot2::aes(target_label, source_label, fill = weight)) +
      ggplot2::geom_tile() + ggplot2::facet_wrap(~stage, scales = "free") +
      sb_neural_weight_scale() +
      ggplot2::labs(title = "Neural connection weights", x = "Receiving unit", y = "Sending unit", fill = "Weight") +
      sb_viz_theme() + ggplot2::theme(panel.grid = ggplot2::element_blank())
  }
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = max(6, 2 + .5 * max(x$parameters$sizes)))
  if (x$type != "neural_topology") attr(p, "stormbreaker_notes") <- c(attr(p, "stormbreaker_notes"),
    "Weight colors use a symmetric pseudo-log scale (sigma = 1); legend values are original weights. No clipping or rescaling of model coefficients.")
  p
}

sb_neural_weight_scale <- function() {
  ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#D77935",
    transform = scales::pseudo_log_trans(sigma = 1),
    breaks = function(limits) {
      maximum <- max(abs(limits))
      if (maximum == 0) return(0)
      if (maximum <= 10) return(pretty(limits, n = 5))
      powers <- 10^seq(max(0, floor(log10(maximum)) - 2), floor(log10(maximum)))
      c(-rev(powers), 0, powers)
    }, guide = ggplot2::guide_colourbar(display = "rectangles"))
}
