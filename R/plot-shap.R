#' Plot existing additive-SHAP explanations
#'
#' @param x A successful [sb_additive_shap()] or [sb_tree_shap()] or [sb_model_shap()] result.
#'   The saved exact/approximate method and score scale are preserved.
#' @param top_n Maximum features shown, selected by mean absolute contribution
#'   for importance/beeswarm and by absolute individual contribution for the
#'   waterfall. Ties retain feature input order.
#' @param groupby Optional group vector in explained-observation order. A named
#'   vector is instead aligned by observation IDs. Groups must be nonmissing.
#'   Uses the same saved model, SHAP values and background for every group;
#'   top_n is ranked globally, then mean absolute SHAP is computed within groups.
#'
#' @return A ggplot object containing plotted contributions and descriptive
#'   feature labels, without the original patient identifiers. These graphics
#'   can still disclose individual predictions or feature values; review them
#'   before sharing. Optional ggplot2 is required.
#'
#' @details Plots preserve the saved score scale and empirical background
#'   reference. Additive explanations are not probability effects; tree SHAP
#'   can explicitly explain class-1 probability on its declared score scale.
#'   Importance is mean absolute contribution over the explained observations,
#'   not a causal importance estimate. The beeswarm uses deterministic binned
#'   offsets and does not change the random state. Without supplied feature
#'   values all points have one colour, rather than invented high/low values.
#'   Optional values are aligned by explicit IDs and colour-normalized within
#'   each feature; constant features map to the colour midpoint.
#'
#'   Waterfalls retain a baseline, the displayed feature terms, one sum of any
#'   remaining terms, and the full reconstructed score. The displayed terms
#'   therefore add up even when top_n omits features. Observation numbers, not
#'   patient IDs, label individual plots. No model is fitted or explained again.
#'   Importance bars display numeric mean absolute contributions. Waterfalls
#'   use directional arrows (red increases, blue decreases), neutral endpoint
#'   markers and a separate numeric column so narrow terms remain readable.
#'   Beeswarm colour instead encodes within-feature values, not SHAP signs.
#'
#' @export
sb_plot_shap_importance <- function(x, top_n = 20L, groupby = NULL) {
  m <- sb_viz_shap_matrix(x)
  top_n <- sb_viz_positive_integer(top_n, "top_n")
  importance <- colMeans(abs(m))
  chosen <- utils::head(order(-importance, seq_along(importance)), top_n)
  data <- data.frame(
    feature = colnames(m)[chosen], mean_absolute_shap = unname(importance[chosen]),
    n_observations = nrow(m), stringsAsFactors = FALSE
  )
  features <- data$feature
  if (!is.null(groupby)) {
    if (!is.atomic(groupby) || !is.null(dim(groupby)) || length(groupby) != nrow(m) || anyNA(groupby)) {
      sb_abort("groupby must supply one nonmissing group per explained observation.", "stormbreaker_bad_visualization_argument")
    }
    ids <- names(groupby)
    if (!is.null(ids)) {
      if (anyNA(ids) || anyDuplicated(ids) || !setequal(ids, rownames(m))) {
        sb_abort("Named groupby IDs must match explained observations exactly.", "stormbreaker_bad_visualization_argument")
      }
      groupby <- groupby[match(rownames(m), ids)]
    }
    groups <- as.character(groupby)
    if (any(!nzchar(trimws(groups)))) sb_abort("groupby labels cannot be empty.", "stormbreaker_bad_visualization_argument")
    levels <- if (is.factor(groupby)) levels(droplevels(groupby)) else unique(groups)
    data <- do.call(rbind, lapply(levels, function(g) {
      rows <- groups == g
      data.frame(feature = features,
        mean_absolute_shap = unname(colMeans(abs(m[rows, chosen, drop = FALSE]))),
        n_observations = sum(rows), group = g, stringsAsFactors = FALSE)
    }))
    data$group <- factor(data$group, levels = levels)
  }
  data$feature_axis <- factor(data$feature, levels = rev(features))
  data$value_label <- sb_shap_number(data$mean_absolute_shap)
  feature_axis <- mean_absolute_shap <- value_label <- group <- NULL
  p <- ggplot2::ggplot(data, ggplot2::aes(x = mean_absolute_shap, y = feature_axis))
  if (is.null(groupby)) {
    p <- p + ggplot2::geom_col(fill = "#2878B5", width = 0.62)
  } else {
    p <- p + ggplot2::aes(fill = group, group = group) +
      ggplot2::geom_col(width = 0.72, position = ggplot2::position_dodge(width = 0.8, orientation = "y")) +
      ggplot2::labs(fill = "Group")
  }
  p +
    ggplot2::geom_text(ggplot2::aes(label = value_label),
      position = if (is.null(groupby)) "identity" else ggplot2::position_dodge(width = 0.8, orientation = "y"),
      hjust = -0.2,
      family = "Arial", size = 3.5, colour = "#253746"
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(
      title = "SHAP importance", x = paste0("Mean |SHAP| (", sb_viz_score_label(x$score_scale), ")"), y = NULL
    ) +
    sb_shap_theme()
}

#' @rdname sb_plot_shap_importance
#' @param feature_values Optional data frame containing the explained feature
#'   values and an explicit ID column. Rows may be reordered but IDs must match
#'   the explained observations exactly; only fitted feature columns are read.
#' @param id_col ID-column name, required with feature_values.
#' @export
sb_plot_shap_beeswarm <- function(x, top_n = 20L, feature_values = NULL, id_col = NULL) {
  m <- sb_viz_shap_matrix(x)
  top_n <- sb_viz_positive_integer(top_n, "top_n")
  importance <- colMeans(abs(m))
  chosen <- utils::head(order(-importance, seq_along(importance)), top_n)
  values <- NULL
  if (!is.null(feature_values)) {
    valid <- is.data.frame(feature_values) && is.character(id_col) && length(id_col) == 1L &&
      !is.na(id_col) && id_col %in% names(feature_values) && all(colnames(m) %in% names(feature_values))
    if (!valid) {
      sb_abort(
        "Provide feature_values with an explicit id_col and all explained feature columns.",
        "stormbreaker_bad_visualization_argument"
      )
    }
    ids <- as.character(feature_values[[id_col]])
    if (anyNA(ids) || anyDuplicated(ids) || !setequal(ids, rownames(m)) || length(ids) != nrow(m)) {
      sb_abort("Feature-value IDs must match explained observations exactly.", "stormbreaker_bad_visualization_argument")
    }
    values <- feature_values[match(rownames(m), ids), colnames(m), drop = FALSE]
    if (!all(vapply(values, function(value) is.numeric(value) && all(is.finite(value)), logical(1)))) {
      sb_abort("Explained feature values must be finite numeric columns.", "stormbreaker_bad_visualization_argument")
    }
  } else if (!is.null(id_col)) sb_abort("id_col is only used with feature_values.", "stormbreaker_bad_visualization_argument")
  rows <- lapply(seq_along(chosen), function(i) {
    j <- chosen[[i]]
    contribution <- m[, j]
    span <- diff(range(contribution))
    bins <- if (span == 0) {
      rep(1L, length(contribution))
    } else {
      floor(30 * (contribution - min(contribution)) / span)
    }
    max_bin <- max(tabulate(match(bins, unique(bins))))
    offset <- numeric(length(contribution))
    for (bin in unique(bins)) {
      idx <- which(bins == bin)
      idx <- idx[order(contribution[idx], idx)]
      sequence <- seq_along(idx) - 1L
      offset[idx] <- ceiling(sequence / 2) * ifelse(sequence %% 2 == 1, 1, -1) * 0.7 / max(2, max_bin)
    }
    tab <- data.frame(
      observation = seq_len(nrow(m)), feature = colnames(m)[j],
      contribution = unname(contribution), y = length(chosen) - i + 1 + offset,
      feature_value = NA_real_, colour_value = NA_real_, stringsAsFactors = FALSE
    )
    if (!is.null(values)) {
      v <- values[[j]]
      width <- diff(range(v))
      tab$feature_value <- v
      tab$colour_value <- if (width == 0) rep(0.5, length(v)) else (v - min(v)) / width
    }
    tab
  })
  data <- sb_view_bind(rows)
  contribution <- y <- colour_value <- NULL
  p <- ggplot2::ggplot(data, ggplot2::aes(x = contribution, y = y)) +
    ggplot2::geom_hline(yintercept = seq_along(chosen), colour = "#E3E8EC", linetype = 3, linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 0, colour = "#9AA7B0", linewidth = 0.45)
  if (is.null(values)) {
    p <- p + ggplot2::geom_point(colour = "#2878B5", size = 1.9, alpha = 0.8)
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = colour_value), size = 1.9, alpha = 0.8) +
      ggplot2::scale_colour_gradientn(
        colours = c("#2878B5", "#AD92B5", "#D45064"), limits = c(0, 1),
        breaks = c(0, 1), labels = c("Low", "High"), name = "Feature value\n(within gene)",
        guide = ggplot2::guide_colourbar(
          barheight = grid::unit(38, "mm"), barwidth = grid::unit(3, "mm"),
          theme = ggplot2::theme(legend.text = ggplot2::element_text(margin = ggplot2::margin(l = 8)))
        )
      )
  }
  p + ggplot2::scale_y_continuous(breaks = seq_along(chosen), labels = rev(colnames(m)[chosen])) +
    ggplot2::labs(
      title = "SHAP summary", x = paste0("SHAP (", sb_viz_score_label(x$score_scale), ")"), y = NULL,
      subtitle = paste(nrow(m), "observations | one dot per sample per gene | left: lower score; right: higher score"),
      caption = sb_viz_wrap(sb_viz_shap_caption(x), width = 105)
    ) + sb_shap_theme()
}

#' @rdname sb_plot_shap_importance
#' @param row One-based explained-observation index for the waterfall.
#' @export
sb_plot_shap_waterfall <- function(x, row = 1L, top_n = 10L) {
  m <- sb_viz_shap_matrix(x)
  row <- sb_viz_positive_integer(row, "row")
  top_n <- sb_viz_positive_integer(top_n, "top_n")
  if (row > nrow(m)) sb_abort("row exceeds the number of explained observations.", "stormbreaker_bad_visualization_argument")
  terms <- m[row, ]
  chosen <- utils::head(order(-abs(terms), seq_along(terms)), top_n)
  other <- setdiff(seq_along(terms), chosen)
  contributions <- unname(terms[chosen])
  feature <- names(terms)[chosen]
  kind <- rep("feature", length(chosen))
  if (length(other)) {
    contributions <- c(contributions, sum(terms[other]))
    feature <- c(feature, paste0("Other features (", length(other), ")"))
    kind <- c(kind, "remainder")
  }
  end <- x$reference_value + cumsum(contributions)
  start <- c(x$reference_value, utils::head(end, -1))
  score <- x$reference_value + sum(terms)
  data <- data.frame(
    feature = c("Background reference", feature, "Final model score"),
    kind = c("reference", kind, "prediction"),
    contribution = c(x$reference_value, contributions, score),
    start = c(0, start, 0), end = c(x$reference_value, end, score), stringsAsFactors = FALSE
  )
  data$y <- rev(seq_len(nrow(data)))
  data$xmin <- pmin(data$start, data$end)
  data$xmax <- pmax(data$start, data$end)
  data$direction <- ifelse(data$kind == "reference", "Reference", ifelse(data$kind == "prediction", "Final score",
    ifelse(data$contribution >= 0, "Increase", "Decrease")
  ))
  data$value_label <- sb_shap_number(data$contribution,
    signed = !data$kind %in% c("reference", "prediction")
  )
  data$label <- data$value_label
  extent <- range(c(x$reference_value, start, end, score))
  span <- diff(extent)
  if (span == 0) span <- max(1, abs(score) * 0.1)
  data$label_x <- extent[2] + span * 0.12
  terms_data <- data[!data$kind %in% c("reference", "prediction"), , drop = FALSE]
  endpoints <- data[data$kind %in% c("reference", "prediction"), , drop = FALSE]
  polygons <- lapply(seq_len(nrow(terms_data)), function(i) {
    d <- terms_data[i, ]
    tip <- min(abs(d$end - d$start) * 0.4, span * 0.016)
    neck <- d$end - sign(d$contribution) * tip
    data.frame(
      px = c(d$start, neck, d$end, neck, d$start),
      py = d$y + c(-0.27, -0.27, 0, 0.27, 0.27),
      group = i, direction = d$direction
    )
  })
  polygons <- sb_view_bind(polygons)
  links <- data[-nrow(data), , drop = FALSE]
  links$next_y <- utils::tail(data$y, -1)
  xmin <- xmax <- y <- direction <- end <- next_y <- label <- NULL
  px <- py <- group <- label_x <- value_label <- NULL
  ggplot2::ggplot(data) +
    ggplot2::geom_segment(ggplot2::aes(y = y, yend = y),
      x = extent[1] - span * 0.06, xend = extent[2] + span * 0.07,
      colour = "#E3E8EC", linetype = 3, linewidth = 0.3
    ) +
    ggplot2::geom_polygon(
      data = polygons,
      ggplot2::aes(x = px, y = py, group = group, fill = direction)
    ) +
    ggplot2::geom_point(
      data = endpoints, ggplot2::aes(x = end, y = y),
      shape = 18, size = 3, colour = "#536777"
    ) +
    ggplot2::geom_segment(
      data = terms_data[terms_data$contribution == 0, , drop = FALSE],
      ggplot2::aes(x = end, xend = end, y = y - 0.2, yend = y + 0.2), colour = "#9AA7B0"
    ) +
    ggplot2::geom_segment(
      data = links, ggplot2::aes(x = end, xend = end, y = y - 0.3, yend = next_y + 0.3),
      colour = "#9AA7B0", linetype = 3, linewidth = 0.4
    ) +
    ggplot2::geom_text(ggplot2::aes(x = label_x, y = y, label = value_label),
      hjust = 0, size = 3.5, family = "Arial", colour = "#253746"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(extent[1] - span * 0.06, extent[2] + span * 0.32),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_continuous(breaks = data$y, labels = data$feature) +
    ggplot2::scale_fill_manual(
      values = c(Increase = "#D45064", Decrease = "#2878B5"),
      breaks = c("Increase", "Decrease"), labels = c("Increases score", "Decreases score"), drop = FALSE
    ) +
    ggplot2::labs(
      title = "Individual prediction",
      subtitle = paste("Observation", row, "| contributions ranked by magnitude | numbers: signed contribution or endpoint score"),
      x = paste0("Score (", sb_viz_score_label(x$score_scale), ")"), y = NULL, fill = NULL,
      caption = sb_viz_wrap(paste(sb_viz_shap_caption(x), "Omitted terms are summed, not discarded."), width = 105)
    ) +
    sb_shap_theme() +
    ggplot2::theme(legend.position = "top", legend.justification = "left")
}

sb_shap_number <- function(x, signed = FALSE) {
  label <- vapply(x, function(value) format(signif(value, 3), trim = TRUE), character(1))
  ifelse(rep_len(signed, length(x)) & x >= 0, paste0("+", label), label)
}

sb_shap_theme <- function() {
  ggplot2::theme_classic(base_size = 11, base_family = "Arial") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold", colour = "#172B3A", margin = ggplot2::margin(b = 8)),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "#253746"),
      axis.text.y = ggplot2::element_text(size = 11),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 12)),
      axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(colour = "#9AA7B0", linewidth = 0.4),
      axis.ticks.x = ggplot2::element_line(colour = "#9AA7B0", linewidth = 0.4),
      legend.title = ggplot2::element_text(size = 10), legend.text = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(16, 22, 12, 12), plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

sb_viz_shap_matrix <- function(x) {
  sb_require_curve_plotting()
  m <- x$contributions
  valid <- inherits(x, c("stormbreaker_additive_shap", "stormbreaker_tree_shap", "stormbreaker_model_shap")) && identical(x$status, "success") &&
    is.matrix(m) && is.numeric(m) && nrow(m) > 0 && ncol(m) > 0 && all(is.finite(m)) &&
    !is.null(colnames(m)) && !is.null(rownames(m)) && !anyNA(colnames(m)) && !anyNA(rownames(m)) &&
    !anyDuplicated(colnames(m)) && !anyDuplicated(rownames(m)) &&
    is.numeric(x$reference_value) && length(x$reference_value) == 1L && is.finite(x$reference_value) &&
    is.character(x$score_scale) && length(x$score_scale) == 1L && !is.na(x$score_scale)
  if (!valid) sb_abort("Supply a valid successful SHAP result with finite named contributions.", "stormbreaker_bad_shap_visualization")
  if (!is.data.frame(x$scores) || !identical(as.character(x$scores$patient_id), rownames(m)) ||
    !isTRUE(all.equal(unname(x$scores$score), unname(x$reference_value + rowSums(m)), tolerance = 1e-10))) {
    sb_abort("Saved SHAP scores do not reconstruct from the reference and contributions.", "stormbreaker_bad_shap_visualization")
  }
  m
}

sb_viz_shap_caption <- function(x) {
  if (inherits(x, c("stormbreaker_tree_shap", "stormbreaker_model_shap"))) return(paste0(
    "Approximate permutation interventional SHAP; ", x$permutations,
    " orders, background n = ", x$background$n, "; scale = ", x$score_scale, ". Not exact TreeSHAP."))
  paste0(
    "Exact empirical interventional contributions on the saved score scale, not probability effects. Background n = ",
    x$background$n, "; reference = ", format(x$reference_value, digits = 4), ". Descriptive, not causal."
  )
}
