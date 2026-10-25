#' Visualize saved feature-selection stability
#'
#' @param x A CV result/comparison or supported composite CV result; a uniquely
#'   named list may identify result sets. External evaluations have no training
#'   fold selection records and are not accepted.
#' @param source Which saved sets to show: `"selector"` or
#'   `"glmnet_active_set"`. These are different notions of selection.
#' @param top_n Maximum features shown in each frequency panel. Selection for
#'   display uses decreasing saved frequency, with input-order ties.
#' @param cohort_label Optional label for a single object.
#'
#' @return A ggplot object with saved frequency or Jaccard data and a
#'   `stormbreaker_selection_status` attribute. The UpSet variant also attaches
#'   all exact intersections as `stormbreaker_intersections`, including features
#'   omitted from the displayed top intersections.
#'
#' @details Frequencies use the saved successful-fold denominator, not all
#'   requested folds. Incomplete selections are explicitly labelled and do not
#'   become complete merely because the plot draws. Jaccard shows only saved
#'   unique fold pairs; no diagonal, missing pair or zero is invented.
#'   A failed model fit does not erase separately successful selector records.
#'   Repetitions, cohorts and models remain separate panels. These are
#'   descriptive stability displays, not biological importance or causal
#'   evidence. An unscreened or fixed feature panel is not reclassified as a
#'   learned feature-selection result.
#'
#' @export
sb_plot_selection_frequency <- function(x, source = c("selector", "glmnet_active_set"),
                                         top_n = 20L, cohort_label = NULL) {
  source <- match.arg(source)
  top_n <- sb_viz_positive_integer(top_n, "top_n")
  view <- sb_viz_selection_view(x, source, cohort_label)
  rows <- lapply(view$entries, function(entry) {
    if (!entry$usable) return(NULL)
    tab <- entry$summary$feature_frequency
    if (!is.data.frame(tab) || !nrow(tab)) return(NULL)
    ranking <- order(-tab$selection_frequency, seq_len(nrow(tab)), na.last = TRUE)
    chosen <- utils::head(ranking[is.finite(tab$selection_frequency[ranking])], top_n)
    tab <- tab[sort(chosen), , drop = FALSE]
    if (!nrow(tab)) return(NULL)
    tab$panel <- entry$panel
    tab$n_successful_folds <- entry$summary$n_successful_folds
    tab$n_folds <- entry$summary$n_folds
    tab$feature_key <- paste(entry$key, tab$feature, sep = "::")
    tab$display_order <- match(tab$feature, entry$summary$feature_frequency$feature[chosen])
    tab
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(sb_viz_empty_selection(view, "Feature selection frequency"))
  data <- sb_view_bind(rows)
  levels <- unlist(lapply(rows, function(row) row$feature_key[order(row$display_order)]), use.names = FALSE)
  data$feature_key <- factor(data$feature_key, levels = rev(levels))
  labels <- stats::setNames(data$feature, as.character(data$feature_key))
  selection_frequency <- feature_key <- panel <- NULL
  p <- ggplot2::ggplot(data, ggplot2::aes(x = selection_frequency, y = feature_key)) +
    ggplot2::geom_col(fill = "#24547A", width = 0.7) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_discrete(labels = labels) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y", labeller = ggplot2::as_labeller(sb_view_wrap_label)) +
    ggplot2::labs(title = "Feature selection frequency", x = "Fraction of successful selection folds", y = NULL,
      caption = paste0("Top ", top_n, " per panel by saved frequency; ties retain input order. Not biological importance. ",
        sum(!view$status$usable), " unavailable entry(s).")) + sb_viz_theme()
  attr(p, "stormbreaker_selection_status") <- view$status
  p
}

#' @rdname sb_plot_selection_frequency
#' @export
sb_plot_selection_jaccard <- function(x, source = c("selector", "glmnet_active_set"), cohort_label = NULL) {
  source <- match.arg(source)
  view <- sb_viz_selection_view(x, source, cohort_label)
  rows <- lapply(view$entries, function(entry) {
    if (!entry$usable) return(NULL)
    tab <- entry$summary$pairwise
    if (!is.data.frame(tab) || !nrow(tab)) return(NULL)
    tab$panel <- entry$panel
    tab$fill_value <- ifelse(tab$status == "success" & is.finite(tab$jaccard), tab$jaccard, NA_real_)
    tab$annotation <- ifelse(is.finite(tab$fill_value), formatC(tab$jaccard, digits = 2, format = "f"), "NA")
    tab
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(sb_viz_empty_selection(view, "Fold-pair selection Jaccard"))
  data <- sb_view_bind(rows)
  fold_i <- fold_j <- fill_value <- annotation <- panel <- NULL
  p <- ggplot2::ggplot(data, ggplot2::aes(x = factor(fold_j), y = factor(fold_i), fill = fill_value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_label(ggplot2::aes(label = annotation), size = 3, fill = "white", alpha = 0.85) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), na.value = "grey90", name = "Jaccard") +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free", labeller = ggplot2::as_labeller(sb_view_wrap_label)) +
    ggplot2::labs(title = "Fold-pair selection Jaccard", x = "Fold j", y = "Fold i",
      caption = "Saved unique fold pairs only; no diagonal or missing pair is invented. Grey is unavailable, not zero.") +
    sb_viz_theme() + ggplot2::theme(panel.grid = ggplot2::element_blank())
  attr(p, "stormbreaker_selection_status") <- view$status
  p
}

#' @rdname sb_plot_selection_frequency
#' @param max_intersections Maximum exact intersections displayed per UpSet
#'   panel, ordered by decreasing size with stable ties. All intersections
#'   remain in the plot attribute. Bars share a scale across panels and their
#'   actual feature counts are labelled.
#' @export
sb_plot_selection_upset <- function(x, source = c("selector", "glmnet_active_set"),
                                     max_intersections = 20L, cohort_label = NULL) {
  source <- match.arg(source)
  max_intersections <- sb_viz_positive_integer(max_intersections, "max_intersections")
  view <- sb_viz_selection_view(x, source, cohort_label)
  pieces <- lapply(view$entries, function(entry) {
    if (!entry$usable) return(NULL)
    records <- entry$records
    if (!is.list(records) || !length(records)) return(NULL)
    successful <- vapply(records, function(record) identical(record$status, "success"), logical(1))
    records <- records[successful]
    if (!length(records)) return(NULL)
    folds <- vapply(records, function(record) as.integer(record$fold), integer(1))
    if (anyNA(folds) || anyDuplicated(folds)) sb_abort("Saved selection folds are missing or duplicated.", "stormbreaker_bad_visualization_input")
    ord <- order(folds)
    records <- records[ord]
    folds <- folds[ord]
    sets <- lapply(records, function(record) unique(as.character(record$selected_features)))
    features <- unique(unlist(sets, use.names = FALSE))
    if (!length(features)) return(NULL)
    membership <- vapply(features, function(feature) paste(folds[vapply(sets, function(set) feature %in% set, logical(1))], collapse = ","), character(1))
    patterns <- unique(membership)
    groups <- lapply(patterns, function(pattern) features[membership == pattern])
    counts <- lengths(groups)
    ord <- order(-counts, seq_along(patterns))
    tab <- data.frame(panel = entry$panel, membership = patterns[ord], count = as.integer(counts[ord]),
      intersection = seq_along(ord), stringsAsFactors = FALSE)
    tab$features <- I(groups[ord])
    tab$shown <- tab$intersection <= max_intersections
    list(intersections = tab, folds = folds)
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(sb_viz_empty_selection(view, "Selection-set intersections"))
  intersections <- sb_view_bind(lapply(pieces, `[[`, "intersections"))
  displayed <- intersections[intersections$shown, , drop = FALSE]
  all_folds <- sort(unique(unlist(lapply(pieces, `[[`, "folds"))))
  n_sets <- length(all_folds)
  baseline <- n_sets + 1.5
  height <- max(3, n_sets * 0.65)
  displayed$bar_bottom <- baseline
  displayed$bar_top <- baseline + height * displayed$count / max(displayed$count)
  dots <- sb_view_bind(lapply(pieces, function(piece) {
    tab <- piece$intersections[piece$intersections$shown, , drop = FALSE]
    sb_view_bind(lapply(seq_len(nrow(tab)), function(i) {
      active <- as.integer(strsplit(tab$membership[[i]], ",", fixed = TRUE)[[1]])
      data.frame(panel = tab$panel[[i]], intersection = tab$intersection[[i]], fold = piece$folds,
        y = match(piece$folds, all_folds), active = piece$folds %in% active)
    }))
  }))
  links <- sb_view_bind(lapply(seq_len(nrow(displayed)), function(i) {
    active <- as.integer(strsplit(displayed$membership[[i]], ",", fixed = TRUE)[[1]])
    data.frame(panel = displayed$panel[[i]], intersection = displayed$intersection[[i]],
      bottom = min(match(active, all_folds)), top = max(match(active, all_folds)))
  }))
  intersection <- y <- active <- panel <- bar_bottom <- bar_top <- count <- bottom <- top <- NULL
  p <- ggplot2::ggplot(dots, ggplot2::aes(x = intersection, y = y)) +
    ggplot2::geom_segment(data = displayed, ggplot2::aes(x = intersection, xend = intersection,
      y = bar_bottom, yend = bar_top), inherit.aes = FALSE, linewidth = 5, colour = "#24547A") +
    ggplot2::geom_text(data = displayed, ggplot2::aes(x = intersection, y = bar_top, label = count),
      inherit.aes = FALSE, vjust = -0.5, size = 3) +
    ggplot2::geom_segment(data = links, ggplot2::aes(x = intersection, xend = intersection, y = bottom, yend = top),
      inherit.aes = FALSE, colour = "#24547A") +
    ggplot2::geom_point(ggplot2::aes(colour = active), size = 2.7) +
    ggplot2::scale_colour_manual(values = c("FALSE" = "grey85", "TRUE" = "#24547A"), guide = "none") +
    ggplot2::scale_y_continuous(breaks = seq_along(all_folds), labels = paste("Fold", all_folds),
      limits = c(0.5, baseline + height + 1)) +
    ggplot2::scale_x_continuous(breaks = sort(unique(displayed$intersection))) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_x", labeller = ggplot2::as_labeller(sb_view_wrap_label)) +
    ggplot2::labs(title = "Selection-set intersections", x = "Exact intersection (size order)", y = NULL,
      caption = paste0("Exact membership among successful folds only; bars label feature counts. ",
        "Never-selected features are outside the union. ", sum(!intersections$shown),
        " intersection(s) omitted from display, retained in attributes.")) +
    sb_viz_theme() + ggplot2::theme(panel.grid = ggplot2::element_blank())
  attr(p, "stormbreaker_selection_status") <- view$status
  attr(p, "stormbreaker_intersections") <- intersections
  p
}

sb_viz_selection_view <- function(x, source, cohort_label) {
  sb_require_curve_plotting()
  leaves <- sb_viz_leaves(x, cohort_label)
  entries <- lapply(leaves, function(leaf) {
    value <- leaf$value
    if (inherits(value, "stormbreaker_external_evaluation")) sb_abort("External evaluations do not contain training-fold selection stability.",
      "stormbreaker_bad_visualization_task")
    summary <- if (source == "selector") value$selection_stability else value$glmnet_active_set_stability
    records <- if (source == "selector") value$selections else lapply(value$coefficients, function(record) {
      list(fold = record$fold, status = record$status, selected_features = record$nonzero_features)
    })
    if (is.null(summary)) summary <- list(status = "not_saved", n_successful_folds = 0L, n_folds = NA_integer_)
    usable <- !summary$status %in% c("not_saved", "not_applicable") && isTRUE(summary$n_successful_folds > 0)
    panel <- paste(sb_viz_panel(leaf), sb_viz_pipeline_id(value), source,
      paste0(summary$n_successful_folds, "/", summary$n_folds, " successful folds | ", summary$status), sep = "\n")
    list(summary = summary, records = records, usable = usable, panel = panel, key = leaf$context$view_id)
  })
  status <- sb_view_bind(lapply(entries, function(entry) data.frame(panel = entry$panel,
    status = entry$summary$status, n_successful_folds = entry$summary$n_successful_folds,
    n_folds = entry$summary$n_folds, usable = entry$usable)))
  list(entries = entries, status = status)
}

sb_viz_empty_selection <- function(view, title) {
  p <- sb_empty_display_plot(title, "No usable saved selection records for this view.",
    "Inspect stormbreaker_selection_status; missing or non-applicable sets are not invented.")
  p$data <- data.frame()
  attr(p, "stormbreaker_selection_status") <- view$status
  p
}
