sb_viz_result_input <- function(x) {
  if (inherits(x, "stormbreaker_run")) return(x$comparison)
  if (inherits(x, "stormbreaker_completed_run")) {
    if (!length(x$evaluations)) sb_abort(
      "No saved external evaluations; use the original run for development CV views.",
      "stormbreaker_bad_visualization_input")
    return(x$evaluations)
  }
  x
}

sb_viz_leaves <- function(x, cohort_label = NULL) {
  x <- sb_viz_result_input(x)
  label <- sb_result_cohort_label(cohort_label)
  if (is.list(x) && !is.object(x)) {
    nms <- names(x)
    if (!length(x) || is.null(nms) || anyNA(nms) || any(!nzchar(trimws(nms))) || anyDuplicated(nms)) {
      sb_abort("Supply a non-empty list with unique, non-blank cohort labels.",
        "stormbreaker_bad_visualization_input")
    }
    if (!is.null(cohort_label)) sb_abort("A named list already supplies cohort labels; leave cohort_label NULL.",
      "stormbreaker_bad_visualization_input")
    return(unlist(lapply(seq_along(x), function(i) {
      leaves <- sb_view_leaves(x[[i]], path = paste0("dataset_", i))
      lapply(leaves, function(leaf) { leaf$dataset <- nms[[i]]; leaf })
    }), recursive = FALSE))
  }
  context <- attr(x, "stormbreaker_view_context")
  leaves <- if (is.null(context)) sb_view_leaves(x) else list(list(value = x, context = context))
  lapply(leaves, function(leaf) {
    leaf$dataset <- if (is.na(label)) "Unnamed cohort" else label
    leaf
  })
}

sb_viz_panel <- function(leaf) {
  paste(leaf$dataset, sb_view_boundary(leaf$context), sb_view_label(leaf$context), sep = "\n")
}

sb_viz_pipeline_id <- function(value) {
  id <- if (inherits(value, "stormbreaker_external_evaluation")) value$pipeline_id else value$pipeline$id
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) {
    sb_abort("A saved visualization entry has no pipeline ID.", "stormbreaker_bad_visualization_input")
  }
  id
}

sb_viz_score_label <- function(scale) {
  switch(scale, class_1_log_odds = "log odds", log_odds = "log odds",
    zero_reference_log_hazard_ratio = "log HR",
    negative_log_time_location = "-log time",
    class_1_probability = "probability", scale)
}

sb_viz_theme <- function() {
  ggplot2::theme_minimal(base_size = 11, base_family = "sans") +
    ggplot2::theme(plot.caption = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 12, face = "bold", colour = "#222222",
        margin = ggplot2::margin(b = 8)),
      axis.text = ggplot2::element_text(colour = "#333333"),
      axis.title = ggplot2::element_text(colour = "#222222"),
      strip.text = ggplot2::element_text(hjust = 0, colour = "#333333"),
      panel.grid = ggplot2::element_line(colour = "#ECECEC", linewidth = .25),
      panel.spacing = grid::unit(10, "pt"),
      plot.margin = ggplot2::margin(9, 12, 9, 9),
      panel.grid.minor = ggplot2::element_blank())
}

sb_viz_wrap <- function(x, width = 105L) {
  paste(strwrap(x, width = width), collapse = "\n")
}

sb_viz_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != floor(x)) {
    sb_abort(paste0("`", name, "` must be one positive integer."),
      "stormbreaker_bad_visualization_argument")
  }
  as.integer(x)
}
