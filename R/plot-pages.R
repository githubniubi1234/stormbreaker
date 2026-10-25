#' Paginate algorithm plots without dropping observations
#' @param plots A ggplot or named list of ggplots, such as sb_plots(fit)$plots.
#' @param panels_per_page Maximum number of query or feature facets on each page.
#' @param observations_per_page Maximum number of labelled query observations
#'   on a plot without query facets.
#' @return Named list of ggplots. All query rows and corresponding layer rows
#'   are retained across pages. Explanatory notes remain in the
#'   stormbreaker_notes attribute, separate from the visible panel.
#' @export
sb_plot_pages <- function(plots, panels_per_page = 6L, observations_per_page = 25L) {
  if (inherits(plots, "ggplot")) plots <- list(plot = plots)
  for (n in list(panels_per_page, observations_per_page)) {
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1 || n != floor(n))
      stop("Page capacities must be positive integers.", call. = FALSE)
  }
  if (length(plots) && is.null(names(plots))) names(plots) <- paste0("plot", seq_along(plots))
  out <- list()
  for (name in names(plots)) {
    p <- plots[[name]]
    d <- p$data
    facet_vars <- unlist(lapply(p$facet$params$facets,
      function(q) all.vars(rlang::get_expr(q))))
    variable <- if (is.data.frame(d) && "query" %in% names(d)) "query" else
      if (is.data.frame(d) && "feature" %in% names(d) && identical(unname(facet_vars), "feature")) "feature" else
      if (is.data.frame(d) && "feature" %in% names(d) && length(unique(d$feature)) > 25L &&
          "feature" %in% all.vars(rlang::get_expr(p$mapping$y))) "feature" else NULL
    if (!inherits(p, "ggplot") || is.null(variable)) {
      out[[name]] <- p
      next
    }
    facet_query <- variable %in% facet_vars
    cap <- if (variable == "feature" && facet_query && identical(p$facet$params$ncol, 1)) 4L else
      if (facet_query) panels_per_page else observations_per_page
    ids <- unique(d[[variable]])
    batches <- split(ids, ceiling(seq_along(ids) / cap))
    for (i in seq_along(batches)) {
      keep <- batches[[i]]
      subset_rows <- function(z) {
        if (is.data.frame(z) && variable %in% names(z)) {
          z <- z[z[[variable]] %in% keep, , drop = FALSE]
          if (is.factor(z[[variable]])) z[[variable]] <- droplevels(z[[variable]])
        }
        z
      }
      page <- p
      page$data <- subset_rows(d)
      page$layers <- lapply(p$layers, function(layer)
        ggplot2::ggproto(NULL, layer, data = subset_rows(layer$data)))
      if (facet_query) {
        free <- p$facet$params$free
        scales <- if (isTRUE(free$x) && isTRUE(free$y)) "free" else
          if (isTRUE(free$x)) "free_x" else if (isTRUE(free$y)) "free_y" else "fixed"
        page <- page + ggplot2::facet_wrap(stats::reformulate(variable),
          ncol = if (variable == "feature" && !is.null(p$facet$params$ncol)) p$facet$params$ncol else 2, scales = scales,
          labeller = if (variable == "query") ggplot2::as_labeller(function(z) paste("Sample", z)) else ggplot2::label_value)
      }
      attr(page, "stormbreaker_notes") <- attr(p, "stormbreaker_notes")
      key <- if (length(batches) == 1L) name else paste0(name, " / ", i)
      out[[key]] <- page
    }
  }
  out
}
