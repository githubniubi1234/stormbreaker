#' Summarize exact intersections of named gene sets
#' @param sets Named list of character gene vectors. Duplicate genes within a
#'   set are counted once; empty sets are retained. Identifiers must already
#'   use the same naming system. No gene annotation or model fitting occurs.
#' @return List with sets, membership, exclusive intersections (including gene
#'   members), and common genes present in every supplied set.
#' @export
sb_gene_overlap <- function(sets) {
  if (!is.list(sets) || length(sets) < 2L || is.null(names(sets)) ||
      anyNA(names(sets)) || anyDuplicated(names(sets)) || any(!nzchar(trimws(names(sets))))) {
    stop("Supply at least two uniquely named gene sets.")
  }
  if (!all(vapply(sets, function(s) is.character(s) && !anyNA(s) &&
      all(nzchar(trimws(s))), logical(1)))) stop("Gene sets must contain nonmissing, nonempty character IDs.")
  sets <- lapply(sets, unique)
  genes <- unique(unlist(sets, use.names = FALSE))
  membership <- matrix(FALSE, length(genes), length(sets), dimnames = list(genes, names(sets)))
  for (j in seq_along(sets)) membership[, j] <- genes %in% sets[[j]]
  keys <- if (length(genes)) apply(membership, 1L, paste0, collapse = ":") else character()
  members <- split(genes, keys)
  combinations <- lapply(members, function(g) names(sets)[membership[g[1L], ]])
  intersections <- data.frame(intersection = vapply(combinations, paste, character(1), collapse = " & "),
    n = lengths(members), row.names = NULL)
  intersections$sets <- unname(combinations)
  intersections$genes <- unname(members)
  intersections <- intersections[order(-intersections$n, intersections$intersection), , drop = FALSE]
  list(sets = sets, membership = membership, intersections = intersections,
    common = Reduce(intersect, sets))
}

#' Plot gene-set overlap as a Venn diagram or intersection bars
#' @param x Result of sb_gene_overlap.
#' @param type auto uses Venn for two to four sets, otherwise bars.
#' @return A plot with exact membership in its gene_overlap attribute. Venn
#'   regions are schematic, not proportional to set sizes. Bars show exclusive
#'   intersections: genes in the named sets and absent from every other set.
#' @export
sb_plot_gene_overlap <- function(x, type = c("auto", "venn", "bar")) {
  type <- match.arg(type)
  x <- sb_gene_overlap(x$sets)
  if (type == "auto") type <- if (length(x$sets) <= 4L) "venn" else "bar"
  if (type == "venn") {
    if (length(x$sets) > 4L) stop("Use bar for more than four sets.")
    if (!requireNamespace("ggvenn", quietly = TRUE)) stop("Install ggvenn for Venn diagrams, or use type='bar'.")
    p <- ggvenn::ggvenn(x$sets, show_percentage = FALSE,
      fill_color = c("#6497B1", "#E69F80", "#8FBF9F", "#B4A7D6"), stroke_size = .5)
  } else {
    sb_require_curve_plotting()
    d <- x$intersections
    d$intersection <- factor(d$intersection, levels = rev(d$intersection))
    intersection <- n <- NULL
    p <- ggplot2::ggplot(d, ggplot2::aes(n, intersection)) +
      ggplot2::geom_col(fill = "#2878A0", width = .65) +
      ggplot2::labs(x = "Genes", y = NULL) + sb_viz_theme()
  }
  attr(p, "gene_overlap") <- x
  p
}
