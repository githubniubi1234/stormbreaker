#' Export a gallery of saved plots to offline HTML
#'
#' @param plots Non-empty, uniquely named list of ggplot objects. Names become
#'   section headings and alternative image text.
#' @param file New `.html` path in an existing directory. Existing files are
#'   refused; there is no overwrite option.
#' @param title Non-empty report title.
#' @param width,height Positive finite plot dimensions in inches, either one
#'   value for all plots or one value per plot in list order.
#' @param groups Optional character vector, one topic per plot in list order.
#'   Topics appear in first-occurrence order as expandable sections. NULL keeps
#'   the flat gallery. Grouping organizes existing plots; it does not compute
#'   missing clinical or algorithm diagnostics.
#' @param expanded Topics initially expanded. NULL expands only the first topic;
#'   character(0) collapses all. Supply all unique groups for complete detail.
#' @param compact Remove plot titles, subtitles and captions from the image and
#'   retain their text in a separate report-body section. Axes, legends, data
#'   and statistical layers stay.
#' @param captions Optional character vector of additional figure legends, one
#'   per plot. Use this for diagnostic interpretation stored separately from
#'   plot labels. All text is escaped, not interpreted as HTML.
#' @param availability Optional data frame of view availability and reasons,
#'   displayed as a selectable table before the gallery.
#' @return Invisibly, the normalized absolute output path. Requires optional
#'   ggplot2 and an R build with SVG support; no Pandoc or internet is needed.
#' @details Embeds rendered SVG images, not R objects or patient-level tables.
#'   Plot titles, labels and visible points remain present, including any
#'   individual values in SHAP plots. This is not an anonymization service:
#'   review all supplied plots before sharing. No metrics are recomputed and
#'   no model is fitted. Rendered images are selectable only as images; inspect
#'   the original ggplot data locally for numeric values. Device state is
#'   restored even if drawing fails. A failed output write may leave a partial
#'   file, which is not automatically removed or overwritten.
#'   Grouped sections use native HTML disclosure controls and work without
#'   JavaScript. Every supplied image remains embedded, including collapsed
#'   details; collapsing is not privacy protection or file-size reduction.
#'   Expand the sections you want before printing. Extra legends are selectable
#'   text in the report body. Supply `compact = FALSE` to keep original labels.
#' @seealso [sb_render_report()], [sb_plot_roc()], [sb_plot_meta_forest()],
#'   [sb_plot_selection_upset()], [sb_plot_shap_waterfall()]
#' @export
sb_render_visual_report <- function(plots, file,
                                    title = "Saved-result visualization gallery",
                                    width = 11, height = 7, groups = NULL,
                                    expanded = NULL, compact = TRUE, captions = NULL,
                                    availability = NULL) {
  text_ok <- function(x) {
    is.character(x) && length(x) == 1L &&
      !is.na(x) && nzchar(trimws(x))
  }
  if (!is.list(plots) || !length(plots) || is.null(names(plots)) ||
    anyNA(names(plots)) || anyDuplicated(names(plots)) ||
    any(!nzchar(trimws(names(plots)))) ||
    !all(vapply(plots, inherits, logical(1), what = "ggplot"))) {
    sb_abort(
      "Supply a non-empty uniquely named list of ggplot objects.",
      "stormbreaker_bad_report_input"
    )
  }
  if (!text_ok(file) || !grepl("\\.html$", file, ignore.case = TRUE) ||
    !dir.exists(dirname(file)) || !text_ok(title)) {
    sb_abort(
      "Supply a title and a new .html path in an existing directory.",
      "stormbreaker_bad_report_input"
    )
  }
  dimensions_ok <- function(x) {
    is.numeric(x) && length(x) %in% c(1L, length(plots)) &&
      all(is.finite(x) & x > 0)
  }
  if (!dimensions_ok(width) || !dimensions_ok(height)) {
    sb_abort(
      "Plot dimensions must be positive finite numbers, length one or number of plots.",
      "stormbreaker_bad_report_input"
    )
  }
  if (!is.logical(compact) || length(compact) != 1L || is.na(compact)) sb_abort(
    "compact must be TRUE or FALSE.", "stormbreaker_bad_report_input")
  vector_ok <- function(x) is.character(x) && length(x) == length(plots) && !anyNA(x)
  if (!is.null(groups) && (!vector_ok(groups) || any(!nzchar(trimws(groups))))) sb_abort(
    "Supply one non-empty group name per plot.", "stormbreaker_bad_report_input")
  if (!is.null(captions) && !vector_ok(captions)) sb_abort(
    "Supply one caption per plot; use an empty string when none is needed.", "stormbreaker_bad_report_input")
  topics <- unique(groups)
  if (!is.null(availability) && (!is.data.frame(availability) || !ncol(availability) ||
      !all(vapply(availability, is.atomic, logical(1))))) sb_abort(
    "Availability must be a data frame with atomic columns.", "stormbreaker_bad_report_input")
  if (!is.null(expanded) && (!is.character(expanded) || anyNA(expanded) ||
      anyDuplicated(expanded) || !all(expanded %in% topics))) sb_abort(
    "Expanded topics must be unique supplied group names.", "stormbreaker_bad_report_input")
  if (is.null(expanded)) expanded <- utils::head(topics, 1L)
  refuse_existing <- function() {
    if (file.exists(file)) {
      sb_abort(
        "The report path already exists; choose a new filename.",
        "stormbreaker_report_exists"
      )
    }
  }
  refuse_existing()
  sb_require_curve_plotting()
  width <- rep(width, length.out = length(plots))
  height <- rep(height, length.out = length(plots))
  notes <- vector("list", length(plots))
  draw <- function(i) {
    path <- tempfile(fileext = ".svg")
    grDevices::svg(path, width = width[i], height = height[i])
    device <- grDevices::dev.cur()
    on.exit(
      {
        if (device %in% grDevices::dev.list()) grDevices::dev.off(device)
      },
      add = TRUE
    )
    plot <- plots[[i]]
    legend <- captions[i]
    if (compact) {
      original_labels <- plot$labels[c("title", "subtitle", "caption")]
      label_text <- vapply(original_labels, function(z) paste(as.character(z), collapse = " "), character(1))
      legend <- unique(c(label_text[label_text != names(plots)[i]], legend))
      plot <- plot + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
    }
    legend <- legend[!is.na(legend) & nzchar(legend)]
    notes[[i]] <<- legend
    print(plot)
    grDevices::dev.off(device)
    svg <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    label <- sb_report_escape(names(plots)[i])
    paste0(
      "<section id='plot-", i, "'><h2>", label,
      "</h2><figure><img alt='", label, "' src='data:image/svg+xml;charset=utf-8,",
      utils::URLencode(svg, reserved = TRUE), "'></figure></section>"
    )
  }
  sections <- vapply(seq_along(plots), draw, character(1))
  nav <- paste0("<li><a href='#plot-", seq_along(plots), "'>",
    sb_report_escape(names(plots)), "</a></li>",
    collapse = ""
  )
  if (length(topics)) {
    sections <- vapply(seq_along(topics), function(i) {
      indices <- which(groups == topics[i])
      paste0("<details id='topic-", i, "'", if (topics[i] %in% expanded) " open" else "",
        "><summary>", sb_report_escape(topics[i]), " <span>(", length(indices), ")</span></summary>",
        paste(sections[indices], collapse = "\n"), "</details>")
    }, character(1))
    nav <- paste0("<li><a href='#topic-", seq_along(topics), "'>",
      sb_report_escape(topics), "</a></li>", collapse = "")
  }
  html <- paste0(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'><title>",
    sb_report_escape(title), "</title><style>",
    "*{box-sizing:border-box}body{margin:0;color:#172b3a;background:#f5f7f9;font:16px/1.55 system-ui,sans-serif}",
    "main{max-width:1280px;margin:auto;padding:28px}header{border-bottom:4px solid #24547a}h1{line-height:1.2}",
    "section{background:white;padding:20px;margin:24px 0}h1,h2,a{overflow-wrap:anywhere}",
    "details{background:white;margin:24px 0;border:1px solid #ccd6df}summary{cursor:pointer;padding:18px;font-size:1.2rem;font-weight:600;overflow-wrap:anywhere}",
    "summary span{font-weight:400;color:#465966}summary:focus-visible{outline:3px solid #bf6b00}details section{margin:0;border-top:1px solid #e2e8ed}",
    "figure{margin:0}img{display:block;width:100%;height:auto}figcaption,footer{font-size:.9rem;color:#465966}figcaption{white-space:pre-line;overflow-wrap:anywhere}",
    "table{border-collapse:collapse;width:100%;font-size:.9rem}td,th{padding:8px;text-align:left;border-bottom:1px solid #ccd6df;overflow-wrap:anywhere}",
    "a{color:#24547a}a:focus-visible{outline:3px solid #bf6b00}footer{border-top:1px solid #ccd6df;padding:18px 0}",
    "@media(max-width:600px){main{padding:12px}section{padding:8px}h1{font-size:1.6rem}}",
    "@media print{body{background:white}main{padding:0;max-width:none}nav{display:none}section{padding:0}figure{break-inside:avoid}}",
    "</style></head><body><main><header><h1>", sb_report_escape(title),
    "</h1></header>",
    "<nav aria-label='Gallery sections'><ol>", nav, "</ol></nav>",
    if (!is.null(availability)) paste0("<details><summary>View availability</summary>", sb_model_report_table(availability), "</details>") else "",
    paste(sections, collapse = "\n"),
    if (any(lengths(notes))) paste0("<section id='report-notes'><h2>Methods and interpretation</h2>",
      paste(vapply(which(lengths(notes) > 0L), function(i) paste0("<h3>",
        sb_report_escape(names(plots)[i]), "</h3><p>",
        paste(sb_report_escape(unique(notes[[i]])), collapse = "<br>"), "</p>"), character(1)), collapse = "\n"),
      "</section>") else "",
    "</main></body></html>"
  )
  refuse_existing()
  writeLines(enc2utf8(html), file, useBytes = TRUE)
  invisible(normalizePath(file, winslash = "/", mustWork = TRUE))
}
