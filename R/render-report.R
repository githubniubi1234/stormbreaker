sb_report_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  for (pair in list(
    c("&", "&amp;"), c("<", "&lt;"), c(">", "&gt;"),
    c('"', "&quot;"), c("'", "&#39;")
  )) {
    x <- gsub(pair[1], pair[2], x, fixed = TRUE)
  }
  x
}

sb_report_table <- function(tab, columns, caption) {
  columns <- unique(c(intersect(c(
    "repeat_id", "heldout_cohort_id", "model_role",
    "validation_scheme", "view_id"
  ), names(tab)), columns))
  tab <- tab[intersect(columns, names(tab))]
  header <- paste0("<th scope='col'>", sb_report_escape(names(tab)), "</th>", collapse = "")
  rows <- vapply(seq_len(nrow(tab)), function(i) {
    cells <- vapply(tab, function(column) {
      value <- column[i]
      if (is.numeric(value)) value <- format(value, digits = 8, trim = TRUE)
      sb_report_escape(value)
    }, character(1))
    paste0("<tr>", paste0("<td>", cells, "</td>", collapse = ""), "</tr>")
  }, character(1))
  paste0(
    "<div class='table-scroll' tabindex='0' role='region' aria-label='",
    sb_report_escape(caption), "'><table><caption>", sb_report_escape(caption),
    "</caption><thead><tr>", header, "</tr></thead><tbody>",
    paste(rows, collapse = ""), "</tbody></table></div>"
  )
}

sb_report_plot <- function(plot, label, panels) {
  plot <- sb_compact_report_plot(plot)
  path <- tempfile(fileext = ".svg")
  grDevices::svg(path, width = 10, height = max(5, 3 * panels))
  device <- grDevices::dev.cur()
  on.exit(
    {
      if (device %in% grDevices::dev.list()) grDevices::dev.off(device)
    },
    add = TRUE
  )
  print(plot)
  grDevices::dev.off(device)
  svg <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  paste0(
    "<figure><img alt='", sb_report_escape(label),
    "' src='data:image/svg+xml;charset=utf-8,", utils::URLencode(svg, reserved = TRUE),
    "'><figcaption>", sb_report_escape(label), "</figcaption></figure>"
  )
}

sb_report_section <- function(x, label, index) {
  external <- inherits(x, c("stormbreaker_external_evaluation", "stormbreaker_external_comparison"))
  results <- if (inherits(x, "stormbreaker_cv_comparison")) {
    x$results
  } else if (
    inherits(x, "stormbreaker_external_comparison")) {
    x$evaluations
  } else {
    list(x)
  }
  if (!length(results)) sb_abort("Report sections cannot be empty.", "stormbreaker_bad_report_input")
  runs <- sb_result_table(x, "runs", cohort_label = label)
  types <- unique(runs$task_type)
  if (length(types) != 1L || !types %in% c("binary", "survival")) {
    sb_abort("Each section must contain one supported task type.", "stormbreaker_bad_report_input")
  }
  metrics <- if (types == "binary") c("roc_auc", "average_precision", "brier_score", "log_loss") else "harrell_c"
  if (types == "survival") {
    requested <- any(vapply(results, function(result) {
      descriptor <- if (external) result$endpoint else result$task
      !is.null(descriptor$prediction_horizon) ||
        any(sb_horizon_view_metrics() %in% names(result$metrics))
    }, logical(1)))
    if (requested) metrics <- c(metrics, sb_horizon_view_metrics())
  }
  boundary <- if (external) "External evaluation" else "Development CV"
  context <- attr(x, "stormbreaker_view_context")
  if (!is.null(context)) boundary <- sb_view_boundary(context)
  successful <- vapply(results, function(result) identical(result$status, "success"), logical(1))
  successful_ids <- vapply(results[successful], function(result) {
    if (external) result$pipeline_id else result$pipeline$id
  }, character(1))
  runs <- runs[runs$pipeline_id %in% successful_ids & runs$status == "success", , drop = FALSE]
  body <- c(
    paste0(
      "<section id='section-", index, "'><h2>", sb_report_escape(label),
      "</h2><p class='boundary'>", boundary, " &middot; ", types, "</p>"
    ),
    paste0(
      "<p>", sum(successful), " of ", length(results),
      " pipelines completed", if (any(!successful)) paste0("; ", sum(!successful), " failed/incomplete"), ".</p>"
    ),
    if (!any(successful)) "<p><strong>No usable results in this section. This is not a successful analysis.</strong></p>",
    sb_report_table(runs, c(
      "pipeline_id", "task_type", "status", "fold",
      "n", "patient_match_status"
    ), "Successful run status")
  )
  for (metric in metrics) {
    tab <- sb_result_table(x, metric = metric, cohort_label = label)
    tab <- tab[sb_performance_display_rows(tab), , drop = FALSE]
    p <- sb_plot_performance(x, metric, cohort_label = label)
    body <- c(
      body, paste0("<details open><summary>", metric, "</summary>"),
      sb_report_table(tab, c(
        "pipeline_id", "estimate", "status", "run_status",
        "scope", "aggregation", "prediction_horizon",
        "time_unit"
      ), paste("Saved", metric)),
      sb_report_plot(p, paste(boundary, label, metric), max(1, ceiling(nrow(tab) / 8))), "</details>"
    )
  }
  for (kind in c("calibration", "dca")) {
    heading <- if (kind == "calibration") "Calibration" else "Decision curves"
    p <- if (kind == "calibration") {
      sb_plot_calibration(x, cohort_label = label)
    } else {
      sb_plot_dca(x, cohort_label = label)
    }
    status <- attr(p, "stormbreaker_curve_status")
    status_counts <- as.data.frame(table(status$metric_status), stringsAsFactors = FALSE)
    names(status_counts) <- c("saved_metric_status", "pipelines")
    displayed <- unique(as.character(p$data$pipeline_id))
    status <- status[status$pipeline_id %in% displayed, , drop = FALSE]
    body <- c(
      body, paste0("<details open><summary>", heading, "</summary>"),
      sb_report_table(status_counts, names(status_counts), paste(heading, "availability counts")),
      sb_report_table(status, c(
        "pipeline_id", "metric_status", "run_status", "scope",
        "prediction_horizon", "time_unit"
      ), paste(heading, "status")),
      sb_report_plot(p, paste(boundary, label, heading), max(1, nrow(status))), "</details>"
    )
  }
  paste(c(body, "</section>"), collapse = "\n")
}

#' Export a self-contained report from saved results
#'
#' Create an offline HTML report without fitting models or recomputing metrics.
#'
#' @param x One result supported by [sb_result_table()], or a uniquely named
#'   list of such results. Whole repeated CV, IECV and stacking CV objects are
#'   expanded into separate saved-entry sections with their context and role.
#'   IECV is labeled development validation, not final independent validation.
#'   Each ordinary list member is a separate report section;
#'   names are caller-supplied cohort labels, not verified cohort identities.
#' @param file A new `.html` file in an existing directory. A companion
#'   `<stem>_failures.csv` is written in the same directory. Both paths must be
#'   new; existing files are refused and there is no overwrite option.
#' @param title One non-empty report title.
#' @return Invisibly, the normalized absolute output path. Requires optional
#'   `ggplot2` and SVG support in R, but no Pandoc or network access.
#' @details The report contains saved scalar metric and status tables,
#'   performance plots, calibration and decision curves. CV and external
#'   sections remain separate. Pipeline order is preserved; no winner,
#'   confidence interval, missing metric or missing curve is inferred.
#'   Survival sections with a saved horizon or fixed-horizon metric also show
#'   `time_auc`, `uno_c`, and `ipcw_brier`, retaining their horizon and unit.
#'   Horizon-free sections do not request these metrics. In mixed sections,
#'   unavailable values remain unavailable and different targets are labeled
#'   separately, never pooled or converted.
#'   Unsupported or malformed inputs are errors, not silently omitted sections.
#'   Default tables and graphs show usable results only. Failure/issue counts
#'   and a companion CSV link remain visible; all-failed sections explicitly
#'   state that no usable results are available. The CSV is header-only when
#'   there are no issues. See [sb_failure_table()] for its record definitions.
#'
#'   Patient-level rows, raw result objects and free-text error messages are
#'   never embedded. Caller-supplied titles, cohort labels and pipeline labels
#'   remain visible and must be reviewed before sharing. This is not a guarantee
#'   of anonymization: plotted calibration curves describe the supplied cohort.
#'   Saved numeric summaries are formatted to eight significant digits.
#'   The companion CSV includes saved error text, which may be sensitive;
#'   review it separately before sharing. Patient-level rows are not exported.
#'   If an I/O error interrupts the two-file write, a completed CSV may remain;
#'   it is not deleted or overwritten automatically. Inspect the reported paths.
#' @seealso [sb_result_table()], [sb_plot_performance()], [sb_plot_calibration()],
#'   [sb_plot_dca()]
#' @export
sb_render_report <- function(x, file, title = "Saved prediction results") {
  valid_text <- function(value) {
    is.character(value) && length(value) == 1L &&
      !is.na(value) && nzchar(trimws(value))
  }
  if (!valid_text(file) || !grepl("\\.html$", file, ignore.case = TRUE) ||
    !dir.exists(dirname(file)) || !valid_text(title)) {
    sb_abort(
      "Supply a non-empty title and a new .html path in an existing directory.",
      "stormbreaker_bad_report_input"
    )
  }
  if (file.exists(file)) {
    sb_abort(
      "The report path already exists; choose a new filename.",
      "stormbreaker_report_exists"
    )
  }
  failure_file <- sb_new_output_path(sub("\\.html$", "_failures.csv", file,
    ignore.case = TRUE
  ), "csv")
  supported <- function(value) {
    inherits(value, c(
      "stormbreaker_binary_cv_result",
      "stormbreaker_cv_result", "stormbreaker_cv_comparison", "stormbreaker_external_evaluation",
      "stormbreaker_external_comparison", "stormbreaker_repeated_cv_comparison",
      "stormbreaker_internal_external_cv", "stormbreaker_stack_cv_result"
    ))
  }
  if (supported(x)) {
    sections <- list(Result = x)
  } else {
    sections <- x
    if (!is.list(sections) || !length(sections) || is.null(names(sections)) ||
      anyDuplicated(names(sections)) || anyNA(names(sections)) ||
      any(!nzchar(trimws(names(sections)))) || !all(vapply(sections, supported, logical(1)))) {
      sb_abort(
        "Supply a supported result or a non-empty uniquely named list of results.",
        "stormbreaker_bad_report_input"
      )
    }
  }
  failures <- sb_failure_table(sections)
  expanded <- list()
  for (i in seq_along(sections)) {
    if (!sb_is_composite_view(sections[[i]])) {
      expanded <- c(expanded, stats::setNames(list(sections[[i]]), names(sections)[i]))
    } else {
      leaves <- sb_view_leaves(sections[[i]])
      for (leaf in leaves) {
        value <- leaf$value
        attr(value, "stormbreaker_view_context") <- leaf$context
        label <- paste(names(sections)[i], sb_view_label(leaf$context),
          leaf$context$view_id,
          sep = " | "
        )
        expanded <- c(expanded, stats::setNames(list(value), label))
      }
    }
  }
  sections <- expanded
  sb_require_curve_plotting()
  body <- vapply(seq_along(sections), function(i) {
    sb_report_section(sections[[i]], names(sections)[i], i)
  }, character(1))
  nav <- paste0("<li><a href='#section-", seq_along(sections), "'>",
    sb_report_escape(names(sections)), "</a></li>",
    collapse = ""
  )
  html <- paste0(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'><title>",
    sb_report_escape(title), "</title><style>",
    "*{box-sizing:border-box}body{margin:0;background:#f5f7f9;color:#172b3a;font:16px/1.55 system-ui,sans-serif}",
    "main{max-width:1180px;margin:auto;padding:28px}header{border-bottom:4px solid #24547a;padding-bottom:20px}",
    "h1{font-size:2rem;line-height:1.2}h2{margin-bottom:.3em}section{margin:32px 0;padding:22px;background:white}",
    ".boundary{color:#24547a;font-weight:650}details{margin:18px 0;border-top:1px solid #ccd6df;padding-top:12px}",
    "summary{cursor:pointer;font-weight:650;font-size:1.2rem}a{color:#24547a}a:focus-visible,summary:focus-visible,.table-scroll:focus-visible{outline:3px solid #bf6b00}",
    ".table-scroll{overflow-x:auto;max-width:100%}table{border-collapse:collapse;width:100%;font-size:.88rem}",
    "caption{text-align:left;font-weight:bold;padding:10px 0}td,th{padding:8px;text-align:left;border-bottom:1px solid #dbe2e8;overflow-wrap:anywhere}",
    "th{background:#edf2f5}figure{margin:18px 0}img{display:block;width:100%;height:auto}figcaption{font-size:.85rem;color:#465966}",
    "footer{font-size:.85rem;border-top:1px solid #ccd6df;padding:18px 0}@media(max-width:600px){main{padding:12px}section{padding:12px}h1{font-size:1.6rem}td,th{padding:5px}}",
    "@media print{body{background:white}main{max-width:none;padding:0}nav{display:none}details::details-content{content-visibility:visible!important}section{padding:0}figure{break-inside:avoid}.table-scroll{overflow:visible}}",
    "</style></head><body><main><header><h1>", sb_report_escape(title),
    "</h1>",
    "<p><a href='", sb_report_escape(utils::URLencode(basename(failure_file), reserved = TRUE)),
    "' download>Download failures CSV</a>: ", nrow(failures),
    " issue records.</p></header>",
    "<nav aria-label='Report sections'><ol>", nav, "</ol></nav>", paste(body, collapse = "\n"),
    "</main></body></html>"
  )
  if (file.exists(file)) {
    sb_abort(
      "The report path now exists; choose a new filename.",
      "stormbreaker_report_exists"
    )
  }
  failure_file <- sb_new_output_path(failure_file, "csv")
  utils::write.csv(failures, failure_file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  writeLines(enc2utf8(html), file, useBytes = TRUE)
  invisible(normalizePath(file, winslash = "/", mustWork = TRUE))
}
