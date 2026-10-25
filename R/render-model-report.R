sb_model_report_image <- function(plot, label, width = 11, height = 7) {
  path <- tempfile(fileext = ".svg")
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(path, width = width, height = height)
  } else {
    grDevices::svg(path, width = width, height = height)
  }
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
    "<img alt='", sb_report_escape(label),
    "' src='data:image/svg+xml;charset=utf-8,", utils::URLencode(svg, reserved = TRUE), "'>"
  )
}

sb_model_report_table <- function(data) {
  cells <- lapply(data, function(column) {
    if (is.numeric(column)) column <- format(signif(column, 6), trim = TRUE)
    sb_report_escape(as.character(column))
  })
  rows <- vapply(seq_len(nrow(data)), function(i) {
    paste0(
      "<tr>",
      paste0("<td>", vapply(cells, `[`, character(1), i), "</td>", collapse = ""), "</tr>"
    )
  }, character(1))
  paste0("<table><thead><tr>", paste0("<th scope='col'>", sb_report_escape(names(data)),
    "</th>",
    collapse = ""
  ), "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table>")
}

#' Render an offline interactive explanation of a saved linear model
#'
#' @param fit A successful supported linear final fit.
#' @param new_data Feature data for every sample to include; outcome columns
#'   are ignored. All rows receive individual plots, without sampling.
#' @param id_col Unique identifier column in new_data.
#' @param background_data Original-training feature reference, held fixed for
#'   every SHAP view and coefficient-SD calculation. The caller is responsible
#'   for supplying the original training reference, not explained test patients.
#' @param file New HTML filename in an existing directory. Existing files refused.
#' @param development Optional original survival development task, needed only
#'   if the saved model has no full baseline. See [sb_add_survival_curve()].
#' @param title Report title.
#' @param top_n Maximum terms in global SHAP and waterfall displays. Remaining
#'   waterfall terms are aggregated by the existing plotting implementation.
#' @param selection_results Optional saved development CV results for gene
#'   selection rankings. No selection counts are inferred from the final fit.
#' @param compact Remove repeated plot titles, subtitles and long captions.
#' @return Invisibly, the absolute HTML path. This single offline file embeds
#'   R-rendered SVG plots and selectable coefficient/contribution tables.
#' @details Reuses [sb_additive_shap()], [sb_plot_shap_beeswarm()],
#'   [sb_plot_shap_importance()] and [sb_plot_shap_waterfall()]; it does not
#'   introduce a second SHAP algorithm. One sample selector switches both the
#'   waterfall and survival view. SHAP is on log-odds/log-hazard score scales,
#'   never probability effects. The baseline is zero-score, not average-patient
#'   survival. No external service or resource is contacted. Patient IDs are
#'   replaced by row numbers in the report, but plotted values and predictions
#'   remain sensitive; this is not anonymization. Rendering all patient views
#'   can take time and produce a large file. Use the separate plot functions
#'   when an interactive report is not required.
#'   Optional svglite produces compact text-based SVG; otherwise the R SVG
#'   device is used. Both preserve the same plotted data and require no network.
#' @export
sb_render_model_report <- function(fit, new_data, id_col, background_data, file,
                                   development = NULL,
                                   title = "Saved model: parameters and explanations",
                                   top_n = 20L, selection_results = NULL, compact = TRUE) {
  text_ok <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x))
  if (!text_ok(file) || !grepl("\\.html$", file, ignore.case = TRUE) ||
    !dir.exists(dirname(file)) || !text_ok(title)) {
    sb_abort(
      "Supply a title and a new .html path in an existing directory.", "stormbreaker_bad_report_input"
    )
  }
  if (file.exists(file)) {
    sb_abort(
      "The report path already exists; choose a new filename.",
      "stormbreaker_report_exists"
    )
  }
  top_n <- sb_viz_positive_integer(top_n, "top_n")
  if (!is.logical(compact) || length(compact) != 1L || is.na(compact)) {
    sb_abort("compact must be TRUE or FALSE.", "stormbreaker_bad_report_input")
  }
  fit <- sb_validate_final_fit_for_prediction(fit)
  coefficients <- sb_model_coefficients(fit, background_data)
  shap <- sb_additive_shap(fit, new_data, id_col, background_data)
  curves <- NULL
  if (fit$task_type == "survival") {
    if (is.null(fit$survival_curve)) {
      if (is.null(development)) {
        sb_abort(
          "Supply the original development task or an augmented survival model.",
          "stormbreaker_missing_survival_curve"
        )
      }
      fit <- sb_add_survival_curve(fit, development)
    }
    curves <- sb_survival_curves(fit, new_data, id_col)
  }
  image <- function(plot, label, width = 11, height = 7) {
    if (compact) plot <- sb_compact_report_plot(plot)
    sb_model_report_image(plot, label, width, height)
  }
  genes <- if (is.null(selection_results)) "" else sb_model_report_gene_section(selection_results, image)
  shap_height <- max(4.2, 1.7 + 0.30 * min(top_n, ncol(shap$contributions)))
  waterfall_height <- max(4.6, shap_height + 0.6)
  raw <- image(sb_plot_coefficients(coefficients), "Raw signed coefficients")
  scaled <- image(sb_plot_coefficients(coefficients, scale = "reference_sd"), "Coefficients per training SD")
  swarm <- image(sb_plot_shap_beeswarm(shap, top_n, new_data, id_col), "Global SHAP beeswarm", height = shap_height)
  importance <- image(sb_plot_shap_importance(shap, top_n), "Global mean absolute SHAP importance", height = shap_height)
  base <- if (!is.null(curves)) {
    paste0(
      "<details><summary>Baseline survival alone</summary>",
      image(sb_plot_survival_curves(curves, rows = integer()), "Original-training zero-score baseline"),
      "</details>"
    )
  } else {
    ""
  }
  n <- nrow(shap$contributions)
  samples <- vapply(seq_len(n), function(i) {
    summary <- paste0(
      "<p><strong>Observation ", i, "</strong> | Model score: ",
      format(signif(shap$scores$score[i], 6)), " | Background score: ",
      format(signif(shap$reference_value, 6)), "</p>"
    )
    if (is.null(curves)) {
      summary <- paste0(
        summary, "<p>Class-1 probability: ",
        sprintf("%.2f%%", 100 * stats::plogis(shap$scores$score[i])), "</p>"
      )
    } else {
      probability <- curves$survival[curves$baseline$time == curves$prediction_horizon, i]
      summary <- paste0(
        summary, "<p>Survival at ", curves$prediction_horizon, " ",
        sb_report_escape(curves$time_unit), ": <strong>", sprintf("%.2f%%", probability * 100), "</strong></p>"
      )
    }
    contributions <- data.frame(
      Feature = colnames(shap$contributions),
      Contribution = unname(shap$contributions[i, ])
    )
    paste0(
      summary,
      image(sb_plot_shap_waterfall(shap, i, top_n), paste("SHAP waterfall for observation", i), height = waterfall_height),
      "<details><summary>SHAP values</summary>",
      sb_model_report_table(contributions), "</details>",
      if (!is.null(curves)) {
        image(
          sb_plot_survival_curves(curves, rows = i),
          paste("Predicted survival for observation", i)
        )
      } else {
        ""
      }
    )
  }, character(1))
  options <- paste0("<option value='", seq_len(n), "'>Observation ", seq_len(n), "</option>", collapse = "")
  templates <- paste0("<template id='sample-", seq_len(n), "'>", samples, "</template>", collapse = "\n")
  html <- paste0(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'><title>", sb_report_escape(title),
    "</title><style>*{box-sizing:border-box}body{margin:0;background:#f5f7f9;color:#172b3a;font:16px/1.6 system-ui,sans-serif}",
    "main{max-width:1240px;margin:auto;padding:28px}header{border-bottom:4px solid #24547a}h1{line-height:1.2}",
    "h1,h2,p,a,td,th{overflow-wrap:anywhere}section{background:white;margin:26px 0;padding:22px}",
    "img{display:block;width:100%;height:auto}nav{display:flex;gap:20px;flex-wrap:wrap;margin:20px 0}",
    "a{color:#24547a}select,button{font:inherit;padding:7px 12px;max-width:100%;background:white;border:1px solid #8496a2;border-radius:4px}",
    "label{display:inline-block;margin-right:10px;font-weight:600}select:focus-visible,summary:focus-visible,a:focus-visible{outline:3px solid #bf6b00;outline-offset:3px}",
    "details{margin:16px 0}summary{cursor:pointer;color:#24547a}table{border-collapse:collapse;width:100%;font-size:.9rem;table-layout:fixed}",
    "th,td{padding:8px;border-bottom:1px solid #dce3e8;text-align:left}th{background:#edf2f5}",
    "[hidden]{display:none!important}.note,footer{color:#465966;font-size:.92rem}footer{border-top:1px solid #ccd6df;padding:20px 0}",
    "@media(max-width:600px){main{padding:12px}section{padding:10px}h1{font-size:1.6rem}th,td{padding:5px;font-size:.8rem}}",
    "@media print{body{background:white}main{max-width:none;padding:0}.controls,nav{display:none}section{break-inside:avoid;padding:0}details>table{display:table}}",
    "</style></head><body><main><header><h1>", sb_report_escape(title), "</h1>",
    "<p>", n, " samples | ", nrow(coefficients$table), " model features | ",
    sb_report_escape(sb_display_learner(fit$pipeline$learner, fit$pipeline$learner_params)),
    "</p><p>SHAP: ", sb_report_escape(sb_report_score_label(shap$score_scale)),
    " | Training reference: ", shap$background$n, "</p></header>",
    "<nav aria-label='Report sections'>", if (nzchar(genes)) "<a href='#genes'>Gene selection</a>",
    "<a href='#parameters'>Coefficients</a><a href='#global'>Global SHAP</a><a href='#individual'>Individual prediction</a></nav>",
    genes,
    "<section id='parameters'><h2>Model coefficients</h2>",
    if (!is.null(coefficients$intercept)) paste0("<p>Intercept: ", format(signif(coefficients$intercept, 6)), "</p>"),
    "<div class='controls'>",
    "<label for='coefficient-view'>Coefficient units</label><select id='coefficient-view'><option value='raw'>Original input units</option><option value='sd'>Per training standard deviation</option></select></div>",
    "<div id='coef-raw'>", raw, "</div><div id='coef-sd' hidden>", scaled, "</div>",
    "<details><summary>All coefficients</summary>", sb_model_report_table(coefficients$table),
    "</details></section>",
    "<section id='global'><h2>Global SHAP</h2><div class='controls'>",
    "<label for='shap-view'>Global view</label><select id='shap-view'><option value='swarm'>Beeswarm</option><option value='importance'>Mean absolute contribution</option></select></div>",
    "<div id='shap-swarm'>", swarm, "</div><div id='shap-importance' hidden>", importance,
    "</div></section>",
    "<section id='individual'><h2>Individual prediction</h2><div class='controls'><label for='sample-select'>Observation</label>",
    "<select id='sample-select'>", options, "</select></div>",
    "<p id='sample-status' role='status' aria-live='polite'>Showing observation 1 of ", n, "</p><div id='sample-panel'>",
    samples[1], "</div>", base, "</section>", templates,
    "<noscript><p>Enable JavaScript to switch views.</p></noscript></main>",
    "<script>'use strict';const el=id=>document.getElementById(id);",
    "el('coefficient-view').addEventListener('change',e=>{el('coef-raw').hidden=e.target.value!=='raw';el('coef-sd').hidden=e.target.value!=='sd';});",
    "el('shap-view').addEventListener('change',e=>{el('shap-swarm').hidden=e.target.value!=='swarm';el('shap-importance').hidden=e.target.value!=='importance';});",
    if (nzchar(genes)) "el('gene-view').addEventListener('change',e=>{el('gene-selector').hidden=e.target.value!=='selector';el('gene-glmnet_active_set').hidden=e.target.value!=='glmnet_active_set';});",
    "el('sample-select').addEventListener('change',e=>{const t=el('sample-'+e.target.value);el('sample-panel').replaceChildren(t.content.cloneNode(true));",
    "el('sample-status').textContent='Showing observation '+e.target.value+' of ", n, "';});</script></body></html>"
  )
  if (file.exists(file)) {
    sb_abort(
      "The report path already exists; choose a new filename.",
      "stormbreaker_report_exists"
    )
  }
  writeLines(enc2utf8(html), file, useBytes = TRUE)
  invisible(normalizePath(file, winslash = "/", mustWork = TRUE))
}

sb_report_score_label <- function(x) {
  gsub("class_1_log_odds", "log odds (class 1)",
    gsub("zero_reference_log_hazard_ratio", "log hazard ratio", x, fixed = TRUE),
    fixed = TRUE
  )
}

sb_compact_report_plot <- function(plot) {
  xlabel <- plot$labels$x
  plot <- plot + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
  if (is.character(xlabel)) plot <- plot + ggplot2::labs(x = sb_report_score_label(xlabel))
  plot
}

sb_model_report_gene_section <- function(results, image) {
  views <- vapply(c("selector", "glmnet_active_set"), function(source) {
    ranking <- sb_gene_selection_ranking(results, source = source)
    tab <- ranking$ranking
    display <- data.frame(
      Rank = tab$rank, Gene = tab$feature,
      Selected = tab$selected_count, Eligible = tab$eligible_count,
      Frequency = ifelse(is.finite(tab$selection_frequency), sprintf("%.1f%%", 100 * tab$selection_frequency), "NA")
    )
    paste0(
      "<div id='gene-", source, "'", if (source != "selector") " hidden", ">",
      if (ranking$status == "incomplete") "<p><strong>Incomplete selection records</strong></p>",
      if (ranking$status == "not_available") {
        "<p>No saved learned-selection counts.</p>"
      } else {
        image(sb_plot_gene_selection_ranking(ranking), paste("Gene selection ranking", source), height = 7)
      },
      "<details><summary>All genes (", nrow(tab), ")</summary>",
      sb_model_report_table(display), "</details></div>"
    )
  }, character(1))
  paste0(
    "<section id='genes'><h2>Gene selection ranking</h2><div class='controls'>",
    "<label for='gene-view'>Count</label><select id='gene-view'><option value='selector'>Selector selections</option>",
    "<option value='glmnet_active_set'>LASSO / Elastic Net nonzero</option></select></div>",
    "<p>Successful CV folds; shared selector copies count once.</p>", paste(views, collapse = ""), "</section>"
  )
}
