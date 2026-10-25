#' Assemble clinical views from saved evaluation results
#' @param x Saved CV or external evaluation results, or a supported collection.
#'   Results from [bi_run()] and [surv_run()] display development CV;
#'   [sb_complete_run()] results display saved external evaluations by cohort.
#' @param thresholds Optional explicit binary probability thresholds. No optimum
#'   is selected; saved decision curves keep their own original thresholds.
#' @param confusion_threshold Optional single threshold for the confusion matrix.
#'   It must be included in thresholds.
#' @param native Optional [sb_native_diagnostics()] result or uniquely named list
#'   of extracted algorithm diagnostics.
#'   These are already fitted-model summaries; no extraction refit is triggered.
#' @param shap Optional saved [sb_additive_shap()] result. Adds global importance
#'   and beeswarm views; never recomputes contributions or creates a surrogate.
#' @param feature_values,id_col Optional feature frame and explicit ID column for
#'   SHAP colouring and dependence. IDs must match explained observations exactly.
#' @param dependence_features Explicit explained feature names for dependence
#'   plots. NULL adds none; no automatic mass per-feature export.
#' @param nomogram Optional exact [sb_nomogram()] result. No new model is fitted.
#' @param detail Default `"brief"` requests ROC, calibration and saved decision
#'   curves only. `"full"` additionally requests precision-recall and probability
#'   distributions. Explicit diagnostics and explanations remain opt-in in both.
#' @param selectors Optional uniquely named list of saved
#'   [sb_selector_diagnostics()] results, shown separately from model evidence.
#' @return A list with named plots, topic groups, separate captions, an availability
#'   table and original diagnostic tables. Unavailable views are recorded, not
#'   replaced with empty figures or fabricated metrics. Partial-context plots
#'   retain their original fold/cohort diagnostics. No model is fitted or selected.
#' @details ROC, precision-recall and probability distributions require binary
#'   predictions; survival inputs retain their saved calibration and decision
#'   curves. Unsupported binary views are recorded as unavailable. Missing DCA
#'   is not recomputed. SHAP and nomograms must be supplied from their exact-model
#'   extractors; no patient card or subgroup module is added. Individual SHAP
#'   values remain sensitive even though original identifiers are not plotted.
#'   Technical diagnostics are opt-in: call an extractor separately and supply
#'   its result through native or selectors. Report assembly never discovers
#'   diagnostics from a model, including when all report topics are expanded.
#'   Figures use readable display labels; original IDs and status remain in
#'   their data and diagnostics. Model settings are kept separately in models.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' development <- bi_example(40, seed = 42)
#' run <- bi_run(development, outcome = "outcome", id = "patient_id",
#'   predictors = c("gene_a", "gene_b"), models = "Logistic", folds = 2,
#'   seed = 42, plots = FALSE)
#' views <- clinical_plots(run)
#' names(views$plots)
#' views$plots$ROC
#' @export
sb_clinician_views <- function(x, thresholds = NULL, confusion_threshold = NULL, native = list(),
                               shap = NULL, feature_values = NULL, id_col = NULL,
                               dependence_features = NULL, nomogram = NULL, selectors = list(),
                               detail = c("brief", "full")) {
  detail <- match.arg(detail)
  sb_require_curve_plotting()
  if (!is.list(selectors) || (length(selectors) && (is.null(names(selectors)) ||
      anyNA(names(selectors)) || anyDuplicated(names(selectors)) || any(!nzchar(trimws(names(selectors)))) ||
      !all(vapply(selectors, inherits, logical(1), what = "stormbreaker_selector_diagnostic"))))) sb_abort(
    "Supply a uniquely named list of extracted selector diagnostics.", "stormbreaker_bad_report_input")
  native_status <- NULL
  if (inherits(native, "stormbreaker_native_diagnostics")) {
    native_status <- native$status
    native <- native$diagnostics
  }
  if (!is.list(native) || (length(native) && (is.null(names(native)) || anyNA(names(native)) ||
      anyDuplicated(names(native)) || any(!nzchar(trimws(names(native))))))) sb_abort(
    "Native diagnostics must be a uniquely named list.", "stormbreaker_bad_report_input")
  if (length(native) && !all(vapply(native, inherits, logical(1), what = "stormbreaker_algorithm_diagnostic"))) sb_abort(
    "Supply extracted algorithm diagnostics, not fits or arbitrary objects.", "stormbreaker_bad_report_input")
  if (!is.null(thresholds) && (!is.numeric(thresholds) || !length(thresholds) ||
      any(!is.finite(thresholds)) || any(thresholds < 0 | thresholds > 1))) sb_abort(
    "Supply explicit thresholds in [0, 1].", "stormbreaker_bad_report_input")
  if (!is.null(confusion_threshold) && (!is.numeric(confusion_threshold) ||
      length(confusion_threshold) != 1L || !is.finite(confusion_threshold) ||
      !confusion_threshold %in% thresholds)) sb_abort(
    "The confusion threshold must be one of the supplied thresholds.", "stormbreaker_bad_report_input")
  # Validate the result container before treating per-view unavailability as data.
  models <- sb_clinician_models(sb_viz_leaves(x, NULL))
  if (!is.null(nomogram) && !inherits(nomogram, "stormbreaker_nomogram")) sb_abort(
    "Supply an exact sb_nomogram result.", "stormbreaker_bad_report_input")
  if (!is.null(shap)) sb_viz_shap_matrix(shap)
  if (!is.null(dependence_features) && (is.null(shap) || !is.character(dependence_features) ||
      !length(dependence_features) || anyNA(dependence_features) || anyDuplicated(dependence_features) ||
      !all(dependence_features %in% colnames(shap$contributions)) || is.null(feature_values) || is.null(id_col))) sb_abort(
    "Dependence plots need saved SHAP, unique explained feature names, feature values and an ID column.",
    "stormbreaker_bad_report_input")
  plots <- groups <- captions <- list()
  widths <- heights <- numeric()
  diagnostics <- list()
  status <- list()
  add <- function(label, group, draw, caption = "", height = 7) {
    value <- tryCatch(draw(), error = identity)
    if (inherits(value, "error")) {
      unavailable <- inherits(value, c("stormbreaker_bad_visualization_task", "stormbreaker_unavailable_diagnostic",
        "stormbreaker_missing_diagnostic_artifacts", "stormbreaker_unsupported_diagnostic"))
      status[[length(status) + 1L]] <<- data.frame(view = label, group = group,
        status = if (unavailable) "unavailable" else "error", reason = conditionMessage(value))
      diagnostics[[label]] <<- list(error_class = class(value), message = conditionMessage(value))
      return(invisible(NULL))
    }
    attrs <- attributes(value)
    notes <- unique(c(caption, attr(value, "stormbreaker_notes"), attr(value, "stormbreaker_figure_notes"),
      value$labels$subtitle, value$labels$caption))
    caption <- paste(notes[!is.na(notes) & nzchar(notes)], collapse = "\n")
    caption <- gsub("0 unavailable curve group(s). ", "", caption, fixed = TRUE)
    lines <- strsplit(caption, "\n", fixed = TRUE)[[1L]]
    caption <- paste(lines[!grepl("^0 saved coordinate\\(s\\) not drawn;", lines)], collapse = "\n")
    value <- value + ggplot2::labs(caption = NULL, subtitle = NULL)
    value <- sb_clinician_clean_plot(value, models)
    details <- attrs[grep("^stormbreaker_.*(status|diagnostics)$", names(attrs))]
    diagnostics[[label]] <<- details
    usable <- is.data.frame(value$data) && nrow(value$data) > 0L
    if (usable && "draw" %in% names(value$data)) usable <- any(value$data$draw %in% TRUE)
    issues <- unique(unlist(lapply(details, function(d) {
      if (!is.data.frame(d)) return(character())
      fields <- intersect(c("curve_status", "metric_status", "run_status", "status"), names(d))
      values <- as.character(unlist(d[fields], use.names = FALSE))
      values[!is.na(values) & values != "success"]
    }), use.names = FALSE))
    reason <- if (length(issues)) paste(issues, collapse = "; ") else
      if (usable) "" else "No usable coordinates; see original diagnostic tables."
    status[[length(status) + 1L]] <<- data.frame(view = label, group = group,
      status = if (!usable) "unavailable" else if (length(issues)) "partial" else "available", reason = reason)
    if (usable) {
      plots[[label]] <<- value
      groups[[label]] <<- group
      captions[[label]] <<- caption
      dimensions <- attr(value, "stormbreaker_suggested_size")
      widths[[label]] <<- if (is.null(dimensions)) 11 else dimensions[["width"]]
      heights[[label]] <<- if (is.null(dimensions)) height else dimensions[["height"]]
    }
  }
  add("ROC", "Performance", function() sb_plot_roc(x))
  if (detail == "full") add("Precision-recall", "Performance", function() sb_plot_pr(x))
  add("Calibration", "Performance", function() sb_plot_calibration(x))
  add("Decision curve", "Clinical decisions", function() sb_plot_dca(x))
  if (detail == "full") add("Probability distribution", "Clinical decisions", function() sb_plot_probability_distribution(sb_probability_distribution(x)))
  if (!is.null(thresholds)) {
    add("Threshold performance", "Clinical decisions", function() sb_plot_thresholds(sb_threshold_metrics(x, thresholds)))
    if (!is.null(confusion_threshold)) add("Confusion matrix", "Clinical decisions",
      function() sb_plot_confusion(sb_threshold_metrics(x, thresholds), confusion_threshold))
  }
  for (label in names(native)) {
    q <- native[[label]]
    add(paste0("Model: ", label), "Model diagnostics", function() sb_plot_algorithm_diagnostic(q),
      paste(c(q$source, q$interpretation), collapse = "\n"))
  }
  for (label in names(selectors)) {
    q <- selectors[[label]]
    add(paste0("Selector: ", label), "Feature selection", function() sb_plot_selector_diagnostic(q),
      q$source, height = if (q$type %in% c("vsurf_importance", "boruta_summary", "boruta_hits"))
        max(5, 2 + .3 * nrow(q$table)) else 7)
  }
  if (!is.null(native_status)) {
    diagnostics[["Native extraction"]] <- native_status
    missing <- native_status[native_status$status != "available", , drop = FALSE]
    if (nrow(missing)) status[[length(status) + 1L]] <- data.frame(
      view = paste0("Model: ", missing$view), group = "Model diagnostics",
      status = missing$status, reason = missing$reason)
  }
  if (!is.null(shap)) {
    height <- max(4.2, 1.8 + .32 * min(20L, ncol(shap$contributions)))
    add("SHAP importance", "Model explanations", function() sb_plot_shap_importance(shap), height = height)
    add("SHAP beeswarm", "Model explanations",
      function() sb_plot_shap_beeswarm(shap, feature_values = feature_values, id_col = id_col), height = height)
    for (feature in dependence_features) add(paste0("SHAP dependence: ", feature), "Model explanations",
      function() sb_plot_shap_dependence(shap, feature, feature_values, id_col), height = 5)
  }
  if (!is.null(nomogram)) add("Original-model nomogram", "Model explanations",
    function() sb_plot_nomogram(nomogram), height = max(5, 2 + .45 * nrow(nomogram$features)))
  structure(list(plots = plots, groups = unlist(groups, use.names = FALSE),
    captions = unlist(captions, use.names = FALSE), status = do.call(rbind, status),
    diagnostics = diagnostics, models = models, width = unname(widths), height = unname(heights)), class = "stormbreaker_clinician_views")
}
