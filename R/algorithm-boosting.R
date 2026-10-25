sb_mboost_diagnostic_ids <- function() {
  c("gamboost_binary", "blackboost_binary",
    sb_mboost_cox_learner_ids(), sb_mboost_aft_learner_ids(),
    sb_mboost_ranking_learner_ids())
}

#' Inspect saved boosting loss and component selection
#' @param fit A final supported mboost, GBM, XGBoost, LightGBM, CatBoost or
#'   CoxBoost, AdaBoost, LogitBoost or hinge-boosting model. Saved state is
#'   sufficient except binary GBM training loss, CoxBoost paths and hinge
#'   iteration diagnostics require diagnostic capture.
#' @param type training_loss or component_usage for mboost; training_loss,
#'   cv_curve or importance for GBM where recorded; importance for XGBoost;
#'   importance/tree_size for LightGBM and importance/interaction for CatBoost.
#'   For CoxBoost use coefficient_path, nonzero_path or cv_curve.
#'   AdaBoost/LogitBoost support split_usage; AdaBoost also supports tree_weights
#'   and tree_size, while LogitBoost supports stump_sequence.
#' @param max_pairs Maximum native CatBoost interaction pairs displayed (20 by
#'   default). Only the plot is limited; the complete reported table is retained.
#' @return An algorithm diagnostic. training_loss retains the native objective
#'   at iterations zero through the original stopping iteration. Loss scale and
#'   sign are unchanged, including the negative concordance objective. Values
#'   from different families or samples are not directly comparable. This is
#'   in-bag training loss, not an external validation or early-stopping curve.
#'
#'   component_usage retains every original component, including unused terms
#'   and an intercept when present, and a per-iteration selection table. Counts
#'   describe how often a component was updated, not its coefficient magnitude,
#'   predictive importance or across-resample stability. Internal feature aliases
#'   are mapped back to their saved input names. BlackBoost updates a single
#'   tree base learner; it is not given a misleading variable-selection plot.
#'   Neither diagnostic changes stopping, fits, coefficients or predictions.
#'
#'   GBM keeps final-fit training errors separate from its original CV search.
#'   XGBoost returns native gain/cover/frequency shares for tree boosters, or
#'   signed coefficients for linear boosters (bias excluded). AFT coefficients
#'   are on the log-time scale, not the package's reversed risk scale. Missing
#'   zero-use feature rows are restored as zeros, with the native table retained.
#'   Original pipelines without evaluation logs do not get invented curves.
#'   LightGBM retains unscaled gain/cover/split counts and counts actual saved
#'   trees, which can be fewer than configured rounds. CatBoost interaction is
#'   a model-internal strength score, not SHAP interaction or a causal effect.
#'   Native zero-based feature indices are mapped to the saved input names;
#'   unreported pairs are not assigned fabricated zero scores.
#'   A single input has no distinct feature pairs and returns an empty table
#'   without calling the native pair extractor. Native extraction failures for
#'   other models remain explicit unavailable diagnostics with their cause.
#'   CoxBoost supports coefficient_path and nonzero_path with diagnostic
#'   capture, and cv_curve from saved tuning for its CV variant without capture.
#'   Coefficients retain the original prediction adapter's scale conversion;
#'   the CV curve is conditional on the previously selected penalty.
#'   AdaBoost and LogitBoost expose split_usage without capture. AdaBoost also
#'   exposes tree_weights and tree_size for actual fitted trees, including
#'   native convergence before the configured maximum. LogitBoost exposes its
#'   stump_sequence with exact native thresholds and signed log-odds votes.
#'   These are model-specific split rules, not clinical decision thresholds.
#'   Hinge boosting supports training_loss, component_usage and coefficient_path
#'   from newly captured original iterations. Native losses start after update
#'   one; no initial loss is invented. Coefficient paths include the initial
#'   margin equation and its centering-adjusted intercept, not probability or
#'   log-odds coefficients. Historical compact fits cannot supply these paths.
#' @export
sb_boosting_diagnostics <- function(fit, type = "training_loss", max_pairs = 20L) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  learner <- fit$pipeline$learner
  if (learner == "bst_hinge_binary") return(sb_hinge_boosting_diagnostic(fit, type))
  if (learner %in% sb_classical_binary_boosting_learner_ids()) {
    return(sb_classical_boosting_diagnostic(fit, type))
  }
  if (learner %in% c("coxboost", "coxboost_cv")) {
    return(sb_coxboost_diagnostic(fit, type))
  }
  if (learner %in% c(sb_lightgbm_binary_learner_ids(), sb_catboost_diagnostic_ids())) {
    return(sb_light_cat_diagnostic(fit, type, max_pairs))
  }
  if (learner %in% c(sb_gbm_diagnostic_ids(), sb_xgboost_learner_ids())) {
    return(sb_gbm_xgboost_diagnostic(fit, type))
  }
  if (!learner %in% sb_mboost_diagnostic_ids()) sb_abort(
    "This learner has no supported native boosting diagnostic adapter.",
    "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("training_loss", "component_usage")) sb_abort(
    "Choose training_loss or component_usage.", "stormbreaker_unsupported_diagnostic")
  if (!requireNamespace("mboost", quietly = TRUE)) sb_abort(
    "Package mboost is required to read its saved diagnostic objects.",
    "stormbreaker_missing_dependency")
  state <- fit$learner_state
  engine <- unserialize(state$engine_raw)
  if (!inherits(engine, "mboost")) sb_abort("The original mboost engine is missing.",
    "stormbreaker_invalid_diagnostic")
  stop <- mboost::mstop(engine)
  if (!is.numeric(stop) || length(stop) != 1L || !is.finite(stop) ||
      stop < 1 || stop != floor(stop)) sb_abort("The original stopping iteration is invalid.",
    "stormbreaker_invalid_diagnostic")
  selections <- NULL
  if (type == "training_loss") {
    loss <- mboost::risk(engine)
    if (!is.numeric(loss) || length(loss) != stop + 1L || any(!is.finite(loss))) sb_abort(
      "The native in-bag loss sequence is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(iteration = seq.int(0L, stop), loss = unname(loss))
    measure <- paste("Native training objective:", engine$family@name)
  } else {
    if (startsWith(learner, "blackboost_")) sb_abort(
      "BlackBoost updates one tree base learner; component counts are not a variable-selection diagnostic.",
      "stormbreaker_unavailable_diagnostic")
    native_variables <- stats::variable.names(engine, usedonly = FALSE)
    selected <- mboost::selected(engine)
    if (!length(native_variables) || anyNA(native_variables) ||
        !is.numeric(selected) || length(selected) != stop || any(!is.finite(selected)) ||
        any(selected < 1 | selected > length(native_variables)) || any(selected != floor(selected))) sb_abort(
      "The original component selection sequence is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
    variables <- unname(native_variables)
    if (!is.null(state$engine_features)) {
      index <- match(variables, state$engine_features)
      matched <- !is.na(index)
      variables[matched] <- state$features[index[matched]]
    }
    if (!all(variables %in% c(state$features, "(Intercept)"))) sb_abort(
      "A component cannot be mapped to its original input feature.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(component_id = seq_along(variables), feature = variables,
      native_component = names(native_variables),
      updates = tabulate(selected, nbins = length(variables)))
    selections <- data.frame(iteration = seq_along(selected), component_id = as.integer(selected),
      feature = variables[selected])
    measure <- "Native component update count"
  }
  structure(list(type = paste0("boosting_", type), table = tab, selections = selections,
    learner = learner, pipeline_id = fit$pipeline$id, stopping_iteration = as.integer(stop),
    measure = measure, source = "original saved mboost engine; no additional fitting",
    validation_scope = "in-bag development training, not held-out evaluation"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_boosting_diagnostic <- function(x) {
  iteration <- loss <- updates <- feature <- NULL
  if (x$type == "boosting_training_loss") {
    return(ggplot2::ggplot(x$table, ggplot2::aes(iteration, loss)) +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::geom_point(data = x$table[c(1, nrow(x$table)), ], colour = "#2878A0", size = 1.5) +
      ggplot2::labs(title = "Boosting training loss", subtitle = x$learner,
        x = "Boosting iteration (zero = initial model)", y = "Native in-bag objective",
        caption = paste0(x$measure, ".\n",
          "Original scale and sign; not validation performance or a new stopping rule.")) + sb_viz_theme())
  }
  tab <- x$table
  tab$feature <- factor(tab$feature, levels = tab$feature[order(tab$updates)])
  p <- ggplot2::ggplot(tab, ggplot2::aes(updates, feature)) +
    ggplot2::geom_col(fill = "#2878A0", width = .65) +
    ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
    ggplot2::labs(title = "Boosting component updates", subtitle = x$learner,
      x = "Number of iterations selecting this component", y = NULL,
      caption = "All original components retained, including zero counts and any intercept.\nUpdate frequency is not feature importance, effect magnitude or resampling stability.") + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 10, height = max(5, 2.5 + .28 * nrow(tab)))
  p
}
