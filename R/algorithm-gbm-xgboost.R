sb_gbm_diagnostic_ids <- function() c("gbm_binary", "gbm_cox", "gbm_cox_cv")

sb_gbm_xgboost_diagnostic <- function(fit, type) {
  learner <- fit$pipeline$learner
  state <- fit$learner_state
  xgb <- learner %in% sb_xgboost_learner_ids()
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("training_loss", "cv_curve", "importance")) sb_abort(
    "This boosting engine supports training_loss, cv_curve or importance where originally available.",
    "stormbreaker_unsupported_diagnostic")
  if (xgb && type != "importance") sb_abort(
    "These XGBoost pipelines did not record an evaluation history; no loss or CV curve is reconstructed.",
    "stormbreaker_unavailable_diagnostic")
  selected <- NULL
  native <- NULL
  if (xgb) {
    model <- sb_load_xgboost_state_model(state, learner)
    native <- as.data.frame(xgboost::xgb.importance(feature_names = state$features, model = model$booster))
    linear <- identical(state$params$booster, "gblinear")
    columns <- if (linear) "Weight" else c("Gain", "Cover", "Frequency")
    if (nrow(native) && (!all(c("Feature", columns) %in% names(native)) ||
        anyDuplicated(native$Feature) || !all(native$Feature %in% state$features) ||
        any(!is.finite(as.matrix(native[columns]))))) sb_abort(
      "Native XGBoost feature scores are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(feature = state$features)
    for (column in columns) {
      values <- rep(0, length(state$features))
      if (nrow(native)) values[match(native$Feature, state$features)] <- native[[column]]
      tab[[column]] <- values
    }
    tab$value <- if (linear) tab$Weight else tab$Gain
    measure <- if (linear) "Original linear coefficient" else "Native split-gain share"
    scope <- if (linear) {
      scale <- if (learner %in% sb_xgboost_aft_learner_ids()) "log survival time" else
        if (learner %in% sb_xgboost_cox_learner_ids()) "log hazard" else "log odds for class 1"
      paste0("Coefficient scale: ", scale, ". Bias/base score excluded; not a complete prediction equation.")
    } else "Native gain, cover and split-frequency shares retained; not permutation importance or SHAP."
  } else if (type == "cv_curve") {
    if (!identical(learner, "gbm_cox_cv")) sb_abort(
      "This fixed GBM pipeline did not run inner CV.", "stormbreaker_unavailable_diagnostic")
    if (identical(state$tuning$search_status, "single_candidate_no_search")) sb_abort(
      "Only one tree count was allowed; no tree-count search or CV curve was computed.",
      "stormbreaker_unavailable_diagnostic")
    loss <- state$tuning$cv_error
    selected <- state$tuning$selected_trees
    if (!is.numeric(loss) || !length(loss) || any(!is.finite(loss)) ||
        length(selected) != 1L || !selected %in% seq_along(loss)) sb_abort(
      "The original GBM tuning sequence is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(iteration = seq_along(loss), loss = loss)
    measure <- "Native Cox inner-CV error"
    scope <- "Original inner-CV search, not validation of prior feature selection or untouched-cohort performance."
  } else {
    if (!requireNamespace("gbm", quietly = TRUE)) sb_abort("Package gbm is required.",
      "stormbreaker_missing_dependency")
    if (learner == "gbm_binary") {
      if (type == "training_loss") {
        if (is.null(fit$diagnostic_artifacts)) sb_abort(
          "Binary GBM training errors require capture_diagnostics when fitting.",
          "stormbreaker_missing_diagnostic_artifacts")
        engine <- fit$diagnostic_artifacts$learner$engine
      } else engine <- unserialize(state$engine_raw)
    } else engine <- unserialize(state$model_raw)
    if (!inherits(engine, "gbm")) sb_abort("The native GBM model is missing.",
      "stormbreaker_invalid_diagnostic")
    if (type == "training_loss") {
      loss <- engine$train.error
      if (!is.numeric(loss) || length(loss) != engine$n.trees || !length(loss) || any(!is.finite(loss))) sb_abort(
        "The original GBM training errors are incomplete or invalid.", "stormbreaker_invalid_diagnostic")
      tab <- data.frame(iteration = seq_along(loss), loss = loss)
      measure <- paste("Native training error:", engine$distribution$name)
      scope <- "Training sequence of the final fit; not the separate tuning search or held-out performance."
    } else {
      value <- gbm::relative.influence(engine, n.trees = engine$n.trees, scale. = FALSE, sort. = FALSE)
      if (length(value) != length(engine$var.names) || any(!is.finite(value)) ||
          !setequal(engine$var.names, state$features)) sb_abort(
        "Native GBM relative influence is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
      native <- data.frame(feature = engine$var.names, value = unname(value))
      tab <- native[match(state$features, native$feature), , drop = FALSE]
      rownames(tab) <- NULL
      measure <- "Native relative influence (unscaled)"
      scope <- "Original split-based influence; not permutation importance, a causal effect or SHAP."
    }
  }
  structure(list(type = paste0("boosted_tree_", type), table = tab, native_table = native,
    learner = learner, pipeline_id = fit$pipeline$id, selected_iteration = selected,
    measure = measure, interpretation = scope, source = "original fitted boosting model; no refitting"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_gbm_xgboost_diagnostic <- function(x) {
  value <- feature <- iteration <- loss <- NULL
  if (x$type == "boosted_tree_importance") {
    tab <- x$table
    tab$feature <- factor(tab$feature, levels = tab$feature[order(abs(tab$value))])
    p <- ggplot2::ggplot(tab, ggplot2::aes(value, feature)) +
      ggplot2::geom_col(fill = "#2878A0", width = .65) +
      ggplot2::geom_vline(xintercept = 0, colour = "#909090") +
      ggplot2::labs(title = x$measure, subtitle = x$learner, x = x$measure, y = NULL,
        caption = paste0(x$interpretation, "\nAll input features retained; zero-use features are not silently dropped.")) + sb_viz_theme()
    attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = max(5, 2.5 + .28 * nrow(tab)))
    return(p)
  }
  p <- ggplot2::ggplot(x$table, ggplot2::aes(iteration, loss)) +
    ggplot2::geom_point(colour = "#2878A0", size = 1) +
    ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
    ggplot2::labs(title = x$measure, subtitle = x$learner, x = "Number of trees", y = "Native error",
      caption = paste0(x$interpretation, "\nNo refitting, error-bar invention or alteration of the selected tree count.")) + sb_viz_theme()
  if (nrow(x$table) > 1L) p <- p + ggplot2::geom_line(colour = "#2878A0")
  if (!is.null(x$selected_iteration)) p <- p + ggplot2::geom_vline(
    xintercept = x$selected_iteration, linetype = 2, colour = "#D77935")
  p
}
