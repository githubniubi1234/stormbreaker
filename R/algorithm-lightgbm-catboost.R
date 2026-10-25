sb_catboost_diagnostic_ids <- function() {
  c("catboost_binary", "catboost_cox", "catboost_aft_survival")
}

sb_light_cat_diagnostic <- function(fit, type, max_pairs) {
  id <- fit$pipeline$learner
  state <- fit$learner_state
  light <- id %in% sb_lightgbm_binary_learner_ids()
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("importance", "tree_size", "interaction", "training_loss", "cv_curve")) sb_abort(
    "Choose a supported importance, tree_size or interaction diagnostic.",
    "stormbreaker_unsupported_diagnostic")
  if (type %in% c("training_loss", "cv_curve")) sb_abort(
    "No original evaluation history is retained by these pipelines; none is reconstructed.",
    "stormbreaker_unavailable_diagnostic")
  native <- NULL
  requested <- NULL
  display_limit <- NULL
  if (type == "importance") {
    tab <- state$feature_importance
    required <- if (light) c("gain", "cover", "frequency") else "importance"
    if (!is.data.frame(tab) || !all(c("feature", required) %in% names(tab)) ||
        anyDuplicated(tab$feature) || !setequal(tab$feature, state$features) ||
        !all(vapply(tab[required], function(z) is.numeric(z) && all(is.finite(z)), logical(1)))) sb_abort(
      "Saved native feature importance is incomplete or invalid.", "stormbreaker_invalid_diagnostic")
    native <- tab
    tab$value <- if (light) tab$gain else tab$importance
    measure <- if (light) "Native LightGBM gain (unscaled)" else "Native CatBoost FeatureImportance"
    interpretation <- if (light) "Saved gain, cover and split counts; not normalized XGBoost gain shares." else
      "Original CatBoost importance; not gain, SHAP, effect direction or a causal effect."
  } else if (type == "tree_size") {
    if (!light) sb_abort("This tree-size adapter is for LightGBM, not CatBoost.",
      "stormbreaker_unavailable_diagnostic")
    if (!requireNamespace("lightgbm", quietly = TRUE)) sb_abort("Package lightgbm is required.",
      "stormbreaker_missing_dependency")
    if (!requireNamespace("jsonlite", quietly = TRUE)) sb_abort("Package jsonlite is required to read the native LightGBM dump.",
      "stormbreaker_missing_dependency")
    booster <- lightgbm::lgb.load(model_str = state$model_string)
    native <- jsonlite::fromJSON(lightgbm::lgb.dump(booster), simplifyVector = FALSE)$tree_info
    leaves <- vapply(native, function(z) as.numeric(z$num_leaves), numeric(1))
    index <- vapply(native, function(z) as.numeric(z$tree_index), numeric(1))
    if (!length(leaves) || any(!is.finite(leaves)) || any(leaves < 1) ||
        any(leaves != floor(leaves)) || !identical(index, as.numeric(seq_along(leaves) - 1L))) sb_abort(
      "Native LightGBM tree counts or indices are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(tree = index + 1L, native_tree_index = index, terminal_nodes = as.integer(leaves))
    # Do not retain the full dump: the diagnostic needs counts, not every split.
    native <- NULL
    requested <- state$nrounds
    measure <- "Terminal nodes per original LightGBM tree"
    interpretation <- paste0("Actual saved trees: ", nrow(tab), "; configured rounds: ", requested,
      ". Fewer trees can be retained when no further valid split is available.")
  } else {
    if (light) sb_abort("This native interaction adapter is for CatBoost, not LightGBM.",
      "stormbreaker_unavailable_diagnostic")
    display_limit <- sb_check_positive_integer(max_pairs, "max_pairs", minimum = 1L)
    if (!requireNamespace("catboost", quietly = TRUE)) sb_abort("Package catboost is required.",
      "stormbreaker_missing_dependency")
    booster <- unserialize(state$model_raw)
    # With one input there are no distinct feature pairs. The native R wrapper
    # tries to name three columns on an empty result and fails before returning.
    native <- if (length(state$features) == 1L) {
      matrix(numeric(), nrow = 0L, ncol = 3L)
    } else {
      tryCatch(
        catboost::catboost.get_feature_importance(booster, type = "Interaction", thread_count = 1L),
        error = function(e) sb_abort(
          paste0("Native CatBoost interaction extraction failed: ", conditionMessage(e),
            ". No interaction scores were substituted."),
          "stormbreaker_unavailable_diagnostic", parent = e)
      )
    }
    if (is.null(native) || !length(native)) native <- matrix(numeric(), nrow = 0L, ncol = 3L)
    if (!is.matrix(native) || ncol(native) != 3L || !is.numeric(native) || any(!is.finite(native)) ||
        any(native[, 1:2, drop = FALSE] < 0) || any(native[, 1:2, drop = FALSE] >= length(state$features)) ||
        any(native[, 1:2, drop = FALSE] != floor(native[, 1:2, drop = FALSE])) ||
        any(native[, 1] == native[, 2])) sb_abort(
      "Native CatBoost interaction indices or scores are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(feature_1 = state$features[native[, 1] + 1L],
      feature_2 = state$features[native[, 2] + 1L], score = unname(native[, 3]),
      native_index_1 = unname(native[, 1]), native_index_2 = unname(native[, 2]))
    measure <- "Native CatBoost interaction strength"
    interpretation <- if (length(state$features) == 1L) {
      "Only one input feature: no distinct feature pairs exist; the native pair extractor was not called."
    } else {
      "Model-internal interaction strength, not SHAP interaction, clinical synergy or causality."
    }
  }
  structure(list(type = paste0("light_cat_", type), table = tab, native_table = native,
    learner = id, pipeline_id = fit$pipeline$id, measure = measure,
    interpretation = interpretation, requested_iterations = requested, max_pairs = display_limit,
    source = "original saved boosting model; no training, permutation or patient prediction"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_light_cat_diagnostic <- function(x) {
  if (x$type == "light_cat_importance") {
    z <- x
    z$type <- "boosted_tree_importance"
    return(sb_plot_gbm_xgboost_diagnostic(z))
  }
  if (x$type == "light_cat_tree_size") {
    z <- x
    z$type <- "forest_tree_size"
    return(sb_plot_forest_diagnostic(z) + ggplot2::labs(
      title = "LightGBM tree-size distribution", subtitle = x$learner,
      caption = paste0(x$interpretation, "\nAll actual trees counted; no artificial padding to the configured iteration count.")))
  }
  tab <- x$table
  tab <- tab[utils::head(order(tab$score, decreasing = TRUE), x$max_pairs), , drop = FALSE]
  pair_labels <- stats::setNames(paste(tab$feature_1, tab$feature_2, sep = " + "), seq_len(nrow(tab)))
  tab$pair <- factor(seq_len(nrow(tab)), levels = rev(seq_len(nrow(tab))))
  score <- pair <- NULL
  p <- ggplot2::ggplot(tab, ggplot2::aes(score, pair)) +
    ggplot2::geom_col(fill = "#2878A0", width = .65) +
    ggplot2::scale_y_discrete(labels = pair_labels) +
    ggplot2::labs(title = x$measure,
      subtitle = paste0(x$learner, ": showing ", nrow(tab), " of ", nrow(x$table), " reported pairs"),
      x = "Native interaction score", y = NULL,
      caption = paste0(x$interpretation, "\nComplete reported pairs remain in the table; unreported pairs are not assigned zero scores.")) + sb_viz_theme()
  if (!nrow(tab)) {
    p <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
      label = "No distinct feature pairs are available to display.") + ggplot2::theme_void() +
      ggplot2::labs(title = x$measure, subtitle = x$learner, caption = x$interpretation)
  }
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = max(5, 2.5 + .28 * nrow(tab)))
  p
}
