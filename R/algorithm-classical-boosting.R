sb_classical_boosting_diagnostic <- function(fit, type) {
  id <- fit$pipeline$learner
  state <- fit$learner_state
  ada <- id == "adaboost_binary"
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("split_usage", if (ada) c("tree_weights", "tree_size") else "stump_sequence")) {
    sb_abort("Use split_usage for either learner; tree_weights/tree_size for AdaBoost or stump_sequence for LogitBoost. Original evaluation histories are not retained.",
      "stormbreaker_unavailable_diagnostic")
  }
  engine <- sb_validate_classical_binary_boosting_engine(unserialize(state$model_raw), id, state$iterations)
  n <- if (ada) length(engine$trees) else nrow(engine$Stump)
  native <- NULL
  if (type == "split_usage") {
    tab <- state$feature_importance
    if (!is.data.frame(tab) || !identical(tab$feature, state$features) ||
        !is.numeric(tab$importance) || any(!is.finite(tab$importance)) ||
        any(tab$importance < 0) || any(tab$importance != floor(tab$importance))) sb_abort(
      "Original boosting split counts are invalid.", "stormbreaker_invalid_diagnostic")
    native <- tab
    tab <- data.frame(feature = tab$feature, count = tab$importance)
    measure <- "Original split usage"
    interpretation <- "Unweighted primary split count; not gain, permutation importance, effect direction or selection stability."
  } else if (type == "tree_weights") {
    native <- engine$weights
    tab <- data.frame(iteration = seq_len(n), weight = as.numeric(native))
    measure <- "Original AdaBoost tree weights"
    interpretation <- "Saved ensemble voting weights, not tree accuracy, probabilities or a training-loss curve."
  } else if (type == "tree_size") {
    tab <- data.frame(tree = seq_len(n), terminal_nodes = vapply(engine$trees,
      function(tree) sum(tree$frame$var == "<leaf>"), integer(1)))
    measure <- "Original AdaBoost tree-size distribution"
    interpretation <- "Every original fitted tree counted; one tree is not the complete ensemble."
  } else {
    native <- engine$Stump
    index <- as.integer(native[, 1])
    if (any(index < 1L) || any(index > length(state$features)) ||
        any(index != native[, 1]) || any(!native[, 3] %in% c(-1, 1))) sb_abort(
      "Native LogitBoost stump indices or directions are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- data.frame(iteration = seq_len(n), feature = state$features[index],
      native_feature_index = index, threshold = as.numeric(native[, 2]),
      contribution_at_or_below = as.numeric(native[, 3]),
      contribution_above = -as.numeric(native[, 3]))
    measure <- "Original LogitBoost stump sequence"
    interpretation <- "Colour: contribution to class-1 log-odds when value <= the saved threshold; above it, the sign reverses.\nExact feature-specific thresholds remain in the table; these are model splits, not clinical cutoffs."
  }
  structure(list(type = paste0("classical_boosting_", type), table = tab, native_table = native,
    learner = id, pipeline_id = fit$pipeline$id, fitted_iterations = n,
    requested_iterations = state$iterations, features = state$features,
    measure = measure, interpretation = interpretation,
    source = "original saved boosting model; no refit or patient prediction"),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_classical_boosting_diagnostic <- function(x) {
  iteration <- weight <- feature <- count <- direction <- NULL
  tab <- x$table
  if (x$type == "classical_boosting_tree_size") {
    z <- x
    z$type <- "forest_tree_size"
    p <- sb_plot_forest_diagnostic(z) +
      ggplot2::labs(x = "Terminal nodes per original AdaBoost tree")
  } else if (x$type == "classical_boosting_tree_weights") {
    p <- ggplot2::ggplot(tab, ggplot2::aes(iteration, weight)) +
      ggplot2::geom_segment(ggplot2::aes(xend = iteration, yend = 0), colour = "#2878A0") +
      ggplot2::geom_point(colour = "#2878A0") +
      ggplot2::labs(x = "Original tree iteration", y = "Native voting weight") + sb_viz_theme()
  } else if (x$type == "classical_boosting_split_usage") {
    tab$feature <- factor(tab$feature, levels = tab$feature[order(tab$count)])
    p <- ggplot2::ggplot(tab, ggplot2::aes(count, feature)) +
      ggplot2::geom_col(fill = "#2878A0", width = .65) +
      ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
      ggplot2::labs(x = "Primary split count", y = NULL) + sb_viz_theme()
  } else {
    tab$feature <- factor(tab$feature, levels = rev(x$features))
    tab$direction <- factor(tab$contribution_at_or_below, levels = c(-1, 1), labels = c("-1", "+1"))
    p <- ggplot2::ggplot(tab, ggplot2::aes(iteration, feature, colour = direction)) +
      ggplot2::geom_point(size = 2) + ggplot2::scale_y_discrete(drop = FALSE) +
      ggplot2::scale_colour_manual(values = c("-1" = "#2878A0", "+1" = "#D77935"), drop = FALSE) +
      ggplot2::labs(x = "Original stump iteration", y = NULL, colour = "At/below threshold") + sb_viz_theme()
  }
  p + ggplot2::labs(title = x$measure,
    subtitle = paste0(x$learner, "; fitted iterations: ", x$fitted_iterations,
      if (x$learner == "adaboost_binary") "; configured maximum: " else "; configured iterations: ",
      x$requested_iterations), caption = x$interpretation)
}
