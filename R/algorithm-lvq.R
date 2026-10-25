#' Inspect the original learned LVQ prototypes
#' @param fit A successful LVQ1, OLVQ1, LVQ2.1 or LVQ3 final fit.
#' @param type prototypes or prototype_distances.
#' @return An algorithm diagnostic using compact saved prediction state.
#'   Prototype profiles retain the model's training-centered/scaled coordinates;
#'   the table additionally contains their exact original-unit coordinates.
#'   Distances are squared Euclidean distances in that same prediction metric,
#'   without clustering, projection, refitting or restandardizing prototypes.
#'   Native prototype order and encoded classes are retained. Class counts are
#'   prototype counts, not sample counts or estimated prevalence. Learned
#'   prototypes are not individual patients or class means. No training-history
#'   or calibration claim is made: this learner maps the difference between
#'   nearest-class squared distances through a fixed logistic function.
#' @export
sb_lvq_diagnostics <- function(fit, type = "prototypes") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!fit$pipeline$learner %in% sb_lvq_binary_learner_ids()) sb_abort(
    "This diagnostic requires an LVQ learner.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("prototypes", "prototype_distances")) sb_abort(
    "Choose prototypes or prototype_distances.", "stormbreaker_unsupported_diagnostic")
  state <- fit$learner_state
  e <- unserialize(state$engine_raw)
  z <- e$prototypes
  features <- state$features
  if (!is.matrix(z) || nrow(z) < 2L || ncol(z) != length(features) ||
      any(!is.finite(z)) || length(e$prototype_class) != nrow(z) ||
      !identical(sort(unique(e$prototype_class)), c("0", "1")) ||
      length(e$center) != ncol(z) || length(e$scale) != ncol(z) ||
      any(!is.finite(c(e$center, e$scale))) || any(e$scale <= 0)) sb_abort(
    "Saved LVQ codebook or training scaling is invalid.", "stormbreaker_bad_learner_state")
  labels <- paste0("P", seq_len(nrow(z)), " (", e$prototype_class, ")")
  original <- sweep(sweep(z, 2L, e$scale, "*"), 2L, e$center, "+")
  if (type == "prototypes") {
    tab <- data.frame(prototype = rep(labels, ncol(z)),
      prototype_index = rep(seq_len(nrow(z)), ncol(z)),
      class = rep(e$prototype_class, ncol(z)),
      feature = rep(features, each = nrow(z)), value = as.vector(z),
      original_value = as.vector(original))
    interpretation <- paste("Original learned prototypes in the fitted training-standardized coordinates;",
      "original-unit coordinates remain in the table. Parentheses identify encoded classes.",
      "Not individual patients, class means, feature importance or SHAP.")
  } else {
    distances <- as.matrix(stats::dist(z))^2
    tab <- data.frame(prototype = rep(labels, nrow(z)), prototype2 = rep(labels, each = nrow(z)),
      class = rep(e$prototype_class, nrow(z)), class2 = rep(e$prototype_class, each = nrow(z)),
      value = as.vector(distances))
    interpretation <- paste("Squared Euclidean distances between original learned prototypes,",
      "using the model's training scaling. Parentheses identify encoded classes.",
      "Not patient distances, calibrated risk, a class separation test or validation performance.")
  }
  structure(list(type = paste0("lvq_", type), table = tab, features = features, prototypes = labels,
    prototype_counts = data.frame(class = c("0", "1"), count = as.integer(table(factor(e$prototype_class, levels = c("0", "1"))))),
    center = stats::setNames(e$center, features), scale = stats::setNames(e$scale, features),
    iterations = e$iterations, learner = fit$pipeline$learner, pipeline_id = fit$pipeline$id,
    source = "original saved LVQ codebook and training scaling", interpretation = interpretation),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_lvq_diagnostic <- function(x) {
  feature <- prototype <- prototype2 <- value <- NULL
  tab <- x$table
  tab$prototype <- factor(tab$prototype, levels = rev(x$prototypes))
  if (x$type == "lvq_prototypes") {
    tab$feature <- factor(tab$feature, levels = x$features)
    p <- ggplot2::ggplot(tab, ggplot2::aes(feature, prototype, fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "#2878A0", mid = "white", high = "#D77935",
        guide = ggplot2::guide_colourbar(display = "rectangles")) +
      ggplot2::labs(title = "Learned LVQ prototypes", x = NULL, y = "Prototype (class)",
        fill = "Standardized\ncoordinate")
    size <- c(width = max(11, 3 + .45 * length(x$features)), height = max(5, 2 + .35 * length(x$prototypes)))
  } else {
    tab$prototype2 <- factor(tab$prototype2, levels = x$prototypes)
    p <- ggplot2::ggplot(tab, ggplot2::aes(prototype2, prototype, fill = value)) +
      ggplot2::geom_tile(colour = "#D9DFE3", linewidth = 0.3) + ggplot2::coord_fixed() +
      ggplot2::scale_fill_gradient(low = "white", high = "#2878A0",
        guide = ggplot2::guide_colourbar(display = "rectangles")) +
      ggplot2::labs(title = "LVQ prototype distances", x = "Prototype (class)", y = "Prototype (class)",
        fill = "Squared\ndistance")
    size <- c(width = max(8, 3 + .4 * length(x$prototypes)), height = max(6, 2 + .4 * length(x$prototypes)))
  }
  p <- p + sb_viz_theme() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  attr(p, "stormbreaker_suggested_size") <- size
  p
}
