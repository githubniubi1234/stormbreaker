#' Inspect original randomGLM ensemble members
#' @param fit A successful randomGLM final fit.
#' @param type inclusion_frequency, bag_coefficients or bag_size.
#' @return An algorithm diagnostic from compact saved ensemble members. Every
#'   fitted feature and bag is retained, including intercept-only bags. Inclusion
#'   means membership in a bag's selected feature list, not nonzero coefficient.
#'   Coefficients are individual bag log-odds slopes in original input units;
#'   omitted features have zero contribution and a separate included flag.
#'   Bag frequency is not cross-validation stability, population prevalence,
#'   uncertainty or importance. The ensemble averages probabilities, not slopes:
#'   these coefficients must not be averaged into a surrogate nomogram.
#'   No training data, refitting, bag resampling or new model is required.
#' @export
sb_randomglm_diagnostics <- function(fit, type = "inclusion_frequency") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (fit$pipeline$learner != "random_glm_binary") sb_abort(
    "This diagnostic requires randomGLM.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("inclusion_frequency", "bag_coefficients", "bag_size")) sb_abort(
    "Choose inclusion_frequency, bag_coefficients or bag_size.", "stormbreaker_unsupported_diagnostic")
  state <- fit$learner_state
  n <- length(state$intercepts)
  features <- state$features
  if (!n || length(state$bag_features) != n || length(state$bag_coefficients) != n) sb_abort(
    "Saved randomGLM ensemble is incomplete.", "stormbreaker_bad_learner_state")
  coefficient <- matrix(0, nrow = n, ncol = length(features))
  included <- matrix(FALSE, nrow = n, ncol = length(features))
  for (i in seq_len(n)) {
    index <- state$bag_features[[i]]
    beta <- state$bag_coefficients[[i]]
    if (length(index) != length(beta) || anyNA(index) || anyDuplicated(index) ||
        any(index < 1L | index > length(features)) || any(!is.finite(beta))) sb_abort(
      "Saved randomGLM bag terms are invalid.", "stormbreaker_bad_learner_state")
    coefficient[i, index] <- beta
    included[i, index] <- TRUE
  }
  if (type == "inclusion_frequency") {
    tab <- data.frame(feature = features, bags = colSums(included),
      fraction = colMeans(included), total_bags = n)
    interpretation <- "Feature inclusion in original ensemble bags; not CV selection stability, uncertainty, causal effect or permutation importance."
  } else if (type == "bag_coefficients") {
    tab <- data.frame(bag = rep(seq_len(n), length(features)), feature = rep(features, each = n),
      coefficient = as.vector(coefficient), included = as.vector(included))
    interpretation <- paste("Original per-bag log-odds coefficients in input units; omitted features have zero contribution and remain flagged.",
      "Across-bag spread is not a confidence interval. The ensemble averages probabilities, not these coefficients.")
  } else {
    tab <- data.frame(bag = seq_len(n), selected_features = rowSums(included))
    interpretation <- "Number of selected features in each original ensemble bag, including intercept-only bags; not learning progress or patient counts."
  }
  structure(list(type = paste0("randomglm_", type), table = tab, features = features,
    intercepts = state$intercepts, total_bags = n, learner = fit$pipeline$learner,
    pipeline_id = fit$pipeline$id, source = "original saved randomGLM ensemble members",
    interpretation = interpretation), class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_randomglm_diagnostic <- function(x) {
  feature <- fraction <- coefficient <- included <- selected_features <- NULL
  tab <- x$table
  if (x$type != "randomglm_bag_size") tab$feature <- factor(tab$feature, levels = rev(x$features))
  if (x$type == "randomglm_inclusion_frequency") {
    p <- ggplot2::ggplot(tab, ggplot2::aes(fraction, feature)) +
      ggplot2::geom_col(fill = "#2878A0", width = .65) +
      ggplot2::scale_x_continuous(limits = c(0, 1), labels = function(z) paste0(round(100 * z), "%")) +
      ggplot2::labs(title = "randomGLM feature inclusion", x = "Ensemble bags including feature", y = NULL)
  } else if (x$type == "randomglm_bag_coefficients") {
    p <- ggplot2::ggplot(tab, ggplot2::aes(coefficient, feature, colour = included)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey75", linewidth = .4) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0, height = .18, seed = 1),
        size = 1.3, alpha = .65) +
      ggplot2::scale_colour_manual(values = c("FALSE" = "#9AA6AF", "TRUE" = "#2878A0"),
        breaks = c("TRUE", "FALSE"), labels = c("Included", "Omitted")) +
      ggplot2::labs(title = "randomGLM bag coefficients", x = "Bag coefficient (log odds per input unit)",
        y = NULL, colour = NULL)
  } else {
    p <- ggplot2::ggplot(tab, ggplot2::aes(selected_features)) +
      ggplot2::geom_bar(fill = "#2878A0", width = .65) +
      ggplot2::scale_x_continuous(breaks = sort(unique(tab$selected_features))) +
      ggplot2::labs(title = "randomGLM bag sizes", x = "Selected predictors", y = "Number of bags")
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11,
    height = if (x$type == "randomglm_bag_size") 5 else max(5, 2 + .3 * length(x$features)))
  p
}
