#' Inspect fitted Naive Bayes class-conditional distributions
#' @param fit A final Gaussian or kernel Naive Bayes model.
#' @param type conditional_density or log_likelihood_ratio (Gaussian only).
#' @param features Predictor names to display; NULL retains all fitted predictors.
#' @param grid_points Gaussian grid size per feature. Kernel models retain their
#'   original density grids and ignore this argument.
#' @return An algorithm diagnostic using compact saved state, without capture or
#'   refitting. Gaussian curves use the original smoothed variances and a
#'   model-implied grid spanning both class means plus/minus four fitted SDs.
#'   Kernel curves retain every native density-grid point, with no new density
#'   estimate or bandwidth choice. Density is not class probability or importance;
#'   axes differ across features because their units differ. Gaussian likelihood
#'   ratios are log f(x|class 1) minus log f(x|class 0), excluding prior log odds.
#'   Prior log odds plus all feature contributions yield the original posterior
#'   log odds under the fitted conditional-independence model.
#'   The result's `source` and `interpretation` fields provide separate figure
#'   legend text; plotting does not put these explanations inside the panels.
#' @examples
#' d <- sb_example_binary(n = 60L, seed = 108L)
#' task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
#' fit <- sb_fit_final(task, sb_pipeline(learner = "gaussian_naive_bayes"))
#' diagnostic <- sb_naive_bayes_diagnostics(fit)
#' diagnostic$interpretation
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   sb_plot_algorithm_diagnostic(diagnostic)
#' }
#' @export
sb_naive_bayes_diagnostics <- function(fit, type = "conditional_density", features = NULL, grid_points = 101L) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  id <- fit$pipeline$learner
  if (!id %in% c("gaussian_naive_bayes", "kernel_naive_bayes_binary")) sb_abort(
    "This diagnostic requires Gaussian or kernel Naive Bayes.", "stormbreaker_unsupported_diagnostic")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("conditional_density", "log_likelihood_ratio")) sb_abort(
    "Choose conditional_density or log_likelihood_ratio.", "stormbreaker_unsupported_diagnostic")
  all_features <- fit$learner_state$features
  if (is.null(features)) features <- all_features
  if (!is.character(features) || !length(features) || anyNA(features) || anyDuplicated(features) ||
      !all(features %in% all_features)) sb_abort("Choose unique fitted predictor names.", "stormbreaker_invalid_diagnostic")
  prior_log_odds <- NULL
  if (id == "gaussian_naive_bayes") {
    points <- sb_check_positive_integer(grid_points, "grid_points", minimum = 2L)
    m <- fit$learner_state$payload
    prior_log_odds <- log(m$priors[["1"]]) - log(m$priors[["0"]])
    tables <- lapply(features, function(feature) {
      mu <- m$means[c("0", "1"), feature]
      sd <- sqrt(m$variances[c("0", "1"), feature])
      if (any(!is.finite(mu)) || any(!is.finite(sd)) || any(sd <= 0)) sb_abort(
        "Original Gaussian parameters are invalid.", "stormbreaker_invalid_diagnostic")
      x <- seq(min(mu - 4 * sd), max(mu + 4 * sd), length.out = points)
      if (type == "log_likelihood_ratio") return(data.frame(feature = feature, x = x,
        value = stats::dnorm(x, mu[2L], sd[2L], log = TRUE) - stats::dnorm(x, mu[1L], sd[1L], log = TRUE)))
      rbind(data.frame(feature = feature, class = "0", x = x, value = stats::dnorm(x, mu[1L], sd[1L])),
        data.frame(feature = feature, class = "1", x = x, value = stats::dnorm(x, mu[2L], sd[2L])))
    })
    native <- list(means = m$means, variances = m$variances, priors = m$priors, epsilon = m$epsilon)
    source <- "original Gaussian means and smoothed variances; model-implied +/-4 SD grid"
  } else {
    if (type == "log_likelihood_ratio") sb_abort(
      "Kernel likelihood contributions require native prediction-floor semantics; this view exposes only the original densities.",
      "stormbreaker_unavailable_diagnostic")
    e <- unserialize(fit$learner_state$engine_raw)
    if (!inherits(e, "NaiveBayes") || !isTRUE(e$usekernel) || !identical(e$levels, c("0", "1"))) sb_abort(
      "Original kernel Naive Bayes model is invalid.", "stormbreaker_invalid_diagnostic")
    tables <- lapply(features, function(feature) do.call(rbind, lapply(c("0", "1"), function(class) {
      density <- e$tables[[feature]][[class]]
      if (!inherits(density, "density") || length(density$x) != length(density$y) ||
          length(density$x) < 2L || any(diff(density$x) <= 0)) sb_abort(
        "Original kernel density grid is invalid.", "stormbreaker_invalid_diagnostic")
      data.frame(feature = feature, class = class, x = density$x, value = density$y)
    })))
    native <- e$tables[features]
    source <- "original saved class-conditional kernel-density grids; no bandwidth change or refitting"
  }
  tab <- do.call(rbind, tables)
  rownames(tab) <- NULL
  if (any(!is.finite(tab$x)) || any(!is.finite(tab$value)) ||
      (type == "conditional_density" && any(tab$value < 0))) sb_abort(
    "Class-conditional curve contains invalid values.", "stormbreaker_invalid_diagnostic")
  structure(list(type = paste0("naive_bayes_", type), table = tab, native_table = native,
    features = features, prior_log_odds = prior_log_odds, learner = id, pipeline_id = fit$pipeline$id,
    source = source, interpretation = if (type == "conditional_density")
      "Fitted class-conditional densities, not posterior class probabilities or variable importance.\nFeature axes use their own units; do not compare density peak heights across features." else
      "Original log likelihood contribution toward class 1; excludes prior log odds and other features.\nNot SHAP, a causal effect, or a validated clinical threshold; original variance smoothing retained."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_naive_bayes_diagnostic <- function(x) {
  feature <- value <- class <- NULL
  tab <- x$table
  tab$feature <- factor(tab$feature, levels = x$features)
  if (x$type == "naive_bayes_conditional_density") {
    p <- ggplot2::ggplot(tab, ggplot2::aes(x, value, colour = class)) +
      ggplot2::geom_line() + ggplot2::scale_colour_manual(values = c("0" = "#2878A0", "1" = "#D77935")) +
      ggplot2::scale_y_continuous(labels = function(z) sprintf("%.1e", z)) +
      ggplot2::labs(title = "Fitted class-conditional distributions", y = "Density", colour = "Encoded class")
  } else {
    p <- ggplot2::ggplot(tab, ggplot2::aes(x, value)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey65", linetype = 2) +
      ggplot2::geom_line(colour = "#2878A0") +
      ggplot2::labs(title = "Gaussian Naive Bayes log likelihood contributions", y = "Log likelihood ratio: class 1 / class 0")
  }
  p + ggplot2::facet_wrap(~feature,
    scales = if (x$type == "naive_bayes_conditional_density") "free" else "free_x",
    ncol = min(3L, length(x$features))) +
    ggplot2::labs(x = "Feature value (original units)") + sb_viz_theme()
}
