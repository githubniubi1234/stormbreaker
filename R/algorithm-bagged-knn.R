#' Inspect the original bagged survival KNN ensemble
#' @param fit Successful bagged survival KNN final fit.
#' @param type bootstrap or bag_curves.
#' @param newdata Explicit query rows required for bag_curves.
#' @return Original bootstrap multiplicities and per-bag coverage, or native
#'   per-bag survival curves and their ensemble mean. Bag-to-bag variation is
#'   not a confidence interval. No refitting or resampling occurs. Bootstrap
#'   row indices are sensitive; repeated rows remain repeated in prediction.
#' @export
sb_bagged_knn_diagnostics <- function(fit, type = "bootstrap", newdata = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (s$learner != "bagged_knn_survival" || !is.character(type) ||
      length(type) != 1L || is.na(type) || !type %in% c("bootstrap", "bag_curves")) sb_abort(
    "Choose bootstrap or bag_curves for bagged survival KNN.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  bags <- methods::slot(e, "base_learners")
  train <- methods::slot(e, "train_data")
  if (type == "bootstrap") {
    tab <- do.call(rbind, lapply(seq_along(bags), function(i) {
      index <- methods::slot(bags[[i]], "bootstrap_sample")
      data.frame(bag = i, training_row = seq_len(nrow(train)),
        multiplicity = tabulate(index, nbins = nrow(train)))
    }))
    summary <- do.call(rbind, lapply(split(tab, tab$bag), function(z)
      data.frame(bag = z$bag[1L], sampled_rows = sum(z$multiplicity),
        unique_rows = sum(z$multiplicity > 0), out_of_bag_rows = sum(z$multiplicity == 0))))
  } else {
    if (is.null(newdata)) sb_abort("Supply explicit queries for bag-specific survival curves.",
      "stormbreaker_unavailable_diagnostic")
    x <- sb_select_final_feature_frame(newdata, s$features, "Query data", "query fitted features")
    if (!nrow(x)) sb_abort("Query data must contain at least one row.", "stormbreaker_invalid_diagnostic")
    x <- data.matrix(sb_bagged_knn_internal_data(x)$data)
    time <- methods::slot(e, "timepoints")
    native_predict <- get("predict", asNamespace("bnnSurvival"))
    values <- lapply(bags, function(b) native_predict(b, train, x, time,
      methods::slot(e, "metric"), methods::slot(e, "weighting_function"), methods::slot(e, "k")))
    if (any(vapply(values, function(v) !is.matrix(v) ||
        !identical(dim(v), c(nrow(x), length(time))) || any(!is.finite(v)) ||
        any(v < 0 | v > 1), logical(1)))) sb_abort(
      "Native bag-specific survival curves are invalid.", "stormbreaker_invalid_diagnostic")
    tab <- do.call(rbind, lapply(seq_along(values), function(i)
      data.frame(bag = i, query = rep(seq_len(nrow(x)), each = length(time)),
        time = rep(time, nrow(x)), survival = as.numeric(t(values[[i]])))))
    summary <- stats::aggregate(survival ~ query + time, tab, mean)
  }
  structure(list(type = paste0("bagged_knn_", type), table = tab, summary = summary,
    learner = s$learner, pipeline_id = fit$pipeline$id, sensitive_training_data = type == "bootstrap",
    source = "Original bnnSurvival bootstrap indices and native base-learner predictions",
    interpretation = paste("All original bags retained. No new bootstrap or fitted model.",
      "Grey curves show individual bags; blue is the original ensemble mean, not a confidence band.",
      "Bootstrap coverage is not OOB predictive performance. Explicit queries are not validation.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_bagged_knn_diagnostic <- function(x) {
  bag <- unique_rows <- time <- survival <- query <- NULL
  if (x$type == "bagged_knn_bootstrap") {
    p <- ggplot2::ggplot(x$summary, ggplot2::aes(bag, unique_rows)) +
      ggplot2::geom_point(colour = "#2878A0") +
      ggplot2::labs(title = "Bootstrap coverage", x = "Bag", y = "Unique training rows")
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(time, survival, group = bag)) +
      ggplot2::geom_step(colour = "#A5ADB4", alpha = 0.35, linewidth = 0.3) +
      ggplot2::geom_step(data = x$summary, ggplot2::aes(group = query), colour = "#2878A0", linewidth = 0.9) +
      ggplot2::facet_wrap(~query, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "Bag-specific survival predictions", x = "Time", y = "Survival probability")
  }
  p + sb_viz_theme()
}
