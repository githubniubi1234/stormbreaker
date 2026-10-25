sb_validate_multiclass_probabilities <- function(probability, classes) {
  if (!is.matrix(probability) || !is.numeric(probability) || !nrow(probability) ||
      !identical(colnames(probability), classes) || any(!is.finite(probability)) ||
      any(probability < 0 | probability > 1) || any(abs(rowSums(probability) - 1) > 1e-6)) sb_abort(
    "Supply finite class probabilities summing to one, in declared class order.",
    "stormbreaker_bad_multiclass_probabilities")
  invisible(probability)
}

#' Evaluate aligned multiclass probabilities
#' @param truth Observed labels, in the same row order as probability.
#' @param probability Numeric matrix with one named column per class.
#' @return A list with overall and per-class metrics and a confusion table.
#'   Balanced accuracy is unweighted mean class recall. Macro F1 is unweighted
#'   mean of class F1; zero predicted counts yield precision NA and F1 zero if
#'   the class is present. If any class is absent, macro recall/F1/AUC are NA
#'   rather than silently excluding it. AUC is one-versus-rest, macro averaged.
#'   Brier score is mean sum of squared class errors (range 0 to 2). Log loss
#'   clips probabilities only for the logarithm at 1e-15. Ties in predicted
#'   class use the first declared class. These are row-weighted point estimates.
#' @export
sb_multiclass_metrics <- function(truth, probability) {
  classes <- colnames(probability)
  if (is.null(classes) || length(classes) < 3L || anyNA(classes) ||
      anyDuplicated(classes) || any(!nzchar(classes))) sb_abort(
    "Probability columns must identify at least three distinct classes.",
    "stormbreaker_bad_multiclass_probabilities")
  sb_validate_multiclass_probabilities(probability, classes)
  if (length(truth) != nrow(probability) || anyNA(truth) ||
      any(!as.character(truth) %in% classes)) sb_abort(
    "Truth must contain known classes, aligned with probability rows.",
    "stormbreaker_bad_multiclass_outcome")
  truth <- factor(truth, levels = classes)
  predicted <- factor(classes[max.col(probability, ties.method = "first")], levels = classes)
  counts <- table(truth = truth, predicted = predicted)
  tp <- diag(counts)
  actual <- rowSums(counts)
  guessed <- colSums(counts)
  per_class <- data.frame(class = classes, n = as.integer(actual),
    precision = ifelse(guessed > 0, tp / guessed, NA_real_),
    recall = ifelse(actual > 0, tp / actual, NA_real_),
    f1 = ifelse(actual > 0, 2 * tp / (actual + guessed), NA_real_),
    auc = vapply(seq_along(classes), function(k)
      sb_binary_auc_stratum(as.integer(truth == classes[k]), probability[, k])$estimate,
      numeric(1)), row.names = NULL)
  target <- outer(as.character(truth), classes, `==`)
  list(overall = data.frame(accuracy = mean(truth == predicted),
    balanced_accuracy = mean(per_class$recall), macro_f1 = mean(per_class$f1),
    macro_auc = mean(per_class$auc),
    log_loss = -mean(log(pmax(probability[cbind(seq_along(truth), as.integer(truth))], 1e-15))),
    brier = mean(rowSums((probability - target)^2))),
    per_class = per_class, confusion = counts, predicted = predicted)
}

#' Compare native multiclass models with shared outer folds
#' @param task A [sb_multiclass_task()].
#' @param methods Supported classifier names; defaults to Lasso, random forest
#'   and radial SVM, a fixed small starting set rather than a data-selected winner. See
#'   sb_method_catalog("multiclass"). Standalone selectors (including Boruta,
#'   varSelRF and VSURF) are not classifiers and are not accepted here.
#' @param folds,seed Outer stratified/grouped fold count and random seed.
#' @param inner_folds,inner_seed Inner tuning settings, used only within each
#'   outer training set. Glmnet standardization is fitted within resampling.
#' @return A `stormbreaker_multiclass_comparison` with a comparison table,
#'   shared fold assignment, row-aligned OOF probabilities, per-class metrics
#'   and each outer fit's tuning/feature results. A failed fold fails that
#'   method's overall evaluation; partial predictions are not scored as complete.
#'   No winner, clinical action threshold, DCA or confusion plot is generated.
#'   Candidate genes must be prespecified, not selected using held-out outcomes.
#' @export
sb_compare_multiclass <- function(task, methods = c("lasso", "random_forest", "svm_radial"), folds = 5L, seed = 1L,
                                   inner_folds = 5L, inner_seed = 1701L) {
  task <- sb_validate_multiclass_task(task)
  methods <- sb_method_names(methods)
  catalog <- sb_method_catalog("multiclass")
  catalog$supported[catalog$method %in% sb_multiclass_standalone_selectors()] <- FALSE
  if (anyDuplicated(methods) || any(!methods %in% catalog$method[catalog$supported])) sb_abort(
    "Choose distinct supported multiclass methods.", "stormbreaker_bad_method")
  assignment <- sb_multiclass_folds(task, folds, seed)
  results <- lapply(methods, function(method) {
    probability <- matrix(NA_real_, task$n, length(task$classes),
      dimnames = list(as.character(task$data[[task$id_col]]), task$classes))
    fold_results <- vector("list", folds)
    for (fold in seq_len(folds)) {
      train <- which(assignment != fold)
      test <- which(assignment == fold)
      training <- sb_multiclass_task(task$data[train, , drop = FALSE], task$id_col,
        task$outcome_col, task$predictors, task$classes, task$group_col)
      fold_results[[fold]] <- tryCatch({
        fit <- sb_fit_multiclass(training, method, inner_folds, inner_seed,
          seed = seed + fold)
        probability[test, ] <- sb_predict_multiclass(fit, task$data[test, , drop = FALSE])
        list(status = "success", fold = fold, training_rows = train, assessment_rows = test,
          inner_fold_id = fit$fold_id, parameters = fit$parameters,
          warnings = fit$warnings,
          features = sb_multiclass_evidence(fit))
      }, error = function(e) list(status = "failed", fold = fold,
        error_message = conditionMessage(e)))
    }
    failed <- vapply(fold_results, function(x) x$status == "failed", logical(1))
    list(method = method, status = if (any(failed)) "failed" else "success",
      error_message = if (any(failed)) paste(vapply(fold_results[failed],
        function(x) paste0("Fold ", x$fold, ": ", x$error_message), character(1)), collapse = "; ") else "",
      probability = probability, folds = fold_results,
      metrics = if (any(failed)) NULL else sb_multiclass_metrics(task$data[[task$outcome_col]], probability))
  })
  names(results) <- methods
  comparison <- do.call(rbind, lapply(results, function(result) {
    row <- data.frame(method = result$method, status = result$status, reason = result$error_message)
    metrics <- if (result$status == "success") result$metrics$overall else data.frame(
      accuracy = NA_real_, balanced_accuracy = NA_real_, macro_f1 = NA_real_,
      macro_auc = NA_real_, log_loss = NA_real_, brier = NA_real_)
    cbind(row, metrics)
  }))
  rownames(comparison) <- NULL
  structure(list(comparison = comparison, results = results, assignment = assignment,
    ids = as.character(task$data[[task$id_col]]), truth = task$data[[task$outcome_col]],
    classes = task$classes, predictors = task$predictors, seed = seed,
    interpretation = "Shared outer-fold evaluation; multinomial tuning within training folds"),
    class = "stormbreaker_multiclass_comparison")
}

#' Explicitly plot a multiclass confusion matrix
#' @param x A result from [sb_multiclass_metrics()].
#' @return A ggplot with exactly one count per true/predicted class pair.
#'   No model fitting occurs. This plot is never generated by default.
#' @export
sb_plot_multiclass_confusion <- function(x) {
  if (!is.list(x) || !inherits(x$confusion, "table")) sb_abort(
    "Supply sb_multiclass_metrics() output.", "stormbreaker_bad_visualization_argument")
  truth <- predicted <- Freq <- NULL
  ggplot2::ggplot(as.data.frame(x$confusion), ggplot2::aes(predicted, truth, fill = Freq)) +
    ggplot2::geom_tile(colour = "white") + ggplot2::geom_text(ggplot2::aes(label = Freq)) +
    ggplot2::scale_fill_gradient(low = "#eff6fa", high = "#73aac5", guide = "none") +
    ggplot2::coord_equal() + sb_viz_theme() +
    ggplot2::labs(title = "Confusion matrix", x = "Predicted class", y = "True class")
}
