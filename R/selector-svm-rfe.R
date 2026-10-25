sb_svm_rfe_path <- function(x, outcome, assessment = NULL) {
  active <- names(x)
  subsets <- predictions <- vector("list", ncol(x))
  while (length(active)) {
    k <- length(active)
    train <- as.matrix(x[active])
    center <- colMeans(train)
    spread <- apply(train, 2L, stats::sd)
    spread[!is.finite(spread) | spread == 0] <- 1
    train <- sweep(sweep(train, 2L, center), 2L, spread, "/")
    model <- e1071::svm(train, factor(outcome, levels = c(0, 1)),
      kernel = "linear", cost = 1, scale = FALSE, probability = FALSE)
    subsets[[k]] <- active
    if (!is.null(assessment)) {
      test <- sweep(sweep(as.matrix(assessment[active]), 2L, center), 2L, spread, "/")
      predictions[[k]] <- as.integer(as.character(stats::predict(model, test)))
    }
    weights <- as.numeric(crossprod(model$coefs, model$SV))^2
    active <- active[-which.min(weights)]
  }
  list(subsets = subsets, predictions = predictions)
}

sb_fit_svm_rfe <- function(x, outcome, group_ids = NULL, folds = 5L, seed = 1L,
                           prepare_fold = NULL) {
  x <- sb_validate_predictor_frame(x, "x")
  if (length(outcome) != nrow(x) || anyNA(outcome) ||
      !all(outcome %in% c(0, 1)) || length(unique(outcome)) != 2L) {
    sb_abort("SVM-RFE requires both binary outcome classes.", "stormbreaker_bad_binary_outcome")
  }
  groups <- if (is.null(group_ids)) as.character(seq_len(nrow(x))) else
    sb_validate_group_ids(group_ids, nrow(x))
  fold_id <- sb_make_group_stratified_fold_id(outcome, groups, folds, seed,
    "stormbreaker_infeasible_inner_folds")
  correct <- vector("list", folds)
  fold_subsets <- vector("list", folds)
  for (fold in seq_len(folds)) {
    heldout <- fold_id == fold
    data <- if (is.null(prepare_fold)) list(train = x[!heldout, , drop = FALSE],
      test = x[heldout, , drop = FALSE]) else prepare_fold(!heldout)
    path <- sb_svm_rfe_path(data$train, outcome[!heldout], data$test)
    correct[[fold]] <- vapply(path$predictions,
      function(pred) sum(pred == outcome[heldout]), integer(1))
    fold_subsets[[fold]] <- path$subsets
  }
  candidates <- seq_len(min(ncol(x), lengths(correct)))
  correct <- do.call(rbind, lapply(correct, function(z) z[candidates]))
  accuracy <- colSums(correct) / nrow(x)
  best <- which.max(accuracy)
  final <- sb_svm_rfe_path(x, outcome)
  structure(list(id = "svm_rfe_binary", source_features = names(x),
    selected_features = final$subsets[[best]],
    evidence = list(curve = data.frame(n_features = candidates,
      accuracy = accuracy, error = 1 - accuracy, selected = candidates == best),
      fold_id = fold_id, group_ids = group_ids, fold_subsets = fold_subsets,
      subsets = final$subsets, correct = correct)), class = "stormbreaker_selector_fit")
}
