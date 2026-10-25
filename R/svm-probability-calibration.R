# The sigmoid optimizer below adapts LIBSVM's sigmoid_train algorithm.
# Copyright (c) 2000-2023 Chih-Chung Chang and Chih-Jen Lin.
# See inst/COPYRIGHTS for the complete BSD-3-Clause notice and disclaimer.

sb_svm_probability_learner_ids <- function() {
  c("svm_linear_binary", "svm_radial_binary",
    "svm_radial_weighted_binary", "svm_polynomial_binary",
    "svm_sigmoid_binary")
}

sb_make_svm_calibration_folds <- function(outcome, seed, group_ids = NULL) {
  n <- length(outcome)
  class_order <- unlist(lapply(unique(outcome), function(value) {
    which(outcome == value)
  }), use.names = FALSE)
  # LIBSVM first groups by class and then shuffles with one uniform per row.
  # Keep that ordering, including the subproblem training order, for singleton
  # groups. R's sample.int() does not reproduce the native shuffle.
  shuffled <- withr::with_seed(seed, {
    perm <- seq_len(n)
    for (i in seq_len(n)) {
      j <- i + floor(stats::runif(1L) * (n - i + 1L))
      tmp <- perm[i]
      perm[i] <- perm[j]
      perm[j] <- tmp
    }
    class_order[perm]
  })
  if (!is.null(group_ids)) {
    group_ids <- sb_validate_group_ids(group_ids, n)
  }
  grouped <- !is.null(group_ids) && anyDuplicated(group_ids) > 0L
  if (grouped) {
    fold_id <- sb_make_group_stratified_fold_id(
      outcome, group_ids, 5L, seed, "stormbreaker_infeasible_inner_folds",
      require_both_classes = FALSE
    )
  } else {
    fold_id <- integer(n)
    for (fold in seq_len(5L)) {
      begin <- floor((fold - 1L) * n / 5L)
      end <- floor(fold * n / 5L)
      if (end > begin) fold_id[shuffled[seq.int(begin + 1L, end)]] <- fold
    }
  }
  list(fold_id = fold_id, order = shuffled, class_order = class_order,
       source = if (grouped) "group_stratified" else "libsvm_row_shuffle")
}

sb_fit_libsvm_sigmoid <- function(decision, positive) {
  prior1 <- sum(positive)
  prior0 <- length(positive) - prior1
  target <- ifelse(positive, (prior1 + 1) / (prior1 + 2), 1 / (prior0 + 2))
  loss <- function(A, B) {
    z <- decision * A + B
    sum(ifelse(z >= 0, target * z, (target - 1) * z) +
          log1p(exp(-abs(z))))
  }
  A <- 0
  B <- log((prior0 + 1) / (prior1 + 1))
  value <- loss(A, B)
  termination <- "max_iterations"
  for (iteration in seq_len(100L)) {
    z <- decision * A + B
    small <- exp(-abs(z))
    p <- ifelse(z >= 0, small / (1 + small), 1 / (1 + small))
    q <- ifelse(z >= 0, 1 / (1 + small), small / (1 + small))
    d2 <- p * q
    h11 <- 1e-12 + sum(decision^2 * d2)
    h22 <- 1e-12 + sum(d2)
    h21 <- sum(decision * d2)
    g1 <- sum(decision * (target - p))
    g2 <- sum(target - p)
    if (abs(g1) < 1e-5 && abs(g2) < 1e-5) {
      termination <- "gradient_tolerance"
      break
    }
    determinant <- h11 * h22 - h21^2
    dA <- -(h22 * g1 - h21 * g2) / determinant
    dB <- -(-h21 * g1 + h11 * g2) / determinant
    directional <- g1 * dA + g2 * dB
    step <- 1
    while (step >= 1e-10) {
      newA <- A + step * dA
      newB <- B + step * dB
      newvalue <- loss(newA, newB)
      if (newvalue < value + 0.0001 * step * directional) {
        A <- newA
        B <- newB
        value <- newvalue
        break
      }
      step <- step / 2
    }
    if (step < 1e-10) {
      termination <- "line_search_failed"
      break
    }
  }
  list(A = A, B = B, termination = termination, iterations = iteration,
       negative_log_likelihood = value,
       positive_target = (prior1 + 1) / (prior1 + 2),
       negative_target = 1 / (prior0 + 2))
}

sb_calibrate_svm_probability <- function(x, outcome, learner, params,
                                          patient_ids, group_ids = NULL,
                                          selector = NULL) {
  inputs <- sb_validate_classical_binary_inputs(x, outcome, learner)
  sb_require_classical_binary_dependency(learner)
  x <- inputs$x
  outcome <- inputs$outcome
  patient_ids <- as.character(patient_ids)
  if (length(patient_ids) != nrow(x) || anyNA(patient_ids) ||
      any(!nzchar(trimws(patient_ids))) || anyDuplicated(patient_ids)) {
    sb_abort("SVM calibration requires a unique analysis ID for each row.",
             "stormbreaker_bad_runtime_patient_ids")
  }
  if (!is.null(group_ids)) {
    group_ids <- sb_validate_group_ids(group_ids, nrow(x))
  }
  folds <- sb_make_svm_calibration_folds(outcome, params$seed, group_ids)
  positive <- outcome[[1L]]
  decision <- rep(NA_real_, nrow(x))
  records <- vector("list", 5L)
  warnings <- character()
  fit_params <- params
  fit_params$probability <- FALSE
  for (fold in seq_len(5L)) {
    train <- folds$order[folds$fold_id[folds$order] != fold]
    assess <- folds$order[folds$fold_id[folds$order] == fold]
    # Native binary subproblems use +1/-1, with +1 first after class grouping.
    train <- c(train[outcome[train] == positive],
               train[outcome[train] != positive])
    features <- colnames(x)
    center <- scale <- class_weights <- NULL
    fit_warnings <- character()
    if (!length(assess)) {
      fit_status <- "empty_assessment"
    } else if (length(unique(outcome[train])) < 2L) {
      # Preserve LIBSVM's defined constant margin for single-class subproblems.
      decision[assess] <- if (!length(train)) 0 else {
        if (outcome[train[[1L]]] == positive) 1 else -1
      }
      features <- character()
      fit_status <- "single_class_constant"
    } else {
      train_x <- x[train, , drop = FALSE]
      assess_x <- x[assess, , drop = FALSE]
      model <- withCallingHandlers({
        if (!is.null(selector)) {
          selected <- sb_fit_runtime_selector(selector, train_x,
              list(task_type = "binary", outcome = outcome[train], group_ids = group_ids[train]))
          train_x <- selector$apply(selected, train_x)
          assess_x <- selector$apply(selected, assess_x)
        }
        sb_fit_classical_binary(train_x, outcome[train], learner,
                                       fit_params)
      }, warning = function(warning) {
        fit_warnings <<- c(fit_warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      })
      features <- model$features
      predicted <- stats::predict(model$engine, as.matrix(assess_x),
                                   decision.values = TRUE)
      native_positive <- as.character(model$engine$levels[
        model$engine$labels[[1L]]])
      direction <- if (native_positive == as.character(positive)) 1 else -1
      decision[assess] <- direction * as.numeric(attr(predicted,
                                                      "decision.values"))
      center <- model$engine$x.scale[["scaled:center"]]
      scale <- model$engine$x.scale[["scaled:scale"]]
      if (identical(params$class_weights, "inverse_frequency_normalized")) {
        counts <- table(factor(outcome[train], levels = c(0L, 1L)))
        class_weights <- stats::setNames(length(train) / (2 * as.numeric(counts)),
                                         names(counts))
      }
      fit_status <- "success"
    }
    records[[fold]] <- list(
      fold = fold, analysis_patient_ids = patient_ids[train],
      assessment_patient_ids = patient_ids[assess], selected_features = features,
      status = fit_status, scale_center = center, scale_scale = scale,
      class_weights = class_weights, warnings = fit_warnings
    )
    warnings <- c(warnings, if (length(fit_warnings)) {
      paste0("Calibration fold ", fold, ": ", fit_warnings)
    })
  }
  if (any(!is.finite(decision))) {
    sb_abort("SVM calibration did not produce finite out-of-fold margins.",
             "stormbreaker_svm_calibration_prediction")
  }
  # Retain LIBSVM's class-grouped order for its optimizer, not assessment order.
  sigmoid <- sb_fit_libsvm_sigmoid(decision[folds$class_order],
                                   outcome[folds$class_order] == positive)
  if (!identical(sigmoid$termination, "gradient_tolerance")) {
    warnings <- c(warnings, paste0("LIBSVM sigmoid: ", sigmoid$termination))
  }
  list(
    lambda = NA_real_, alpha = NA_real_, selection_rule = "libsvm_platt",
    measure = "smoothed_binomial_likelihood", inner_folds = 5L,
    seed = params$seed, fold_id = stats::setNames(folds$fold_id, patient_ids),
    alignment = "decision_value",
    inner_selection_method = if (is.null(selector)) "none" else {
      "svm_fold_local_probability_calibration"
    },
    inner_selections = records, curve = NULL, warnings = warnings,
    standardize = params$standardize,
    calibration = c(list(
      method = "libsvm_platt_oof", positive_class = as.character(positive),
      fold_source = folds$source,
      group_ids = stats::setNames(if (is.null(group_ids)) patient_ids else {
        group_ids
      }, patient_ids),
      out_of_fold = data.frame(patient_id = patient_ids, outcome = outcome,
        fold = folds$fold_id, decision = decision, stringsAsFactors = FALSE)
    ), sigmoid)
  )
}

sb_svm_calibrated_runtime_learner <- function(learner, calibration) {
  id <- learner$id
  params <- learner$params
  learner$fit_context <- calibration
  learner$fit <- function(x, outcome) {
    fit_params <- params
    fit_params$probability <- FALSE
    model <- sb_fit_classical_binary(x, outcome, id, fit_params)
    positive <- as.character(model$engine$levels[model$engine$labels[[1L]]])
    model$engine$probA <- calibration$A
    model$engine$probB <- if (positive == calibration$positive_class) {
      calibration$B
    } else -calibration$B
    model$engine$compprob <- TRUE
    model$params <- params
    model
  }
  learner
}
