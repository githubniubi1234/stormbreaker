sb_binary_auc_stratum <- function(outcome, probability) {
  outcome <- as.integer(outcome)
  n_positive <- as.integer(sum(outcome == 1L))
  n_negative <- as.integer(sum(outcome == 0L))
  comparable_pairs <- unname(n_positive * n_negative)
  result <- list(
    status = "not_evaluable",
    estimate = NA_real_,
    n = as.integer(length(outcome)),
    n_positive = n_positive,
    n_negative = n_negative,
    correctly_ordered_pairs = 0,
    tied_pairs = 0,
    comparable_pairs = comparable_pairs
  )
  if (n_positive == 0L || n_negative == 0L) {
    return(result)
  }

  ranks <- rank(probability, ties.method = "average")
  pair_credit <- sum(ranks[outcome == 1L]) -
    n_positive * (n_positive + 1) / 2
  ordered <- order(probability)
  runs <- rle(probability[ordered])
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  ordered_outcome <- outcome[ordered]
  tied_pairs <- sum(vapply(seq_along(starts), function(index) {
    values <- ordered_outcome[starts[[index]]:ends[[index]]]
    sum(values == 1L) * sum(values == 0L)
  }, numeric(1)))

  result$status <- "success"
  result$estimate <- unname(pair_credit / comparable_pairs)
  result$correctly_ordered_pairs <- unname(pair_credit - 0.5 * tied_pairs)
  result$tied_pairs <- unname(tied_pairs)
  result
}

sb_binary_auc_metric <- function(predictions) {
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  per_fold <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      c("outcome", "probability", "status"),
      drop = FALSE
    ]
    audit <- data.frame(
      fold = fold_id,
      status = "prediction_failed",
      n = nrow(fold_predictions),
      n_positive = as.integer(sum(fold_predictions$outcome == 1L)),
      n_negative = as.integer(sum(fold_predictions$outcome == 0L)),
      estimate = NA_real_,
      correctly_ordered_pairs = 0,
      tied_pairs = 0,
      comparable_pairs = 0,
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    prediction_failed <- anyNA(fold_predictions$status) ||
      any(fold_predictions$status != "success") ||
      anyNA(fold_predictions$probability) ||
      any(!is.finite(fold_predictions$probability))
    if (prediction_failed) {
      audit$error_message <- paste0(
        "OOF probabilities for this fold are failed or non-finite."
      )
      return(audit)
    }
    stratum <- sb_binary_auc_stratum(
      outcome = fold_predictions$outcome,
      probability = fold_predictions$probability
    )
    audit$n <- stratum$n
    audit$n_positive <- stratum$n_positive
    audit$n_negative <- stratum$n_negative
    audit$correctly_ordered_pairs <- stratum$correctly_ordered_pairs
    audit$tied_pairs <- stratum$tied_pairs
    audit$comparable_pairs <- stratum$comparable_pairs
    if (!identical(stratum$status, "success")) {
      audit$status <- "not_evaluable"
      audit$error_message <- "Both outcome classes are required in this fold."
      return(audit)
    }

    audit$status <- "success"
    audit$estimate <- stratum$estimate
    audit
  })
  per_fold <- do.call(rbind, per_fold)
  rownames(per_fold) <- NULL
  evaluable <- per_fold$status == "success"
  totals <- c(
    correctly_ordered_pairs = 0,
    tied_pairs = 0,
    comparable_pairs = 0
  )
  if (any(evaluable)) {
    totals <- colSums(per_fold[
      evaluable,
      names(totals),
      drop = FALSE
    ])
  }
  complete <- length(fold_ids) > 0L && all(evaluable)
  estimate <- if (complete) {
    unname(
      (totals[["correctly_ordered_pairs"]] +
         0.5 * totals[["tied_pairs"]]) /
        totals[["comparable_pairs"]]
    )
  } else {
    NA_real_
  }

  list(
    metric = "roc_auc",
    scope = "within_outer_assessment_fold",
    aggregation = "summed_positive_negative_pairs",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    status = if (complete) "success" else "incomplete",
    estimate = estimate,
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    correctly_ordered_pairs = unname(
      totals[["correctly_ordered_pairs"]]
    ),
    tied_pairs = unname(totals[["tied_pairs"]]),
    comparable_pairs = unname(totals[["comparable_pairs"]]),
    per_fold = per_fold
  )
}

sb_empty_binary_precision_recall_curve <- function() {
  data.frame(
    fold = integer(),
    threshold = numeric(),
    n_at_threshold = integer(),
    positive_at_threshold = integer(),
    negative_at_threshold = integer(),
    true_positive = integer(),
    false_positive = integer(),
    precision = numeric(),
    recall = numeric(),
    recall_increment = numeric(),
    contribution = numeric(),
    stringsAsFactors = FALSE
  )
}

sb_binary_precision_recall_curve <- function(
    outcome,
    probability,
    fold_id) {
  ordered <- order(probability, decreasing = TRUE)
  ordered_probability <- probability[ordered]
  ordered_outcome <- as.integer(outcome[ordered])
  score_runs <- rle(ordered_probability)
  ends <- cumsum(score_runs$lengths)
  starts <- ends - score_runs$lengths + 1L
  positive_at_threshold <- vapply(seq_along(starts), function(index) {
    sum(ordered_outcome[starts[[index]]:ends[[index]]] == 1L)
  }, integer(1))
  negative_at_threshold <- score_runs$lengths - positive_at_threshold
  true_positive <- cumsum(positive_at_threshold)
  false_positive <- cumsum(negative_at_threshold)
  n_positive <- sum(ordered_outcome == 1L)
  precision <- true_positive / (true_positive + false_positive)
  recall <- true_positive / n_positive
  recall_increment <- c(recall[[1L]], diff(recall))

  data.frame(
    fold = rep(as.integer(fold_id), length(score_runs$values)),
    threshold = as.numeric(score_runs$values),
    n_at_threshold = as.integer(score_runs$lengths),
    positive_at_threshold = as.integer(positive_at_threshold),
    negative_at_threshold = as.integer(negative_at_threshold),
    true_positive = as.integer(true_positive),
    false_positive = as.integer(false_positive),
    precision = unname(precision),
    recall = unname(recall),
    recall_increment = unname(recall_increment),
    contribution = unname(recall_increment * precision),
    stringsAsFactors = FALSE
  )
}

sb_binary_average_precision_metric <- function(predictions) {
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  fold_results <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      c("outcome", "probability", "status"),
      drop = FALSE
    ]
    audit <- data.frame(
      fold = fold_id,
      status = "prediction_failed",
      n = as.integer(nrow(fold_predictions)),
      n_scored = 0L,
      n_positive = as.integer(sum(fold_predictions$outcome == 1L)),
      n_negative = as.integer(sum(fold_predictions$outcome == 0L)),
      n_thresholds = 0L,
      positive_weighted_sum = 0,
      estimate = NA_real_,
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    curve <- sb_empty_binary_precision_recall_curve()
    prediction_failed <- anyNA(fold_predictions$status) ||
      any(fold_predictions$status != "success") ||
      anyNA(fold_predictions$outcome) ||
      any(!fold_predictions$outcome %in% c(0L, 1L)) ||
      anyNA(fold_predictions$probability) ||
      any(!is.finite(fold_predictions$probability)) ||
      any(
        fold_predictions$probability < 0 |
          fold_predictions$probability > 1
      )
    if (prediction_failed) {
      audit$error_message <- paste0(
        "OOF outcomes or probabilities for this fold are failed or invalid."
      )
      return(list(audit = audit, curve = curve))
    }
    if (audit$n_positive == 0L || audit$n_negative == 0L) {
      audit$status <- "not_evaluable"
      audit$error_message <- "Both outcome classes are required in this fold."
      return(list(audit = audit, curve = curve))
    }

    curve <- sb_binary_precision_recall_curve(
      as.integer(fold_predictions$outcome),
      fold_predictions$probability,
      fold_id
    )
    positive_weighted_sum <- sum(
      curve$positive_at_threshold * curve$precision
    )
    audit$status <- "success"
    audit$n_scored <- audit$n
    audit$n_thresholds <- as.integer(nrow(curve))
    audit$positive_weighted_sum <- unname(positive_weighted_sum)
    audit$estimate <- unname(positive_weighted_sum / audit$n_positive)
    list(audit = audit, curve = curve)
  })
  per_fold <- do.call(rbind, lapply(fold_results, `[[`, "audit"))
  rownames(per_fold) <- NULL
  curves <- lapply(fold_results, `[[`, "curve")
  names(curves) <- as.character(fold_ids)
  evaluable <- per_fold$status == "success"
  complete <- length(fold_ids) > 0L && all(evaluable)
  positive_weighted_sum <- unname(sum(
    per_fold$positive_weighted_sum[evaluable]
  ))
  n_positive <- as.integer(sum(per_fold$n_positive))

  list(
    metric = "average_precision",
    scope = "within_outer_assessment_fold",
    aggregation = "positive_weighted_mean_of_fold_average_precision",
    definition = "sum((recall_i - recall_previous) * precision_i)",
    interpolation = "none",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    status = if (complete) "success" else "incomplete",
    estimate = if (complete) {
      unname(positive_weighted_sum / n_positive)
    } else {
      NA_real_
    },
    n = as.integer(sum(per_fold$n)),
    n_scored = as.integer(sum(per_fold$n_scored[evaluable])),
    n_positive = n_positive,
    n_negative = as.integer(sum(per_fold$n_negative)),
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    positive_weighted_sum = positive_weighted_sum,
    per_fold = per_fold,
    curves = curves
  )
}

sb_binary_brier_reference <- function(predictions, fold_ids, reference_input) {
  n <- as.integer(nrow(predictions))
  if (is.null(reference_input)) {
    return(list(
      type = "outer_analysis_prevalence",
      source = "outer_analysis_outcome_only",
      status = "not_requested",
      estimate = NA_real_,
      n = n,
      n_scored = 0L,
      n_folds = as.integer(length(fold_ids)),
      n_evaluable_folds = 0L,
      squared_error_sum = 0,
      per_fold = NULL
    ))
  }

  valid_schema <- is.data.frame(reference_input) &&
    all(c("fold", "analysis_prevalence") %in% names(reference_input))
  per_fold <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      c("outcome"),
      drop = FALSE
    ]
    audit <- data.frame(
      fold = fold_id,
      status = "reference_failed",
      analysis_prevalence = NA_real_,
      n = as.integer(nrow(fold_predictions)),
      n_scored = 0L,
      squared_error_sum = 0,
      estimate = NA_real_,
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    matched <- if (valid_schema) {
      which(as.integer(reference_input$fold) == fold_id)
    } else {
      integer()
    }
    if (length(matched) != 1L) {
      audit$error_message <- paste0(
        "Exactly one outer-analysis prevalence is required for this fold."
      )
      return(audit)
    }
    probability <- reference_input$analysis_prevalence[[matched]]
    valid_probability <- is.numeric(probability) && length(probability) == 1L &&
      !is.na(probability) && is.finite(probability) &&
      probability >= 0 && probability <= 1
    valid_outcome <- !anyNA(fold_predictions$outcome) &&
      all(fold_predictions$outcome %in% c(0L, 1L))
    if (!valid_probability || !valid_outcome) {
      audit$error_message <- paste0(
        "The analysis prevalence or assessment outcomes are invalid."
      )
      return(audit)
    }

    squared_error <- (
      as.integer(fold_predictions$outcome) - probability
    )^2
    audit$status <- "success"
    audit$analysis_prevalence <- probability
    audit$n_scored <- audit$n
    audit$squared_error_sum <- unname(sum(squared_error))
    audit$estimate <- unname(mean(squared_error))
    audit
  })
  per_fold <- do.call(rbind, per_fold)
  rownames(per_fold) <- NULL
  evaluable <- per_fold$status == "success"
  complete <- length(fold_ids) > 0L && all(evaluable)
  n_scored <- as.integer(sum(per_fold$n_scored[evaluable]))
  if (complete) {
    matched <- match(
      as.integer(predictions$fold),
      as.integer(reference_input$fold)
    )
    probability <- reference_input$analysis_prevalence[matched]
    squared_error <- (
      as.integer(predictions$outcome) - probability
    )^2
    squared_error_sum <- unname(sum(squared_error))
    estimate <- unname(mean(squared_error))
  } else {
    squared_error_sum <- unname(sum(
      per_fold$squared_error_sum[evaluable]
    ))
    estimate <- NA_real_
  }

  list(
    type = "outer_analysis_prevalence",
    source = "outer_analysis_outcome_only",
    status = if (complete) "success" else "incomplete",
    estimate = estimate,
    n = n,
    n_scored = n_scored,
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    squared_error_sum = squared_error_sum,
    per_fold = per_fold
  )
}

sb_binary_brier_metric <- function(predictions, reference_input = NULL) {
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  per_fold <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      c("outcome", "probability", "status"),
      drop = FALSE
    ]
    audit <- data.frame(
      fold = fold_id,
      status = "prediction_failed",
      n = as.integer(nrow(fold_predictions)),
      n_scored = 0L,
      squared_error_sum = 0,
      estimate = NA_real_,
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    prediction_failed <- anyNA(fold_predictions$status) ||
      any(fold_predictions$status != "success") ||
      anyNA(fold_predictions$outcome) ||
      any(!fold_predictions$outcome %in% c(0L, 1L)) ||
      anyNA(fold_predictions$probability) ||
      any(!is.finite(fold_predictions$probability)) ||
      any(
        fold_predictions$probability < 0 |
          fold_predictions$probability > 1
      )
    if (prediction_failed) {
      audit$error_message <- paste0(
        "OOF outcomes or probabilities for this fold are failed or invalid."
      )
      return(audit)
    }

    squared_error <- (
      as.integer(fold_predictions$outcome) - fold_predictions$probability
    )^2
    audit$status <- "success"
    audit$n_scored <- audit$n
    audit$squared_error_sum <- unname(sum(squared_error))
    audit$estimate <- unname(mean(squared_error))
    audit
  })
  per_fold <- do.call(rbind, per_fold)
  rownames(per_fold) <- NULL
  evaluable <- per_fold$status == "success"
  complete <- length(fold_ids) > 0L && all(evaluable)
  n <- as.integer(sum(per_fold$n))
  n_scored <- as.integer(sum(per_fold$n_scored[evaluable]))
  if (complete) {
    squared_error <- (
      as.integer(predictions$outcome) - predictions$probability
    )^2
    squared_error_sum <- unname(sum(squared_error))
    estimate <- unname(mean(squared_error))
  } else {
    squared_error_sum <- unname(sum(
      per_fold$squared_error_sum[evaluable]
    ))
    estimate <- NA_real_
  }
  reference <- sb_binary_brier_reference(
    predictions,
    fold_ids,
    reference_input
  )
  skill_complete <- complete && identical(reference$status, "success") &&
    is.finite(reference$squared_error_sum) &&
    reference$squared_error_sum > 0
  skill_score <- list(
    metric = "brier_skill_score",
    reference = reference$type,
    formula = "1 - model_squared_error / reference_squared_error",
    status = if (skill_complete) "success" else "incomplete",
    estimate = if (skill_complete) {
      unname(1 - squared_error_sum / reference$squared_error_sum)
    } else {
      NA_real_
    }
  )

  list(
    metric = "brier_score",
    scope = "patient_level_out_of_fold",
    aggregation = "summed_squared_error_over_patients",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    status = if (complete) "success" else "incomplete",
    estimate = estimate,
    n = n,
    n_scored = n_scored,
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    squared_error_sum = squared_error_sum,
    per_fold = per_fold,
    reference = reference,
    skill_score = skill_score
  )
}

sb_binary_case_log_loss <- function(outcome, probability) {
  outcome <- as.integer(outcome)
  positive <- outcome == 1L
  case_loss <- numeric(length(outcome))
  case_loss[positive] <- -log(probability[positive])
  case_loss[!positive] <- -log1p(-probability[!positive])
  unname(case_loss)
}

sb_binary_log_loss_reference <- function(predictions, fold_ids, reference_input) {
  n <- as.integer(nrow(predictions))
  if (is.null(reference_input)) {
    return(list(
      type = "outer_analysis_prevalence",
      source = "outer_analysis_outcome_only",
      probability_clipping = "none",
      status = "not_requested",
      estimate = NA_real_,
      n = n,
      n_scored = 0L,
      n_finite_loss = 0L,
      n_infinite_loss = 0L,
      n_folds = as.integer(length(fold_ids)),
      n_evaluable_folds = 0L,
      loss_sum = 0,
      per_fold = NULL
    ))
  }

  valid_schema <- is.data.frame(reference_input) &&
    all(c("fold", "analysis_prevalence") %in% names(reference_input))
  per_fold <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      "outcome",
      drop = FALSE
    ]
    audit <- data.frame(
      fold = fold_id,
      status = "reference_failed",
      analysis_prevalence = NA_real_,
      n = as.integer(nrow(fold_predictions)),
      n_scored = 0L,
      n_finite_loss = 0L,
      n_infinite_loss = 0L,
      loss_sum = 0,
      estimate = NA_real_,
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    matched <- if (valid_schema) {
      which(as.integer(reference_input$fold) == fold_id)
    } else {
      integer()
    }
    if (length(matched) != 1L) {
      audit$error_message <- paste0(
        "Exactly one outer-analysis prevalence is required for this fold."
      )
      return(audit)
    }
    probability <- reference_input$analysis_prevalence[[matched]]
    valid_probability <- is.numeric(probability) && length(probability) == 1L &&
      !is.na(probability) && is.finite(probability) &&
      probability >= 0 && probability <= 1
    valid_outcome <- !anyNA(fold_predictions$outcome) &&
      all(fold_predictions$outcome %in% c(0L, 1L))
    if (!valid_probability || !valid_outcome) {
      audit$error_message <- paste0(
        "The analysis prevalence or assessment outcomes are invalid."
      )
      return(audit)
    }

    case_loss <- sb_binary_case_log_loss(
      fold_predictions$outcome,
      rep(probability, nrow(fold_predictions))
    )
    audit$status <- "success"
    audit$analysis_prevalence <- probability
    audit$n_scored <- audit$n
    audit$n_finite_loss <- as.integer(sum(is.finite(case_loss)))
    audit$n_infinite_loss <- as.integer(sum(is.infinite(case_loss)))
    audit$loss_sum <- unname(sum(case_loss))
    audit$estimate <- unname(mean(case_loss))
    audit
  })
  per_fold <- do.call(rbind, per_fold)
  rownames(per_fold) <- NULL
  evaluable <- per_fold$status == "success"
  complete <- length(fold_ids) > 0L && all(evaluable)
  n_scored <- as.integer(sum(per_fold$n_scored[evaluable]))
  loss_sum <- unname(sum(per_fold$loss_sum[evaluable]))

  list(
    type = "outer_analysis_prevalence",
    source = "outer_analysis_outcome_only",
    probability_clipping = "none",
    status = if (complete) "success" else "incomplete",
    estimate = if (complete) unname(loss_sum / n) else NA_real_,
    n = n,
    n_scored = n_scored,
    n_finite_loss = as.integer(sum(per_fold$n_finite_loss[evaluable])),
    n_infinite_loss = as.integer(sum(per_fold$n_infinite_loss[evaluable])),
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    loss_sum = loss_sum,
    per_fold = per_fold
  )
}

sb_binary_log_loss_metric <- function(predictions, reference_input = NULL) {
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  per_fold <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      c("outcome", "probability", "status"),
      drop = FALSE
    ]
    audit <- data.frame(
      fold = fold_id,
      status = "prediction_failed",
      n = as.integer(nrow(fold_predictions)),
      n_scored = 0L,
      n_finite_loss = 0L,
      n_infinite_loss = 0L,
      loss_sum = 0,
      estimate = NA_real_,
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    prediction_failed <- anyNA(fold_predictions$status) ||
      any(fold_predictions$status != "success") ||
      anyNA(fold_predictions$outcome) ||
      any(!fold_predictions$outcome %in% c(0L, 1L)) ||
      anyNA(fold_predictions$probability) ||
      any(!is.finite(fold_predictions$probability)) ||
      any(
        fold_predictions$probability < 0 |
          fold_predictions$probability > 1
      )
    if (prediction_failed) {
      audit$error_message <- paste0(
        "OOF outcomes or probabilities for this fold are failed or invalid."
      )
      return(audit)
    }

    case_loss <- sb_binary_case_log_loss(
      fold_predictions$outcome,
      fold_predictions$probability
    )
    audit$status <- "success"
    audit$n_scored <- audit$n
    audit$n_finite_loss <- as.integer(sum(is.finite(case_loss)))
    audit$n_infinite_loss <- as.integer(sum(is.infinite(case_loss)))
    audit$loss_sum <- unname(sum(case_loss))
    audit$estimate <- unname(mean(case_loss))
    audit
  })
  per_fold <- do.call(rbind, per_fold)
  rownames(per_fold) <- NULL
  evaluable <- per_fold$status == "success"
  complete <- length(fold_ids) > 0L && all(evaluable)
  n <- as.integer(sum(per_fold$n))
  n_scored <- as.integer(sum(per_fold$n_scored[evaluable]))
  loss_sum <- unname(sum(per_fold$loss_sum[evaluable]))
  reference <- sb_binary_log_loss_reference(
    predictions,
    fold_ids,
    reference_input
  )
  skill_complete <- complete && identical(reference$status, "success") &&
    is.finite(reference$loss_sum) && reference$loss_sum > 0
  skill_score <- list(
    metric = "log_loss_skill_score",
    reference = reference$type,
    formula = "1 - model_log_loss / reference_log_loss",
    status = if (skill_complete) "success" else "incomplete",
    estimate = if (skill_complete) {
      unname(1 - loss_sum / reference$loss_sum)
    } else {
      NA_real_
    }
  )

  list(
    metric = "log_loss",
    scope = "patient_level_out_of_fold",
    aggregation = "summed_case_loss_over_patients",
    definition = "negative_binary_log_likelihood",
    probability_clipping = "none",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    metric_direction = "lower_is_better",
    status = if (complete) "success" else "incomplete",
    estimate = if (complete) unname(loss_sum / n) else NA_real_,
    n = n,
    n_scored = n_scored,
    n_finite_loss = as.integer(sum(per_fold$n_finite_loss[evaluable])),
    n_infinite_loss = as.integer(sum(per_fold$n_infinite_loss[evaluable])),
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    loss_sum = loss_sum,
    per_fold = per_fold,
    reference = reference,
    skill_score = skill_score
  )
}

sb_binary_empty_calibration_fit <- function(fit_type, terms) {
  list(
    fit_type = fit_type,
    status = "not_fitted",
    coefficients = data.frame(
      term = terms,
      estimate = rep(NA_real_, length(terms)),
      stringsAsFactors = FALSE
    ),
    converged = NA,
    iterations = 0L,
    boundary = NA,
    deviance = NA_real_,
    warnings = character(),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_binary_calibration_has_separation <- function(outcome, logit_probability) {
  if (length(unique(logit_probability)) < 2L) {
    return(FALSE)
  }
  positive <- logit_probability[outcome == 1L]
  negative <- logit_probability[outcome == 0L]
  min(positive) >= max(negative) || max(positive) <= min(negative)
}

sb_binary_fit_calibration <- function(outcome, logit_probability, fit_type) {
  intercept_only <- identical(fit_type, "calibration_in_the_large")
  terms <- if (intercept_only) {
    "(Intercept)"
  } else {
    c("(Intercept)", "logit_probability")
  }
  audit <- sb_binary_empty_calibration_fit(fit_type, terms)
  if (!intercept_only && sb_binary_calibration_has_separation(
    outcome,
    logit_probability
  )) {
    audit$status <- "failed"
    audit$error_class <- "stormbreaker_calibration_separation"
    audit$error_message <- paste0(
      "Complete or quasi-complete separation means that no finite ",
      "calibration slope maximum-likelihood estimate exists."
    )
    return(audit)
  }
  design <- if (intercept_only) {
    matrix(
      1,
      nrow = length(outcome),
      ncol = 1L,
      dimnames = list(NULL, "(Intercept)")
    )
  } else {
    cbind(
      "(Intercept)" = 1,
      logit_probability = logit_probability
    )
  }
  offset <- if (intercept_only) {
    logit_probability
  } else {
    rep(0, length(outcome))
  }
  captured_warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      stats::glm.fit(
        x = design,
        y = as.integer(outcome),
        family = stats::binomial(),
        offset = offset,
        control = stats::glm.control()
      ),
      warning = function(warning) {
        captured_warnings <<- c(
          captured_warnings,
          conditionMessage(warning)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) error
  )
  audit$warnings <- unique(captured_warnings)
  if (inherits(fit, "condition")) {
    audit$status <- "failed"
    audit$error_class <- class(fit)[[1L]]
    audit$error_message <- conditionMessage(fit)
    return(audit)
  }

  audit$coefficients$estimate <- unname(fit$coefficients[terms])
  audit$converged <- isTRUE(fit$converged)
  audit$iterations <- as.integer(fit$iter)
  audit$boundary <- isTRUE(fit$boundary)
  audit$deviance <- unname(fit$deviance)
  if (!audit$converged) {
    audit$status <- "failed"
    audit$error_class <-
      "stormbreaker_binary_calibration_nonconvergence"
    audit$error_message <- paste0(
      "The ", fit_type, " calibration fit did not converge."
    )
    return(audit)
  }
  if (any(!is.finite(audit$coefficients$estimate))) {
    audit$status <- "failed"
    audit$error_class <-
      "stormbreaker_binary_calibration_nonfinite_coefficients"
    audit$error_message <- paste0(
      "The ", fit_type,
      " calibration fit returned non-finite coefficients."
    )
    return(audit)
  }
  audit$status <- "success"
  audit
}

sb_binary_calibration_metric <- function(predictions) {
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  per_fold <- lapply(fold_ids, function(fold_id) {
    fold_predictions <- predictions[
      predictions$fold == fold_id,
      c("outcome", "probability", "status"),
      drop = FALSE
    ]
    prediction_failed <- anyNA(fold_predictions$status) ||
      any(fold_predictions$status != "success") ||
      anyNA(fold_predictions$outcome) ||
      any(!fold_predictions$outcome %in% c(0L, 1L)) ||
      anyNA(fold_predictions$probability) ||
      any(!is.finite(fold_predictions$probability)) ||
      any(
        fold_predictions$probability < 0 |
          fold_predictions$probability > 1
      )
    data.frame(
      fold = fold_id,
      status = if (prediction_failed) {
        "prediction_failed"
      } else {
        "success"
      },
      n = as.integer(nrow(fold_predictions)),
      n_scored = if (prediction_failed) {
        0L
      } else {
        as.integer(nrow(fold_predictions))
      },
      n_positive = as.integer(sum(
        fold_predictions$outcome == 1L,
        na.rm = TRUE
      )),
      n_negative = as.integer(sum(
        fold_predictions$outcome == 0L,
        na.rm = TRUE
      )),
      n_boundary_probability = if (prediction_failed) {
        0L
      } else {
        as.integer(sum(fold_predictions$probability %in% c(0, 1)))
      },
      error_message = if (prediction_failed) {
        "OOF outcomes or probabilities for this fold are failed or invalid."
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  })
  per_fold <- do.call(rbind, per_fold)
  rownames(per_fold) <- NULL
  evaluable_folds <- per_fold$status == "success"
  complete_predictions <- length(fold_ids) > 0L && all(evaluable_folds)
  n <- as.integer(sum(per_fold$n))
  n_scored <- as.integer(sum(per_fold$n_scored[evaluable_folds]))
  n_positive <- as.integer(sum(per_fold$n_positive))
  n_negative <- as.integer(sum(per_fold$n_negative))
  n_boundary_probability <- as.integer(sum(
    per_fold$n_boundary_probability[evaluable_folds]
  ))
  intercept_fit <- sb_binary_empty_calibration_fit(
    "calibration_in_the_large",
    "(Intercept)"
  )
  slope_fit <- sb_binary_empty_calibration_fit(
    "logistic_recalibration",
    c("(Intercept)", "logit_probability")
  )
  result <- list(
    metric = "calibration_intercept_slope",
    scope = "pooled_patient_level_out_of_fold",
    purpose = "descriptive_evaluation_only",
    probability_clipping = "none",
    positive_class = 1L,
    prediction_direction = "higher_is_more_likely_positive",
    status = "incomplete",
    input_status = if (complete_predictions) {
      "success"
    } else {
      "prediction_failed"
    },
    calibration_intercept = NA_real_,
    calibration_slope = NA_real_,
    intercept_target = 0,
    slope_target = 1,
    n = n,
    n_scored = n_scored,
    n_positive = n_positive,
    n_negative = n_negative,
    n_boundary_probability = n_boundary_probability,
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable_folds)),
    per_fold = per_fold,
    intercept_fit = intercept_fit,
    slope_fit = slope_fit,
    error_class = NA_character_,
    error_message = NA_character_
  )
  if (!complete_predictions) {
    result$error_message <- paste0(
      "Complete successful OOF probabilities are required for calibration."
    )
    return(result)
  }
  if (n_positive == 0L || n_negative == 0L) {
    result$input_status <- "not_evaluable"
    result$error_message <- paste0(
      "Both outcome classes are required for calibration."
    )
    return(result)
  }
  if (n_boundary_probability > 0L) {
    result$input_status <- "not_evaluable"
    result$error_class <-
      "stormbreaker_binary_calibration_boundary_probability"
    result$error_message <- paste0(
      "Exact boundary probabilities cannot enter a finite-logit calibration ",
      "fit without probability clipping, which is not performed."
    )
    return(result)
  }

  logit_probability <- stats::qlogis(predictions$probability)
  intercept_fit <- sb_binary_fit_calibration(
    predictions$outcome,
    logit_probability,
    "calibration_in_the_large"
  )
  slope_fit <- sb_binary_fit_calibration(
    predictions$outcome,
    logit_probability,
    "logistic_recalibration"
  )
  result$intercept_fit <- intercept_fit
  result$slope_fit <- slope_fit
  if (identical(intercept_fit$status, "success")) {
    result$calibration_intercept <-
      intercept_fit$coefficients$estimate[[1L]]
  }
  if (identical(slope_fit$status, "success")) {
    result$calibration_slope <- slope_fit$coefficients$estimate[[2L]]
  }
  fits_complete <- identical(intercept_fit$status, "success") &&
    identical(slope_fit$status, "success")
  if (fits_complete) {
    result$status <- "success"
    return(result)
  }
  failed_fit <- if (!identical(intercept_fit$status, "success")) {
    intercept_fit
  } else {
    slope_fit
  }
  result$error_class <- failed_fit$error_class
  result$error_message <- paste(
    stats::na.omit(c(intercept_fit$error_message, slope_fit$error_message)),
    collapse = " "
  )
  result
}
