sb_external_paired_metric_spec <- function(task_type, has_horizon) {
  if (identical(task_type, "binary")) {
    return(data.frame(
      metric = c(
        "roc_auc", "average_precision", "brier_score", "log_loss"
      ),
      direction = c(
        "higher_is_better", "higher_is_better",
        "lower_is_better", "lower_is_better"
      ),
      stringsAsFactors = FALSE
    ))
  }
  metric <- "harrell_c"
  direction <- "higher_is_better"
  if (has_horizon) {
    metric <- c(metric, "ipcw_brier_score", "time_auc", "uno_c")
    direction <- c(
      direction,
      "lower_is_better", "higher_is_better", "higher_is_better"
    )
  }
  data.frame(metric = metric, direction = direction, stringsAsFactors = FALSE)
}

sb_external_paired_metric_rows <- function(scored, task_type, metric_spec) {
  rows <- sb_external_bootstrap_metric_rows(scored, task_type)
  positions <- match(metric_spec$metric, rows$metric)
  if (anyNA(positions)) {
    sb_abort(
      "A required paired-comparison metric is unavailable from the runtime.",
      "stormbreaker_missing_paired_metric"
    )
  }
  rows[positions, , drop = FALSE]
}

sb_external_paired_candidate_contract <- function(evaluation, reference) {
  candidate <- sb_external_bootstrap_contract(evaluation)
  same_design <- identical(candidate$task_type, reference$task_type) &&
    identical(candidate$n, reference$n) &&
    identical(candidate$has_horizon, reference$has_horizon) &&
    identical(candidate$prediction_horizon, reference$prediction_horizon) &&
    identical(candidate$time_unit, reference$time_unit)
  endpoint_columns <- if (identical(reference$task_type, "binary")) {
    c("patient_id", "outcome")
  } else {
    c("patient_id", "time", "event")
  }
  same_endpoint <- same_design && identical(
    candidate$scored[endpoint_columns],
    reference$scored[endpoint_columns]
  )
  if (!same_endpoint) {
    sb_abort(
      paste0(
        "Every paired external evaluation must contain the same patients ",
        "and endpoint rows in the same endpoint-defined order."
      ),
      "stormbreaker_unpaired_external_evaluation"
    )
  }
  candidate
}

sb_external_paired_empty_differences <- function(replicates = FALSE) {
  result <- data.frame(
    reference_pipeline_id = character(),
    candidate_pipeline_id = character(),
    contrast = character(),
    metric = character(),
    metric_direction = character(),
    reference_estimate = numeric(),
    candidate_estimate = numeric(),
    estimate = numeric(),
    status = character(),
    error_class = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
  if (replicates) {
    result <- data.frame(
      replicate = integer(),
      result,
      stringsAsFactors = FALSE
    )
  }
  result
}

sb_external_paired_empty_intervals <- function() {
  data.frame(
    reference_pipeline_id = character(),
    candidate_pipeline_id = character(),
    contrast = character(),
    metric = character(),
    metric_direction = character(),
    reference_point_estimate = numeric(),
    candidate_point_estimate = numeric(),
    point_estimate = numeric(),
    lower = numeric(),
    upper = numeric(),
    confidence_level = numeric(),
    status = character(),
    interval_available = logical(),
    replicate_completeness = character(),
    n_requested = integer(),
    n_finite_successful = integer(),
    n_unavailable = integer(),
    n_nonfinite = integer(),
    valid_fraction = numeric(),
    error_class = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
}

sb_external_paired_difference_rows <- function(
    source,
    reference_pipeline_id,
    candidate_pipeline_id,
    metric_spec) {
  candidate <- source[
    source$pipeline_id == candidate_pipeline_id,
    ,
    drop = FALSE
  ]
  reference <- source[
    source$pipeline_id == reference_pipeline_id,
    ,
    drop = FALSE
  ]
  has_replicate <- "replicate" %in% names(source)
  candidate_key <- if (has_replicate) {
    paste(candidate$replicate, candidate$metric, sep = "\r")
  } else {
    candidate$metric
  }
  reference_key <- if (has_replicate) {
    paste(reference$replicate, reference$metric, sep = "\r")
  } else {
    reference$metric
  }
  reference <- reference[match(candidate_key, reference_key), , drop = FALSE]

  available <- !is.na(candidate$status) & !is.na(reference$status) &
    candidate$status == "success" & reference$status == "success"
  difference <- candidate$estimate - reference$estimate
  difference[!available] <- NA_real_
  error_class <- rep(NA_character_, nrow(candidate))
  error_message <- rep(NA_character_, nrow(candidate))
  for (index in which(!available)) {
    candidate_failed <- !identical(candidate$status[[index]], "success")
    source_row <- if (candidate_failed) candidate else reference
    error_class[[index]] <- source_row$error_class[[index]]
    error_message[[index]] <- source_row$error_message[[index]]
    if (is.na(error_class[[index]]) || !nzchar(error_class[[index]])) {
      error_class[[index]] <- "stormbreaker_unavailable_paired_metric"
    }
    if (is.na(error_message[[index]]) || !nzchar(error_message[[index]])) {
      error_message[[index]] <- paste0(
        "Both models require an available `", candidate$metric[[index]],
        "` estimate in the same patient draw."
      )
    }
  }
  result <- data.frame(
    reference_pipeline_id = rep(reference_pipeline_id, nrow(candidate)),
    candidate_pipeline_id = rep(candidate_pipeline_id, nrow(candidate)),
    contrast = rep("candidate_minus_reference", nrow(candidate)),
    metric = candidate$metric,
    metric_direction = metric_spec$direction[match(
      candidate$metric,
      metric_spec$metric
    )],
    reference_estimate = reference$estimate,
    candidate_estimate = candidate$estimate,
    estimate = difference,
    status = ifelse(available, "success", "not_evaluable"),
    error_class = error_class,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
  if (has_replicate) {
    result <- data.frame(
      replicate = candidate$replicate,
      result,
      stringsAsFactors = FALSE
    )
  }
  result
}

sb_external_paired_interval_rows <- function(
    point_differences,
    differences,
    confidence_level,
    reps) {
  interval_input <- point_differences[c(
    "metric", "estimate", "status", "error_class", "error_message"
  )]
  replicate_input <- differences[c(
    "replicate", "metric", "estimate", "status",
    "error_class", "error_message"
  )]
  intervals <- sb_external_percentile_intervals(
    interval_input,
    replicate_input,
    confidence_level,
    reps
  )
  positions <- match(intervals$metric, point_differences$metric)
  data.frame(
    reference_pipeline_id = point_differences$reference_pipeline_id[positions],
    candidate_pipeline_id = point_differences$candidate_pipeline_id[positions],
    contrast = point_differences$contrast[positions],
    metric = intervals$metric,
    metric_direction = point_differences$metric_direction[positions],
    reference_point_estimate =
      point_differences$reference_estimate[positions],
    candidate_point_estimate =
      point_differences$candidate_estimate[positions],
    intervals[names(intervals) != "metric"],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

sb_external_paired_model_rows <- function(
    contract,
    metric_spec,
    pipeline_id,
    replicate = NULL,
    indices = NULL) {
  scored <- if (is.null(indices)) {
    contract$scored
  } else {
    sampled <- contract$scored[indices, , drop = FALSE]
    rownames(sampled) <- NULL
    sampled
  }
  rows <- sb_external_paired_metric_rows(
    scored,
    contract$task_type,
    metric_spec
  )
  result <- data.frame(
    pipeline_id = rep(pipeline_id, nrow(rows)),
    rows,
    stringsAsFactors = FALSE
  )
  if (!is.null(replicate)) {
    result <- data.frame(
      replicate = rep.int(as.integer(replicate), nrow(result)),
      result,
      stringsAsFactors = FALSE
    )
  }
  result
}

#' Paired external-patient bootstrap differences from one locked comparison
#'
#' `sb_bootstrap_external_comparison()` estimates paired differences between
#' every eligible candidate prediction in a completed [sb_compare_external()]
#' result and one explicitly named reference prediction. In every bootstrap
#' replicate, one patient-index vector is drawn and applied unchanged to the
#' reference and all candidates. Each difference is therefore calculated as
#' `candidate - reference` within the same resampled external patients.
#'
#' The reference must be chosen before inspecting external outcomes. The
#' function does not choose it, rank candidates, select a winner, or report a
#' p value. Binary comparisons are limited to ROC AUC, non-interpolated average
#' precision, Brier score, and unclipped log loss. Survival comparisons use
#' Harrell C and, when the shared endpoint has a locked horizon, IPCW Brier,
#' cumulative/dynamic time-dependent AUC, and truncated Uno C. The result labels
#' whether higher or lower values are preferable; this does not convert the
#' intervals into a model-selection rule.
#'
#' Every metric is recomputed from fixed, endpoint-aligned predictions in each
#' ordinary patient bootstrap draw. For horizon survival metrics, the external
#' marginal censoring distribution is re-estimated inside that draw. Percentile
#' intervals use finite successful paired differences and retain unavailable
#' and non-finite counts. They describe independent external-patient sampling
#' conditional on the locked models and predictions. They do not include
#' development fitting, feature selection, tuning, model-selection uncertainty,
#' multiplicity adjustment, transportability, or clinical-effectiveness proof.
#'
#' Calibration and decision-curve contrasts are deliberately outside this
#' function: a raw difference is not a universal calibration-quality measure,
#' and model-to-model net-benefit comparison requires common predeclared
#' thresholds and separate reporting semantics. Use [sb_bootstrap_external()]
#' for each model's existing calibration and model-versus-strategy uncertainty.
#'
#' A failed candidate remains in `candidate_status` and cannot erase valid
#' candidate-reference comparisons. A failed reference cannot define shared
#' contrasts and is rejected.
#'
#' @param comparison A `stormbreaker_external_comparison` returned by
#'   [sb_compare_external()] with at least two distinct predictions.
#' @param reference_pipeline_id The exact pipeline ID of the reference
#'   prediction. It must be declared independently of external outcomes.
#' @param reps Number of ordinary patient-level bootstrap replicates. Must be an
#'   integer of at least 2. Defaults to 1000.
#' @param confidence_level Percentile interval confidence level strictly between
#'   0 and 1. Defaults to 0.95.
#' @param seed Positive integer seed. The caller's random state is restored.
#'
#' @return A `stormbreaker_external_bootstrap_comparison` containing the
#'   reference ID, explicit prediction target, candidate eligibility, shared
#'   resample audit, model-level point and replicate metrics, paired point and
#'   replicate differences, and type-7 percentile intervals for
#'   `candidate - reference`.
#' @export
#'
#' @examples
#' development <- sb_example_binary(n = 50, seed = 1)
#' task <- sb_binary_task(
#'   development,
#'   "patient_id",
#'   "outcome",
#'   c("gene_a", "gene_b", "gene_c")
#' )
#' pipelines <- list(
#'   sb_pipeline(learner = "logistic"),
#'   sb_pipeline(learner = "gaussian_naive_bayes")
#' )
#' fits <- lapply(pipelines, function(pipeline) {
#'   sb_fit_final(task, pipeline)
#' })
#' external <- sb_example_binary(n = 30, seed = 2)
#' external$patient_id <- paste0("external_", external$patient_id)
#' predictions <- lapply(fits, function(fit) {
#'   sb_predict_final(
#'     fit,
#'     external[c("patient_id", "gene_a", "gene_b", "gene_c")],
#'     "patient_id"
#'   )
#' })
#' endpoint <- sb_binary_endpoint(external, "patient_id", "outcome")
#' comparison <- sb_compare_external(predictions, endpoint)
#' paired <- sb_bootstrap_external_comparison(
#'   comparison,
#'   reference_pipeline_id = predictions[[1L]]$pipeline_id,
#'   reps = 20,
#'   seed = 3
#' )
#' paired$intervals
sb_bootstrap_external_comparison <- function(
    comparison,
    reference_pipeline_id,
    reps = 1000L,
    confidence_level = 0.95,
    seed = 1L) {
  valid_comparison <- inherits(
    comparison,
    "stormbreaker_external_comparison"
  ) && is.list(comparison) && is.list(comparison$evaluations) &&
    length(comparison$evaluations) >= 2L
  if (!valid_comparison) {
    sb_abort(
      paste0(
        "`comparison` must be returned by `sb_compare_external()` and ",
        "contain at least two predictions."
      ),
      "stormbreaker_bad_external_comparison"
    )
  }
  pipeline_ids <- vapply(comparison$evaluations, function(evaluation) {
    pipeline_id <- evaluation[["pipeline_id"]]
    if (sb_external_scalar_character(pipeline_id)) {
      pipeline_id
    } else {
      NA_character_
    }
  }, character(1L))
  valid_reference <- sb_external_scalar_character(reference_pipeline_id) &&
    reference_pipeline_id %in% pipeline_ids
  if (!valid_reference) {
    sb_abort(
      "`reference_pipeline_id` must identify one compared prediction.",
      "stormbreaker_bad_reference_pipeline"
    )
  }
  reps <- sb_check_positive_integer(reps, "reps", minimum = 2L)
  seed <- sb_check_positive_integer(seed, "seed")
  valid_confidence <- is.numeric(confidence_level) &&
    length(confidence_level) == 1L && !is.na(confidence_level) &&
    is.finite(confidence_level) && confidence_level > 0 &&
    confidence_level < 1
  if (!valid_confidence) {
    sb_abort(
      "`confidence_level` must be one finite number strictly between 0 and 1.",
      "stormbreaker_bad_confidence_level"
    )
  }
  confidence_level <- as.numeric(confidence_level)

  reference_position <- match(reference_pipeline_id, pipeline_ids)
  reference <- sb_external_bootstrap_contract(
    comparison$evaluations[[reference_position]]
  )
  metric_spec <- sb_external_paired_metric_spec(
    reference$task_type,
    reference$has_horizon
  )
  candidate_ids <- pipeline_ids[-reference_position]
  candidate_evaluations <- comparison$evaluations[-reference_position]
  candidate_attempts <- lapply(candidate_evaluations, function(evaluation) {
    sb_external_attempt(
      sb_external_paired_candidate_contract(evaluation, reference)
    )
  })
  candidate_status <- do.call(rbind, lapply(
    seq_along(candidate_ids),
    function(index) {
      error <- candidate_attempts[[index]]$error
      data.frame(
        candidate_pipeline_id = candidate_ids[[index]],
        status = if (is.null(error)) "success" else "failed",
        error_class = if (is.null(error)) {
          NA_character_
        } else {
          class(error)[[1L]]
        },
        error_message = if (is.null(error)) {
          NA_character_
        } else {
          conditionMessage(error)
        },
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(candidate_status) <- NULL
  eligible <- which(candidate_status$status == "success")
  contracts <- list(reference)
  names(contracts) <- reference_pipeline_id
  if (length(eligible) > 0L) {
    eligible_contracts <- lapply(
      candidate_attempts[eligible],
      `[[`,
      "value"
    )
    names(eligible_contracts) <- candidate_ids[eligible]
    contracts <- c(contracts, eligible_contracts)
  }

  model_point_estimates <- do.call(rbind, lapply(names(contracts), function(id) {
    sb_external_paired_model_rows(
      contracts[[id]],
      metric_spec,
      id
    )
  }))
  rownames(model_point_estimates) <- NULL
  index_plan <- withr::with_seed(seed, {
    lapply(seq_len(reps), function(replicate) {
      sample.int(reference$n, reference$n, replace = TRUE)
    })
  })
  sampled <- lapply(seq_len(reps), function(replicate) {
      indices <- index_plan[[replicate]]
      reference_sample <- reference$scored[indices, , drop = FALSE]
      rownames(reference_sample) <- NULL
      list(
        rows = do.call(rbind, lapply(names(contracts), function(id) {
          sb_external_paired_model_rows(
            contracts[[id]],
            metric_spec,
            id,
            replicate,
            indices
          )
        })),
        audit = sb_external_bootstrap_audit_row(
          reference_sample,
          reference$task_type,
          replicate
        )
      )
  })
  model_replicates <- do.call(rbind, lapply(sampled, `[[`, "rows"))
  resample_audit <- do.call(rbind, lapply(sampled, `[[`, "audit"))
  rownames(model_replicates) <- NULL
  rownames(resample_audit) <- NULL

  point_pieces <- lapply(candidate_ids[eligible], function(candidate_id) {
    sb_external_paired_difference_rows(
      model_point_estimates,
      reference_pipeline_id,
      candidate_id,
      metric_spec
    )
  })
  difference_pieces <- lapply(candidate_ids[eligible], function(candidate_id) {
    sb_external_paired_difference_rows(
      model_replicates,
      reference_pipeline_id,
      candidate_id,
      metric_spec
    )
  })
  if (length(point_pieces) == 0L) {
    point_differences <- sb_external_paired_empty_differences()
    differences <- sb_external_paired_empty_differences(replicates = TRUE)
    intervals <- sb_external_paired_empty_intervals()
  } else {
    point_differences <- do.call(rbind, point_pieces)
    differences <- do.call(rbind, difference_pieces)
    rownames(point_differences) <- NULL
    rownames(differences) <- NULL
    interval_pieces <- lapply(candidate_ids[eligible], function(candidate_id) {
      point_rows <- point_differences[
        point_differences$candidate_pipeline_id == candidate_id,
        ,
        drop = FALSE
      ]
      difference_rows <- differences[
        differences$candidate_pipeline_id == candidate_id,
        ,
        drop = FALSE
      ]
      sb_external_paired_interval_rows(
        point_rows,
        difference_rows,
        confidence_level,
        reps
      )
    })
    intervals <- do.call(rbind, interval_pieces)
    rownames(intervals) <- NULL
  }

  complete <- all(candidate_status$status == "success") &&
    nrow(intervals) > 0L && all(intervals$status == "success")
  result <- list(
    status = if (complete) "success" else "incomplete",
    task_type = reference$task_type,
    reference_pipeline_id = reference_pipeline_id,
    candidate_status = candidate_status,
    n = reference$n,
    prediction_target = reference$prediction_target,
    reps = reps,
    confidence_level = confidence_level,
    seed = seed,
    metrics = metric_spec$metric,
    pairing = list(
      resampling_unit = "patient",
      draw_sharing = "same_patient_indices_for_all_models",
      contrast = "candidate_minus_reference"
    ),
    interval_method = sb_external_bootstrap_interval_method(confidence_level),
    resampling = sb_external_bootstrap_resampling(
      reference$task_type,
      reference$has_horizon
    ),
    model_point_estimates = model_point_estimates,
    point_differences = point_differences,
    model_replicates = model_replicates,
    differences = differences,
    intervals = intervals,
    resample_audit = resample_audit,
    provenance = sb_external_bootstrap_provenance()
  )
  structure(result, class = "stormbreaker_external_bootstrap_comparison")
}

#' @export
print.stormbreaker_external_bootstrap_comparison <- function(x, ...) {
  cat("<stormbreaker_external_bootstrap_comparison:", x$status, ">\n")
  cat("paired fixed-prediction external patient bootstrap\n")
  cat("Reference:", x$reference_pipeline_id, "| Task:", x$task_type, "\n")
  cat(
    "Patients:", x$n,
    "| Replicates:", x$reps,
    "| Confidence:", format(x$confidence_level), "\n"
  )
  if (nrow(x$intervals) > 0L) {
    print(
      x$intervals[c(
        "candidate_pipeline_id", "metric", "point_estimate",
        "lower", "upper", "status"
      )],
      row.names = FALSE,
      digits = 3
    )
  }
  failed <- sum(x$candidate_status$status == "failed")
  if (failed > 0L) {
    cat("Failed candidates retained:", failed, "\n")
  }
  cat(
    paste0(
      "No ranking, winner, or p value is produced; the reference must be ",
      "declared before external outcomes are inspected.\n"
    )
  )
  invisible(x)
}
