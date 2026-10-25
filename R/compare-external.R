sb_external_comparison_value <- function(evaluation, path, default) {
  value <- evaluation
  for (name in path) {
    if (!is.list(value) || is.null(value[[name]])) {
      return(default)
    }
    value <- value[[name]]
  }
  if (length(value) != 1L) {
    return(default)
  }
  value[[1L]]
}

sb_external_binary_comparison_row <- function(evaluation) {
  value <- function(metric, field, default) {
    sb_external_comparison_value(
      evaluation,
      c("metrics", metric, field),
      default
    )
  }
  data.frame(
    pipeline_id = evaluation$pipeline_id,
    learner = evaluation$learner,
    status = evaluation$status,
    patient_match_status = evaluation$patient_match$status,
    roc_auc_status = value("roc_auc", "status", NA_character_),
    roc_auc = value("roc_auc", "estimate", NA_real_),
    comparable_pairs = value("roc_auc", "comparable_pairs", NA_real_),
    average_precision_status = value(
      "average_precision", "status", NA_character_
    ),
    average_precision = value(
      "average_precision", "estimate", NA_real_
    ),
    brier_status = value("brier_score", "status", NA_character_),
    brier_score = value("brier_score", "estimate", NA_real_),
    log_loss_status = value("log_loss", "status", NA_character_),
    log_loss = value("log_loss", "estimate", NA_real_),
    mean_calibration_status = value(
      "mean_calibration", "status", NA_character_
    ),
    mean_calibration_ratio = value(
      "mean_calibration", "mean_calibration_ratio", NA_real_
    ),
    calibration_status = value(
      "calibration", "status", NA_character_
    ),
    calibration_intercept = value(
      "calibration", "calibration_intercept", NA_real_
    ),
    calibration_slope = value(
      "calibration", "calibration_slope", NA_real_
    ),
    moderate_calibration_status = value(
      "moderate_calibration", "status", NA_character_
    ),
    moderate_calibration_ici = value(
      "moderate_calibration", "ici", NA_real_
    ),
    moderate_calibration_e50 = value(
      "moderate_calibration", "e50", NA_real_
    ),
    moderate_calibration_e90 = value(
      "moderate_calibration", "e90", NA_real_
    ),
    n = evaluation$n,
    error_stage = evaluation$error_stage,
    error_class = evaluation$error_class,
    error_message = evaluation$error_message,
    stringsAsFactors = FALSE
  )
}

sb_external_survival_comparison_row <- function(evaluation, has_horizon) {
  value <- function(metric, field, default) {
    sb_external_comparison_value(
      evaluation,
      c("metrics", metric, field),
      default
    )
  }
  row <- data.frame(
    pipeline_id = evaluation$pipeline_id,
    learner = evaluation$learner,
    status = evaluation$status,
    patient_match_status = evaluation$patient_match$status,
    harrell_c_status = value("harrell_c", "status", NA_character_),
    harrell_c = value("harrell_c", "estimate", NA_real_),
    comparable_pairs = value(
      "harrell_c", "comparable_pairs", NA_real_
    ),
    n = evaluation$n,
    error_stage = evaluation$error_stage,
    error_class = evaluation$error_class,
    error_message = evaluation$error_message,
    stringsAsFactors = FALSE
  )
  if (!has_horizon) {
    return(row)
  }

  weak_value <- function(component, field, default) {
    sb_external_comparison_value(
      evaluation,
      c("metrics", "weak_calibration", component, field),
      default
    )
  }
  data.frame(
    pipeline_id = row$pipeline_id,
    learner = row$learner,
    status = row$status,
    patient_match_status = row$patient_match_status,
    harrell_c_status = row$harrell_c_status,
    harrell_c = row$harrell_c,
    comparable_pairs = row$comparable_pairs,
    ipcw_brier_status = value(
      "ipcw_brier", "status", NA_character_
    ),
    ipcw_brier = value("ipcw_brier", "estimate", NA_real_),
    time_auc_status = value("time_auc", "status", NA_character_),
    time_auc = value("time_auc", "estimate", NA_real_),
    time_auc_pair_denominator = value(
      "time_auc", "pair_denominator", NA_real_
    ),
    uno_c_status = value("uno_c", "status", NA_character_),
    uno_c = value("uno_c", "estimate", NA_real_),
    uno_c_weighted_pair_denominator = value(
      "uno_c", "weighted_pair_denominator", NA_real_
    ),
    calibration_status = value(
      "calibration", "status", NA_character_
    ),
    calibration_ici = value("calibration", "ici", NA_real_),
    calibration_e50 = value("calibration", "e50", NA_real_),
    calibration_e90 = value("calibration", "e90", NA_real_),
    mean_calibration_status = weak_value(
      "mean_calibration", "status", NA_character_
    ),
    mean_calibration_ratio = sb_external_comparison_value(
      evaluation,
      c("metrics", "weak_calibration", "mean_calibration_ratio"),
      NA_real_
    ),
    calibration_slope_status = weak_value(
      "slope_fit", "status", NA_character_
    ),
    calibration_slope = sb_external_comparison_value(
      evaluation,
      c("metrics", "weak_calibration", "calibration_slope"),
      NA_real_
    ),
    n = row$n,
    error_stage = row$error_stage,
    error_class = row$error_class,
    error_message = row$error_message,
    stringsAsFactors = FALSE
  )
}

sb_external_prediction_target_state <- function(prediction, endpoint) {
  if (!is.list(prediction)) {
    return(list(known = FALSE))
  }
  task_type <- prediction[["task_type"]]
  if (identical(task_type, "binary")) {
    return(NULL)
  }
  if (!identical(task_type, "survival")) {
    return(list(known = FALSE))
  }

  contract_attempt <- sb_external_attempt(
    sb_external_prediction_contract(prediction)
  )
  if (is.null(contract_attempt$error)) {
    contract <- contract_attempt$value
    return(list(
      known = TRUE,
      has_horizon = contract$has_horizon,
      prediction_horizon = contract$prediction_horizon,
      time_unit = if (contract$has_horizon) {
        contract$time_unit
      } else {
        endpoint$time_unit
      }
    ))
  }

  signature <- prediction[["signature"]]
  if (inherits(signature, "stormbreaker_linear_signature")) {
    signature_attempt <- sb_external_attempt(
      sb_validate_linear_signature(signature)
    )
    if (is.null(signature_attempt$error) &&
        identical(signature_attempt$value$task_type, "survival")) {
      signature <- signature_attempt$value
      has_horizon <- "prediction_horizon" %in% names(signature)
      return(list(
        known = TRUE,
        has_horizon = has_horizon,
        prediction_horizon = if (has_horizon) {
          signature$prediction_horizon
        } else {
          NULL
        },
        time_unit = if (has_horizon) {
          signature$time_unit
        } else {
          endpoint$time_unit
        }
      ))
    }
  }

  rows <- prediction[["predictions"]]
  if (!is.data.frame(rows)) {
    return(list(known = FALSE))
  }
  horizon_fields <- c(
    "prediction_horizon", "time_unit", "survival_probability"
  )
  horizon_counts <- vapply(horizon_fields, function(name) {
    sum(names(rows) == name)
  }, integer(1L))
  if (all(horizon_counts == 0L) && sum(names(rows) == "risk") == 1L) {
    return(list(
      known = TRUE,
      has_horizon = FALSE,
      prediction_horizon = NULL,
      time_unit = endpoint$time_unit
    ))
  }
  if (!all(horizon_counts == 1L) || nrow(rows) == 0L) {
    return(list(known = FALSE))
  }

  horizon <- rows[["prediction_horizon"]]
  time_unit <- rows[["time_unit"]]
  valid_horizon <- is.numeric(horizon) && length(horizon) == nrow(rows) &&
    !anyNA(horizon) && all(is.finite(horizon)) && all(horizon > 0) &&
    length(unique(as.numeric(horizon))) == 1L
  valid_unit <- is.character(time_unit) && length(time_unit) == nrow(rows) &&
    !anyNA(time_unit) && length(unique(time_unit)) == 1L &&
    unique(time_unit) %in% c("day", "month", "year")
  if (!valid_horizon || !valid_unit) {
    return(list(known = FALSE))
  }
  list(
    known = TRUE,
    has_horizon = TRUE,
    prediction_horizon = as.numeric(horizon[[1L]]),
    time_unit = time_unit[[1L]]
  )
}

sb_external_comparison_target <- function(predictions, endpoint) {
  if (identical(endpoint$task_type, "binary")) {
    return(list(
      target_status = "complete",
      has_horizon = FALSE,
      prediction_horizon = NULL,
      time_unit = NULL,
      horizon_source = "not_applicable"
    ))
  }
  if (!is.null(endpoint$prediction_horizon)) {
    return(list(
      target_status = "complete",
      has_horizon = TRUE,
      prediction_horizon = endpoint$prediction_horizon,
      time_unit = endpoint$time_unit,
      horizon_source = "endpoint"
    ))
  }

  states <- lapply(predictions, function(prediction) {
    sb_external_prediction_target_state(prediction, endpoint)
  })
  states <- Filter(Negate(is.null), states)
  known <- Filter(function(state) isTRUE(state$known), states)
  if (length(known) == 0L) {
    return(list(
      target_status = "unknown",
      has_horizon = NA,
      prediction_horizon = NULL,
      time_unit = endpoint$time_unit,
      horizon_source = "unknown"
    ))
  }

  reference <- known[[1L]]
  same_target <- vapply(known, function(target) {
    same_presence <- identical(target$has_horizon, reference$has_horizon)
    same_presence && (!reference$has_horizon || (
      identical(target$time_unit, reference$time_unit) &&
        identical(
          as.numeric(target$prediction_horizon),
          as.numeric(reference$prediction_horizon)
        )
    ))
  }, logical(1L))
  if (!all(same_target)) {
    sb_abort(
      paste0(
        "Survival predictions compared without an endpoint horizon must all ",
        "be risk-only or share one locked time unit and horizon."
      ),
      "stormbreaker_external_target_mismatch"
    )
  }
  reference$known <- NULL
  reference$target_status <- if (length(known) == length(states)) {
    "complete"
  } else {
    "partial"
  }
  reference$horizon_source <- if (reference$has_horizon) {
    "prediction"
  } else {
    "none"
  }
  reference
}

#' Compare locked predictions on one external endpoint
#'
#' `sb_compare_external()` evaluates a non-empty set of already computed final
#' predictions against the same outcome-only external endpoint. It calls
#' [sb_evaluate_external()] once per prediction, so exact patient-ID alignment,
#' horizon semantics, metric definitions, and failure behavior are unchanged.
#' The function has no access to fitted models or external predictor values and
#' does not refit, recalibrate, or replace a learner.
#'
#' The returned table copies scalar metrics from each complete evaluation.
#' Every full `stormbreaker_external_evaluation` remains available in
#' `evaluations`; a failed or incomplete evaluation stays visible and does not
#' erase the other predictions. Optional decision curves are retained as named
#' complete objects rather than flattened or integrated.
#'
#' For survival, the comparison retains one effective `prediction_target`. A
#' horizon declared by the endpoint is authoritative. If the endpoint omits a
#' horizon, a horizon already locked into the predictions is retained instead;
#' all otherwise evaluable survival predictions must then be either risk-only
#' or share that same time unit and horizon. This prevents unlike time targets
#' from being presented as one comparison while preserving the endpoint's
#' original declaration separately in `endpoint`. A failed prediction can
#' still contribute target metadata when its locked horizon and time unit are
#' recoverable. Otherwise the target is marked `partial` or `unknown` rather
#' than being mislabeled as risk-only.
#'
#' This is a descriptive comparison of predictions that must have been locked
#' before consulting the external outcomes. It does not rank pipelines, select
#' a winner or threshold, estimate paired uncertainty, or prove
#' transportability. In particular, the external endpoint must not be used to
#' choose which fitted pipeline becomes the final model.
#'
#' @param predictions A non-empty list of distinct final-prediction objects
#'   returned by [sb_predict_final()] or [sb_predict_signature()]. Each object
#'   keeps its own failure status, but every available pipeline ID must be
#'   distinct.
#' @param endpoint One outcome-only object returned by [sb_binary_endpoint()] or
#'   [sb_surv_endpoint()]. The complete patient-ID set is reused for every
#'   prediction.
#'
#' @return A `stormbreaker_external_comparison` containing a descriptive
#'   `comparison` table, every complete named external `evaluation`, the task
#'   type, effective `prediction_target`, endpoint summary, and patient count.
#'   When present, named complete `decision_curves` are copied from the
#'   evaluations.
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
#' sb_compare_external(predictions, endpoint)
sb_compare_external <- function(predictions, endpoint) {
  valid_collection <- is.list(predictions) && length(predictions) > 0L &&
    all(vapply(
      predictions,
      function(prediction) {
        inherits(prediction, "stormbreaker_final_prediction")
      },
      logical(1L)
    ))
  if (!valid_collection) {
    sb_abort(
      paste0(
        "`predictions` must be a non-empty list returned by ",
        "`sb_predict_final()` or `sb_predict_signature()`."
      ),
      "stormbreaker_bad_prediction_collection"
    )
  }

  pipeline_ids <- vapply(predictions, function(prediction) {
    pipeline_id <- prediction[["pipeline_id"]]
    if (!is.character(pipeline_id) || length(pipeline_id) != 1L ||
        is.na(pipeline_id) || !nzchar(pipeline_id)) {
      return(NA_character_)
    }
    pipeline_id
  }, character(1L))
  if (anyNA(pipeline_ids)) {
    sb_abort(
      "Every final prediction must retain one non-empty pipeline ID.",
      "stormbreaker_bad_prediction_collection"
    )
  }
  if (anyDuplicated(pipeline_ids) > 0L) {
    sb_abort(
      "Every compared external prediction must have a distinct pipeline ID.",
      "stormbreaker_duplicate_pipeline_id"
    )
  }

  endpoint_contract <- sb_external_endpoint_contract(endpoint)
  prediction_target <- sb_external_comparison_target(
    predictions,
    endpoint_contract
  )
  evaluations <- lapply(predictions, function(prediction) {
    sb_evaluate_external(prediction, endpoint, reverse = FALSE)
  })
  names(evaluations) <- pipeline_ids
  has_horizon <- isTRUE(prediction_target$has_horizon)
  rows <- lapply(evaluations, function(evaluation) {
    if (identical(endpoint_contract$task_type, "binary")) {
      sb_external_binary_comparison_row(evaluation)
    } else {
      sb_external_survival_comparison_row(evaluation, has_horizon)
    }
  })
  comparison <- do.call(rbind, rows)
  rownames(comparison) <- NULL
  comparison$prediction_variation <- vapply(evaluations, function(evaluation) {
    sb_prediction_variation(evaluation$scored_predictions, endpoint_contract$task_type)$status
  }, character(1))

  complete <- all(vapply(
    evaluations,
    function(evaluation) identical(evaluation$status, "success"),
    logical(1L)
  ))
  result <- list(
    status = if (complete) "success" else "incomplete",
    task_type = endpoint_contract$task_type,
    comparison = comparison,
    evaluations = evaluations,
    n_patients = as.integer(length(endpoint_contract$patient_ids)),
    prediction_target = prediction_target,
    endpoint = endpoint_contract$summary
  )
  has_decision_curve <- any(vapply(
    evaluations,
    function(evaluation) {
      !is.null(evaluation$metrics[["decision_curve"]])
    },
    logical(1L)
  ))
  if (has_decision_curve) {
    decision_curves <- lapply(evaluations, function(evaluation) {
      evaluation$metrics[["decision_curve"]]
    })
    names(decision_curves) <- pipeline_ids
    result$decision_curves <- decision_curves
  }
  structure(result, class = "stormbreaker_external_comparison")
}

#' @export
print.stormbreaker_external_comparison <- function(x, ...) {
  cat("<stormbreaker_external_comparison:", x$status, ">\n")
  cat(
    "Pipelines:", nrow(x$comparison),
    "| Patients:", x$n_patients,
    "| Task:", x$task_type, "\n"
  )
  metric_columns <- if (identical(x$task_type, "binary")) {
    c("roc_auc", "average_precision", "brier_score", "log_loss")
  } else {
    columns <- "harrell_c"
    if ("ipcw_brier" %in% names(x$comparison)) {
      columns <- c(columns, "ipcw_brier", "time_auc", "uno_c")
    }
    columns
  }
  print(
    x$comparison[c("pipeline_id", "status", metric_columns)],
    row.names = FALSE,
    digits = 3
  )
  cat(
    paste0(
      "No pipeline is ranked or selected; external outcomes must not ",
      "choose the final model.\n"
    )
  )
  invisible(x)
}
