sb_external_attempt <- function(code) {
  tryCatch(
    list(value = force(code), error = NULL),
    error = function(error) list(value = NULL, error = error)
  )
}

sb_external_scalar_character <- function(value) {
  is.character(value) && length(value) == 1L &&
    !is.na(value) && nzchar(value)
}

sb_external_column <- function(data, name) {
  positions <- which(names(data) == name)
  if (length(positions) != 1L) {
    sb_abort(
      paste0("External prediction column `", name, "` is missing or duplicated."),
      "stormbreaker_bad_external_prediction"
    )
  }
  data[[positions]]
}

sb_external_prediction_contract <- function(prediction) {
  if (!inherits(prediction, "stormbreaker_final_prediction") ||
      !is.list(prediction)) {
    sb_abort(
      paste0(
        "`prediction` must be returned by `sb_predict_final()`, ",
        "`sb_predict_stack_final()`, or `sb_predict_signature()`."
      ),
      "stormbreaker_bad_external_prediction"
    )
  }
  if (identical(prediction[["status"]], "failed")) {
    sb_abort(
      "A failed final prediction cannot be externally evaluated.",
      "stormbreaker_failed_external_prediction"
    )
  }
  task_type <- prediction[["task_type"]]
  pipeline_id <- prediction[["pipeline_id"]]
  learner <- prediction[["learner"]]
  id_col <- prediction[["id_col"]]
  rows <- prediction[["predictions"]]
  valid_metadata <- sb_external_scalar_character(task_type) &&
    task_type %in% c("binary", "survival") &&
    sb_external_scalar_character(pipeline_id) &&
    sb_external_scalar_character(learner) &&
    sb_external_scalar_character(id_col) &&
    identical(prediction[["status"]], "success") &&
    is.data.frame(rows) && nrow(rows) >= 1L
  if (!valid_metadata) {
    sb_abort(
      "The final prediction metadata are incomplete or inconsistent.",
      "stormbreaker_bad_external_prediction"
    )
  }

  patient_ids <- as.character(sb_external_column(rows, "patient_id"))
  row_pipeline <- sb_external_column(rows, "pipeline_id")
  row_learner <- sb_external_column(rows, "learner")
  row_status <- sb_external_column(rows, "status")
  valid_rows <- length(patient_ids) == nrow(rows) &&
    !anyNA(patient_ids) && all(nzchar(trimws(patient_ids))) &&
    !anyDuplicated(patient_ids) &&
    is.character(row_pipeline) &&
    length(row_pipeline) == nrow(rows) &&
    !anyNA(row_pipeline) && all(row_pipeline == pipeline_id) &&
    is.character(row_learner) && length(row_learner) == nrow(rows) &&
    !anyNA(row_learner) && all(row_learner == learner) &&
    is.character(row_status) && length(row_status) == nrow(rows) &&
    !anyNA(row_status) && all(row_status == "success")
  if (!valid_rows) {
    sb_abort(
      "External prediction rows have invalid IDs, provenance, or status.",
      "stormbreaker_bad_external_prediction"
    )
  }

  if (identical(task_type, "binary")) {
    positive_class <- sb_external_column(rows, "positive_class")
    probability <- sb_external_column(rows, "probability")
    valid_probability <- is.numeric(probability) &&
      length(probability) == nrow(rows) &&
      !anyNA(probability) && all(is.finite(probability)) &&
      all(probability >= 0 & probability <= 1)
    valid_positive <- is.numeric(positive_class) &&
      length(positive_class) == nrow(rows) &&
      !anyNA(positive_class) && all(positive_class == 1)
    if (!valid_probability || !valid_positive) {
      sb_abort(
        "External binary predictions must be finite class-one probabilities.",
        "stormbreaker_bad_external_prediction"
      )
    }
    decision_thresholds <- if (
      "decision_thresholds" %in% names(prediction)
    ) {
      sb_validate_decision_thresholds(
        prediction[["decision_thresholds"]]
      )
    } else {
      NULL
    }
    return(list(
      task_type = task_type,
      pipeline_id = pipeline_id,
      learner = learner,
      id_col = id_col,
      patient_ids = patient_ids,
      probability = unname(probability),
      decision_thresholds = decision_thresholds,
      has_horizon = FALSE,
      prediction_horizon = NULL,
      time_unit = NULL
    ))
  }

  risk <- sb_external_column(rows, "risk")
  valid_risk <- is.numeric(risk) && length(risk) == nrow(rows) &&
    !anyNA(risk) && all(is.finite(risk))
  if (!valid_risk) {
    sb_abort(
      "External survival risks must be finite numeric values.",
      "stormbreaker_bad_external_prediction"
    )
  }
  decision_thresholds <- if (
    "decision_thresholds" %in% names(prediction)
  ) {
    sb_validate_decision_thresholds(prediction[["decision_thresholds"]])
  } else {
    NULL
  }
  horizon_fields <- c(
    "prediction_horizon", "time_unit", "survival_probability"
  )
  horizon_present <- vapply(
    horizon_fields,
    function(name) sum(names(rows) == name) == 1L,
    logical(1)
  )
  horizon_named <- vapply(
    horizon_fields,
    function(name) sum(names(rows) == name),
    integer(1)
  )
  if (any(horizon_named > 1L) ||
      (any(horizon_present) && !all(horizon_present))) {
    sb_abort(
      "External survival horizon prediction columns are incomplete.",
      "stormbreaker_bad_external_prediction"
    )
  }
  if ("log_cumulative_hazard" %in% names(rows) && !all(horizon_present)) {
    sb_abort("A stable survival scale requires complete horizon predictions.",
      "stormbreaker_bad_external_prediction")
  }
  if (!any(horizon_present)) {
    if (!is.null(decision_thresholds)) {
      sb_abort(
        paste0(
          "Survival decision thresholds require locked-horizon ",
          "probabilities."
        ),
        "stormbreaker_bad_external_prediction"
      )
    }
    return(list(
      task_type = task_type,
      pipeline_id = pipeline_id,
      learner = learner,
      id_col = id_col,
      patient_ids = patient_ids,
      risk = unname(risk),
      has_horizon = FALSE,
      prediction_horizon = NULL,
      time_unit = NULL,
      survival_probability = NULL,
      decision_thresholds = NULL
    ))
  }

  horizon <- rows[["prediction_horizon"]]
  time_unit <- rows[["time_unit"]]
  probability <- rows[["survival_probability"]]
  valid_horizon <- is.numeric(horizon) && length(horizon) == nrow(rows) &&
    !anyNA(horizon) && all(is.finite(horizon)) && all(horizon > 0) &&
    length(unique(as.numeric(horizon))) == 1L
  valid_unit <- is.character(time_unit) && length(time_unit) == nrow(rows) &&
    !anyNA(time_unit) && length(unique(time_unit)) == 1L &&
    unique(time_unit) %in% c("day", "month", "year")
  valid_probability <- is.numeric(probability) &&
    length(probability) == nrow(rows) &&
    !anyNA(probability) && all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid_horizon || !valid_unit || !valid_probability) {
    sb_abort(
      paste0(
        "External survival horizon, unit, or patient probabilities are ",
        "invalid."
      ),
      "stormbreaker_bad_external_prediction"
    )
  }
  list(
    task_type = task_type,
    pipeline_id = pipeline_id,
    learner = learner,
    id_col = id_col,
    patient_ids = patient_ids,
    risk = unname(risk),
    has_horizon = TRUE,
    prediction_horizon = as.numeric(horizon[[1L]]),
    time_unit = time_unit[[1L]],
    survival_probability = unname(probability),
    log_cumulative_hazard = sb_external_stable_scale(rows, probability),
    decision_thresholds = decision_thresholds
  )
}

sb_external_endpoint_contract <- function(endpoint) {
  if (inherits(endpoint, "stormbreaker_binary_endpoint")) {
    endpoint <- sb_validate_binary_endpoint_contract(endpoint)
    patient_ids <- as.character(endpoint$data[[endpoint$id_col]])
    return(list(
      task_type = "binary",
      patient_ids = patient_ids,
      outcome = as.integer(endpoint$data[[endpoint$outcome_col]]),
      time = NULL,
      event = NULL,
      time_unit = NULL,
      prediction_horizon = NULL,
      summary = list(
        id_col = endpoint$id_col,
        outcome_col = endpoint$outcome_col,
        n = endpoint$n,
        n_positive = endpoint$n_positive,
        n_negative = endpoint$n_negative,
        positive_class = endpoint$positive_class,
        negative_class = endpoint$negative_class,
        independent_unit = endpoint$independent_unit
      )
    ))
  }
  if (inherits(endpoint, "stormbreaker_surv_endpoint")) {
    endpoint <- sb_validate_surv_endpoint_contract(endpoint)
    horizon <- if ("prediction_horizon" %in% names(endpoint)) {
      endpoint[["prediction_horizon"]]
    } else {
      NULL
    }
    summary <- list(
      id_col = endpoint$id_col,
      time_col = endpoint$time_col,
      event_col = endpoint$event_col,
      time_unit = endpoint$time_unit,
      n = endpoint$n,
      n_events = endpoint$n_events,
      n_censored = endpoint$n_censored,
      independent_unit = endpoint$independent_unit
    )
    if (!is.null(horizon)) {
      summary <- append(
        summary,
        list(prediction_horizon = horizon),
        after = match("time_unit", names(summary))
      )
    }
    return(list(
      task_type = "survival",
      patient_ids = as.character(endpoint$data[[endpoint$id_col]]),
      outcome = NULL,
      time = endpoint$data[[endpoint$time_col]],
      event = as.integer(endpoint$data[[endpoint$event_col]]),
      time_unit = endpoint$time_unit,
      prediction_horizon = horizon,
      summary = summary
    ))
  }
  sb_abort(
    paste0(
      "`endpoint` must be created by `sb_binary_endpoint()` or ",
      "`sb_surv_endpoint()`."
    ),
    "stormbreaker_bad_endpoint"
  )
}

sb_external_horizon_semantics <- function(prediction, endpoint) {
  if (!identical(prediction$task_type, "survival")) {
    return(invisible(TRUE))
  }
  endpoint_has_horizon <- !is.null(endpoint$prediction_horizon)
  if (!prediction$has_horizon) {
    if (endpoint_has_horizon) {
      sb_abort(
        paste0(
          "A risk-only prediction cannot evaluate a declared survival ",
          "horizon."
        ),
        "stormbreaker_external_horizon_mismatch"
      )
    }
    return(invisible(TRUE))
  }
  if (!identical(prediction$time_unit, endpoint$time_unit)) {
    sb_abort(
      "Prediction and endpoint survival time units do not match.",
      "stormbreaker_external_time_unit_mismatch"
    )
  }
  if (endpoint_has_horizon && !identical(
    as.numeric(prediction$prediction_horizon),
    as.numeric(endpoint$prediction_horizon)
  )) {
    sb_abort(
      "Prediction and endpoint survival horizons do not match.",
      "stormbreaker_external_horizon_mismatch"
    )
  }
  invisible(TRUE)
}

sb_align_external_patients <- function(prediction, endpoint) {
  exact <- length(prediction$patient_ids) == length(endpoint$patient_ids) &&
    setequal(prediction$patient_ids, endpoint$patient_ids)
  if (!exact) {
    sb_abort(
      paste0(
        "Prediction and endpoint patient IDs must form exactly the same set."
      ),
      "stormbreaker_external_patient_mismatch"
    )
  }
  matched <- match(endpoint$patient_ids, prediction$patient_ids)
  if (identical(prediction$task_type, "binary")) {
    scored <- data.frame(
      patient_id = endpoint$patient_ids,
      outcome = endpoint$outcome,
      probability = prediction$probability[matched],
      status = rep("success", length(matched)),
      stringsAsFactors = FALSE
    )
  } else {
    scored <- data.frame(
      patient_id = endpoint$patient_ids,
      time = endpoint$time,
      event = endpoint$event,
      risk = prediction$risk[matched],
      stringsAsFactors = FALSE
    )
    if (prediction$has_horizon) {
      scored$prediction_horizon <- rep(
        prediction$prediction_horizon,
        length(matched)
      )
      scored$time_unit <- rep(prediction$time_unit, length(matched))
      scored$survival_probability <-
        prediction$survival_probability[matched]
      if (!is.null(prediction$log_cumulative_hazard)) {
        scored$log_cumulative_hazard <- prediction$log_cumulative_hazard[matched]
      }
    }
    scored$status <- rep("success", length(matched))
  }
  list(
    scored = scored,
    patient_match = list(
      status = "success",
      prediction_n = as.integer(length(prediction$patient_ids)),
      endpoint_n = as.integer(length(endpoint$patient_ids)),
      exact_set = TRUE,
      output_order = "external_endpoint"
    )
  )
}

sb_external_evaluation_boundary <- function() {
  list(
    prediction_input = "precomputed_state_only_prediction",
    model_access = "none",
    outcome_use = "metric_evaluation_only",
    patient_alignment = "exact_character_id_match"
  )
}

sb_external_evaluation_provenance <- function() {
  list(
    package = "stormbreaker",
    package_version = as.character(utils::packageVersion("stormbreaker")),
    r_version = R.version.string,
    survival_version = as.character(utils::packageVersion("survival"))
  )
}

sb_external_safe_metadata <- function(prediction) {
  result <- list(
    task_type = NA_character_,
    pipeline_id = NA_character_,
    learner = NA_character_,
    prediction_n = 0L
  )
  if (!is.list(prediction)) {
    return(result)
  }
  for (name in c("task_type", "pipeline_id", "learner")) {
    value <- prediction[[name]]
    if (sb_external_scalar_character(value)) {
      result[[name]] <- value
    }
  }
  rows <- prediction[["predictions"]]
  if (is.data.frame(rows)) {
    result$prediction_n <- as.integer(nrow(rows))
  }
  result
}

sb_external_evaluation_failure <- function(
    prediction,
    stage,
    error,
    endpoint_n = 0L,
    patient_match = NULL) {
  metadata <- sb_external_safe_metadata(prediction)
  if (is.null(patient_match)) {
    patient_match <- list(
      status = "not_checked",
      prediction_n = metadata$prediction_n,
      endpoint_n = as.integer(endpoint_n),
      exact_set = NA,
      output_order = "external_endpoint"
    )
  }
  result <- list(
    status = "failed",
    task_type = metadata$task_type,
    pipeline_id = metadata$pipeline_id,
    learner = metadata$learner,
    n = as.integer(endpoint_n),
    evaluation_boundary = sb_external_evaluation_boundary(),
    patient_match = patient_match,
    scored_predictions = data.frame(),
    metrics = list(),
    endpoint = list(),
    provenance = sb_external_evaluation_provenance(),
    error_stage = stage,
    error_class = class(error)[[1L]],
    error_message = conditionMessage(error)
  )
  class(result) <- "stormbreaker_external_evaluation"
  result
}

#' Evaluate precomputed final predictions against an external endpoint
#'
#' `sb_evaluate_external()` accepts a completed endpoint-free prediction from
#' [sb_predict_final()] or [sb_predict_signature()] and a separate outcome-only
#' object from
#' [sb_binary_endpoint()] or [sb_surv_endpoint()]. Patients must match exactly
#' by ID. Outcomes are attached only after prediction.
#' This function has no access to a fit or pipeline
#' runtime and never refits, re-predicts, recalibrates, derives a threshold, or
#' selects a model.
#' The returned `scored_predictions` table follows the external endpoint's row
#' order after exact character-ID matching, which may differ from prediction
#' input order. This choice is recorded in `patient_match$output_order`.
#'
#' Binary evaluation reports external-cohort ROC AUC, non-interpolated average
#' precision, Brier score, unclipped log loss, observed-to-expected mean
#' calibration, descriptive weak-calibration intercept/slope, and a separate
#' moderate-calibration curve. Binary O:E is the observed event proportion
#' divided by mean predicted event probability and targets one. It is distinct
#' from the calibration intercept and cannot diagnose miscalibration within the
#' prediction range. A zero expected event count or a numerically non-finite
#' ratio leaves O:E unavailable without clipping or continuity correction. The
#' binary curve is a three-knot restricted natural cubic spline logistic model
#' of outcome on
#' `qlogis(probability)` and reports ICI, E50, and E90. Exact probabilities zero
#' or one are not clipped and leave only this finite-logit curve unavailable.
#' Complete or quasi-complete separation likewise leaves the curve unavailable,
#' including when logistic IRLS reports numerical convergence. The curve is
#' descriptive and is never applied as recalibration. When the development task
#' declared decision thresholds, binary evaluation also reports net benefit at
#' those unchanged thresholds using the action rule `probability >= threshold`,
#' with treat-all and treat-none references. The same TP/FP/TN/FN counts report
#' sensitivity, specificity, positive predictive value, and negative predictive
#' value. Zero denominators remain `NA`; predictive values retain the external
#' cohort prevalence. The function does not optimize, smooth, integrate, or
#' recommend a threshold. Survival
#' evaluation reports pooled Harrell C. If the saved
#' prediction contains a horizon, it also reports an IPCW Brier score, a
#' cumulative/dynamic time-dependent AUC, and truncated Uno IPCW concordance
#' using the external endpoint's marginal reverse Kaplan-Meier censoring
#' distribution. For time-dependent AUC, the
#' marker is `1 - S(t)`, observed events at or before the locked horizon are
#' cases, patients followed strictly beyond it are controls, and patients
#' censored at or before it are excluded. Cases receive inverse censoring
#' weights immediately before their event times and exact marker ties receive
#' half credit. This relies on the marginal non-informative censoring
#' assumption. Uno C instead uses the saved risk score: an eligible index has an
#' observed event strictly before the horizon, later observed times and
#' censoring at the same recorded time are comparators, and its weight is
#' 1 / G(T-)^2. Events tied at the same time are not compared, while exact risk
#' ties receive half credit. It does not use the survival probability.
#' If the horizon prediction also carries study-declared event-
#' risk thresholds, external evaluation reports a survival decision curve using
#' `1 - survival_probability >= threshold`. It estimates event probability at
#' the horizon by Kaplan-Meier within each threshold-selected external group,
#' multiplies by the selected fraction, and reports treat-all and treat-none
#' references. A nonempty group ending follow-up before the horizon is
#' unavailable unless its Kaplan-Meier survival already equals the absorbing
#' value zero. This point estimator assumes independent non-informative
#' censoring in the relevant risk groups and censoring independent of model
#' predictors. It is distinct from the marginal IPCW estimators and never
#' refits or changes a prediction. The same locked predictions and aligned
#' external outcomes also
#' enter a descriptive moderate-calibration curve: a three-knot restricted
#' natural cubic spline of `log(-log(S(t)))` is fitted in an Efron Cox model
#' using complete right-censored follow-up, and the curve is read at the locked
#' horizon. ICI, E50, and E90 summarize absolute error on the survival-
#' probability scale. This assumes proportional hazards and independent
#' non-informative censoring. The curve is not applied as recalibration,
#' probabilities are not clipped or changed, no universal acceptable cutoff is
#' asserted, and it remains separate from weak calibration. At the same locked
#' horizon, mean calibration is the Kaplan-Meier observed event probability
#' divided by mean predicted event probability (O:E), while the calibration
#' slope is the coefficient from an Efron Cox model after administrative
#' censoring at the horizon on `log(-log(survival_probability))`. Both target 1.
#' Exact boundary or constant predictions can leave only the slope unavailable;
#' unsupported Kaplan-Meier follow-up can leave only the mean ratio unavailable.
#' No probability is clipped or recalibrated, and no survival calibration
#' intercept, p value, or universal cutoff is reported. The slope assumes a
#' proportional-hazards recalibration relation, and both components assume
#' independent non-informative censoring. A calibration failure does not erase
#' otherwise valid external metrics.
#' No event-Kaplan-Meier reference or Brier skill score is fitted.
#'
#' Exact ID matching cannot prove that cohorts are independent or that outcome
#' definitions, time origins, eligibility, predictor measurement, and upstream
#' preprocessing are comparable and outcome-blind. Those study-design
#' conditions remain the researcher's responsibility. Returned values are
#' point estimates for the supplied endpoint. Time-dependent AUC measures
#' discrimination, Brier score measures overall prediction error, and the curve
#' measures calibration; none of these metrics by themselves
#' proves transportability, clinical utility, significance, or a decision
#' threshold. A decision curve is a point-estimate clinical-utility analysis;
#' this function does not derive a threshold, report decision-curve uncertainty,
#' or prove that an action is clinically effective.
#'
#' @references Vickers AJ, Cronin AM, Elkin EB, Gonen M (2008). Extensions to
#'   decision curve analysis. *BMC Medical Informatics and Decision Making*,
#'   8, 53. \doi{10.1186/1472-6947-8-53}.
#'
#' McLernon DJ, Giardiello D, Van Calster B, et al. (2023). Assessing
#'   performance and clinical usefulness in prediction models with survival
#'   outcomes: practical guidance for Cox proportional hazards models.
#'   *Annals of Internal Medicine*, 176(1), 105-114.
#'   \doi{10.7326/M22-0844}.
#'
#' Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ (2011). On the C-statistics
#'   for evaluating overall adequacy of risk prediction procedures with
#'   censored survival data. *Statistics in Medicine*, 30(10), 1105-1117.
#'   \doi{10.1002/sim.4154}.
#'
#' @param prediction A successful `stormbreaker_final_prediction` returned by
#'   [sb_predict_final()], [sb_predict_stack_final()], or
#'   [sb_predict_signature()].
#' @param endpoint An outcome-only object returned by [sb_binary_endpoint()] or
#'   [sb_surv_endpoint()] with exactly the same patient IDs.
#' @param reverse For binary outcomes, `TRUE` (default) uses cohort outcomes to
#'   flip probabilities when AUC < 0.5; `FALSE` keeps direction, and survival
#'   outcomes require `FALSE` (used when omitted).
#'
#' @return A `stormbreaker_external_evaluation` containing ID-aligned endpoint
#'   and prediction values, task-specific point-estimate metrics, the evaluation
#'   boundary audit, and explicit incomplete or failed status.
#' @details
#' For ordinary Cox final predictions whose probabilities round to zero or one,
#' the optional log_cumulative_hazard column preserves log(H0(h)) + x beta from
#' the saved model when finite. External time AUC and complementary-log-log
#' calibration, including bootstrap, use this monotone scale. Probabilities,
#' Brier scores and decision curves are not clipped or changed. Without this
#' saved scale, boundary-probability calibration keeps its existing unavailable
#' status; arbitrary risk scores are never substituted.
#' @export
#'
#' @examples
#' development <- sb_example_binary(n = 30, seed = 1)
#' task <- sb_binary_task(
#'   development, "patient_id", "outcome", c("gene_a", "gene_b")
#' )
#' fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
#' external <- sb_example_binary(n = 10, seed = 2)
#' external$patient_id <- paste0("external_", seq_len(nrow(external)))
#' prediction <- sb_predict_final(
#'   fit,
#'   external[c("patient_id", "gene_a", "gene_b")],
#'   "patient_id"
#' )
#' endpoint <- sb_binary_endpoint(external, "patient_id", "outcome")
#' sb_evaluate_external(prediction, endpoint)
sb_evaluate_external <- function(prediction, endpoint, reverse = TRUE) {
  if (missing(reverse) && inherits(endpoint, "stormbreaker_surv_endpoint")) {
    reverse <- FALSE
  }
  if (!is.logical(reverse) || length(reverse) != 1L || is.na(reverse)) {
    sb_abort("`reverse` must be TRUE or FALSE.", "stormbreaker_bad_reverse")
  }
  if (reverse) {
    return(sb_external_reverse_evaluation(prediction, endpoint))
  }
  prediction_attempt <- sb_external_attempt(
    sb_external_prediction_contract(prediction)
  )
  if (!is.null(prediction_attempt$error)) {
    return(sb_external_evaluation_failure(
      prediction,
      "prediction_contract",
      prediction_attempt$error
    ))
  }
  prediction_contract <- prediction_attempt$value

  endpoint_attempt <- sb_external_attempt(
    sb_external_endpoint_contract(endpoint)
  )
  if (!is.null(endpoint_attempt$error)) {
    return(sb_external_evaluation_failure(
      prediction,
      "endpoint_contract",
      endpoint_attempt$error
    ))
  }
  endpoint_contract <- endpoint_attempt$value
  endpoint_n <- length(endpoint_contract$patient_ids)
  if (!identical(
    prediction_contract$task_type,
    endpoint_contract$task_type
  )) {
    error <- tryCatch(
      sb_abort(
        "Prediction and endpoint task types do not match.",
        "stormbreaker_external_task_mismatch"
      ),
      error = function(error) error
    )
    return(sb_external_evaluation_failure(
      prediction,
      "endpoint_contract",
      error,
      endpoint_n
    ))
  }

  horizon_attempt <- sb_external_attempt(
    sb_external_horizon_semantics(
      prediction_contract,
      endpoint_contract
    )
  )
  if (!is.null(horizon_attempt$error)) {
    return(sb_external_evaluation_failure(
      prediction,
      "horizon_semantics",
      horizon_attempt$error,
      endpoint_n
    ))
  }

  alignment_attempt <- sb_external_attempt(
    sb_align_external_patients(
      prediction_contract,
      endpoint_contract
    )
  )
  if (!is.null(alignment_attempt$error)) {
    patient_match <- list(
      status = "failed",
      prediction_n = as.integer(length(prediction_contract$patient_ids)),
      endpoint_n = as.integer(endpoint_n),
      exact_set = FALSE,
      output_order = "external_endpoint"
    )
    return(sb_external_evaluation_failure(
      prediction,
      "patient_match",
      alignment_attempt$error,
      endpoint_n,
      patient_match
    ))
  }
  aligned <- alignment_attempt$value
  scored <- aligned$scored

  metric_attempt <- sb_external_attempt({
    if (identical(prediction_contract$task_type, "binary")) {
      metrics <- list(
        roc_auc = sb_external_binary_auc_metric(scored),
        average_precision = sb_external_binary_average_precision_metric(
          scored
        ),
        brier_score = sb_external_binary_brier_metric(scored),
        log_loss = sb_external_binary_log_loss_metric(scored),
        mean_calibration = sb_external_binary_mean_calibration_metric(scored),
        calibration = sb_external_binary_calibration_metric(scored)
      )
      if (!is.null(prediction_contract$decision_thresholds)) {
        metrics$decision_curve <- sb_binary_decision_curve_metric(
          scored,
          prediction_contract$decision_thresholds,
          scope = "patient_level_external_cohort"
        )
      }
      metrics
    } else {
      metrics <- list(
        harrell_c = sb_external_harrell_c_metric(scored)
      )
      if (prediction_contract$has_horizon) {
        ipcw_brier <- sb_external_ipcw_brier_metric(
          scored,
          endpoint_contract
        )
        metrics$ipcw_brier <- ipcw_brier
        metrics$time_auc <- sb_external_time_auc_metric(
          scored,
          endpoint_contract,
          ipcw_brier$censoring_curve
        )
        metrics$uno_c <- sb_external_uno_c_metric(
          scored,
          endpoint_contract,
          ipcw_brier$censoring_curve
        )
        if (!is.null(prediction_contract$decision_thresholds)) {
          metrics$decision_curve <- sb_survival_decision_curve_metric(
            scored,
            prediction_contract$decision_thresholds,
            prediction_contract$prediction_horizon,
            prediction_contract$time_unit,
            scope = "patient_level_external_cohort"
          )
        }
      }
      metrics
    }
  })
  if (!is.null(metric_attempt$error)) {
    return(sb_external_evaluation_failure(
      prediction,
      "metric_evaluation",
      metric_attempt$error,
      endpoint_n,
      aligned$patient_match
    ))
  }
  metrics <- metric_attempt$value
  if (identical(prediction_contract$task_type, "binary")) {
    metrics$moderate_calibration <- tryCatch(
      sb_binary_moderate_calibration_metric(
        scored,
        scope = "patient_level_external_cohort"
      ),
      error = function(error) {
        sb_binary_moderate_calibration_error_metric(
          error,
          scope = "patient_level_external_cohort",
          n = nrow(scored)
        )
      }
    )
  }
  if (identical(prediction_contract$task_type, "survival") &&
      prediction_contract$has_horizon) {
    metrics$calibration <- tryCatch(
      sb_survival_calibration_metric(
        patient_id = scored$patient_id,
        time = scored$time,
        event = scored$event,
        survival_probability = scored$survival_probability,
        log_cumulative_hazard = scored$log_cumulative_hazard,
        prediction_status = scored$status,
        prediction_horizon = prediction_contract$prediction_horizon,
        time_unit = endpoint_contract$time_unit,
        scope = "patient_level_external_cohort_at_locked_horizon"
      ),
      error = function(error) {
        sb_survival_calibration_error_metric(
          error = error,
          prediction_horizon = prediction_contract$prediction_horizon,
          time_unit = endpoint_contract$time_unit,
          scope = "patient_level_external_cohort_at_locked_horizon",
          n = nrow(scored)
        )
      }
    )
    metrics$weak_calibration <- tryCatch(
      sb_survival_weak_calibration_metric(
        patient_id = scored$patient_id,
        time = scored$time,
        event = scored$event,
        survival_probability = scored$survival_probability,
        log_cumulative_hazard = scored$log_cumulative_hazard,
        prediction_status = scored$status,
        prediction_horizon = prediction_contract$prediction_horizon,
        time_unit = endpoint_contract$time_unit,
        scope = "patient_level_external_cohort_at_locked_horizon"
      ),
      error = function(error) {
        sb_survival_weak_calibration_error_metric(
          error = error,
          prediction_horizon = prediction_contract$prediction_horizon,
          time_unit = endpoint_contract$time_unit,
          scope = "patient_level_external_cohort_at_locked_horizon",
          n = nrow(scored)
        )
      }
    )
  }
  required_complete <- if (identical(
    prediction_contract$task_type,
    "binary"
  )) {
    all(vapply(
      metrics[c("roc_auc", "average_precision", "brier_score", "log_loss")],
      function(metric) identical(metric$status, "success"),
      logical(1)
    ))
  } else {
    complete <- identical(metrics$harrell_c$status, "success")
    if (prediction_contract$has_horizon) {
      complete <- complete && identical(
        metrics$ipcw_brier$status,
        "success"
      ) && identical(metrics$time_auc$status, "success") &&
        identical(metrics$uno_c$status, "success")
      if (!is.null(prediction_contract$decision_thresholds)) {
        complete <- complete && identical(
          metrics$decision_curve$status,
          "success"
        )
      }
    }
    complete
  }

  result <- list(
    status = if (required_complete) "success" else "incomplete",
    task_type = prediction_contract$task_type,
    pipeline_id = prediction_contract$pipeline_id,
    learner = prediction_contract$learner,
    n = as.integer(nrow(scored)),
    evaluation_boundary = sb_external_evaluation_boundary(),
    patient_match = aligned$patient_match,
    scored_predictions = scored,
    metrics = metrics,
    endpoint = endpoint_contract$summary,
    provenance = sb_external_evaluation_provenance(),
    error_stage = NA_character_,
    error_class = NA_character_,
    error_message = NA_character_
  )
  class(result) <- "stormbreaker_external_evaluation"
  result
}

#' @export
print.stormbreaker_external_evaluation <- function(x, ...) {
  cat("<stormbreaker_external_evaluation:", x$status, ">\n")
  if (!is.null(x$direction_correction)) {
    cat("Post-hoc direction correction using evaluation labels; reversed:",
        x$direction_correction$reversed, "\n")
  }
  cat("precomputed-prediction external evaluation\n")
  cat("Pipeline:", x$pipeline_id, "| Task:", x$task_type, "\n")
  cat("Patients:", x$n, "\n")
  if (identical(x$status, "failed")) {
    cat("Failed stage:", x$error_stage, "|", x$error_class, "\n")
    return(invisible(x))
  }
  if (identical(x$task_type, "binary")) {
    cat("External ROC AUC:", format(x$metrics$roc_auc$estimate), "\n")
    cat("External Brier score:", format(x$metrics$brier_score$estimate), "\n")
    mean_calibration <- x$metrics$mean_calibration
    if (identical(mean_calibration$status, "success")) {
      cat(
        "External mean calibration O:E:",
        format(mean_calibration$mean_calibration_ratio), "\n"
      )
    } else {
      cat("External mean calibration O:E: unavailable\n")
    }
    moderate_calibration <- x$metrics$moderate_calibration
    if (identical(moderate_calibration$status, "success")) {
      cat(
        "External binary moderate calibration ICI:",
        format(moderate_calibration$ici),
        "| E50:", format(moderate_calibration$e50),
        "| E90:", format(moderate_calibration$e90), "\n"
      )
    } else {
      cat("External binary moderate calibration: unavailable\n")
    }
    if ("decision_curve" %in% names(x$metrics)) {
      decision_curve <- x$metrics$decision_curve
      if (identical(decision_curve$status, "success")) {
        cat(
          "External binary decision curve:",
          decision_curve$n_thresholds,
          "declared thresholds\n"
        )
      } else {
        cat("External binary decision curve: unavailable\n")
      }
    }
  } else {
    cat("External Harrell C:", format(x$metrics$harrell_c$estimate), "\n")
    if ("ipcw_brier" %in% names(x$metrics)) {
      cat(
        "External IPCW Brier:",
        format(x$metrics$ipcw_brier$estimate), "\n"
      )
      cat(
        "External time-dependent AUC:",
        format(x$metrics$time_auc$estimate), "\n"
      )
      cat(
        "External Uno C:",
        format(x$metrics$uno_c$estimate), "\n"
      )
      calibration <- x$metrics$calibration
      if (identical(calibration$status, "success")) {
        cat(
          "External survival calibration ICI:",
          format(calibration$ici),
          "| E50:", format(calibration$e50),
          "| E90:", format(calibration$e90), "\n"
        )
      } else {
        cat("External survival calibration: unavailable\n")
      }
      weak_calibration <- x$metrics$weak_calibration
      if (identical(weak_calibration$mean_calibration$status, "success")) {
        cat(
          "External mean calibration O:E:",
          format(weak_calibration$mean_calibration_ratio), "\n"
        )
      } else {
        cat("External mean calibration O:E: unavailable\n")
      }
      if (identical(weak_calibration$slope_fit$status, "success")) {
        cat(
          "External survival calibration slope:",
          format(weak_calibration$calibration_slope), "\n"
        )
      } else {
        cat("External survival calibration slope: unavailable\n")
      }
      if ("decision_curve" %in% names(x$metrics)) {
        decision_curve <- x$metrics$decision_curve
        if (identical(decision_curve$status, "success")) {
          label <- if (decision_curve$n_thresholds == 1L) {
            "declared threshold"
          } else {
            "declared thresholds"
          }
          cat(
            "External survival decision curve:",
            decision_curve$n_thresholds,
            label, "\n"
          )
        } else {
          cat("External survival decision curve: unavailable\n")
        }
      }
    }
  }
  invisible(x)
}
