sb_external_bootstrap_empty_point_estimates <- function() {
  data.frame(
    metric = character(),
    estimate = numeric(),
    status = character(),
    error_class = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
}

sb_external_bootstrap_empty_intervals <- function() {
  data.frame(
    metric = character(),
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

sb_external_bootstrap_empty_replicates <- function() {
  data.frame(
    replicate = integer(),
    metric = character(),
    estimate = numeric(),
    status = character(),
    error_class = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
}

sb_external_bootstrap_empty_audit <- function() {
  data.frame(replicate = integer(), stringsAsFactors = FALSE)
}

sb_external_bootstrap_safe_integer <- function(value) {
  if (is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value == trunc(value) &&
      value >= -.Machine$integer.max && value <= .Machine$integer.max) {
    return(as.integer(value))
  }
  NA_integer_
}

sb_external_bootstrap_safe_number <- function(value) {
  if (is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value)) {
    return(as.numeric(value))
  }
  NA_real_
}

sb_external_bootstrap_safe_metadata <- function(evaluation) {
  result <- list(
    task_type = NA_character_,
    pipeline_id = NA_character_,
    learner = NA_character_,
    n = 0L
  )
  if (!is.list(evaluation)) {
    return(result)
  }
  for (name in c("task_type", "pipeline_id", "learner")) {
    value <- evaluation[[name]]
    if (sb_external_scalar_character(value)) {
      result[[name]] <- value
    }
  }
  scored <- evaluation[["scored_predictions"]]
  if (is.data.frame(scored)) {
    result$n <- as.integer(nrow(scored))
  }
  result
}

sb_external_bootstrap_unknown_target <- function(task_type) {
  if (identical(task_type, "binary")) {
    return(list(
      target_status = "complete",
      has_horizon = FALSE,
      prediction_horizon = NULL,
      time_unit = NULL,
      horizon_source = "not_applicable"
    ))
  }
  list(
    target_status = "unknown",
    has_horizon = NA,
    prediction_horizon = NULL,
    time_unit = NULL,
    horizon_source = "unknown"
  )
}

sb_external_bootstrap_provenance <- function() {
  kinds <- RNGkind()
  list(
    package = "stormbreaker",
    package_version = as.character(utils::packageVersion("stormbreaker")),
    r_version = R.version.string,
    survival_version = as.character(utils::packageVersion("survival")),
    rng_kind = unname(kinds[[1L]]),
    normal_kind = unname(kinds[[2L]]),
    sample_kind = unname(kinds[[3L]])
  )
}

sb_external_bootstrap_resampling <- function(task_type, has_horizon = NULL) {
  censoring <- if (!identical(task_type, "survival")) {
    "not_applicable"
  } else if (is.null(has_horizon)) {
    "unknown_not_validated"
  } else if (!isTRUE(has_horizon)) {
    "not_applicable_risk_only"
  } else {
    "reestimated_external_reverse_kaplan_meier_each_replicate"
  }
  list(
    unit = "patient",
    scheme = "ordinary_nonparametric_with_replacement",
    paired_endpoint_prediction = TRUE,
    prediction = "fixed",
    model_access = "none",
    censoring = censoring
  )
}

sb_external_bootstrap_interval_method <- function(confidence_level) {
  alpha <- 1 - confidence_level
  list(
    name = "percentile",
    quantile_type = 7L,
    probabilities = c(
      lower = alpha / 2,
      upper = 1 - alpha / 2
    )
  )
}

sb_external_bootstrap_failure <- function(
    evaluation,
    reps,
    confidence_level,
    seed,
    stage,
    error,
    has_horizon = NULL,
    prediction_target = NULL) {
  metadata <- sb_external_bootstrap_safe_metadata(evaluation)
  safe_confidence <- sb_external_bootstrap_safe_number(confidence_level)
  if (is.null(prediction_target)) {
    prediction_target <- sb_external_bootstrap_unknown_target(
      metadata$task_type
    )
  }
  result <- list(
    status = "failed",
    task_type = metadata$task_type,
    pipeline_id = metadata$pipeline_id,
    learner = metadata$learner,
    n = metadata$n,
    prediction_target = prediction_target,
    reps = sb_external_bootstrap_safe_integer(reps),
    confidence_level = safe_confidence,
    seed = sb_external_bootstrap_safe_integer(seed),
    interval_method = if (is.finite(safe_confidence) &&
        safe_confidence > 0 && safe_confidence < 1) {
      sb_external_bootstrap_interval_method(safe_confidence)
    } else {
      list(
        name = "percentile",
        quantile_type = 7L,
        probabilities = c(lower = NA_real_, upper = NA_real_)
      )
    },
    resampling = sb_external_bootstrap_resampling(
      metadata$task_type,
      has_horizon
    ),
    point_estimates = sb_external_bootstrap_empty_point_estimates(),
    intervals = sb_external_bootstrap_empty_intervals(),
    replicates = sb_external_bootstrap_empty_replicates(),
    resample_audit = sb_external_bootstrap_empty_audit(),
    provenance = sb_external_bootstrap_provenance(),
    error_stage = stage,
    error_class = class(error)[[1L]],
    error_message = conditionMessage(error)
  )
  class(result) <- "stormbreaker_external_bootstrap"
  result
}

sb_external_bootstrap_column <- function(data, name) {
  positions <- which(names(data) == name)
  if (length(positions) != 1L) {
    sb_abort(
      paste0(
        "External evaluation column `", name, "` is missing or duplicated."
      ),
      "stormbreaker_bad_external_evaluation"
    )
  }
  data[[positions]]
}

sb_external_bootstrap_contract <- function(evaluation) {
  if (!inherits(evaluation, "stormbreaker_external_evaluation") ||
      !is.list(evaluation)) {
    sb_abort(
      "`evaluation` must be returned by `sb_evaluate_external()`.",
      "stormbreaker_bad_external_evaluation"
    )
  }
  if (identical(evaluation[["status"]], "failed")) {
    sb_abort(
      "A failed external evaluation cannot be bootstrapped.",
      "stormbreaker_failed_external_evaluation"
    )
  }
  task_type <- evaluation[["task_type"]]
  pipeline_id <- evaluation[["pipeline_id"]]
  learner <- evaluation[["learner"]]
  scored <- evaluation[["scored_predictions"]]
  valid_metadata <- evaluation[["status"]] %in% c("success", "incomplete") &&
    length(evaluation[["status"]]) == 1L &&
    sb_external_scalar_character(task_type) &&
    task_type %in% c("binary", "survival") &&
    sb_external_scalar_character(pipeline_id) &&
    sb_external_scalar_character(learner) &&
    is.data.frame(scored) && nrow(scored) >= 1L &&
    is.numeric(evaluation[["n"]]) && length(evaluation[["n"]]) == 1L &&
    !is.na(evaluation[["n"]]) && is.finite(evaluation[["n"]]) &&
    evaluation[["n"]] == nrow(scored)
  if (!isTRUE(valid_metadata)) {
    sb_abort(
      "The external evaluation metadata are incomplete or inconsistent.",
      "stormbreaker_bad_external_evaluation"
    )
  }

  patient_id <- as.character(
    sb_external_bootstrap_column(scored, "patient_id")
  )
  row_status <- sb_external_bootstrap_column(scored, "status")
  valid_rows <- length(patient_id) == nrow(scored) &&
    !anyNA(patient_id) && all(nzchar(trimws(patient_id))) &&
    !anyDuplicated(patient_id) && is.character(row_status) &&
    length(row_status) == nrow(scored) && !anyNA(row_status) &&
    all(row_status == "success")
  if (!valid_rows) {
    sb_abort(
      "External evaluation rows have invalid patient IDs or status.",
      "stormbreaker_bad_external_evaluation"
    )
  }

  if (identical(task_type, "binary")) {
    outcome <- sb_external_bootstrap_column(scored, "outcome")
    probability <- sb_external_bootstrap_column(scored, "probability")
    decision_thresholds <-
      sb_external_bootstrap_decision_thresholds(evaluation)
    valid_outcome <- (is.numeric(outcome) || is.logical(outcome)) &&
      length(outcome) == nrow(scored) && !anyNA(outcome) &&
      (!is.numeric(outcome) || all(is.finite(outcome))) &&
      all(outcome %in% c(0, 1))
    valid_probability <- is.numeric(probability) &&
      length(probability) == nrow(scored) && !anyNA(probability) &&
      all(is.finite(probability)) &&
      all(probability >= 0 & probability <= 1)
    if (!valid_outcome || !valid_probability) {
      sb_abort(
        "External binary endpoint or probability rows are invalid.",
        "stormbreaker_bad_external_evaluation"
      )
    }
    canonical <- data.frame(
      patient_id = patient_id,
      outcome = as.integer(outcome),
      probability = unname(probability),
      status = rep("success", nrow(scored)),
      stringsAsFactors = FALSE
    )
    return(list(
      task_type = task_type,
      pipeline_id = pipeline_id,
      learner = learner,
      n = as.integer(nrow(canonical)),
      prediction_target = list(
        target_status = "complete",
        has_horizon = FALSE,
        prediction_horizon = NULL,
        time_unit = NULL,
        horizon_source = "not_applicable"
      ),
      has_horizon = FALSE,
      prediction_horizon = NULL,
      time_unit = NULL,
      decision_thresholds = decision_thresholds,
      scored = canonical
    ))
  }

  endpoint_summary <- evaluation[["endpoint"]]
  endpoint_time_unit <- if (is.list(endpoint_summary)) {
    endpoint_summary[["time_unit"]]
  } else {
    NULL
  }
  endpoint_horizon <- if (is.list(endpoint_summary)) {
    endpoint_summary[["prediction_horizon"]]
  } else {
    NULL
  }
  valid_endpoint_unit <- sb_external_scalar_character(endpoint_time_unit) &&
    endpoint_time_unit %in% c("day", "month", "year")
  valid_endpoint_horizon <- is.null(endpoint_horizon) || (
    is.numeric(endpoint_horizon) && length(endpoint_horizon) == 1L &&
      !is.na(endpoint_horizon) && is.finite(endpoint_horizon) &&
      endpoint_horizon > 0
  )
  if (!valid_endpoint_unit || !valid_endpoint_horizon) {
    sb_abort(
      "External survival endpoint target metadata are incomplete or invalid.",
      "stormbreaker_bad_external_evaluation"
    )
  }

  time <- sb_external_bootstrap_column(scored, "time")
  event <- sb_external_bootstrap_column(scored, "event")
  risk <- sb_external_bootstrap_column(scored, "risk")
  valid_time <- is.numeric(time) && length(time) == nrow(scored) &&
    !anyNA(time) && all(is.finite(time)) && all(time >= 0)
  valid_event <- (is.numeric(event) || is.logical(event)) &&
    length(event) == nrow(scored) && !anyNA(event) &&
    (!is.numeric(event) || all(is.finite(event))) &&
    all(event %in% c(0, 1))
  valid_risk <- is.numeric(risk) && length(risk) == nrow(scored) &&
    !anyNA(risk) && all(is.finite(risk))
  if (!valid_time || !valid_event || !valid_risk) {
    sb_abort(
      "External survival endpoint or risk rows are invalid.",
      "stormbreaker_bad_external_evaluation"
    )
  }

  horizon_names <- c(
    "prediction_horizon", "time_unit", "survival_probability"
  )
  horizon_counts <- vapply(
    horizon_names,
    function(name) sum(names(scored) == name),
    integer(1)
  )
  if (any(horizon_counts > 1L) ||
      (any(horizon_counts == 1L) && !all(horizon_counts == 1L))) {
    sb_abort(
      "External survival horizon fields are incomplete or duplicated.",
      "stormbreaker_bad_external_evaluation"
    )
  }
  canonical <- data.frame(
    patient_id = patient_id,
    time = unname(time),
    event = as.integer(event),
    risk = unname(risk),
    stringsAsFactors = FALSE
  )
  has_horizon <- all(horizon_counts == 1L)
  prediction_horizon <- NULL
  time_unit <- NULL
  if (has_horizon) {
    horizon <- scored[["prediction_horizon"]]
    unit <- scored[["time_unit"]]
    probability <- scored[["survival_probability"]]
    valid_horizon <- is.numeric(horizon) &&
      length(horizon) == nrow(scored) && !anyNA(horizon) &&
      all(is.finite(horizon)) && all(horizon > 0) &&
      length(unique(as.numeric(horizon))) == 1L
    valid_unit <- is.character(unit) && length(unit) == nrow(scored) &&
      !anyNA(unit) && length(unique(unit)) == 1L &&
      unique(unit) %in% c("day", "month", "year")
    valid_probability <- is.numeric(probability) &&
      length(probability) == nrow(scored) && !anyNA(probability) &&
      all(is.finite(probability)) &&
      all(probability >= 0 & probability <= 1)
    if (!valid_horizon || !valid_unit || !valid_probability) {
      sb_abort(
        "External survival horizon metadata or probabilities are invalid.",
        "stormbreaker_bad_external_evaluation"
      )
    }
    if (!identical(unit[[1L]], endpoint_time_unit) ||
        (!is.null(endpoint_horizon) && !identical(
          as.numeric(horizon[[1L]]),
          as.numeric(endpoint_horizon)
        ))) {
      sb_abort(
        paste0(
          "External survival scored predictions and endpoint do not retain ",
          "the same time target."
        ),
        "stormbreaker_bad_external_evaluation"
      )
    }
    prediction_horizon <- as.numeric(horizon[[1L]])
    time_unit <- unit[[1L]]
    canonical$prediction_horizon <- rep(
      prediction_horizon,
      nrow(canonical)
    )
    canonical$time_unit <- rep(time_unit, nrow(canonical))
    canonical$survival_probability <- unname(probability)
    stable <- sb_external_stable_scale(scored, probability)
    if (!is.null(stable)) canonical$log_cumulative_hazard <- stable
  } else if (!is.null(endpoint_horizon)) {
    sb_abort(
      paste0(
        "A risk-only external evaluation cannot retain an endpoint prediction ",
        "horizon."
      ),
      "stormbreaker_bad_external_evaluation"
    )
  }
  canonical$status <- rep("success", nrow(canonical))
  prediction_target <- list(
    target_status = "complete",
    has_horizon = has_horizon,
    prediction_horizon = prediction_horizon,
    time_unit = endpoint_time_unit,
    horizon_source = if (!has_horizon) {
      "none"
    } else if (!is.null(endpoint_horizon)) {
      "endpoint"
    } else {
      "prediction"
    }
  )
  list(
    task_type = task_type,
    pipeline_id = pipeline_id,
    learner = learner,
    n = as.integer(nrow(canonical)),
    prediction_target = prediction_target,
    has_horizon = has_horizon,
    prediction_horizon = prediction_horizon,
    time_unit = time_unit,
    decision_thresholds = if (has_horizon) {
      sb_external_bootstrap_decision_thresholds(evaluation)
    } else {
      NULL
    },
    scored = canonical
  )
}

sb_external_bootstrap_metric_row <- function(
    metric,
    estimate,
    status,
    error_class = NA_character_,
    error_message = NA_character_) {
  data.frame(
    metric = metric,
    estimate = as.numeric(estimate),
    status = as.character(status),
    error_class = as.character(error_class),
    error_message = as.character(error_message),
    stringsAsFactors = FALSE
  )
}

sb_external_bootstrap_metric_error <- function(metric, error) {
  sb_external_bootstrap_metric_row(
    metric,
    NA_real_,
    "metric_failed",
    class(error)[[1L]],
    conditionMessage(error)
  )
}

sb_external_bootstrap_metric_result <- function(metric_name, metric) {
  error_class <- metric[["error_class"]]
  error_message <- metric[["error_message"]]
  if (!sb_external_scalar_character(error_class)) {
    error_class <- NA_character_
  }
  if (!sb_external_scalar_character(error_message)) {
    error_message <- NA_character_
  }
  sb_external_bootstrap_metric_row(
    metric_name,
    metric[["estimate"]],
    metric[["status"]],
    error_class,
    error_message
  )
}

sb_external_bootstrap_metric_value <- function(metric_name, code) {
  tryCatch({
    sb_external_bootstrap_metric_result(metric_name, force(code))
  }, error = function(error) {
    sb_external_bootstrap_metric_error(metric_name, error)
  })
}

sb_external_bootstrap_calibration_rows <- function(scored) {
  calibration <- tryCatch(
    sb_external_binary_calibration_metric(scored),
    error = function(error) error
  )
  metric_names <- c("calibration_intercept", "calibration_slope")
  if (inherits(calibration, "condition")) {
    rows <- lapply(
      metric_names,
      sb_external_bootstrap_metric_error,
      error = calibration
    )
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    return(result)
  }
  if (!identical(calibration$input_status, "success")) {
    status <- calibration$input_status
    if (!sb_external_scalar_character(status)) {
      status <- "not_evaluable"
    }
    rows <- lapply(metric_names, function(metric_name) {
      sb_external_bootstrap_metric_row(
        metric_name,
        NA_real_,
        status,
        calibration$error_class,
        calibration$error_message
      )
    })
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    return(result)
  }

  fits <- list(calibration$intercept_fit, calibration$slope_fit)
  estimates <- c(
    calibration$calibration_intercept,
    calibration$calibration_slope
  )
  rows <- lapply(seq_along(metric_names), function(index) {
    fit <- fits[[index]]
    estimate <- if (identical(fit$status, "success")) {
      estimates[[index]]
    } else {
      NA_real_
    }
    sb_external_bootstrap_metric_row(
      metric_names[[index]],
      estimate,
      fit$status,
      fit$error_class,
      fit$error_message
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

sb_external_bootstrap_binary_moderate_calibration_rows <- function(scored) {
  metric_names <- c(
    "moderate_calibration_ici",
    "moderate_calibration_e50",
    "moderate_calibration_e90"
  )
  calibration <- tryCatch(
    sb_binary_moderate_calibration_metric(
      scored,
      scope = "patient_level_external_bootstrap_sample"
    ),
    error = function(error) error
  )
  if (inherits(calibration, "condition")) {
    result <- do.call(rbind, lapply(
      metric_names,
      sb_external_bootstrap_metric_error,
      error = calibration
    ))
    rownames(result) <- NULL
    return(result)
  }

  successful <- identical(calibration$status, "success")
  estimates <- if (successful) {
    c(calibration$ici, calibration$e50, calibration$e90)
  } else {
    rep(NA_real_, length(metric_names))
  }
  status <- if (sb_external_scalar_character(calibration$status)) {
    calibration$status
  } else {
    "metric_failed"
  }
  result <- do.call(rbind, lapply(seq_along(metric_names), function(index) {
    sb_external_bootstrap_metric_row(
      metric = metric_names[[index]],
      estimate = estimates[[index]],
      status = status,
      error_class = calibration$error_class,
      error_message = calibration$error_message
    )
  }))
  rownames(result) <- NULL
  result
}

sb_external_bootstrap_survival_calibration_rows <- function(scored) {
  metric_names <- c("calibration_ici", "calibration_e50", "calibration_e90")
  calibration <- tryCatch(
    sb_survival_calibration_metric(
      patient_id = scored$patient_id,
      time = scored$time,
      event = scored$event,
      survival_probability = scored$survival_probability,
      log_cumulative_hazard = scored$log_cumulative_hazard,
      prediction_status = scored$status,
      prediction_horizon = unique(scored$prediction_horizon)[[1L]],
      time_unit = unique(scored$time_unit)[[1L]],
      scope = "patient_level_external_cohort_at_locked_horizon"
    ),
    error = function(error) error
  )
  if (inherits(calibration, "condition")) {
    result <- do.call(rbind, lapply(
      metric_names,
      sb_external_bootstrap_metric_error,
      error = calibration
    ))
    rownames(result) <- NULL
    return(result)
  }

  successful <- identical(calibration$status, "success")
  estimates <- if (successful) {
    c(calibration$ici, calibration$e50, calibration$e90)
  } else {
    rep(NA_real_, length(metric_names))
  }
  status <- if (sb_external_scalar_character(calibration$status)) {
    calibration$status
  } else {
    "metric_failed"
  }
  result <- do.call(rbind, lapply(seq_along(metric_names), function(index) {
    sb_external_bootstrap_metric_row(
      metric = metric_names[[index]],
      estimate = estimates[[index]],
      status = status,
      error_class = calibration$error_class,
      error_message = calibration$error_message
    )
  }))
  rownames(result) <- NULL
  result
}

sb_external_bootstrap_survival_weak_calibration_rows <- function(scored) {
  metric_names <- c("mean_calibration_ratio", "calibration_slope")
  calibration <- tryCatch(
    sb_survival_weak_calibration_metric(
      patient_id = scored$patient_id,
      time = scored$time,
      event = scored$event,
      survival_probability = scored$survival_probability,
      log_cumulative_hazard = scored$log_cumulative_hazard,
      prediction_status = scored$status,
      prediction_horizon = unique(scored$prediction_horizon)[[1L]],
      time_unit = unique(scored$time_unit)[[1L]],
      scope = "patient_level_external_cohort_at_locked_horizon"
    ),
    error = function(error) error
  )
  if (inherits(calibration, "condition")) {
    result <- do.call(rbind, lapply(
      metric_names,
      sb_external_bootstrap_metric_error,
      error = calibration
    ))
    rownames(result) <- NULL
    return(result)
  }

  components <- list(calibration$mean_calibration, calibration$slope_fit)
  estimates <- c(
    calibration$mean_calibration_ratio,
    calibration$calibration_slope
  )
  result <- do.call(rbind, lapply(seq_along(metric_names), function(index) {
    component <- components[[index]]
    status <- if (sb_external_scalar_character(component$status)) {
      component$status
    } else {
      "metric_failed"
    }
    estimate <- if (identical(status, "success")) {
      estimates[[index]]
    } else {
      NA_real_
    }
    error_class <- if (sb_external_scalar_character(component$error_class)) {
      component$error_class
    } else {
      NA_character_
    }
    error_message <- if (sb_external_scalar_character(
        component$error_message
      )) {
      component$error_message
    } else {
      NA_character_
    }
    sb_external_bootstrap_metric_row(
      metric = metric_names[[index]],
      estimate = estimate,
      status = status,
      error_class = error_class,
      error_message = error_message
    )
  }))
  rownames(result) <- NULL
  result
}

sb_external_bootstrap_metric_rows <- function(scored, task_type) {
  if (identical(task_type, "binary")) {
    rows <- list(
      sb_external_bootstrap_metric_value(
        "roc_auc",
        sb_external_binary_auc_metric(scored)
      ),
      sb_external_bootstrap_metric_value(
        "average_precision",
        sb_external_binary_average_precision_metric(scored)
      ),
      sb_external_bootstrap_metric_value(
        "brier_score",
        sb_external_binary_brier_metric(scored)
      ),
      sb_external_bootstrap_metric_value(
        "log_loss",
        sb_external_binary_log_loss_metric(scored)
      ),
      sb_external_bootstrap_metric_value(
        "mean_calibration_ratio",
        sb_external_binary_mean_calibration_metric(scored)
      ),
      sb_external_bootstrap_calibration_rows(scored),
      sb_external_bootstrap_binary_moderate_calibration_rows(scored)
    )
  } else {
    rows <- list(sb_external_bootstrap_metric_value(
      "harrell_c",
      sb_external_harrell_c_metric(scored)
    ))
    if (all(c(
      "prediction_horizon", "time_unit", "survival_probability"
    ) %in% names(scored))) {
      endpoint <- list(time_unit = unique(scored$time_unit)[[1L]])
      brier <- tryCatch(
        sb_external_ipcw_brier_metric(scored, endpoint),
        error = function(error) error
      )
      if (inherits(brier, "condition")) {
        rows[[length(rows) + 1L]] <- sb_external_bootstrap_metric_error(
          "ipcw_brier_score",
          brier
        )
        censoring_curve <- sb_empty_censoring_curve()
      } else {
        rows[[length(rows) + 1L]] <- sb_external_bootstrap_metric_result(
          "ipcw_brier_score",
          brier
        )
        censoring_curve <- brier$censoring_curve
      }
      rows[[length(rows) + 1L]] <- sb_external_bootstrap_metric_value(
        "time_auc",
        sb_external_time_auc_metric(scored, endpoint, censoring_curve)
      )
      rows[[length(rows) + 1L]] <- sb_external_bootstrap_metric_value(
        "uno_c",
        sb_external_uno_c_metric(scored, endpoint, censoring_curve)
      )
      rows[[length(rows) + 1L]] <-
        sb_external_bootstrap_survival_calibration_rows(scored)
      rows[[length(rows) + 1L]] <-
        sb_external_bootstrap_survival_weak_calibration_rows(scored)
    }
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

sb_external_bootstrap_audit_row <- function(scored, task_type, replicate) {
  common <- data.frame(
    replicate = as.integer(replicate),
    n_draws = as.integer(nrow(scored)),
    n_unique_patients = as.integer(length(unique(scored$patient_id))),
    stringsAsFactors = FALSE
  )
  if (identical(task_type, "binary")) {
    common$n_positive <- as.integer(sum(scored$outcome == 1L))
    common$n_negative <- as.integer(sum(scored$outcome == 0L))
  } else {
    common$n_events <- as.integer(sum(scored$event == 1L))
    common$n_censored <- as.integer(sum(scored$event == 0L))
  }
  common
}

sb_external_bootstrap_replicate <- function(
    scored,
    indices,
    task_type,
    replicate,
    decision_thresholds = NULL) {
  sampled <- scored[indices, , drop = FALSE]
  rownames(sampled) <- NULL
  rows <- sb_external_bootstrap_metric_rows(sampled, task_type)
  rows <- data.frame(
    replicate = rep.int(as.integer(replicate), nrow(rows)),
    rows,
    stringsAsFactors = FALSE
  )
  decision_curve <- NULL
  if (!is.null(decision_thresholds)) {
    decision_curve <- if (identical(task_type, "binary")) {
      sb_external_bootstrap_decision_curve_replicate(
        sampled,
        decision_thresholds,
        replicate
      )
    } else {
      sb_external_bootstrap_survival_decision_curve_replicate(
        sampled,
        decision_thresholds,
        replicate
      )
    }
  }
  list(
    metrics = rows,
    audit = sb_external_bootstrap_audit_row(
      sampled,
      task_type,
      replicate
    ),
    decision_curve = decision_curve
  )
}

sb_external_bootstrap_resamples <- function(
    scored,
    reps,
    seed,
    task_type,
    decision_thresholds = NULL) {
  pieces <- withr::with_seed(seed, {
    lapply(seq_len(reps), function(replicate) {
      indices <- sample.int(
        nrow(scored),
        nrow(scored),
        replace = TRUE
      )
      sb_external_bootstrap_replicate(
        scored,
        indices,
        task_type,
        replicate,
        decision_thresholds
      )
    })
  })
  replicates <- do.call(rbind, lapply(pieces, `[[`, "metrics"))
  audit <- do.call(rbind, lapply(pieces, `[[`, "audit"))
  decision_curve <- NULL
  if (!is.null(decision_thresholds)) {
    decision_curve <- do.call(
      rbind,
      lapply(pieces, `[[`, "decision_curve")
    )
    rownames(decision_curve) <- NULL
  }
  rownames(replicates) <- NULL
  rownames(audit) <- NULL
  list(
    replicates = replicates,
    audit = audit,
    decision_curve = decision_curve
  )
}

sb_external_percentile_intervals <- function(
    point,
    replicates,
    confidence_level,
    reps) {
  probabilities <- sb_external_bootstrap_interval_method(
    confidence_level
  )$probabilities
  rows <- lapply(seq_len(nrow(point)), function(index) {
    point_row <- point[index, , drop = FALSE]
    metric_rows <- replicates[
      replicates$metric == point_row$metric,
      ,
      drop = FALSE
    ]
    successful <- !is.na(metric_rows$status) &
      metric_rows$status == "success"
    finite_successful <- successful &
      !is.na(metric_rows$estimate) & is.finite(metric_rows$estimate)
    nonfinite <- successful & !finite_successful
    unavailable <- !successful
    n_finite <- as.integer(sum(finite_successful))
    n_nonfinite <- as.integer(sum(nonfinite))
    n_unavailable <- as.integer(sum(unavailable))
    completeness <- if (n_finite == reps) {
      "complete"
    } else if (n_finite == 0L) {
      "none"
    } else {
      "partial"
    }
    point_finite <- identical(point_row$status[[1L]], "success") &&
      !is.na(point_row$estimate[[1L]]) &&
      is.finite(point_row$estimate[[1L]])
    interval_available <- point_finite && n_finite >= 2L
    bounds <- c(lower = NA_real_, upper = NA_real_)
    if (interval_available) {
      bounds <- stats::quantile(
        metric_rows$estimate[finite_successful],
        probs = probabilities,
        names = FALSE,
        type = 7
      )
    }
    row_status <- if (interval_available && completeness == "complete") {
      "success"
    } else {
      "incomplete"
    }
    error_class <- NA_character_
    error_message <- NA_character_
    if (!point_finite) {
      error_class <- point_row$error_class[[1L]]
      error_message <- point_row$error_message[[1L]]
      if (!sb_external_scalar_character(error_class)) {
        error_class <- "stormbreaker_nonfinite_or_unavailable_point_estimate"
      }
      if (!sb_external_scalar_character(error_message)) {
        error_message <- paste0(
          "A finite successful point estimate is required for `",
          point_row$metric[[1L]], "`."
        )
      }
    } else if (n_finite < 2L) {
      error_class <- "stormbreaker_insufficient_bootstrap_replicates"
      error_message <- paste0(
        "At least two finite successful bootstrap replicates are required for `",
        point_row$metric[[1L]], "`."
      )
    } else if (completeness != "complete") {
      error_class <- "stormbreaker_incomplete_bootstrap_replicates"
      error_message <- paste0(
        "The displayed interval for `", point_row$metric[[1L]],
        "` uses only its finite successful replicates."
      )
    }
    data.frame(
      metric = point_row$metric[[1L]],
      point_estimate = point_row$estimate[[1L]],
      lower = unname(bounds[[1L]]),
      upper = unname(bounds[[2L]]),
      confidence_level = confidence_level,
      status = row_status,
      interval_available = interval_available,
      replicate_completeness = completeness,
      n_requested = as.integer(reps),
      n_finite_successful = n_finite,
      n_unavailable = n_unavailable,
      n_nonfinite = n_nonfinite,
      valid_fraction = unname(n_finite / reps),
      error_class = error_class,
      error_message = error_message,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Bootstrap uncertainty for a completed external evaluation
#'
#' `sb_bootstrap_external()` applies an ordinary nonparametric patient
#' bootstrap to the paired endpoint and precomputed prediction rows already
#' stored by [sb_evaluate_external()]. The fitted model and every supplied
#' prediction remain fixed. The result therefore describes external-patient
#' sampling uncertainty conditional on that locked model and those predictions;
#' it does not include uncertainty from development fitting, feature selection,
#' tuning, or model selection.
#' For survival, `prediction_target` separately retains the endpoint time unit,
#' whether a prediction horizon is present, its value when present, and whether
#' that effective horizon came from the endpoint or the locked prediction. A
#' risk-only evaluation is recorded as horizon-free rather than assigned a
#' surrogate time point.
#'
#' Binary replicates recompute ROC AUC, non-interpolated average precision,
#' Brier score, unclipped log loss, observed-to-expected mean calibration, and
#' separate calibration intercept and slope targets. Mean calibration may be
#' zero when a draw contains no observed events and is unavailable when its
#' expected event count is zero or its ratio is numerically non-finite; no
#' clipping or continuity correction is applied. When the source binary
#' evaluation contains researcher-declared
#' decision thresholds, the same sampled patient table also recomputes model,
#' treat-all, and treat-none net benefit. An optional nested decision-curve
#' result retains the raw counts, threshold operating characteristics, and
#' curves, then reports type-7 percentile
#' intervals for the paired differences of model minus treat-all and model
#' minus treat-none. Model minus treat-all is calculated directly from true
#' negatives and false negatives rather than by subtracting two potentially
#' large, nearly equal net-benefit values. These intervals are pointwise at the
#' declared thresholds, not simultaneous coverage for the curve as a whole.
#' They provide no p value, threshold recommendation, automatic decision, or
#' proof of clinical effectiveness. They are not confidence intervals for the
#' threshold operating characteristics. Scalar target tables remain separate
#' from the nested decision-curve result. Survival replicates recompute pooled
#' Harrell C. For
#' a locked-horizon
#' survival prediction, every replicate also refits only its external marginal
#' reverse Kaplan-Meier censoring curve before recomputing IPCW Brier score,
#' cumulative/dynamic time-dependent AUC, and truncated Uno IPCW concordance.
#' The three metrics share that one curve but retain separate estimability and
#' failure status. For time-dependent AUC,
#' observed events at or before the locked horizon are inverse-censoring-
#' weighted cases, follow-up strictly beyond it defines controls, censoring at
#' or before it is excluded, and exact marker ties receive half credit.
#' Uno C uses the fixed risk score, eligible events strictly before the horizon,
#' squared inverse censoring weights at those event times, and all comparable
#' later observations plus same-time censoring; exact risk ties receive half
#' credit. No replicate refits a model or prediction.
#' When the locked-horizon survival evaluation contains researcher-declared
#' decision thresholds, the same sampled patient table also re-estimates event
#' Kaplan-Meier survival within each threshold-selected group and in the complete
#' sample. The model, saved survival probabilities, horizon, and thresholds stay
#' fixed. Paired model-minus-treat-all and model-minus-treat-none differences
#' are formed within each draw, and type-7 percentile intervals use only finite
#' successful replicates while retaining explicit unavailable counts. These are
#' pointwise conditional external-patient sampling intervals, not simultaneous
#' bands, development-fitting or selection uncertainty, optimism correction,
#' p values, proof of treatment effectiveness, or a rule for choosing a model
#' or threshold. They assume independent patients and independent
#' non-informative censoring within the relevant risk groups. No recalibration,
#' development refit, or new prediction is performed. Every point calculation
#' and patient resample
#' for a binary task also fits one separate three-knot restricted natural cubic
#' spline logistic moderate-calibration curve to `qlogis()` of the fixed class-1
#' probabilities and reports ICI, E50, and E90 as separate targets. Knot
#' locations and the curve are re-estimated once per sample; repeated draws
#' retain their empirical row weight. Exact probabilities zero or one are not
#' clipped and leave these targets unavailable without erasing other binary
#' metrics. A complete or quasi-completely separated sample likewise leaves only
#' these targets unavailable. This curve is never applied as recalibration and
#' is distinct from the weak-calibration intercept and slope.
#' For a horizon survival task, every point calculation and patient resample
#' also fits one three-knot Cox moderate-calibration curve to the fixed survival
#' probabilities and sampled outcomes, then reports ICI, E50, and E90 as
#' separate targets. Its knots and curve are re-estimated once per sample;
#' repeated draws remain distinct row instances. This evaluation fit uses
#' complete right-censored follow-up, assumes proportional hazards and
#' independent non-informative censoring, and is read at the locked horizon. It
#' does not clip or alter predictions, is not applied as recalibration, and is
#' not itself a survival calibration intercept or slope. Separately, every
#' locked-horizon point calculation and patient draw reports mean calibration as
#' Kaplan-Meier observed event probability divided by mean predicted event
#' probability, and fits an Efron Cox calibration slope after administratively
#' censoring the endpoint at the horizon. Both target 1 and retain independent
#' finite/unavailable/non-finite replicate accounting. Exact boundary or
#' constant predictions can leave only the slope unavailable, while unsupported
#' Kaplan-Meier horizon follow-up can leave only the mean ratio unavailable. No
#' probability is clipped or recalibrated, and no survival calibration
#' intercept, p value, or universal cutoff is produced. These intervals retain
#' the fixed model and predictions and do not contain development fitting,
#' feature selection, tuning, OOF, or model-selection uncertainty.
#'
#' Scalar intervals are percentile intervals using finite successful replicate
#' estimates and R's type-7 empirical quantiles. All unavailable and non-finite
#' replicates are retained and counted. Bounds based on an incomplete finite
#' subset are reported as provisional with `incomplete` status; no minimum
#' success fraction is imposed.
#'
#' The resampling unit is one independent patient row. These intervals are not
#' valid by construction for clustered specimens, repeated measurements,
#' families, sites, or temporal dependence. They are descriptive for the
#' supplied sampling frame and do not prove external-cohort independence,
#' transportability, calibration, clinical utility, adequate sample size,
#' significance, or a decision threshold. ICI/E50/E90 have no universal
#' acceptable cutoff. Time-dependent AUC and Uno C are discrimination measures
#' under the marginal non-informative censoring assumption.
#'
#' @references
#' Vickers AJ, Cronin AM, Elkin EB, Gonen M (2008). Extensions to decision
#'   curve analysis. *BMC Medical Informatics and Decision Making*, 8, 53.
#'   \doi{10.1186/1472-6947-8-53}.
#'
#' Vickers AJ et al. (2023). Decision curve analysis: confidence intervals and
#'   hypothesis testing for net benefit. *Diagnostic and Prognostic Research*,
#'   7, 11. \doi{10.1186/s41512-023-00148-y}.
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
#' @param evaluation A successful or incomplete
#'   `stormbreaker_external_evaluation` returned by [sb_evaluate_external()].
#' @param reps Number of ordinary patient-level bootstrap replicates. Must be an
#'   integer of at least 2. Defaults to 1000.
#' @param confidence_level Percentile interval confidence level strictly between
#'   0 and 1. Defaults to 0.95.
#' @param seed Positive integer seed. Reproducibility requires the same input,
#'   R version, and RNG kinds. The caller's random state is restored.
#'
#' @return A `stormbreaker_external_bootstrap` containing recomputed point
#'   targets, every target-by-replicate estimate and status, percentile
#'   intervals, the effective `prediction_target`, replicate accounting,
#'   resampling semantics, and provenance.
#'   Binary and locked-horizon survival scalar targets include
#'   observed/expected mean calibration; locked-horizon survival also includes
#'   Uno C and calibration slope as separate rows.
#'   Threshold-enabled binary and locked-horizon survival evaluations
#'   additionally contain a nested `decision_curve` component with raw
#'   replicate curves and paired pointwise net-benefit-difference intervals.
#' @export
#'
#' @examples
#' development <- sb_example_binary(n = 30, seed = 1)
#' task <- sb_binary_task(
#'   development, "patient_id", "outcome", c("gene_a", "gene_b"),
#'   decision_thresholds = c(0.2, 0.5, 0.8)
#' )
#' fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
#' external <- sb_example_binary(n = 12, seed = 2)
#' external$patient_id <- paste0("external_", seq_len(nrow(external)))
#' prediction <- sb_predict_final(
#'   fit,
#'   external[c("patient_id", "gene_a", "gene_b")],
#'   "patient_id"
#' )
#' endpoint <- sb_binary_endpoint(external, "patient_id", "outcome")
#' evaluation <- sb_evaluate_external(prediction, endpoint)
#' uncertainty <- sb_bootstrap_external(evaluation, reps = 20, seed = 3)
#' uncertainty$decision_curve$intervals
sb_bootstrap_external <- function(
    evaluation,
    reps = 1000L,
    confidence_level = 0.95,
    seed = 1L) {
  argument_attempt <- sb_external_attempt({
    checked_reps <- sb_check_positive_integer(reps, "reps", minimum = 2L)
    checked_seed <- sb_check_positive_integer(seed, "seed")
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
    list(
      reps = checked_reps,
      confidence_level = as.numeric(confidence_level),
      seed = checked_seed
    )
  })
  if (!is.null(argument_attempt$error)) {
    return(sb_external_bootstrap_failure(
      evaluation,
      reps,
      confidence_level,
      seed,
      "arguments",
      argument_attempt$error
    ))
  }
  arguments <- argument_attempt$value

  contract_attempt <- sb_external_attempt(
    sb_external_bootstrap_contract(evaluation)
  )
  if (!is.null(contract_attempt$error)) {
    return(sb_external_bootstrap_failure(
      evaluation,
      arguments$reps,
      arguments$confidence_level,
      arguments$seed,
      "evaluation_contract",
      contract_attempt$error
    ))
  }
  contract <- contract_attempt$value
  point <- sb_external_bootstrap_metric_rows(
    contract$scored,
    contract$task_type
  )
  resample_attempt <- sb_external_attempt(
    sb_external_bootstrap_resamples(
      contract$scored,
      arguments$reps,
      arguments$seed,
      contract$task_type,
      contract$decision_thresholds
    )
  )
  if (!is.null(resample_attempt$error)) {
    return(sb_external_bootstrap_failure(
      evaluation,
      arguments$reps,
      arguments$confidence_level,
      arguments$seed,
      "resampling",
      resample_attempt$error,
      has_horizon = contract$has_horizon,
      prediction_target = contract$prediction_target
    ))
  }
  sampled <- resample_attempt$value
  decision_curve <- NULL
  if (!is.null(contract$decision_thresholds)) {
    decision_curve <- if (identical(contract$task_type, "binary")) {
      sb_external_bootstrap_decision_curve_result(
        contract$scored,
        contract$decision_thresholds,
        sampled$decision_curve,
        arguments$confidence_level,
        arguments$reps
      )
    } else {
      sb_external_bootstrap_survival_decision_curve_result(
        contract$scored,
        contract$decision_thresholds,
        sampled$decision_curve,
        arguments$confidence_level,
        arguments$reps
      )
    }
  }
  intervals <- sb_external_percentile_intervals(
    point,
    sampled$replicates,
    arguments$confidence_level,
    arguments$reps
  )
  complete <- nrow(intervals) >= 1L &&
    all(intervals$status == "success")
  result <- list(
    status = if (complete) "success" else "incomplete",
    task_type = contract$task_type,
    pipeline_id = contract$pipeline_id,
    learner = contract$learner,
    n = contract$n,
    prediction_target = contract$prediction_target,
    reps = arguments$reps,
    confidence_level = arguments$confidence_level,
    seed = arguments$seed,
    interval_method = sb_external_bootstrap_interval_method(
      arguments$confidence_level
    ),
    resampling = sb_external_bootstrap_resampling(
      contract$task_type,
      contract$has_horizon
    ),
    point_estimates = point,
    intervals = intervals,
    replicates = sampled$replicates,
    resample_audit = sampled$audit,
    provenance = sb_external_bootstrap_provenance(),
    error_stage = NA_character_,
    error_class = NA_character_,
    error_message = NA_character_
  )
  if (!is.null(decision_curve)) {
    result <- append(
      result,
      list(decision_curve = decision_curve),
      after = match("resample_audit", names(result))
    )
  }
  class(result) <- "stormbreaker_external_bootstrap"
  result
}

#' @export
print.stormbreaker_external_bootstrap <- function(x, ...) {
  cat("<stormbreaker_external_bootstrap:", x$status, ">\n")
  cat("fixed-prediction external patient bootstrap\n")
  cat("Pipeline:", x$pipeline_id, "| Task:", x$task_type, "\n")
  target <- x[["prediction_target"]]
  if (identical(x$task_type, "survival") && is.list(target)) {
    if (identical(target$target_status, "complete") &&
        isTRUE(target$has_horizon)) {
      cat(
        "Target:", format(target$prediction_horizon), target$time_unit,
        "| Horizon source:", target$horizon_source, "\n"
      )
    } else if (identical(target$target_status, "complete") &&
        identical(target$has_horizon, FALSE)) {
      cat(
        "Target: horizon-free risk | Endpoint time unit:",
        target$time_unit, "\n"
      )
    } else {
      cat("Target: unknown\n")
    }
  }
  cat(
    "Patients:", x$n,
    "| Replicates:", x$reps,
    "| Confidence:", format(x$confidence_level), "\n"
  )
  if (identical(x$status, "failed")) {
    cat("Failed stage:", x$error_stage, "|", x$error_class, "\n")
    return(invisible(x))
  }
  for (index in seq_len(nrow(x$intervals))) {
    interval <- x$intervals[index, , drop = FALSE]
    bounds <- if (isTRUE(interval$interval_available)) {
      paste0("[", format(interval$lower), ", ", format(interval$upper), "]")
    } else {
      "unavailable"
    }
    cat(
      interval$metric, ":", format(interval$point_estimate),
      "|", bounds,
      "|", interval$replicate_completeness,
      paste0("(", interval$n_finite_successful, "/", interval$n_requested, ")"),
      "\n"
    )
  }
  if ("decision_curve" %in% names(x)) {
    if (identical(x$task_type, "binary")) {
      cat(
        "External binary DCA paired pointwise intervals:",
        x$decision_curve$n_thresholds,
        "declared thresholds\n"
      )
    } else {
      status <- if (identical(x$decision_curve$status, "success")) {
        ""
      } else {
        paste0(" | ", x$decision_curve$status)
      }
      cat(
        "External survival DCA paired pointwise intervals:",
        x$decision_curve$n_thresholds,
        paste0("declared thresholds", status, "\n")
      )
    }
  }
  invisible(x)
}
