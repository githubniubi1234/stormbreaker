sb_external_meta_provenance <- function() {
  list(
    package = "stormbreaker",
    package_version = as.character(utils::packageVersion("stormbreaker")),
    r_version = R.version.string,
    metafor_version = as.character(utils::packageVersion("metafor"))
  )
}

sb_external_prediction_error_meta_metric <- function(task_type, metric) {
  if (is.null(metric)) {
    return(if (identical(task_type, "binary")) {
      "brier_score"
    } else {
      "ipcw_brier_score"
    })
  }
  if (!sb_external_scalar_character(metric)) {
    sb_abort(
      "`metric` must be NULL or one non-empty metric name.",
      "stormbreaker_bad_external_meta_metric"
    )
  }
  allowed <- if (identical(task_type, "binary")) {
    c("brier_score", "log_loss")
  } else {
    "ipcw_brier_score"
  }
  if (!metric %in% allowed) {
    sb_abort(
      paste0(
        "For a ", task_type, " task, `metric` must be one of: ",
        paste0("`", allowed, "`", collapse = ", "), "."
      ),
      "stormbreaker_bad_external_meta_metric"
    )
  }
  metric
}

sb_external_prediction_error_meta_target <- function(
    comparison,
    cohort_id,
    metric) {
  if (!identical(metric, "ipcw_brier_score")) {
    return(list(prediction_horizon = NULL, time_unit = NULL))
  }
  target <- comparison[["prediction_target"]]
  valid_target <- is.list(target) &&
    identical(target[["target_status"]], "complete") &&
    identical(target[["has_horizon"]], TRUE) &&
    is.numeric(target[["prediction_horizon"]]) &&
    length(target[["prediction_horizon"]]) == 1L &&
    !is.na(target[["prediction_horizon"]]) &&
    is.finite(target[["prediction_horizon"]]) &&
    target[["prediction_horizon"]] > 0 &&
    sb_external_scalar_character(target[["time_unit"]]) &&
    target[["time_unit"]] %in% c("day", "month", "year")
  if (!valid_target) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` requires a complete locked horizon and ",
        "time unit to meta-analyse `", metric, "` differences."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }
  list(
    prediction_horizon = as.numeric(target[["prediction_horizon"]]),
    time_unit = target[["time_unit"]]
  )
}

sb_external_paired_metric_estimates_valid <- function(values, metric) {
  if (length(values) == 0L || anyNA(values) || !all(is.finite(values))) {
    return(FALSE)
  }
  if (metric %in% c(
    "roc_auc", "harrell_c", "time_auc", "uno_c", "brier_score"
  )) {
    return(all(values >= 0 & values <= 1))
  }
  all(values >= 0)
}

sb_external_paired_difference_estimates_valid <- function(values, metric) {
  if (length(values) == 0L || anyNA(values) || !all(is.finite(values))) {
    return(FALSE)
  }
  if (metric %in% c(
    "roc_auc", "harrell_c", "time_auc", "uno_c", "brier_score"
  )) {
    return(all(values >= -1 & values <= 1))
  }
  TRUE
}

sb_external_paired_difference_meta_member <- function(
    comparison,
    cohort_id,
    candidate_pipeline_id,
    metric = NULL,
    metric_family = c("discrimination", "prediction_error")) {
  metric_family <- match.arg(metric_family)
  if (!inherits(
    comparison,
    "stormbreaker_external_bootstrap_comparison"
  ) || !is.list(comparison)) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` is not a result from ",
        "`sb_bootstrap_external_comparison()`."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  task_type <- comparison[["task_type"]]
  reference_pipeline_id <- comparison[["reference_pipeline_id"]]
  n <- comparison[["n"]]
  reps <- comparison[["reps"]]
  valid_metadata <- sb_external_scalar_character(comparison[["status"]]) &&
    comparison[["status"]] %in% c("success", "incomplete") &&
    sb_external_scalar_character(task_type) &&
    task_type %in% c("binary", "survival") &&
    sb_external_scalar_character(reference_pipeline_id) &&
    is.numeric(n) && length(n) == 1L && !is.na(n) && is.finite(n) &&
    n >= 1 && n <= .Machine$integer.max && n == as.integer(n) &&
    is.numeric(reps) && length(reps) == 1L && !is.na(reps) &&
    is.finite(reps) && reps >= 2 && reps <= .Machine$integer.max &&
    reps == as.integer(reps)
  if (!valid_metadata) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` has incomplete paired-bootstrap metadata."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  candidate_status <- comparison[["candidate_status"]]
  valid_candidate_table <- is.data.frame(candidate_status) &&
    all(c("candidate_pipeline_id", "status") %in% names(candidate_status)) &&
    is.character(candidate_status$candidate_pipeline_id) &&
    is.character(candidate_status$status)
  if (!valid_candidate_table) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` does not retain usable candidate-status metadata."
      ),
      "stormbreaker_bad_external_meta"
    )
  }
  candidate_rows <- candidate_status[
    !is.na(candidate_status$candidate_pipeline_id) &
      candidate_status$candidate_pipeline_id == candidate_pipeline_id,
    ,
    drop = FALSE
  ]
  candidate_available <- nrow(candidate_rows) == 1L &&
    identical(candidate_rows$status[[1L]], "success")
  if (!candidate_available) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` requires one successful candidate `",
        candidate_pipeline_id, "` in the paired bootstrap comparison."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  if (identical(metric_family, "discrimination")) {
    metric <- sb_external_discrimination_meta_metric(task_type, metric)
    metric_target <- sb_external_discrimination_meta_target(
      comparison,
      cohort_id,
      metric
    )
    expected_direction <- "higher_is_better"
  } else {
    metric <- sb_external_prediction_error_meta_metric(task_type, metric)
    metric_target <- sb_external_prediction_error_meta_target(
      comparison,
      cohort_id,
      metric
    )
    expected_direction <- "lower_is_better"
  }
  metrics <- comparison[["metrics"]]
  if (!is.character(metrics) || !metric %in% metrics) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` does not retain the requested paired `",
        metric, "` metric."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  point <- comparison[["point_differences"]]
  replicates <- comparison[["differences"]]
  required_point <- c(
    "reference_pipeline_id", "candidate_pipeline_id", "contrast", "metric",
    "metric_direction", "reference_estimate", "candidate_estimate",
    "estimate", "status"
  )
  required_replicates <- c("replicate", required_point)
  valid_tables <- is.data.frame(point) &&
    all(required_point %in% names(point)) &&
    is.character(point$reference_pipeline_id) &&
    is.character(point$candidate_pipeline_id) &&
    is.character(point$contrast) && is.character(point$metric) &&
    is.character(point$metric_direction) &&
    is.numeric(point$reference_estimate) &&
    is.numeric(point$candidate_estimate) && is.numeric(point$estimate) &&
    is.character(point$status) && is.data.frame(replicates) &&
    all(required_replicates %in% names(replicates)) &&
    is.numeric(replicates$replicate) &&
    is.character(replicates$reference_pipeline_id) &&
    is.character(replicates$candidate_pipeline_id) &&
    is.character(replicates$contrast) && is.character(replicates$metric) &&
    is.character(replicates$metric_direction) &&
    is.numeric(replicates$reference_estimate) &&
    is.numeric(replicates$candidate_estimate) &&
    is.numeric(replicates$estimate) && is.character(replicates$status)
  if (!valid_tables) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` does not retain usable paired-difference tables."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  point_rows <- point[
    !is.na(point$reference_pipeline_id) &
      point$reference_pipeline_id == reference_pipeline_id &
      !is.na(point$candidate_pipeline_id) &
      point$candidate_pipeline_id == candidate_pipeline_id &
      !is.na(point$metric) & point$metric == metric,
    ,
    drop = FALSE
  ]
  replicate_rows <- replicates[
    !is.na(replicates$reference_pipeline_id) &
      replicates$reference_pipeline_id == reference_pipeline_id &
      !is.na(replicates$candidate_pipeline_id) &
      replicates$candidate_pipeline_id == candidate_pipeline_id &
      !is.na(replicates$metric) & replicates$metric == metric,
    ,
    drop = FALSE
  ]
  point_difference <-
    point_rows$candidate_estimate - point_rows$reference_estimate
  replicate_differences <-
    replicate_rows$candidate_estimate - replicate_rows$reference_estimate
  point_arithmetic_consistent <- identical(
    unname(point_rows$estimate),
    unname(point_difference)
  )
  replicate_arithmetic_consistent <- identical(
    unname(replicate_rows$estimate),
    unname(replicate_differences)
  )
  expected_replicates <- seq_len(as.integer(reps))
  replicate_ids_complete <- nrow(replicate_rows) == as.integer(reps) &&
    !anyNA(replicate_rows$replicate) &&
    all(is.finite(replicate_rows$replicate)) &&
    all(replicate_rows$replicate >= 1 &
      replicate_rows$replicate <= as.integer(reps)) &&
    all(replicate_rows$replicate == floor(replicate_rows$replicate)) &&
    identical(
      sort(as.integer(replicate_rows$replicate)),
      expected_replicates
    )
  point_available <- nrow(point_rows) == 1L &&
    point_arithmetic_consistent &&
    identical(point_rows$status[[1L]], "success") &&
    identical(point_rows$contrast[[1L]], "candidate_minus_reference") &&
    identical(point_rows$metric_direction[[1L]], expected_direction) &&
    sb_external_paired_metric_estimates_valid(
      c(point_rows$reference_estimate, point_rows$candidate_estimate),
      metric
    ) &&
    sb_external_paired_difference_estimates_valid(
      point_difference,
      metric
    )
  replicates_complete <- replicate_ids_complete &&
    replicate_arithmetic_consistent &&
    !anyNA(replicate_rows$status) &&
    all(replicate_rows$status == "success") &&
    !anyNA(replicate_rows$contrast) &&
    all(replicate_rows$contrast == "candidate_minus_reference") &&
    !anyNA(replicate_rows$metric_direction) &&
    all(replicate_rows$metric_direction == expected_direction) &&
    sb_external_paired_metric_estimates_valid(
      c(
        replicate_rows$reference_estimate,
        replicate_rows$candidate_estimate
      ),
      metric
    ) &&
    sb_external_paired_difference_estimates_valid(
      replicate_differences,
      metric
    )
  if (!point_available || !replicates_complete) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` requires one finite successful paired `",
        metric, "` point difference and all ", as.integer(reps),
        " paired bootstrap differences, numbered exactly once from 1 through ",
        "`reps`. Every saved difference must equal candidate minus reference; ",
        "no failed replicate is dropped and no difference is clipped."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  standard_error <- stats::sd(replicate_differences)
  if (!is.finite(standard_error) || standard_error <= 0) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` has zero or non-finite paired bootstrap uncertainty for `", metric,
        "`."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  list(
    task_type = task_type,
    reference_pipeline_id = reference_pipeline_id,
    candidate_pipeline_id = candidate_pipeline_id,
    metric = metric,
    n = as.integer(n),
    reps = as.integer(reps),
    reference_point_estimate = point_rows$reference_estimate[[1L]],
    candidate_point_estimate = point_rows$candidate_estimate[[1L]],
    point_difference = as.numeric(point_difference[[1L]]),
    standard_error = as.numeric(standard_error),
    prediction_horizon = metric_target$prediction_horizon,
    time_unit = metric_target$time_unit
  )
}

sb_external_paired_difference_meta_core <- function(
    comparisons,
    candidate_pipeline_id,
    confidence_level,
    metric,
    metric_family = c("discrimination", "prediction_error")) {
  metric_family <- match.arg(metric_family)
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
  if (!sb_external_scalar_character(candidate_pipeline_id)) {
    sb_abort(
      "`candidate_pipeline_id` must be one non-empty pipeline ID.",
      "stormbreaker_bad_external_meta"
    )
  }

  valid_input <- is.list(comparisons) && !inherits(
    comparisons,
    "stormbreaker_external_bootstrap_comparison"
  )
  if (!valid_input) {
    sb_abort(
      paste0(
        "`comparisons` must be a named list of results from ",
        "`sb_bootstrap_external_comparison()`."
      ),
      "stormbreaker_bad_external_meta"
    )
  }
  if (length(comparisons) < 2L) {
    sb_abort(
      "At least two independent external cohorts are required.",
      "stormbreaker_insufficient_external_cohorts"
    )
  }

  cohort_ids <- names(comparisons)
  valid_ids <- !is.null(cohort_ids) &&
    length(cohort_ids) == length(comparisons) &&
    !anyNA(cohort_ids) && all(nzchar(trimws(cohort_ids))) &&
    anyDuplicated(cohort_ids) == 0L
  if (!valid_ids) {
    sb_abort(
      "Outer list names must be unique, non-blank cohort IDs.",
      "stormbreaker_bad_external_meta"
    )
  }

  members <- lapply(seq_along(comparisons), function(index) {
    sb_external_paired_difference_meta_member(
      comparisons[[index]],
      cohort_ids[[index]],
      candidate_pipeline_id,
      metric,
      metric_family
    )
  })
  task_types <- vapply(members, `[[`, character(1L), "task_type")
  if (length(unique(task_types)) != 1L) {
    sb_abort(
      "All cohorts must use the same task type.",
      "stormbreaker_external_meta_task_mismatch"
    )
  }
  reference_ids <- vapply(
    members,
    `[[`,
    character(1L),
    "reference_pipeline_id"
  )
  if (length(unique(reference_ids)) != 1L) {
    sb_abort(
      paste0(
        "All cohorts must retain the same reference and candidate pipeline ",
        "IDs. Matching labels do not prove identical fitted states."
      ),
      "stormbreaker_external_meta_model_mismatch"
    )
  }

  metric <- members[[1L]]$metric
  has_horizon <- metric %in% c("time_auc", "uno_c", "ipcw_brier_score")
  metric_target <- list(
    has_horizon = has_horizon,
    prediction_horizon = NULL,
    time_unit = NULL
  )
  if (has_horizon) {
    horizons <- vapply(members, `[[`, numeric(1L), "prediction_horizon")
    time_units <- vapply(members, `[[`, character(1L), "time_unit")
    same_horizon <- all(horizons == horizons[[1L]])
    same_unit <- all(time_units == time_units[[1L]])
    if (!same_horizon || !same_unit) {
      sb_abort(
        paste0(
          "All cohorts must use the same time unit and prediction horizon ",
          "to meta-analyse paired `", metric, "` differences."
        ),
        "stormbreaker_external_meta_target_mismatch"
      )
    }
    metric_target$prediction_horizon <- horizons[[1L]]
    metric_target$time_unit <- time_units[[1L]]
  }

  if (!requireNamespace("metafor", quietly = TRUE)) {
    sb_abort(
      paste0(
        "Package `metafor` is required for paired external difference ",
        "meta-analysis. Install it with `install.packages(\"metafor\")`."
      ),
      "stormbreaker_missing_optional_dependency"
    )
  }

  cohorts <- data.frame(
    cohort_id = cohort_ids,
    n_patients = vapply(members, `[[`, integer(1L), "n"),
    reference_point_estimate = vapply(
      members,
      `[[`,
      numeric(1L),
      "reference_point_estimate"
    ),
    candidate_point_estimate = vapply(
      members,
      `[[`,
      numeric(1L),
      "candidate_point_estimate"
    ),
    point_difference = vapply(
      members,
      `[[`,
      numeric(1L),
      "point_difference"
    ),
    standard_error = vapply(members, `[[`, numeric(1L), "standard_error"),
    n_bootstrap = vapply(members, `[[`, integer(1L), "reps"),
    stringsAsFactors = FALSE
  )
  fit <- metafor::rma.uni(
    yi = cohorts$point_difference,
    sei = cohorts$standard_error,
    method = "REML",
    test = "knha",
    level = 100 * confidence_level,
    slab = cohort_ids
  )
  prediction <- stats::predict(fit)
  summary <- data.frame(
    metric = metric,
    contrast = "candidate_minus_reference",
    k = as.integer(nrow(cohorts)),
    estimate = as.numeric(fit$b),
    lower = as.numeric(fit$ci.lb),
    upper = as.numeric(fit$ci.ub),
    prediction_lower = as.numeric(prediction$pi.lb),
    prediction_upper = as.numeric(prediction$pi.ub),
    tau2 = as.numeric(fit$tau2),
    tau = sqrt(as.numeric(fit$tau2)),
    i2_percent = as.numeric(fit$I2),
    stringsAsFactors = FALSE
  )

  list(
    task_type = task_types[[1L]],
    reference_pipeline_id = reference_ids[[1L]],
    candidate_pipeline_id = candidate_pipeline_id,
    metric = metric,
    metric_target = metric_target,
    confidence_level = as.numeric(confidence_level),
    cohorts = cohorts,
    summary = summary,
    comparisons = comparisons
  )
}
