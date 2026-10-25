sb_external_calibration_meta_target <- function(
    bootstrap,
    cohort_id,
    task_type,
    estimand_label) {
  if (identical(task_type, "binary")) {
    return(list(prediction_horizon = NULL, time_unit = NULL))
  }
  target <- bootstrap[["prediction_target"]]
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
        "time unit to meta-analyse ", estimand_label, "."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }
  list(
    prediction_horizon = as.numeric(target[["prediction_horizon"]]),
    time_unit = target[["time_unit"]]
  )
}

sb_external_mean_calibration_meta_member <- function(
    bootstrap,
    cohort_id) {
  if (!inherits(bootstrap, "stormbreaker_external_bootstrap") ||
      !is.list(bootstrap)) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` is not a result from `sb_bootstrap_external()`."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  task_type <- bootstrap[["task_type"]]
  pipeline_id <- bootstrap[["pipeline_id"]]
  learner <- bootstrap[["learner"]]
  n <- bootstrap[["n"]]
  reps <- bootstrap[["reps"]]
  valid_metadata <- sb_external_scalar_character(task_type) &&
    task_type %in% c("binary", "survival") &&
    sb_external_scalar_character(pipeline_id) &&
    sb_external_scalar_character(learner) &&
    is.numeric(n) && length(n) == 1L && !is.na(n) && is.finite(n) &&
    n >= 1 && n <= .Machine$integer.max && n == as.integer(n) &&
    is.numeric(reps) && length(reps) == 1L && !is.na(reps) &&
    is.finite(reps) && reps >= 2 && reps <= .Machine$integer.max &&
    reps == as.integer(reps)
  if (!valid_metadata) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` has incomplete external-bootstrap metadata."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  metric_target <- sb_external_calibration_meta_target(
    bootstrap,
    cohort_id,
    task_type,
    "survival mean calibration"
  )
  point <- bootstrap[["point_estimates"]]
  replicates <- bootstrap[["replicates"]]
  valid_tables <- is.data.frame(point) &&
    all(c("metric", "estimate", "status") %in% names(point)) &&
    is.character(point$metric) && is.numeric(point$estimate) &&
    is.character(point$status) &&
    is.data.frame(replicates) &&
    all(c("replicate", "metric", "estimate", "status") %in%
      names(replicates)) &&
    is.numeric(replicates$replicate) &&
    is.character(replicates$metric) && is.numeric(replicates$estimate) &&
    is.character(replicates$status)
  if (!valid_tables) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` does not retain usable bootstrap metric tables."
      ),
      "stormbreaker_bad_external_meta"
    )
  }

  metric <- "mean_calibration_ratio"
  point_rows <- point[
    !is.na(point$metric) & point$metric == metric,
    ,
    drop = FALSE
  ]
  replicate_rows <- replicates[
    !is.na(replicates$metric) & replicates$metric == metric,
    ,
    drop = FALSE
  ]
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
    identical(point_rows$status[[1L]], "success") &&
    length(point_rows$estimate) == 1L &&
    !is.na(point_rows$estimate[[1L]]) &&
    is.finite(point_rows$estimate[[1L]]) &&
    point_rows$estimate[[1L]] > 0
  replicates_complete <- replicate_ids_complete &&
    !anyNA(replicate_rows$status) &&
    all(replicate_rows$status == "success") &&
    !anyNA(replicate_rows$estimate) &&
    all(is.finite(replicate_rows$estimate)) &&
    all(replicate_rows$estimate > 0)
  if (!point_available || !replicates_complete) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` requires one finite successful `", metric,
        "` point estimate and all ", as.integer(reps),
        " successful bootstrap estimates, numbered exactly once from 1 ",
        "through `reps`, and strictly greater than zero. No failed replicate ",
        "is dropped and no continuity correction is applied."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  point_estimate <- point_rows$estimate[[1L]]
  log_replicates <- log(replicate_rows$estimate)
  log_se <- stats::sd(log_replicates)
  if (!is.finite(log_se) || log_se <= 0) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` has zero or non-finite bootstrap uncertainty for `", metric,
        "` on the log scale."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  list(
    task_type = task_type,
    pipeline_id = pipeline_id,
    learner = learner,
    n = as.integer(n),
    reps = as.integer(reps),
    point_estimate = as.numeric(point_estimate),
    log_estimate = log(point_estimate),
    log_se = as.numeric(log_se),
    prediction_horizon = metric_target$prediction_horizon,
    time_unit = metric_target$time_unit
  )
}

#' Meta-analyse locked-model external mean calibration across cohorts
#'
#' `sb_meta_external_mean_calibration()` pools one observed-to-expected event
#' ratio (O:E) for the same locked prediction model evaluated in multiple
#' independent external cohorts. It accepts binary results and fixed-horizon
#' right-censored survival results from [sb_bootstrap_external()]. It does not
#' pool calibration slopes, calibration curves, prediction error, or clinical
#' utility.
#'
#' For binary outcomes, O:E is the observed event proportion divided by the
#' mean predicted event probability. For survival outcomes, it is the
#' Kaplan-Meier observed event probability divided by the mean predicted event
#' probability at the locked horizon. Every survival cohort must retain the
#' exact same time unit and horizon. O:E equal to one indicates agreement in
#' mean risk; values below one indicate average overprediction, and values above
#' one indicate average underprediction. O:E alone cannot reveal miscalibration
#' within the prediction range.
#'
#' Each cohort estimate is transformed with `log()`, and its within-cohort
#' standard error is the sample standard deviation of all ordinary
#' patient-bootstrap estimates on that log scale. A random-effects model is
#' fitted with REML and Knapp-Hartung inference through `metafor::rma.uni()`.
#' The pooled estimate, confidence interval, and prediction interval are
#' exponentiated back to the ratio scale. Between-cohort variance and its square
#' root remain on the log O:E scale.
#'
#' The point estimate and every requested target bootstrap replicate must occur
#' exactly once, be successful, finite, and strictly positive. The function
#' neither drops unavailable replicates nor applies a continuity correction.
#' This means that a cohort with zero observed events in the point estimate or
#' any bootstrap draw is not eligible for this log-ratio meta-analysis.
#' Failures of unrelated bootstrap metrics do not invalidate an otherwise
#' complete O:E target.
#'
#' At least two cohorts are required. Heterogeneity and prediction intervals
#' remain especially uncertain with few cohorts, and a pooled O:E near one can
#' hide opposite directions of miscalibration across cohorts. Matching pipeline
#' and learner labels are bookkeeping only: the caller must establish cohort
#' independence, common eligibility and outcome definitions, and an identical
#' fitted model, preprocessing, predictor orientation, and prediction horizon.
#' Bootstrap uncertainty is conditional on each locked model and prediction set
#' and does not recover development-fitting or model-selection uncertainty.
#'
#' @references
#' Debray TPA, Damen JAAG, Snell KIE, et al. (2019). A framework for
#'   meta-analysis of prediction model studies with binary and time-to-event
#'   outcomes. *Statistical Methods in Medical Research*, 28(9), 2768-2786.
#'   \doi{10.1177/0962280218785504}.
#'
#' Snell KIE, Ensor J, Debray TPA, Moons KGM, Riley RD (2018). Meta-analysis of
#'   prediction model performance across multiple studies: Which scale helps
#'   ensure between-study normality for the C-statistic and calibration
#'   measures? *Statistical Methods in Medical Research*, 27(11), 3505-3522.
#'   \doi{10.1177/0962280217705678}.
#'
#' @param bootstraps A named list of at least two
#'   `stormbreaker_external_bootstrap` objects. Names must be unique, non-blank
#'   cohort IDs. All objects must describe the same task, pipeline ID, and
#'   learner label.
#' @param confidence_level Confidence level strictly between zero and one for
#'   the pooled confidence interval and between-cohort prediction interval.
#'   Defaults to 0.95.
#'
#' @return A `stormbreaker_external_mean_calibration_meta` containing cohort
#'   O:E estimates and bootstrap standard errors, the random-effects summary,
#'   explicit time target, method labels, software provenance, cautions, and the
#'   unchanged named source `bootstraps`. No p value, ranking, winner, or
#'   recalibrated model is returned.
#' @export
#'
#' @examples
#' \dontrun{
#' pooled <- sb_meta_external_mean_calibration(list(
#'   hospital_a = bootstrap_a,
#'   hospital_b = bootstrap_b
#' ))
#' pooled$summary
#' }
sb_meta_external_mean_calibration <- function(
    bootstraps,
    confidence_level = 0.95) {
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

  valid_input <- is.list(bootstraps) &&
    !inherits(bootstraps, "stormbreaker_external_bootstrap")
  if (!valid_input) {
    sb_abort(
      paste0(
        "`bootstraps` must be a named list of results from ",
        "`sb_bootstrap_external()`."
      ),
      "stormbreaker_bad_external_meta"
    )
  }
  if (length(bootstraps) < 2L) {
    sb_abort(
      "At least two independent external cohorts are required.",
      "stormbreaker_insufficient_external_cohorts"
    )
  }

  cohort_ids <- names(bootstraps)
  valid_ids <- !is.null(cohort_ids) &&
    length(cohort_ids) == length(bootstraps) &&
    !anyNA(cohort_ids) && all(nzchar(trimws(cohort_ids))) &&
    anyDuplicated(cohort_ids) == 0L
  if (!valid_ids) {
    sb_abort(
      "Outer list names must be unique, non-blank cohort IDs.",
      "stormbreaker_bad_external_meta"
    )
  }

  members <- lapply(seq_along(bootstraps), function(index) {
    sb_external_mean_calibration_meta_member(
      bootstraps[[index]],
      cohort_ids[[index]]
    )
  })
  task_types <- vapply(members, `[[`, character(1L), "task_type")
  if (length(unique(task_types)) != 1L) {
    sb_abort(
      "All cohorts must use the same task type.",
      "stormbreaker_external_meta_task_mismatch"
    )
  }
  pipeline_ids <- vapply(members, `[[`, character(1L), "pipeline_id")
  learners <- vapply(members, `[[`, character(1L), "learner")
  if (length(unique(pipeline_ids)) != 1L ||
      length(unique(learners)) != 1L) {
    sb_abort(
      paste0(
        "All cohorts must retain the same pipeline ID and learner label. ",
        "These labels do not prove identical fitted state or formula."
      ),
      "stormbreaker_external_meta_model_mismatch"
    )
  }

  task_type <- task_types[[1L]]
  metric_target <- list(
    has_horizon = identical(task_type, "survival"),
    prediction_horizon = NULL,
    time_unit = NULL
  )
  if (identical(task_type, "survival")) {
    horizons <- vapply(members, `[[`, numeric(1L), "prediction_horizon")
    time_units <- vapply(members, `[[`, character(1L), "time_unit")
    if (!all(horizons == horizons[[1L]]) ||
        !all(time_units == time_units[[1L]])) {
      sb_abort(
        paste0(
          "All survival cohorts must use the same time unit and prediction ",
          "horizon to meta-analyse mean calibration."
        ),
        "stormbreaker_external_meta_target_mismatch"
      )
    }
    metric_target$prediction_horizon <- horizons[[1L]]
    metric_target$time_unit <- time_units[[1L]]
  }

  estimand <- if (identical(task_type, "binary")) {
    paste0(
      "Exponentiated random-effects mean of log observed-to-expected event ",
      "ratio for one locked binary model across an assumed distribution of ",
      "independent exchangeable cohorts"
    )
  } else {
    paste0(
      "Exponentiated random-effects mean of log observed-to-expected event ",
      "ratio for one locked survival model at ",
      format(metric_target$prediction_horizon), " ", metric_target$time_unit,
      " across an assumed distribution of independent exchangeable cohorts"
    )
  }

  if (!requireNamespace("metafor", quietly = TRUE)) {
    sb_abort(
      paste0(
        "Package `metafor` is required for external mean-calibration ",
        "meta-analysis. Install it with `install.packages(\"metafor\")`."
      ),
      "stormbreaker_missing_optional_dependency"
    )
  }

  cohorts <- data.frame(
    cohort_id = cohort_ids,
    n_patients = vapply(members, `[[`, integer(1L), "n"),
    point_estimate = vapply(members, `[[`, numeric(1L), "point_estimate"),
    log_estimate = vapply(members, `[[`, numeric(1L), "log_estimate"),
    log_se = vapply(members, `[[`, numeric(1L), "log_se"),
    n_bootstrap = vapply(members, `[[`, integer(1L), "reps"),
    stringsAsFactors = FALSE
  )
  fit <- metafor::rma.uni(
    yi = cohorts$log_estimate,
    sei = cohorts$log_se,
    method = "REML",
    test = "knha",
    level = 100 * confidence_level,
    slab = cohort_ids
  )
  prediction <- stats::predict(fit)
  k <- as.integer(nrow(cohorts))
  summary <- data.frame(
    metric = "mean_calibration_ratio",
    k = k,
    estimate = exp(as.numeric(fit$b)),
    lower = exp(as.numeric(fit$ci.lb)),
    upper = exp(as.numeric(fit$ci.ub)),
    prediction_lower = exp(as.numeric(prediction$pi.lb)),
    prediction_upper = exp(as.numeric(prediction$pi.ub)),
    tau2_log = as.numeric(fit$tau2),
    tau_log = sqrt(as.numeric(fit$tau2)),
    i2_percent = as.numeric(fit$I2),
    stringsAsFactors = FALSE
  )
  caution <- paste0(
    if (k < 10L) {
      paste0(
        "Only ", k, " cohorts were supplied; heterogeneity and the ",
        "between-cohort prediction interval are especially uncertain. "
      )
    } else {
      paste0(
        "The prediction interval assumes independent, exchangeable cohorts ",
        "and a correctly specified random-effects model. "
      )
    },
    "A pooled O:E near 1 can hide opposite miscalibration across cohorts or ",
    "within the prediction range."
  )

  structure(
    list(
      status = "success",
      task_type = task_type,
      pipeline_id = pipeline_ids[[1L]],
      learner = learners[[1L]],
      metric = "mean_calibration_ratio",
      metric_target = metric_target,
      confidence_level = as.numeric(confidence_level),
      estimand = estimand,
      interpretation = list(
        target = 1,
        below_one = "average_overprediction",
        above_one = "average_underprediction"
      ),
      method = list(
        model = "random_effects",
        effect_scale = "log_observed_expected_ratio",
        within_cohort_standard_error =
          "sample_sd_of_log_patient_bootstrap_replicates",
        tau2_estimator = "REML",
        inference = "Knapp-Hartung",
        back_transformation = "exponential"
      ),
      provenance = sb_external_meta_provenance(),
      cohorts = cohorts,
      summary = summary,
      prediction_interval_caution = caution,
      bootstraps = bootstraps
    ),
    class = "stormbreaker_external_mean_calibration_meta"
  )
}

#' @export
print.stormbreaker_external_mean_calibration_meta <- function(x, ...) {
  confidence <- 100 * x$confidence_level
  summary <- x$summary[1L, , drop = FALSE]
  cat("<stormbreaker_external_mean_calibration_meta:", x$status, ">\n")
  cat("fixed-model external mean calibration meta-analysis\n")
  cat("Metric: O:E | Cohorts:", summary$k, "| Task:", x$task_type, "\n")
  if (isTRUE(x$metric_target$has_horizon)) {
    cat(
      "Target:", format(x$metric_target$prediction_horizon),
      x$metric_target$time_unit, "\n"
    )
  }
  cat(
    "Summary O:E:", format(summary$estimate),
    paste0("| ", format(confidence), "% CI:"),
    format(summary$lower), "to", format(summary$upper), "\n"
  )
  cat(
    paste0(format(confidence), "% prediction interval:"),
    format(summary$prediction_lower), "to",
    format(summary$prediction_upper), "\n"
  )
  cat(
    "Interpretation: O:E below 1 means overprediction; above 1 means ",
    "underprediction.\n",
    sep = ""
  )
  cat(
    "Scope: mean calibration only; this does not pool calibration slope or ",
    "curves.\n",
    sep = ""
  )
  cat("Caution:", x$prediction_interval_caution, "\n")
  invisible(x)
}
