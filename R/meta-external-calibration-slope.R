sb_external_calibration_slope_meta_member <- function(
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
    "survival calibration slope"
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

  metric <- "calibration_slope"
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
    is.finite(point_rows$estimate[[1L]])
  replicates_complete <- replicate_ids_complete &&
    !anyNA(replicate_rows$status) &&
    all(replicate_rows$status == "success") &&
    !anyNA(replicate_rows$estimate) &&
    all(is.finite(replicate_rows$estimate))
  if (!point_available || !replicates_complete) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` requires one finite successful `", metric,
        "` point estimate and all ", as.integer(reps),
        " successful finite bootstrap estimates, numbered exactly once from ",
        "1 through `reps`. No failed replicate is dropped and no slope is ",
        "transformed or clipped."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  standard_error <- stats::sd(replicate_rows$estimate)
  if (!is.finite(standard_error) || standard_error <= 0) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` has zero or non-finite bootstrap uncertainty for `", metric,
        "` on the original scale."
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
    point_estimate = as.numeric(point_rows$estimate[[1L]]),
    standard_error = as.numeric(standard_error),
    prediction_horizon = metric_target$prediction_horizon,
    time_unit = metric_target$time_unit
  )
}

#' Meta-analyse locked-model external calibration slopes across cohorts
#'
#' `sb_meta_external_calibration_slope()` pools one calibration-slope estimate
#' for the same locked prediction model evaluated in multiple independent
#' external cohorts. It accepts binary results and fixed-horizon right-censored
#' survival results from [sb_bootstrap_external()]. Binary and survival tasks
#' are analysed separately and are never mixed.
#'
#' For binary outcomes, the source slope is from logistic recalibration of the
#' observed outcome on the logit predicted probability with a freely estimated
#' intercept. For survival outcomes, the source slope is the coefficient from
#' an Efron Cox recalibration model after administrative censoring at the locked
#' horizon, using `log(-log(S(h)))` as the transformed prediction. Every
#' survival cohort must retain the exact same time unit and horizon.
#'
#' Calibration slopes are meta-analysed on their original, unbounded scale.
#' The within-cohort standard error is the sample standard deviation of every
#' ordinary patient-bootstrap slope. A random-effects model is fitted with REML
#' and Knapp-Hartung inference through `metafor::rma.uni()`. Finite negative
#' slopes and negative interval bounds are retained rather than clipped. A
#' slope of one is the ideal scaling target under the stated recalibration
#' model; a slope below one commonly indicates predictions that are too
#' extreme, whereas a slope above one commonly indicates insufficient spread.
#' A negative slope indicates inverse ordering under that recalibration model,
#' not merely excessive spread.
#'
#' The point estimate and every requested target bootstrap replicate must occur
#' exactly once and be successful and finite. The function does not drop failed
#' replicates, transform slopes, or substitute another calibration measure.
#' Failures of unrelated bootstrap metrics do not invalidate an otherwise
#' complete calibration-slope target.
#'
#' At least two cohorts are required. Heterogeneity and prediction intervals
#' remain especially uncertain with few cohorts. A pooled slope near one does
#' not establish mean or moderate calibration and can coexist with a poor O:E
#' ratio or calibration curve. This is a univariate synthesis and does not use
#' within-cohort correlations with discrimination or other calibration
#' measures. Matching pipeline and learner labels are bookkeeping only: the
#' caller must establish cohort independence, common eligibility and outcome
#' definitions, and an identical fitted model, preprocessing, predictor
#' orientation, and prediction horizon. Bootstrap uncertainty is conditional
#' on each locked model and prediction set and does not recover development-
#' fitting or model-selection uncertainty.
#'
#' @references
#' Snell KIE, Ensor J, Debray TPA, Moons KGM, Riley RD (2018). Meta-analysis of
#'   prediction model performance across multiple studies: Which scale helps
#'   ensure between-study normality for the C-statistic and calibration
#'   measures? *Statistical Methods in Medical Research*, 27(11), 3505-3522.
#'   \doi{10.1177/0962280217705678}.
#'
#' Debray TPA, Damen JAAG, Snell KIE, et al. (2019). A framework for
#'   meta-analysis of prediction model studies with binary and time-to-event
#'   outcomes. *Statistical Methods in Medical Research*, 28(9), 2768-2786.
#'   \doi{10.1177/0962280218785504}.
#'
#' Snell KIE, Hua H, Debray TPA, et al. (2016). Multivariate meta-analysis of
#'   individual participant data helped externally validate the performance and
#'   implementation of a prediction model. *Journal of Clinical Epidemiology*,
#'   69, 40-50. \doi{10.1016/j.jclinepi.2015.05.009}.
#'
#' Stevens RJ, Poppe KK (2020). Validation of clinical prediction models: what
#'   does the "calibration slope" really measure? *Journal of Clinical
#'   Epidemiology*, 118, 93-99. \doi{10.1016/j.jclinepi.2019.09.016}.
#'
#' @param bootstraps A named list of at least two
#'   `stormbreaker_external_bootstrap` objects. Names must be unique, non-blank
#'   cohort IDs. All objects must describe the same task, pipeline ID, and
#'   learner label.
#' @param confidence_level Confidence level strictly between zero and one for
#'   the pooled confidence interval and between-cohort prediction interval.
#'   Defaults to 0.95.
#'
#' @return A `stormbreaker_external_calibration_slope_meta` containing cohort
#'   slope estimates and bootstrap standard errors, the random-effects summary,
#'   explicit time target, method labels, software provenance, cautions, and the
#'   unchanged named source `bootstraps`. No p value, ranking, winner, or
#'   recalibrated model is returned.
#' @export
#'
#' @examples
#' \dontrun{
#' pooled <- sb_meta_external_calibration_slope(list(
#'   hospital_a = bootstrap_a,
#'   hospital_b = bootstrap_b
#' ))
#' pooled$summary
#' }
sb_meta_external_calibration_slope <- function(
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
    sb_external_calibration_slope_meta_member(
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
          "horizon to meta-analyse calibration slope."
        ),
        "stormbreaker_external_meta_target_mismatch"
      )
    }
    metric_target$prediction_horizon <- horizons[[1L]]
    metric_target$time_unit <- time_units[[1L]]
  }

  estimand <- if (identical(task_type, "binary")) {
    paste0(
      "Random-effects mean calibration slope for one locked binary model ",
      "across an assumed distribution of independent exchangeable cohorts"
    )
  } else {
    paste0(
      "Random-effects mean calibration slope for one locked survival model ",
      "at ", format(metric_target$prediction_horizon), " ",
      metric_target$time_unit,
      " across an assumed distribution of independent exchangeable cohorts"
    )
  }

  if (!requireNamespace("metafor", quietly = TRUE)) {
    sb_abort(
      paste0(
        "Package `metafor` is required for external calibration-slope ",
        "meta-analysis. Install it with `install.packages(\"metafor\")`."
      ),
      "stormbreaker_missing_optional_dependency"
    )
  }

  cohorts <- data.frame(
    cohort_id = cohort_ids,
    n_patients = vapply(members, `[[`, integer(1L), "n"),
    point_estimate = vapply(members, `[[`, numeric(1L), "point_estimate"),
    standard_error = vapply(members, `[[`, numeric(1L), "standard_error"),
    n_bootstrap = vapply(members, `[[`, integer(1L), "reps"),
    stringsAsFactors = FALSE
  )
  fit <- metafor::rma.uni(
    yi = cohorts$point_estimate,
    sei = cohorts$standard_error,
    method = "REML",
    test = "knha",
    level = 100 * confidence_level,
    slab = cohort_ids
  )
  prediction <- stats::predict(fit)
  k <- as.integer(nrow(cohorts))
  summary <- data.frame(
    metric = "calibration_slope",
    k = k,
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
    "A pooled slope near 1 does not establish mean or moderate calibration; ",
    "inspect O:E and calibration curves separately."
  )

  structure(
    list(
      status = "success",
      task_type = task_type,
      pipeline_id = pipeline_ids[[1L]],
      learner = learners[[1L]],
      metric = "calibration_slope",
      metric_target = metric_target,
      confidence_level = as.numeric(confidence_level),
      estimand = estimand,
      interpretation = list(
        target = 1,
        below_one = "predictions_too_extreme_under_recalibration_model",
        above_one = "predictions_too_narrow_under_recalibration_model"
      ),
      method = list(
        model = "random_effects",
        effect_scale = "original_calibration_slope",
        within_cohort_standard_error =
          "sample_sd_of_patient_bootstrap_slope_replicates",
        tau2_estimator = "REML",
        inference = "Knapp-Hartung",
        back_transformation = "none"
      ),
      provenance = sb_external_meta_provenance(),
      cohorts = cohorts,
      summary = summary,
      prediction_interval_caution = caution,
      bootstraps = bootstraps
    ),
    class = "stormbreaker_external_calibration_slope_meta"
  )
}

#' @export
print.stormbreaker_external_calibration_slope_meta <- function(x, ...) {
  confidence <- 100 * x$confidence_level
  summary <- x$summary[1L, , drop = FALSE]
  cat("<stormbreaker_external_calibration_slope_meta:", x$status, ">\n")
  cat("fixed-model original-scale calibration-slope meta-analysis\n")
  cat(
    "Metric: calibration slope | Cohorts:", summary$k,
    "| Task:", x$task_type, "\n"
  )
  if (isTRUE(x$metric_target$has_horizon)) {
    cat(
      "Target:", format(x$metric_target$prediction_horizon),
      x$metric_target$time_unit, "\n"
    )
  }
  cat(
    "Summary slope:", format(summary$estimate),
    paste0("| ", format(confidence), "% CI:"),
    format(summary$lower), "to", format(summary$upper), "\n"
  )
  cat(
    paste0(format(confidence), "% prediction interval:"),
    format(summary$prediction_lower), "to",
    format(summary$prediction_upper), "\n"
  )
  cat(
    "Interpretation: target 1; a slope near 1 does not establish mean or ",
    "moderate calibration.\n",
    sep = ""
  )
  cat(
    "Scope: one univariate performance synthesis; no model ranking or ",
    "recalibration.\n",
    sep = ""
  )
  cat("Caution:", x$prediction_interval_caution, "\n")
  invisible(x)
}
