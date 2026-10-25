sb_external_discrimination_meta_metric <- function(task_type, metric) {
  if (is.null(metric)) {
    return(if (identical(task_type, "binary")) "roc_auc" else "harrell_c")
  }
  if (!sb_external_scalar_character(metric)) {
    sb_abort(
      "`metric` must be NULL or one non-empty metric name.",
      "stormbreaker_bad_external_meta_metric"
    )
  }
  allowed <- if (identical(task_type, "binary")) {
    "roc_auc"
  } else {
    c("harrell_c", "time_auc", "uno_c")
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

sb_external_discrimination_meta_target <- function(
    bootstrap,
    cohort_id,
    metric) {
  if (!metric %in% c("time_auc", "uno_c")) {
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
        "time unit to meta-analyse `", metric, "`."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }
  list(
    prediction_horizon = as.numeric(target[["prediction_horizon"]]),
    time_unit = target[["time_unit"]]
  )
}

sb_external_discrimination_meta_member <- function(
    bootstrap,
    cohort_id,
    metric = NULL) {
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

  metric <- sb_external_discrimination_meta_metric(task_type, metric)
  metric_target <- sb_external_discrimination_meta_target(
    bootstrap,
    cohort_id,
    metric
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

  point_rows <- point[!is.na(point$metric) & point$metric == metric, , drop = FALSE]
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
    point_rows$estimate[[1L]] > 0 && point_rows$estimate[[1L]] < 1
  replicates_complete <- replicate_ids_complete &&
    !anyNA(replicate_rows$status) &&
    all(replicate_rows$status == "success") &&
    !anyNA(replicate_rows$estimate) &&
    all(is.finite(replicate_rows$estimate)) &&
    all(replicate_rows$estimate > 0 & replicate_rows$estimate < 1)
  if (!point_available || !replicates_complete) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id, "` requires one finite successful `", metric,
        "` point estimate and all ", as.integer(reps),
        " successful bootstrap estimates, numbered exactly once from 1 ",
        "through `reps`, and strictly between 0 and 1. No failed replicate ",
        "is dropped and no boundary estimate is clipped."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  point_estimate <- point_rows$estimate[[1L]]
  logit_replicates <- stats::qlogis(replicate_rows$estimate)
  logit_se <- stats::sd(logit_replicates)
  if (!is.finite(logit_se) || logit_se <= 0) {
    sb_abort(
      paste0(
        "Cohort `", cohort_id,
        "` has zero or non-finite bootstrap uncertainty for `", metric,
        "` on the logit scale."
      ),
      "stormbreaker_external_meta_target_unavailable"
    )
  }

  list(
    task_type = task_type,
    pipeline_id = pipeline_id,
    learner = learner,
    metric = metric,
    n = as.integer(n),
    reps = as.integer(reps),
    point_estimate = as.numeric(point_estimate),
    logit_estimate = stats::qlogis(point_estimate),
    logit_se = as.numeric(logit_se),
    prediction_horizon = metric_target$prediction_horizon,
    time_unit = metric_target$time_unit
  )
}

#' Meta-analyse locked-model external discrimination across cohorts
#'
#' `sb_meta_external_discrimination()` performs one deliberately narrow
#' performance meta-analysis for the same locked prediction model evaluated in
#' multiple independent external cohorts. The default is ROC AUC for binary
#' tasks and horizon-free Harrell concordance for survival tasks. A survival
#' caller may instead explicitly request fixed-horizon cumulative/dynamic
#' `time_auc` or truncated `uno_c`; the function never chooses among these
#' estimands automatically. It does not pool calibration, prediction error, or
#' clinical utility.
#'
#' Fixed-horizon `time_auc` and `uno_c` require every source bootstrap to retain
#' a complete `prediction_target` and all cohorts to share the exact same time
#' unit and prediction horizon. A horizon may have been supplied by the endpoint
#' or inherited from the locked prediction; only the effective target must
#' match. Harrell C remains horizon-free even when its source bootstrap also
#' contains horizon-specific metrics. The three survival discrimination
#' definitions are never pooled together.
#'
#' Each input is a completed [sb_bootstrap_external()] result. The cohort point
#' estimate is transformed with `qlogis()`, and its within-cohort standard error
#' is the sample standard deviation of all ordinary patient-bootstrap estimates
#' on that transformed scale. A random-effects model is then fitted with REML
#' and Knapp-Hartung inference through `metafor::rma.uni()`. The pooled estimate,
#' confidence interval, and prediction interval are back-transformed with
#' `plogis()`. Between-cohort variance and its square root remain on the logit
#' discrimination scale.
#'
#' Every requested target bootstrap replicate must occur exactly once, be
#' successful and finite, and lie strictly between zero and one. The function
#' neither drops unavailable replicates nor applies boundary clipping or
#' continuity corrections. Failures of unrelated bootstrap metrics do not
#' invalidate an otherwise complete discrimination target.
#'
#' At least two cohorts are required because one cohort cannot identify
#' between-cohort heterogeneity. Estimates of heterogeneity and prediction
#' intervals remain especially uncertain with few cohorts, and frequentist
#' coverage can be poor in some small or imbalanced meta-analyses. The returned
#' caution is therefore part of the result rather than a model-selection rule.
#'
#' Matching `pipeline_id` and learner labels are necessary bookkeeping only.
#' They do not prove that cohorts are independent, use common eligibility and
#' outcome definitions, or received the identical fitted state, coefficient
#' formula, preprocessing, and predictor orientation. Those design conditions
#' must be established by the caller. The bootstrap uncertainty is conditional
#' on each locked model and prediction set and does not recover development-
#' fitting, feature-selection, tuning, or model-selection uncertainty.
#'
#' @references
#' Snell KIE, Ensor J, Debray TPA, Moons KGM, Riley RD (2018). Meta-analysis of
#'   prediction model performance across multiple studies: Which scale helps
#'   ensure between-study normality for the C-statistic and calibration
#'   measures? *Statistical Methods in Medical Research*, 27(11), 3505-3522.
#'   \doi{10.1177/0962280217705678}.
#'
#' Debray TPA, Damen JAAG, Snell KIE, et al. (2017). A guide to systematic
#'   review and meta-analysis of prediction model performance. *BMJ*, 356,
#'   i6460. \doi{10.1136/bmj.i6460}.
#'
#' Debray TPA, Damen JAAG, Snell KIE, et al. (2019). A framework for
#'   meta-analysis of prediction model studies with binary and time-to-event
#'   outcomes. *Statistical Methods in Medical Research*, 28(9), 2768-2786.
#'   \doi{10.1177/0962280218785504}.
#'
#' Partlett C, Riley RD (2017). Random effects meta-analysis: Coverage
#'   performance of 95% confidence and prediction intervals following REML
#'   estimation. *Statistics in Medicine*, 36(2), 301-317.
#'   \doi{10.1002/sim.7140}.
#'
#' @param bootstraps A named list of at least two
#'   `stormbreaker_external_bootstrap` objects. Names must be unique, non-blank
#'   cohort IDs. All objects must describe the same task, pipeline ID, and
#'   learner label.
#' @param confidence_level Confidence level strictly between zero and one for
#'   the pooled confidence interval and between-cohort prediction interval.
#'   Defaults to 0.95.
#' @param metric Discrimination estimand. `NULL` selects `"roc_auc"` for a
#'   binary task and `"harrell_c"` for a survival task. Binary tasks otherwise
#'   accept only `"roc_auc"`; survival tasks accept `"harrell_c"`,
#'   `"time_auc"`, or `"uno_c"`. Time-specific choices require one common
#'   retained time unit and prediction horizon.
#'
#' @return A `stormbreaker_external_discrimination_meta` containing cohort
#'   point estimates and bootstrap standard errors, the random-effects summary,
#'   explicit `metric_target`, method labels, software provenance,
#'   prediction-interval caution, and the unchanged named source `bootstraps`.
#'   No p value, ranking, or winner is returned.
#' @export
#'
#' @examples
#' \dontrun{
#' pooled <- sb_meta_external_discrimination(list(
#'   hospital_a = bootstrap_a,
#'   hospital_b = bootstrap_b
#' ))
#' pooled$summary
#' sb_meta_external_discrimination(
#'   list(hospital_a = bootstrap_a, hospital_b = bootstrap_b),
#'   metric = "time_auc"
#' )
#' }
sb_meta_external_discrimination <- function(
    bootstraps,
    confidence_level = 0.95,
    metric = NULL) {
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
    sb_external_discrimination_meta_member(
      bootstraps[[index]],
      cohort_ids[[index]],
      metric
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
  if (length(unique(pipeline_ids)) != 1L || length(unique(learners)) != 1L) {
    sb_abort(
      paste0(
        "All cohorts must retain the same pipeline ID and learner label. ",
        "These labels do not prove identical fitted state or formula."
      ),
      "stormbreaker_external_meta_model_mismatch"
    )
  }

  metric <- members[[1L]]$metric
  has_horizon <- metric %in% c("time_auc", "uno_c")
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
          "to meta-analyse `", metric, "`."
        ),
        "stormbreaker_external_meta_target_mismatch"
      )
    }
    metric_target$prediction_horizon <- horizons[[1L]]
    metric_target$time_unit <- time_units[[1L]]
  }
  estimand <- switch(
    metric,
    roc_auc = paste0(
      "Back-transformed random-effects mean of logit ROC AUC for one ",
      "locked model across an assumed distribution of independent ",
      "exchangeable cohorts"
    ),
    harrell_c = paste0(
      "Back-transformed random-effects mean of logit horizon-free Harrell ",
      "concordance for one locked model across an assumed distribution of ",
      "independent exchangeable cohorts"
    ),
    time_auc = paste0(
      "Back-transformed random-effects mean of logit cumulative/dynamic ",
      "time-dependent AUC at ", format(metric_target$prediction_horizon),
      " ", metric_target$time_unit, " across an assumed distribution of ",
      "independent exchangeable cohorts"
    ),
    uno_c = paste0(
      "Back-transformed random-effects mean of logit truncated Uno ",
      "concordance at ",
      format(metric_target$prediction_horizon), " ", metric_target$time_unit,
      " across an assumed distribution of independent exchangeable cohorts"
    )
  )

  if (!requireNamespace("metafor", quietly = TRUE)) {
    sb_abort(
      paste0(
        "Package `metafor` is required for external discrimination ",
        "meta-analysis. Install it with `install.packages(\"metafor\")`."
      ),
      "stormbreaker_missing_optional_dependency"
    )
  }

  cohorts <- data.frame(
    cohort_id = cohort_ids,
    n_patients = vapply(members, `[[`, integer(1L), "n"),
    point_estimate = vapply(members, `[[`, numeric(1L), "point_estimate"),
    logit_estimate = vapply(members, `[[`, numeric(1L), "logit_estimate"),
    logit_se = vapply(members, `[[`, numeric(1L), "logit_se"),
    n_bootstrap = vapply(members, `[[`, integer(1L), "reps"),
    stringsAsFactors = FALSE
  )
  fit <- metafor::rma.uni(
    yi = cohorts$logit_estimate,
    sei = cohorts$logit_se,
    method = "REML",
    test = "knha",
    level = 100 * confidence_level,
    slab = cohort_ids
  )
  prediction <- stats::predict(fit)
  k <- as.integer(nrow(cohorts))
  summary <- data.frame(
    metric = metric,
    k = k,
    estimate = stats::plogis(as.numeric(fit$b)),
    lower = stats::plogis(as.numeric(fit$ci.lb)),
    upper = stats::plogis(as.numeric(fit$ci.ub)),
    prediction_lower = stats::plogis(as.numeric(prediction$pi.lb)),
    prediction_upper = stats::plogis(as.numeric(prediction$pi.ub)),
    tau2_logit = as.numeric(fit$tau2),
    tau_logit = sqrt(as.numeric(fit$tau2)),
    i2_percent = as.numeric(fit$I2),
    stringsAsFactors = FALSE
  )
  caution <- if (k < 10L) {
    paste0(
      "Only ", k, " cohorts were supplied; heterogeneity and the ",
      "between-cohort prediction interval are especially uncertain."
    )
  } else {
    paste0(
      "The prediction interval assumes independent, exchangeable cohorts and ",
      "a correctly specified random-effects model."
    )
  }

  structure(
    list(
      status = "success",
      task_type = task_types[[1L]],
      pipeline_id = pipeline_ids[[1L]],
      learner = learners[[1L]],
      metric = metric,
      metric_target = metric_target,
      confidence_level = as.numeric(confidence_level),
      estimand = estimand,
      method = list(
        model = "random_effects",
        effect_scale = "logit_discrimination_probability",
        within_cohort_standard_error =
          "sample_sd_of_logit_patient_bootstrap_replicates",
        tau2_estimator = "REML",
        inference = "Knapp-Hartung",
        back_transformation = "inverse_logit"
      ),
      provenance = sb_external_meta_provenance(),
      cohorts = cohorts,
      summary = summary,
      prediction_interval_caution = caution,
      bootstraps = bootstraps
    ),
    class = "stormbreaker_external_discrimination_meta"
  )
}

#' @export
print.stormbreaker_external_discrimination_meta <- function(x, ...) {
  confidence <- 100 * x$confidence_level
  summary <- x$summary[1L, , drop = FALSE]
  cat("<stormbreaker_external_discrimination_meta:", x$status, ">\n")
  cat("fixed-model external discrimination meta-analysis\n")
  cat(
    "Metric:", x$metric,
    "| Cohorts:", summary$k,
    "| Task:", x$task_type, "\n"
  )
  if (isTRUE(x$metric_target$has_horizon)) {
    cat(
      "Target:", format(x$metric_target$prediction_horizon),
      x$metric_target$time_unit, "\n"
    )
  }
  cat(
    "Random effects: REML on logit discrimination | ",
    "Knapp-Hartung inference\n",
    sep = ""
  )
  cat(
    "Pooled:", format(summary$estimate),
    paste0("| ", format(confidence), "% CI: ["),
    format(summary$lower), ",", format(summary$upper), "]\n"
  )
  cat(
    paste0(format(confidence), "% prediction interval: ["),
    format(summary$prediction_lower), ",",
    format(summary$prediction_upper), "]\n"
  )
  cat(
    "Tau (logit discrimination):", format(summary$tau_logit),
    "| I-squared:", paste0(format(summary$i2_percent), "%\n")
  )
  cat(x$prediction_interval_caution, "\n")
  cat(
    paste0(
      "Cohort independence and common locked-model identity are caller ",
      "responsibilities; no model ranking or selection is performed.\n"
    )
  )
  invisible(x)
}
