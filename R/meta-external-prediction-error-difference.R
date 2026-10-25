#' Meta-analyse paired external prediction-error differences across cohorts
#'
#' `sb_meta_external_prediction_error_difference()` synthesizes one direct,
#' paired comparison between a predeclared candidate and reference prediction
#' across independent external cohorts. Each source must come from
#' [sb_bootstrap_external_comparison()], which applies the same patient draw to
#' both models. The cohort effect is the raw prediction-error difference
#' `candidate - reference`; its within-cohort standard error is the sample
#' standard deviation of all corresponding paired bootstrap differences.
#'
#' Binary tasks default to Brier score and may explicitly request log loss.
#' Survival tasks use the fixed-horizon IPCW Brier score and therefore require
#' every cohort to retain the exact same prediction horizon and time unit.
#' Lower values are better for every supported metric, so a negative difference
#' means lower prediction error for the candidate.
#'
#' Differences are meta-analysed on their raw scale with REML and Knapp-Hartung
#' inference through `metafor::rma.uni()`. Point estimates and every requested
#' target-metric replicate must be present exactly once, successful, and
#' finite. Each saved difference must exactly equal its candidate estimate
#' minus its reference estimate, and the synthesis uses that recomputed value.
#' Failed replicates are not discarded and confidence or prediction limits are
#' not clipped.
#'
#' A lower prediction error is not by itself evidence of clinical utility, nor
#' does it identify whether a difference arose from calibration,
#' discrimination, or both. The function does not report a p value, rank
#' models, select a winner, make indirect or network comparisons, or use
#' external outcomes to choose the reference, candidate, or metric. Prediction
#' error depends on cohort outcome incidence and case mix; the survival metric
#' additionally relies on the package's stated marginal censoring assumptions.
#' Matching pipeline labels do not prove identical fitted states across
#' cohorts, and cohort independence remains the study team's responsibility.
#'
#' @param comparisons A named list of at least two
#'   `stormbreaker_external_bootstrap_comparison` objects returned by
#'   [sb_bootstrap_external_comparison()]. Outer names are treated as cohort
#'   IDs and must be unique and non-blank.
#' @param candidate_pipeline_id Exact pipeline ID of one candidate retained in
#'   every source comparison. The reference and candidate must have been
#'   declared independently of external outcomes.
#' @param confidence_level Confidence level strictly between 0 and 1. Defaults
#'   to 0.95.
#' @param metric `NULL` for the task-specific default, or one explicit
#'   prediction-error metric. Binary tasks allow `"brier_score"` or
#'   `"log_loss"`; survival tasks allow only `"ipcw_brier_score"`.
#'
#' @return A `stormbreaker_external_prediction_error_difference_meta`
#'   containing the reference and candidate IDs, explicit metric target,
#'   cohort point differences and paired-bootstrap standard errors,
#'   random-effects summary, method labels, software provenance, cautions, and
#'   the unchanged named source `comparisons`. No p value, ranking, winner, or
#'   model is returned.
#' @export
#'
#' @references
#' Gneiting T, Raftery AE (2007). Strictly proper scoring rules, prediction,
#' and estimation. *Journal of the American Statistical Association*,
#' 102(477), 359-378. \doi{10.1198/016214506000001437}
#'
#' Debray TPA, Damen JAAG, Snell KIE, et al. (2017). A guide to systematic
#' review and meta-analysis of prediction model performance. *BMJ*, 356,
#' i6460. \doi{10.1136/bmj.i6460}
#'
#' McLernon DJ, Giardiello D, Van Calster B, et al. (2023). Assessing
#' performance and clinical usefulness in prediction models with survival
#' outcomes: practical guidance for Cox proportional hazards models.
#' *Annals of Internal Medicine*, 176(1), 105-114.
#' \doi{10.7326/M22-0844}
#'
#' @examples
#' \dontrun{
#' pooled_difference <- sb_meta_external_prediction_error_difference(
#'   list(hospital_a = paired_a, hospital_b = paired_b),
#'   candidate_pipeline_id = "candidate_model"
#' )
#' pooled_difference$summary
#' }
sb_meta_external_prediction_error_difference <- function(
    comparisons,
    candidate_pipeline_id,
    confidence_level = 0.95,
    metric = NULL) {
  core <- sb_external_paired_difference_meta_core(
    comparisons = comparisons,
    candidate_pipeline_id = candidate_pipeline_id,
    confidence_level = confidence_level,
    metric = metric,
    metric_family = "prediction_error"
  )
  metric <- core$metric
  metric_target <- core$metric_target
  estimand <- switch(
    metric,
    brier_score = paste0(
      "Random-effects mean paired candidate-minus-reference binary Brier ",
      "score difference across an assumed distribution of independent ",
      "exchangeable cohorts"
    ),
    log_loss = paste0(
      "Random-effects mean paired candidate-minus-reference binary log-loss ",
      "difference across an assumed distribution of independent ",
      "exchangeable cohorts"
    ),
    ipcw_brier_score = paste0(
      "Random-effects mean paired candidate-minus-reference IPCW Brier ",
      "score difference at ", format(metric_target$prediction_horizon), " ",
      metric_target$time_unit,
      " across an assumed distribution of independent exchangeable cohorts"
    )
  )
  k <- core$summary$k[[1L]]
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
    "Raw-scale normal-model intervals are not clipped. A negative difference ",
    "means lower candidate prediction error only, not established clinical ",
    "utility or uniformly better calibration and discrimination."
  )

  structure(
    list(
      status = "success",
      task_type = core$task_type,
      reference_pipeline_id = core$reference_pipeline_id,
      candidate_pipeline_id = core$candidate_pipeline_id,
      contrast = "candidate_minus_reference",
      metric = metric,
      metric_target = metric_target,
      confidence_level = core$confidence_level,
      estimand = estimand,
      interpretation = list(
        null = 0,
        negative = "lower_candidate_prediction_error",
        positive = "lower_reference_prediction_error"
      ),
      method = list(
        model = "random_effects",
        effect_scale = "raw_candidate_minus_reference_prediction_error",
        within_cohort_standard_error =
          "sample_sd_of_paired_patient_bootstrap_differences",
        tau2_estimator = "REML",
        inference = "Knapp-Hartung",
        back_transformation = "none"
      ),
      provenance = sb_external_meta_provenance(),
      cohorts = core$cohorts,
      summary = core$summary,
      prediction_interval_caution = caution,
      comparisons = core$comparisons
    ),
    class = "stormbreaker_external_prediction_error_difference_meta"
  )
}

#' @export
print.stormbreaker_external_prediction_error_difference_meta <- function(
    x,
    ...) {
  confidence <- 100 * x$confidence_level
  summary <- x$summary[1L, , drop = FALSE]
  cat(
    "<stormbreaker_external_prediction_error_difference_meta:",
    x$status,
    ">\n"
  )
  cat("paired external prediction error difference meta-analysis\n")
  cat(
    "Contrast:", x$candidate_pipeline_id, "-", x$reference_pipeline_id,
    "(candidate - reference)\n"
  )
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
    sprintf(
      "Mean difference: %.3f | %.1f%% CI: %.3f to %.3f\n",
      summary$estimate,
      confidence,
      summary$lower,
      summary$upper
    )
  )
  cat(
    sprintf(
      "Between-cohort %.1f%% prediction interval: %.3f to %.3f\n",
      confidence,
      summary$prediction_lower,
      summary$prediction_upper
    )
  )
  cat(
    sprintf(
      "tau^2: %.4f | tau: %.4f | I^2: %.1f%%\n",
      summary$tau2,
      summary$tau,
      summary$i2_percent
    )
  )
  cat(x$prediction_interval_caution, "\n")
  cat(
    paste0(
      "Source comparisons are retained; no ranking or model selection is ",
      "performed. Model identity and cohort independence remain caller ",
      "responsibilities.\n"
    )
  )
  invisible(x)
}
