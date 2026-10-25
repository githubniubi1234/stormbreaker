sb_external_discrimination_difference_meta_member <- function(
    comparison,
    cohort_id,
    candidate_pipeline_id,
    metric = NULL) {
  sb_external_paired_difference_meta_member(
    comparison = comparison,
    cohort_id = cohort_id,
    candidate_pipeline_id = candidate_pipeline_id,
    metric = metric,
    metric_family = "discrimination"
  )
}

#' Meta-analyse paired external discrimination differences across cohorts
#'
#' `sb_meta_external_discrimination_difference()` synthesizes one direct,
#' paired comparison between a predeclared candidate and reference prediction
#' across independent external cohorts. Each source must come from
#' [sb_bootstrap_external_comparison()], which applies the same patient draw to
#' both models. The cohort effect is the raw discrimination difference
#' `candidate - reference`; its within-cohort standard error is the sample
#' standard deviation of all corresponding paired bootstrap differences.
#'
#' Binary tasks use ROC AUC. Survival tasks default to horizon-free Harrell
#' concordance and may explicitly request cumulative/dynamic `time_auc` or
#' truncated `uno_c`. Time-specific metrics require every source comparison to
#' retain the exact same prediction horizon and time unit. The function never
#' mixes discrimination definitions or chooses a metric automatically beyond
#' the documented task-specific default.
#'
#' Differences are meta-analysed on their raw scale with REML and Knapp-Hartung
#' inference through `metafor::rma.uni()`. Point estimates and every requested
#' target-metric replicate must be present exactly once, successful, and
#' finite. Each saved difference must exactly equal its candidate estimate
#' minus its reference estimate, and the synthesis uses that recomputed value.
#' Failed replicates are not discarded, differences are not clipped, and
#' confidence or prediction limits outside the mathematical interval `[-1, 1]`
#' remain visible as a warning about the normal random-effects model.
#'
#' A positive difference indicates higher candidate discrimination only. It is
#' not evidence of better calibration, prediction error, clinical utility, or
#' overall usefulness. The function does not report a p value, rank models,
#' select a winner, make indirect or network comparisons, or use external
#' outcomes to choose the reference or candidate. Matching pipeline labels do
#' not prove identical fitted states across cohorts, and cohort independence
#' remains the study team's responsibility.
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
#'   discrimination metric. Binary tasks allow only `"roc_auc"`; survival tasks
#'   allow `"harrell_c"`, `"time_auc"`, or `"uno_c"`.
#'
#' @return A `stormbreaker_external_discrimination_difference_meta` containing
#'   the reference and candidate IDs, explicit metric target, cohort point
#'   differences and paired-bootstrap standard errors, random-effects summary,
#'   method labels, software provenance, cautions, and the unchanged named
#'   source `comparisons`. No p value, ranking, winner, or model is returned.
#' @export
#'
#' @references
#' Haile SR, Guerra B, Soriano JB, Puhan MA, 3CIA collaboration (2017).
#' Multiple Score Comparison: a network meta-analysis approach to comparison
#' and external validation of prognostic scores. *BMC Medical Research
#' Methodology*, 17, 172. \doi{10.1186/s12874-017-0433-2}
#'
#' Debray TPA, Damen JAAG, Snell KIE, et al. (2017). A guide to systematic
#' review and meta-analysis of prediction model performance. *BMJ*, 356,
#' i6460. \doi{10.1136/bmj.i6460}
#'
#' @examples
#' \dontrun{
#' pooled_difference <- sb_meta_external_discrimination_difference(
#'   list(hospital_a = paired_a, hospital_b = paired_b),
#'   candidate_pipeline_id = "candidate_model"
#' )
#' pooled_difference$summary
#' }
sb_meta_external_discrimination_difference <- function(
    comparisons,
    candidate_pipeline_id,
    confidence_level = 0.95,
    metric = NULL) {
  core <- sb_external_paired_difference_meta_core(
    comparisons = comparisons,
    candidate_pipeline_id = candidate_pipeline_id,
    confidence_level = confidence_level,
    metric = metric,
    metric_family = "discrimination"
  )
  metric <- core$metric
  metric_target <- core$metric_target
  estimand <- switch(
    metric,
    roc_auc = paste0(
      "Random-effects mean paired candidate-minus-reference ROC AUC ",
      "difference across an assumed distribution of independent ",
      "exchangeable cohorts"
    ),
    harrell_c = paste0(
      "Random-effects mean paired candidate-minus-reference horizon-free ",
      "Harrell concordance difference across an assumed distribution of ",
      "independent exchangeable cohorts"
    ),
    time_auc = paste0(
      "Random-effects mean paired candidate-minus-reference cumulative/",
      "dynamic time-dependent AUC difference at ",
      format(metric_target$prediction_horizon), " ", metric_target$time_unit,
      " across an assumed distribution of independent exchangeable cohorts"
    ),
    uno_c = paste0(
      "Random-effects mean paired candidate-minus-reference truncated Uno ",
      "concordance difference at ",
      format(metric_target$prediction_horizon), " ", metric_target$time_unit,
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
    "Raw-scale normal-model intervals are not clipped to [-1, 1]. A positive ",
    "difference means higher candidate discrimination only, not better ",
    "calibration, clinical utility, or overall performance."
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
        positive = "higher_candidate_discrimination",
        negative = "higher_reference_discrimination"
      ),
      method = list(
        model = "random_effects",
        effect_scale = "raw_candidate_minus_reference_discrimination",
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
    class = "stormbreaker_external_discrimination_difference_meta"
  )
}

#' @export
print.stormbreaker_external_discrimination_difference_meta <- function(
    x,
    ...) {
  confidence <- 100 * x$confidence_level
  summary <- x$summary[1L, , drop = FALSE]
  cat(
    "<stormbreaker_external_discrimination_difference_meta:",
    x$status,
    ">\n"
  )
  cat("paired external discrimination difference meta-analysis\n")
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
