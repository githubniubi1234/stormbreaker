#' Select a fitted model using validation cohorts only
#'
#' Select the highest equally weighted mean validation ROC AUC (binary) or
#' Harrell C (survival). The returned object is an existing training-only fit,
#' ready for [sb_predict_final()] and `saveRDS()`. No model is refitted.
#'
#' @param comparison A result from [sb_compare_cohorts()]. At least one cohort
#'   must have the explicit role `"validation"`.
#' @return A `stormbreaker_final_fit` with a `model_selection` audit containing
#'   the validation cohort names, metric, per-cohort scores and full ranking.
#' @details A candidate is eligible only when its final fit and development CV
#'   succeeded, and every validation cohort has a successful finite primary
#'   metric. Ineligible candidates remain in the ranking. Ties preserve the
#'   original pipeline declaration order. Test and external cohorts never
#'   enter selection. Keep the final test cohort out of the comparison call
#'   until the model is saved; select once and assess it on the held-out data.
#'   Validation scores used to choose a model are not unbiased estimates of
#'   the selected model's performance. The function does not choose a decision
#'   threshold, recalibrate predictions or add validation patients to training.
#'   Exactly constant predictions in any validation cohort, or in every
#'   development assessment fold, are ineligible even when AUC/C is 0.5.
#'   Variation is checked within folds, not across their different baselines.
#'   Saved metrics and fitted states are not changed. The selection audit
#'   retains prediction diagnostics and exclusion reasons. Missing prediction
#'   evidence cannot establish eligibility; test/external predictions are not
#'   inspected. This is not a near-constant tolerance or a performance cutoff.
#' @export
#' @examples
#' make_task <- function(prefix, seed) {
#'   d <- sb_example_binary(40, seed)
#'   d$patient_id <- paste0(prefix, d$patient_id)
#'   sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
#' }
#' compared <- sb_compare_cohorts(
#'   make_task("train", 1), list(Validation = make_task("val", 2)),
#'   roles = c(Validation = "validation"),
#'   models = c("logistic", "gaussian_naive_bayes"), folds = 3, view = FALSE
#' )
#' model <- sb_select_model(compared)
#' sb_predict_final(model, sb_example_binary(20, 3), "patient_id")
sb_select_model <- function(comparison) {
  if (!inherits(comparison, "stormbreaker_cohort_comparison")) {
    sb_abort("Supply the result of sb_compare_cohorts().",
      "stormbreaker_bad_cohort_comparison")
  }
  cohorts <- names(comparison$cohort_roles)[comparison$cohort_roles == "validation"]
  if (!length(cohorts)) {
    sb_abort("Model selection requires at least one cohort with role 'validation'.",
      "stormbreaker_no_validation_cohort")
  }
  metric <- if (identical(comparison$task_type, "binary")) "roc_auc" else "harrell_c"
  ids <- vapply(comparison$model_pipelines, function(x) x$id, character(1))
  fit_ids <- vapply(comparison$fits, function(x) x$pipeline$id, character(1))
  fit_index <- match(ids, fit_ids)
  cv <- comparison$development$comparison
  eligible <- !is.na(fit_index) & vapply(fit_index, function(i) {
    !is.na(i) && identical(comparison$fits[[i]]$status, "success")
  }, logical(1)) & cv$status[match(ids, cv$pipeline_id)] == "success"
  exclusion_reason <- rep("", length(ids))
  diagnostics <- list()
  check_variation <- function(results, field, cohort, scope) {
    states <- vapply(ids, function(id) {
      sb_prediction_variation(results[[id]][[field]], comparison$task_type)$status
    }, character(1))
    diagnostics[[length(diagnostics) + 1L]] <<- data.frame(
      pipeline_id = ids, cohort = cohort, prediction_variation = unname(states))
    invalid <- !states %in% c("variable", "partly_constant")
    reason <- ifelse(states == "constant", paste0("constant_", scope, "_predictions"),
      paste0("unavailable_", scope, "_prediction_evidence"))
    if (scope == "validation") reason <- paste0(reason, ":", cohort)
    add <- which(invalid)
    exclusion_reason[add] <<- ifelse(nzchar(exclusion_reason[add]),
      paste(exclusion_reason[add], reason[add], sep = ";"), reason[add])
    eligible <<- eligible & !invalid
  }
  check_variation(comparison$development$results, "predictions", "development_cv", "development")
  scores <- matrix(NA_real_, nrow = length(ids), ncol = length(cohorts),
    dimnames = list(ids, cohorts))
  for (cohort in cohorts) {
    tab <- comparison$cohorts[[cohort]]$comparison
    rows <- match(ids, tab$pipeline_id)
    scores[, cohort] <- tab[[metric]][rows]
    eligible <- eligible & tab$status[rows] == "success" &
      tab[[paste0(metric, "_status")]][rows] == "success" &
      is.finite(scores[, cohort])
    check_variation(comparison$cohorts[[cohort]]$evaluations, "scored_predictions", cohort, "validation")
  }
  eligible[is.na(eligible)] <- FALSE
  exclusion_reason[!eligible & !nzchar(exclusion_reason)] <- "incomplete_fit_or_metrics"
  means <- rowMeans(scores)
  means[!eligible] <- NA_real_
  ranking <- data.frame(pipeline_id = ids, eligible = eligible,
    validation_mean = unname(means), declared_order = seq_along(ids),
    exclusion_reason = exclusion_reason,
    stringsAsFactors = FALSE)
  if (!any(eligible)) {
    sb_abort(paste0("No eligible model has successful development CV, a final fit, ",
      "nonconstant prediction evidence and finite ", metric, " in every validation cohort. ",
      "Exclusions: ", paste(unique(exclusion_reason), collapse = "; "), "."),
      "stormbreaker_no_eligible_model")
  }
  best <- which.max(replace(means, !eligible, -Inf))
  ranking$selected <- seq_along(ids) == best
  model <- comparison$fits[[fit_index[best]]]
  model$model_selection <- list(metric = metric, cohorts = cohorts,
    rule = "maximum_equal_weight_validation_mean; ties_in_declaration_order",
    scores = scores,
    prediction_diagnostics = do.call(rbind, diagnostics),
    ranking = ranking[order(!eligible, -means, ranking$declared_order), , drop = FALSE])
  model
}
