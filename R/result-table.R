#' Extract saved cross-validation or external evaluation tables
#'
#' Read saved metrics, patient predictions, or fold-run records without fitting,
#' evaluating, ranking, or modifying the supplied result.
#'
#' @param x A single result from [sb_run_cv()], a shared-fold comparison from
#'   [sb_compare_cv()], an evaluation from [sb_evaluate_external()], or a
#'   single-cohort comparison from [sb_compare_external()]. Also accepts whole
#'   repeated CV, internal-external CV (IECV), and stacking CV results.
#'   Arbitrary cross-cohort collections are not accepted by this table function.
#' @param component One of `"metrics"`, `"predictions"`, or `"runs"`.
#' @param metric For `component = "metrics"`, one explicitly selected saved
#'   scalar metric: `"roc_auc"`, `"average_precision"`, `"brier_score"`,
#'   `"log_loss"` (binary), or `"harrell_c"`, `"time_auc"`, `"uno_c"`,
#'   `"ipcw_brier"` (survival). The last three are saved fixed-horizon
#'   cumulative/dynamic AUC, truncated Uno concordance, and IPCW Brier score.
#'   Leave this argument `NULL` for prediction and run tables. Nested reference,
#'   skill, calibration and decision-curve objects are not flattened here.
#' @param cohort_label Optional non-empty display label supplied by the caller.
#'   A cohort identity is never inferred. External tables always include this
#'   column (missing when omitted); CV tables include it only when supplied.
#'
#' @return A data frame in supplied pipeline order. Metric tables contain
#'   `pipeline_id`, `task_type`, `metric`, `estimate`, `status`, `run_status`,
#'   `scope`, and `aggregation`. `status` belongs to the selected metric;
#'   `run_status` belongs to the complete CV run or external evaluation.
#'   External metric tables also retain evaluation errors, metric errors and
#'   patient-match status. Prediction and run tables
#'   preserve saved columns and values and append `task_type`. Rows from each
#'   pipeline remain in saved order; output row names are reset.
#'   For external objects, predictions are saved scored rows including outcomes;
#'   runs contain one evaluation summary per pipeline, not fictitious CV folds.
#'   The three fixed-horizon metrics additionally retain `prediction_horizon`,
#'   `time_unit`, and saved metric errors. Other metric tables keep their
#'   existing columns. Different horizons or units are distinct targets; no
#'   conversion or pooling is performed.
#'
#' @details Composite objects add `repeat_id`, `heldout_cohort_id`, `model_role`,
#'   `validation_scheme`, and `view_id` (a structural entry key). These identify
#'   saved repetitions, cohorts, and member/stack outputs without changing
#'   `pipeline_id`. Unused context fields remain missing. Rows stay in stored
#'   order; no averaging, pooling, or confidence interval is inferred. IECV is
#'   method development, not final independent external validation. Standalone
#'   stacking includes outer-CV members and the stack, not inner training folds.
#'   When a stack descriptor omits its horizon or unit, an unambiguous value
#'   saved in its prediction rows is retained. Conflicting or absent row
#'   metadata is not guessed; metric-level metadata still takes precedence.
#'   Failed methods and missing or infinite estimates remain in the
#'   table. A metric can be evaluable even when the complete run failed at
#'   another requested output. No estimate is filled, clipped, recalibrated, or
#'   recomputed. Binary ROC AUC, average precision, and survival Harrell C keep
#'   their saved fold-aware aggregation; they are not replaced by a pooled
#'   calculation. Confidence intervals and a winning model are not inferred.
#'   External metrics retain their saved external-cohort scope and aggregation.
#'   A failed external evaluation with no saved supported metric has an `NA`
#'   estimate and metric status `"not_saved"`; its failure information remains.
#'   An unsupported metric or a missing metric in a nonfailed result is an error.
#'   For the three fixed-horizon metrics, an absent value in a successful
#'   horizon-free result is instead `"not_requested"`; an absent requested
#'   value or an absent value in a failed result is `"not_saved"`. Its estimate
#'   remains `NA`. A horizon stored in the metric takes precedence over the
#'   task/endpoint descriptor; unavailable metadata remains missing.
#'   Failed evaluations with no scored rows contribute no prediction rows;
#'   inspect the run table for their failures. Cohort labels are descriptive,
#'   not evidence of independence. No cohorts are pooled by this function.
#'
#'   Run records retain their existing error information. An error stage is not
#'   invented when the source did not save one. Prediction tables can contain
#'   patient identifiers and outcomes; review them before sharing or exporting.
#'
#' @seealso [sb_plot_performance()]
#' @export
#'
#' @examples
#' d <- sb_example_binary(n = 60, seed = 1)
#' task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
#' result <- sb_run_cv(
#'   task, sb_pipeline(learner = "logistic"), folds = 3, seed = 1
#' )
#' sb_result_table(result, metric = "roc_auc")
#' head(sb_result_table(result, "predictions"))
#' sb_result_table(result, "runs")
sb_result_table <- function(x,
                            component = c("metrics", "predictions", "runs"),
                            metric = NULL, cohort_label = NULL) {
  component <- match.arg(component)
  label <- sb_result_cohort_label(cohort_label)
  if (sb_is_composite_view(x) || !is.null(attr(x, "stormbreaker_view_context"))) {
    return(sb_composite_table(x, component, metric, cohort_label))
  }
  if (inherits(x, c("stormbreaker_external_evaluation",
                    "stormbreaker_external_comparison"))) {
    return(sb_external_result_table(x, component, metric, cohort_label))
  }
  is_cv <- function(value) {
    inherits(value, c("stormbreaker_binary_cv_result", "stormbreaker_cv_result"))
  }
  if (inherits(x, "stormbreaker_cv_comparison")) {
    results <- x$results
  } else if (is_cv(x)) {
    results <- list(x)
  } else {
    sb_abort(
      "`x` must be a CV result/comparison or external evaluation/comparison.",
      "stormbreaker_bad_result_view"
    )
  }
  if (!is.list(results) || length(results) == 0L ||
      !all(vapply(results, is_cv, logical(1)))) {
    sb_abort(
      "`x` must contain at least one supported CV result.",
      "stormbreaker_bad_result_view"
    )
  }
  if (identical(component, "metrics")) {
    supported <- sb_scalar_view_metrics()
    if (!is.character(metric) || length(metric) != 1L ||
        is.na(metric) || !metric %in% supported) {
      sb_abort(
        paste0(
          "Unsupported `metric`; explicitly choose one of: ",
          paste(supported, collapse = ", "), "."
        ),
        "stormbreaker_bad_result_metric"
      )
    }
  } else if (!is.null(metric)) {
    sb_abort(
      "`metric` is only used for `component = \"metrics\"`; otherwise leave it NULL.",
      "stormbreaker_unused_result_metric"
    )
  }

  rows <- lapply(results, function(result) {
    task_type <- if (inherits(result, "stormbreaker_binary_cv_result")) {
      "binary"
    } else {
      "survival"
    }
    if (!identical(component, "metrics")) {
      tab <- result[[component]]
      tab$task_type <- rep(task_type, nrow(tab))
      return(tab)
    }
    if (metric %in% sb_horizon_view_metrics()) {
      return(sb_horizon_metric_row(result, metric, task_type))
    }
    value <- result$metrics[[metric]]
    if (is.null(value)) {
      sb_abort(
        paste0(
          "Metric `", metric, "` is not saved in pipeline `",
          result$pipeline$id, "`. No metric is recomputed by this function."
        ),
        "stormbreaker_result_metric_not_saved"
      )
    }
    data.frame(
      pipeline_id = result$pipeline$id,
      task_type = task_type,
      metric = metric,
      estimate = value$estimate,
      status = value$status,
      run_status = result$status,
      scope = value$scope,
      aggregation = value$aggregation,
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  if (!is.null(cohort_label)) output$cohort_label <- rep(label, nrow(output))
  output
}
