#' Continue a chosen pipeline without repeating model comparison
#' @param run Result from [bi_run()] or [surv_run()].
#' @param pipeline Exact name from names(run$pipelines), chosen using development
#'   evidence, not external outcomes. No winner is selected automatically.
#' @param external Named list of independent cohort data frames. Use the same
#'   ID, predictor and endpoint column names as development data.
#' @param plots Include available algorithm and external ROC plots.
#' @param explain Optional external cohort name to explain. Development data
#'   supply the SHAP background; external outcomes never enter explanation.
#' @param reverse For binary outcomes, `TRUE` (default) uses cohort outcomes to
#'   flip probabilities when AUC < 0.5; `FALSE` keeps direction, and survival
#'   outcomes require `FALSE` (used when omitted).
#' @return A common bundle with comparison, fit, genes, predictions, metrics,
#'   evaluations, plots and issues. Predictions contain patient identifiers.
#'   `predictions` retains original model outputs; `evaluated_predictions`
#'   matches the reported metrics and contains original probabilities and flip
#'   flags for processed evaluations. `direction_correction` records each cohort.
#'   External failures are reported individually; they never trigger refitting.
#'   Survival risk_groups records counts at the locked development median even
#'   when plots is FALSE. Single groups are reported in issues, not repaired by
#'   choosing a new external threshold; prediction metrics remain available.
#'   With survival plots, risk_groups also retains logrank_p, logrank_status
#'   and logrank_error. Failed log-rank tests warn and appear in issues while
#'   usable plots remain available. Without plots the test is not computed.
#'   Printing shows a short summary without patient rows. Saved algorithm plots
#'   are reused unless a new explanation is requested.
#' @examples
#' development <- bi_example(40, seed = 42)
#' validation <- bi_example(30, seed = 17)
#' run <- bi_run(development, outcome = "outcome", id = "patient_id",
#'   predictors = c("gene_a", "gene_b"), models = "Logistic", folds = 2,
#'   seed = 42, plots = FALSE)
#' finished <- finish_run(run, "Logistic",
#'   external = list(Validation = validation), plots = FALSE)
#' result_table(finished$evaluations$Validation, metric = "roc_auc")
#' fit <- finished$fit
#' prediction <- predict_final(fit, validation[c("patient_id", "gene_a", "gene_b")],
#'   "patient_id")
#' head(prediction$predictions)
#' @export
sb_complete_run <- function(run, pipeline, external = list(), plots = TRUE, explain = NULL, reverse = TRUE) {
  if (!is.logical(reverse) || length(reverse) != 1L || is.na(reverse)) {
    sb_abort("`reverse` must be TRUE or FALSE.", "stormbreaker_bad_reverse")
  }
  if (missing(reverse) && identical(run$task_type, "survival")) reverse <- FALSE
  if (!is.character(pipeline) || length(pipeline) != 1L || is.na(pipeline) ||
    !pipeline %in% names(run$pipelines)) {
    stop("Choose one exact name from names(run$pipelines).")
  }
  task <- run$development_task
  if (is.null(task)) stop("This older run does not retain development data; use a current bi_run/surv_run result.")
  if (!is.list(external) || (length(external) && (is.null(names(external)) ||
    anyNA(names(external)) || any(!nzchar(names(external))) || anyDuplicated(names(external))))) {
    stop("Supply external cohorts as a uniquely named list of data frames.")
  }
  if (!is.null(explain) && (length(explain) != 1L || is.na(explain) || !explain %in% names(external))) {
    stop("explain must name one supplied external cohort.")
  }
  fit <- run$details[[pipeline]]$fit
  if (is.null(fit)) {
    fit <- sb_fit_final(task, run$pipelines[[pipeline]],
      capture_diagnostics = isTRUE(plots) || !is.null(explain)
    )
  }
  result <- structure(
    list(
      pipeline = pipeline, task_type = run$task_type,
      comparison = run$comparison$comparison, fit = fit, genes = data.frame(),
      predictions = list(), evaluated_predictions = list(), direction_correction = list(),
      metrics = list(), evaluations = list(), plots = list(),
      explanation = NULL, risk_groups = list(), issues = list()
    ),
    class = c("stormbreaker_completed_run", "list")
  )
  if (!identical(fit$status, "success")) {
    result$issues$fit <- fit
    return(result)
  }
  state <- fit$selection$state
  stages <- if (identical(state$id, "ordered_chain")) state$stages else list(state)
  if (!length(stages) || is.null(stages[[1]]$selected_features)) {
    result$genes <- data.frame(stage = 0L, method = "Input features", gene = task$predictors)
  } else {
    result$genes <- do.call(rbind, lapply(seq_along(stages), function(i) {
      data.frame(stage = i, method = stages[[i]]$id, gene = stages[[i]]$selected_features)
    }))
  }
  if (isTRUE(plots) || !is.null(explain)) {
    if (isTRUE(plots)) result$plots$comparison <- run$plots
    figures <- if (is.null(explain)) run$details[[pipeline]]$figures else NULL
    if (is.null(figures)) {
      figures <- tryCatch(sb_plots(fit,
        new_data = if (!is.null(explain)) external[[explain]][c(task$id_col, task$predictors)] else NULL,
        id_col = task$id_col,
        background_data = if (!is.null(explain)) task$data[c(task$id_col, task$predictors)] else NULL
      ), error = identity)
    }
    if (inherits(figures, "error")) {
      result$issues$algorithm <- conditionMessage(figures)
    } else {
      if (isTRUE(plots)) result$plots$algorithm <- figures$plots
      result$explanation <- figures$shap
      result$issues$algorithm <- figures$unavailable
    }
  }
  reference <- if (run$task_type == "survival" && length(external)) {
    tryCatch(surv_risk_reference(fit, task), error = identity)
  } else {
    NULL
  }
  for (cohort in names(external)) {
    attempt <- tryCatch(
      {
        d <- external[[cohort]]
        prediction <- sb_predict_final(fit, d[c(task$id_col, task$predictors)], task$id_col)
        endpoint <- if (run$task_type == "binary") {
          sb_binary_endpoint(d, task$id_col, task$outcome_col)
        } else {
          sb_surv_endpoint(d, task$id_col, task$time_col, task$event_col, task$time_unit,
            prediction_horizon = task$prediction_horizon
          )
        }
        evaluation <- sb_evaluate_external(prediction, endpoint, reverse = reverse)
        list(prediction = prediction, evaluation = evaluation)
      },
      error = identity
    )
    if (inherits(attempt, "error")) {
      result$issues$external[[cohort]] <- conditionMessage(attempt)
      next
    }
    result$predictions[[cohort]] <- attempt$prediction$predictions
    result$evaluations[[cohort]] <- attempt$evaluation
    result$evaluated_predictions[[cohort]] <- attempt$evaluation$scored_predictions
    result$direction_correction[[cohort]] <- attempt$evaluation$direction_correction
    result$metrics[[cohort]] <- attempt$evaluation$metrics
    if (!identical(attempt$evaluation$status, "success")) {
      result$issues$external[[cohort]] <- attempt$evaluation
    } else if (isTRUE(plots) && run$task_type == "binary") {
      figure <- tryCatch(sb_plot_roc(attempt$evaluation, cohort_label = cohort), error = identity)
      if (inherits(figure, "error")) {
        result$issues$plot[[cohort]] <- conditionMessage(figure)
      } else {
        result$plots$external[[cohort]] <- sb_compact_plot(figure)
      }
    } else if (run$task_type == "survival") {
      if (inherits(reference, "error")) {
        result$issues$risk[[cohort]] <- conditionMessage(reference)
        next
      }
      scores <- attempt$prediction$predictions$risk
      n_high <- sum(scores > reference$cutoff)
      n_low <- length(scores) - n_high
      status <- if (n_high > 0 && n_low > 0) "two_groups" else "single_group"
      result$risk_groups[[cohort]] <- data.frame(
        status = status,
        cutoff = reference$cutoff, n_low = n_low, n_high = n_high,
        logrank_p = NA_real_, logrank_status = "not_computed",
        logrank_error = NA_character_
      )
      if (status == "single_group") {
        result$issues$risk[[cohort]] <- sprintf(
          "single_group: locked development cutoff %.8g gives Low=%d, High=%d; high-versus-low HR and log-rank comparison are unavailable. The cutoff was not changed.",
          reference$cutoff, n_low, n_high
        )
      }
      if (!isTRUE(plots)) next
      risk <- tryCatch(
        {
          ext_task <- sb_surv_task(
            external[[cohort]], task$id_col, task$time_col,
            task$event_col, task$predictors, task$time_unit
          )
          surv_risk_plots(reference, ext_task)
        },
        error = identity
      )
      if (inherits(risk, "error")) {
        result$issues$plot[[cohort]] <- conditionMessage(risk)
      } else {
        result$plots$external[[cohort]] <- risk$plot
        result$risk_groups[[cohort]]$logrank_p <- risk$logrank_p
        result$risk_groups[[cohort]]$logrank_status <- risk$logrank_status
        result$risk_groups[[cohort]]$logrank_error <- risk$logrank_error
        if (identical(risk$logrank_status, "failed")) {
          result$issues$risk[[cohort]] <- paste0("Log-rank test failed: ", risk$logrank_error)
        }
      }
    }
  }
  result
}
