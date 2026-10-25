#' Compare development models across named validation and test cohorts
#'
#' Run shared-fold development CV, fit each declared pipeline on development
#' patients, and apply those same fitted states to every supplied cohort.
#' The result opens as a model-by-cohort matrix in an interactive R session.
#'
#' @param development A task from [sb_binary_task()] or [sb_surv_task()].
#' @param cohorts A non-empty named list of same-task validation, test or
#'   external cohort tasks. Patient/group IDs must be distinct across cohorts.
#' @param roles A named character vector matching `cohorts`, containing
#'   `"validation"`, `"test"` or `"external"`. NULL labels all as external.
#' @param pipelines Explicit [sb_pipeline()] specifications. If supplied,
#'   leave `models` NULL.
#' @param models Model IDs for [sb_model_pipelines()]; NULL uses its full menu
#'   when `pipelines` is also NULL.
#' @param folds,seed Shared development CV settings.
#' @param inner_folds,inner_seed Inner CV settings for glmnet menu models.
#' @param development_name Display name for development CV.
#' @param view Open the results page after computation. Defaults to FALSE;
#'   set TRUE to explicitly request a browser.
#' @param progress Show model, fold, final-fit and evaluation progress.
#'   Defaults to `interactive()`; use TRUE for script logs.
#' @param report_file Optional HTML destination. Saves the page even when
#'   `view = FALSE`. A new path is required.
#' @return A `stormbreaker_cohort_comparison` with development comparison,
#'   final fits, cohort predictions/evaluations, named result sets, roles,
#'   and an optional `report_file`.
#' @details Cohort outcomes are used only for evaluation, after all models
#'   have been fitted on development data. No winner is chosen and validation
#'   or test patients are not added to the fitting data. To assess a model
#'   selected using validation results, keep a separate test cohort.
#' @export
#' @examples
#' make_task <- function(prefix, seed, n) {
#'   d <- sb_example_binary(n, seed)
#'   d$patient_id <- paste0(prefix, d$patient_id)
#'   sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
#' }
#' result <- sb_compare_cohorts(
#'   make_task("dev_", 1, 50),
#'   list(Validation = make_task("val_", 2, 30),
#'        Test = make_task("test_", 3, 30)),
#'   roles = c(Validation = "validation", Test = "test"),
#'   models = c("logistic", "gaussian_naive_bayes"),
#'   folds = 3, view = FALSE
#' )
#' result$cohorts$Test$comparison
sb_compare_cohorts <- function(development, cohorts, roles = NULL,
                               pipelines = NULL, models = NULL,
                               folds = 5L, seed = 1L, inner_folds = 5L,
                               inner_seed = 1701L, development_name = "Development",
                               view = FALSE, report_file = NULL,
                               progress = interactive()) {
  sb_check_progress(progress)
  sb_check_page_request(view, report_file)
  if (!is.character(development_name) || length(development_name) != 1L ||
      is.na(development_name) || !nzchar(trimws(development_name))) {
    sb_abort("Supply one non-empty development_name.", "stormbreaker_bad_cohort_roles")
  }
  if (!is.list(cohorts) || is.object(cohorts) || !length(cohorts) ||
      is.null(names(cohorts)) || development_name %in% names(cohorts)) {
    sb_abort("Supply named cohorts distinct from development_name.",
      "stormbreaker_bad_cohort_roles")
  }
  if (is.null(roles)) roles <- stats::setNames(rep("external", length(cohorts)), names(cohorts))
  if (!is.character(roles) || anyNA(roles) || is.null(names(roles)) ||
      anyDuplicated(names(roles)) || !setequal(names(roles), names(cohorts)) ||
      !all(roles %in% c("validation", "test", "external"))) {
    sb_abort("roles must name every cohort as validation, test or external.",
      "stormbreaker_bad_cohort_roles")
  }
  roles <- roles[names(cohorts)]
  validated <- sb_validate_iecv_cohorts(c(stats::setNames(list(development), development_name), cohorts))
  development <- validated$tasks[[1L]]
  cohorts <- validated$tasks[-1L]
  if (!is.null(pipelines) && !is.null(models)) {
    sb_abort("Supply pipelines or models, not both.", "stormbreaker_bad_model_choice")
  }
  if (is.null(pipelines)) {
    pipelines <- sb_model_pipelines(validated$task_type, models, inner_folds, inner_seed)
  }
  pipelines <- sb_validate_iecv_pipelines(pipelines, validated$task_type)
  started <- proc.time()[["elapsed"]]
  cv <- sb_compare_cv(development, pipelines, folds = folds, seed = seed, progress = progress)
  fits <- sb_with_selector_reuse(lapply(seq_along(pipelines), function(i) {
    pipeline <- pipelines[[i]]
    sb_progress_message(progress, "Final fit ", i, "/", length(pipelines), ": ",
      sb_display_selector(pipeline$selector), " + ",
      sb_display_learner(pipeline$learner, pipeline$learner_params))
    fit <- sb_fit_final(development, pipeline)
    sb_progress_message(progress, "Final fit ", i, "/", length(pipelines), " complete: ", fit$status)
    fit
  }))
  names(fits) <- names(pipelines)
  predictions <- lapply(seq_along(cohorts), function(i) {
    cohort <- cohorts[[i]]
    sb_progress_message(progress, "Predict ", names(cohorts)[i])
    new_data <- cohort$data[c(cohort$id_col, development$predictors)]
    lapply(fits, function(fit) sb_predict_final(fit, new_data, cohort$id_col))
  })
  names(predictions) <- names(cohorts)
  evaluations <- lapply(seq_along(cohorts), function(i) {
    sb_progress_message(progress, "Evaluate ", names(cohorts)[i])
    endpoint <- sb_iecv_endpoint(cohorts[[i]], validated$task_type)
    sb_compare_external(predictions[[i]], endpoint)
  })
  names(evaluations) <- names(cohorts)
  result <- list(
    task_type = validated$task_type, development = cv, fits = fits,
    predictions = predictions, cohorts = evaluations,
    results = c(stats::setNames(list(cv), development_name), evaluations),
    cohort_roles = c(stats::setNames("development_cv", development_name), roles),
    model_pipelines = pipelines,
    elapsed_seconds = unname(proc.time()[["elapsed"]] - started)
  )
  class(result) <- "stormbreaker_cohort_comparison"
  sb_attach_result_page(result, view, report_file)
}

#' @export
print.stormbreaker_cohort_comparison <- function(x, ...) {
  cat("<stormbreaker_cohort_comparison>\n")
  cat("Task:", x$task_type, "| Pipelines:", length(x$model_pipelines),
    "| Evaluation cohorts:", length(x$cohorts), "\n")
  cat("Results:", paste(names(x$results), collapse = ", "), "\n")
  if (!is.null(x$report_file)) cat("Results page:", x$report_file, "\n")
  invisible(x)
}
