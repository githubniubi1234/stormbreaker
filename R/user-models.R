#' List task-specific models for a user-supplied gene set
#'
#' The menu contains every admitted learner with selector `none`, plus explicit
#' Ridge and alpha-0.5 elastic-net configurations for the glmnet learner.
#' Existing glmnet learner IDs denote Lasso (alpha 1). Model variants are not
#' an exhaustive hyperparameter grid. Dependencies are reported, not installed.
#' @param task_type Explicitly `"binary"` or `"survival"`.
#' @return A data frame with model IDs, learner IDs, glmnet alpha, registry
#'   semantics and dependency availability. Missing dependencies stay visible.
#' @examples
#' catalog <- sb_model_catalog("binary")
#' head(catalog[c("model_id", "r_packages", "dependency_status")])
#' @export
sb_model_catalog <- function(task_type) {
  if (!is.character(task_type) || length(task_type) != 1L || is.na(task_type) ||
      !task_type %in% c("binary", "survival")) {
    sb_abort("Choose task_type explicitly: 'binary' or 'survival'.",
      "stormbreaker_bad_task_type")
  }
  registry <- sb_pipeline_registry()
  rows <- registry[registry$task_type == task_type &
    registry$selector == "none" & registry$admission_status == "admitted", , drop = FALSE]
  rows <- rows[!duplicated(rows$learner), , drop = FALSE]
  rows$model_id <- rows$learner
  rows$glmnet_alpha <- ifelse(grepl("^glmnet_", rows$learner), 1, NA_real_)
  g <- rows[rows$learner == if (task_type == "binary") "glmnet_logistic" else "glmnet_cox", , drop = FALSE]
  ridge <- elastic <- g
  ridge$model_id <- paste0("ridge_", if (task_type == "binary") "logistic" else "cox")
  elastic$model_id <- paste0("elastic_net_", if (task_type == "binary") "logistic" else "cox")
  ridge$glmnet_alpha <- 0
  elastic$glmnet_alpha <- 0.5
  rows <- rbind(rows, ridge, elastic)
  rownames(rows) <- NULL
  rows <- rows[c("model_id", "learner", "glmnet_alpha",
    setdiff(names(rows), c("model_id", "learner", "glmnet_alpha")))]
  rows$model_label <- vapply(seq_len(nrow(rows)), function(i) {
    parameters <- if (is.na(rows$glmnet_alpha[i])) list() else list(alpha = rows$glmnet_alpha[i])
    sb_display_learner(rows$learner[i], parameters, rows$model_id[i])
  }, character(1))
  sb_add_dependency_status(rows)
}

#' Build chosen models, or all models, without an extra gene selector
#'
#' Every pipeline uses selector `none`: the task's user-supplied predictors
#' enter fitting without an added screening stage. Intrinsic sparsity, subset
#' selection, component estimation and tree splits remain part of the selected
#' learner; this does not force every coefficient to be nonzero.
#' @param task_type Explicitly `"binary"` or `"survival"`.
#' @param models Character model IDs from [sb_model_catalog()], in desired
#'   order. `NULL` means all catalog entries for this task, including entries
#'   whose optional dependency is missing. No winner is selected.
#' @param inner_folds,inner_seed Inner CV settings for glmnet configurations.
#'   Other learners retain their documented `sb_pipeline()` defaults.
#' @return A named list of [sb_pipeline()] specifications, ready for
#'   [sb_compare_cv()] or [sb_fit_final()].
#' @examples
#' sb_model_pipelines("binary", c("logistic", "gaussian_naive_bayes"))
#' sb_model_pipelines("survival", c("coxph", "survreg_weibull"))
#' @export
sb_model_pipelines <- function(task_type, models = NULL,
                               inner_folds = 5L, inner_seed = 1701L) {
  menu <- sb_model_catalog(task_type)
  if (is.null(models)) models <- menu$model_id
  if (!is.character(models) || !length(models) || anyNA(models) ||
      anyDuplicated(models) || any(!models %in% menu$model_id)) {
    sb_abort("Supply distinct model IDs from sb_model_catalog() for this task.",
      "stormbreaker_bad_model_choice")
  }
  selected <- menu[match(models, menu$model_id), , drop = FALSE]
  pipelines <- lapply(seq_len(nrow(selected)), function(i) {
    args <- list(selector = "none", learner = selected$learner[i])
    if (!is.na(selected$glmnet_alpha[i])) {
      args$glmnet_alpha <- selected$glmnet_alpha[i]
      args$inner_folds <- inner_folds
      args$inner_seed <- inner_seed
    }
    do.call(sb_pipeline, args)
  })
  names(pipelines) <- models
  pipelines
}

#' Compare a small default set or explicitly chosen models
#'
#' Users create a binary or survival task and explicitly supply `predictors`
#' to that task constructor. This wrapper does not infer the task from outcome
#' values, discover genes, add a selector, select a winning model or read
#' external data. It runs the chosen pipelines on the same outer folds.
#' Failed models remain in the underlying comparison rather than being replaced.
#' @param task A task from [sb_binary_task()] or [sb_surv_task()], with the
#'   user's gene set supplied as `predictors`.
#' @param models Defaults to `"recommended"`: logistic/Lasso/random forest for
#'   binary tasks, or Cox/Lasso Cox/ranger for survival. These are fixed,
#'   complementary starting models, not selected using observed performance.
#'   Supply IDs from [sb_model_catalog()] to choose models, or explicitly use
#'   `NULL` for all models. Missing engines remain failures, not substitutions.
#' @param folds,seed Shared outer resampling settings passed to [sb_compare_cv()].
#' @param inner_folds,inner_seed Glmnet inner CV settings; see [sb_model_pipelines()].
#' @param view Open the results matrix after comparison. Defaults to FALSE;
#'   set TRUE to explicitly request a browser.
#' @param report_file Optional new HTML path; saves a page even with `view = FALSE`.
#' @return A `stormbreaker_cv_comparison` with additional `model_pipelines`,
#'   `input_genes`, `requested_models` and total `elapsed_seconds`. Its table
#'   also contains `model_id`. Total time includes menu construction and CV,
#'   but not final fitting or external prediction.
#' @examples
#' d <- sb_example_binary(n = 60, seed = 42)
#' task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
#' result <- sb_compare_models(
#'   task, models = c("logistic", "gaussian_naive_bayes"), folds = 3, seed = 42
#' )
#' result$comparison[c("model_id", "status", "roc_auc")]
#' @export
sb_compare_models <- function(task, models = "recommended", folds = 5L, seed = 1L,
                               inner_folds = 5L, inner_seed = 1701L,
                               view = FALSE, report_file = NULL) {
  sb_check_page_request(view, report_file)
  started <- proc.time()[["elapsed"]]
  type <- if (inherits(task, "stormbreaker_binary_task")) "binary" else
    if (inherits(task, "stormbreaker_surv_task")) "survival" else
      sb_abort("Create the task with sb_binary_task() or sb_surv_task().",
        "stormbreaker_bad_task")
  if (identical(models, "recommended")) models <- if (type == "binary")
    c("logistic", "glmnet_logistic", "random_forest_binary") else
    c("coxph", "glmnet_cox", "ranger_survival")
  pipelines <- sb_model_pipelines(type, models, inner_folds, inner_seed)
  compared <- sb_compare_cv(task, pipelines, folds = folds, seed = seed)
  ids <- vapply(pipelines, function(pipeline) pipeline$id, character(1))
  compared$comparison$model_id <- names(ids)[match(compared$comparison$pipeline_id, ids)]
  compared$model_pipelines <- pipelines
  compared$input_genes <- task$predictors
  compared$requested_models <- names(pipelines)
  compared$elapsed_seconds <- unname(proc.time()[["elapsed"]] - started)
  sb_attach_result_page(compared, view, report_file)
}
