#' List short model names for doctor-facing entrypoints
#' @param task One of binary, survival or multiclass.
#' @return A model menu. Short names map to existing native configurations.
#'   Binary rows declare can_predict, can_extract_features and stage_selector.
#' @export
sb_models <- function(task = c("binary", "survival", "multiclass")) {
  task <- match.arg(task)
  if (task == "multiclass") {
    menu <- sb_method_catalog(task)
    menu$name <- menu$method
    menu$name[menu$method == "random_forest"] <- "RF"
    menu$name[menu$method == "svm_radial"] <- "SVM"
    menu$name[menu$method == "elastic_net"] <- "ElasticNet"
    return(menu)
  }
  menu <- sb_model_catalog(task)
  menu$name <- gsub("_(binary|survival|fixed|logistic|cox)(?=_|$)", "", menu$model_id, perl = TRUE)
  aliases <- c(
    logistic = "Logistic", coxph = "Cox", glmnet_logistic = "Lasso", glmnet_cox = "Lasso",
    ridge_logistic = "Ridge", ridge_cox = "Ridge", elastic_net_logistic = "ElasticNet", elastic_net_cox = "ElasticNet",
    random_forest_binary = "RF", ranger_survival = "RF", gbm_binary = "GBM", gbm_cox = "GBM",
    gaussian_naive_bayes = "NaiveBayes", bernoulli_naive_bayes_binary = "BernoulliNB",
    nearest_centroid_binary = "NearestCentroid", svm_radial_binary = "SVM",
    survivalsvm = "SVM", xgboost_binary = "XGBoost", xgboost_cox = "XGBoost",
    coxboost_cv = "CoxBoost", coxboost = "CoxBoost_fixed", lda_binary = "LDA", qda_binary = "QDA",
    plsrcox = "PLSR", plsrcox_cv = "PLSR_cv"
  )
  hits <- match(menu$model_id, names(aliases))
  menu$name[!is.na(hits)] <- unname(aliases[hits[!is.na(hits)]])
  menu$name <- make.unique(menu$name, sep = "-")
  if (task == "binary") menu <- sb_binary_unified_menu(menu)
  menu
}

sb_doctor_model_ids <- function(models, menu) {
  if (is.null(models)) {
    return(menu$model_id)
  }
  if (!is.character(models) || !length(models) || anyNA(models)) stop("Choose model names from sb_models().")
  if ("plsrcox" %in% menu$model_id) {
    old <- c(plsrcox = "plsrcox", plsrcox_cv = "plsrcox_cv", mboost = "glmboost_cox")
    alias <- match(tolower(models), names(old))
    models[!is.na(alias)] <- unname(old[alias[!is.na(alias)]])
  }
  index <- match(tolower(models), tolower(menu$name))
  index[is.na(index)] <- match(sb_canonical_learner_id(models[is.na(index)]), menu$model_id)
  if (anyNA(index)) stop("Unknown model: ", paste(models[is.na(index)], collapse = ", "), ". See sb_models().")
  if (anyDuplicated(index)) stop("Choose each model once.")
  menu$model_id[index]
}

sb_doctor_pipeline <- function(model, selector, task) {
  if (task == "binary" && is.null(selector)) return(sb_binary_model_route(model))
  ids <- sb_doctor_model_ids(model, sb_models(task))
  if (length(ids) != 1L) stop("Choose one final learner for a pipeline.")
  base <- sb_model_pipelines(task, ids)[[1L]]
  sb_doctor_pipeline_from_base(base, selector, task)
}

sb_doctor_pipeline_from_base <- function(base, selector, task) {
  if (is.null(selector)) {
    return(base)
  }
  if (!inherits(selector, "stormbreaker_selector_chain_spec") || selector$task_type != task) {
    stop("Selector and model tasks must match.")
  }
  stages <- lapply(selector$ids, sb_doctor_stage, learner = base$learner, base = base)
  if (length(stages) == 1L) {
    return(stages[[1L]])
  }
  sb_new_chain_pipeline(base, stages)
}

#' Create an ordered modeling route with existing native defaults
#' @param model For binary tasks, one name or an ordered character vector from
#'   [sb_models()], for example c("Lasso", "RF", "GBM"). Intermediate
#'   methods extract features; the last method predicts. Survival takes one
#'   final model name. No model is fitted by this constructor.
#' @param selector Legacy selector specification, retained for compatibility.
#'   New binary code supplies an ordered model vector instead. Survival may
#'   use [surv_selector()].
#' @return An existing-compatible pipeline; no model is fitted.
#' @export
bi_pipeline <- function(model = "RF", selector = NULL) sb_doctor_pipeline(model, selector, "binary")

#' @rdname bi_pipeline
#' @export
surv_pipeline <- function(model = "Cox", selector = NULL) sb_doctor_pipeline(model, selector, "survival")

sb_doctor_pipeline_set <- function(task, models, stages, selector = NULL) {
  if (task == "binary") return(sb_binary_model_plan(models, stages, selector))
  menu <- sb_models(task)
  ids <- sb_doctor_model_ids(models, menu)
  names_chosen <- menu$name[match(ids, menu$model_id)]
  bases <- sb_model_pipelines(task, ids)
  names(bases) <- names_chosen
  if (!is.numeric(stages) || !length(stages) || anyNA(stages) || !all(stages %in% 1:3)) {
    stop("stages must contain 1, 2 or 3.")
  }
  pipelines <- list()
  skipped <- list()
  add <- function(model, selections = NULL) {
    label <- paste(c(if (!is.null(selections)) selections$methods, model), collapse = " + ")
    p <- tryCatch(sb_doctor_pipeline_from_base(bases[[model]], selections, task), error = identity)
    if (inherits(p, "error")) {
      skipped[[label]] <<- conditionMessage(p)
    } else {
      pipelines[[label]] <<- p
    }
  }
  if (!is.null(selector)) {
    for (m in names_chosen) add(m, selector)
  } else {
    if (1 %in% stages) for (m in names_chosen) add(m)
    if (2 %in% stages) {
      for (s in names(sb_doctor_selector_menu(task))) {
        for (m in names_chosen) add(m, sb_doctor_selector(s, task))
      }
    }
    if (3 %in% stages) {
      presets <- list(c("Boruta", "StepCox", "Cox"), c("Lasso", "StepCox", "CoxBoost"))
      for (p in presets) if (p[3L] %in% names_chosen) add(p[3L], sb_doctor_selector(p[1:2], task))
    }
  }
  list(pipelines = pipelines, unavailable = skipped)
}

sb_doctor_data_task <- function(data, task, id, predictors, outcome = NULL,
                                time = NULL, event = NULL, time_unit = "month",
                                prediction_horizon = NULL, decision_thresholds = NULL) {
  if (inherits(data, c("stormbreaker_binary_task", "stormbreaker_surv_task", "stormbreaker_multiclass_task"))) {
    return(data)
  }
  if (!is.data.frame(data)) stop("Supply a data frame or a task.")
  if (is.null(id)) {
    id <- ".sb_row_id"
    if (id %in% names(data)) stop("Specify the ID column explicitly.")
    data[[id]] <- paste0("row", seq_len(nrow(data)))
  }
  if (is.null(predictors)) predictors <- setdiff(names(data), c(id, outcome, time, event))
  switch(task,
    binary = sb_binary_task(data, id, outcome, predictors,
      decision_thresholds = decision_thresholds),
    survival = sb_surv_task(data, id, time, event, predictors, time_unit,
      prediction_horizon = prediction_horizon, decision_thresholds = decision_thresholds),
    multiclass = sb_multiclass_task(data, id, outcome, predictors)
  )
}

sb_doctor_run <- function(task, type, models, stages, selector, folds, seed, plots, pipelines = NULL) {
  if ((type == "binary" && !inherits(task, "stormbreaker_binary_task")) ||
    (type == "survival" && !inherits(task, "stormbreaker_surv_task"))) {
    stop("Task and entrypoint do not match.")
  }
  if (is.null(pipelines)) {
    plan <- sb_doctor_pipeline_set(type, models, stages, selector)
  } else {
    if (!is.null(models) || !is.null(selector)) stop("Use pipelines or models/selector, not both.")
    if (!is.list(pipelines) || !length(pipelines) || is.null(names(pipelines)) ||
      anyNA(names(pipelines)) || any(!nzchar(names(pipelines))) || anyDuplicated(names(pipelines))) {
      stop("Supply pipelines as a non-empty, uniquely named list.")
    }
    pipelines <- lapply(pipelines, sb_validate_pipeline_contract, task_type = type)
    plan <- list(pipelines = pipelines, unavailable = list())
  }
  if (type == "survival" && anyDuplicated(sb_task_group_ids(task))) {
    for (name in names(plan$pipelines)) {
      p <- plan$pipelines[[name]]
      stage_ids <- if (p$selector == "ordered_chain") {
        vapply(p$selector_params$stages, `[[`, character(1), "selector")
      } else {
        p$selector
      }
      if (p$learner %in% sb_analysis_tuned_survival_learner_ids() || "coxboost_nonzero" %in% stage_ids) {
        plan$unavailable[[name]] <- "Native survival CV does not accept patient groups; not run on repeated patients."
        plan$pipelines[[name]] <- NULL
      }
    }
  }
  if (!length(plan$pipelines)) {
    stop("No requested combination is supported. Choose models/selectors from the task menu.")
  }
  result <- sb_compare_cv(task, unname(plan$pipelines), folds = folds, seed = seed)
  overview <- list()
  plot_errors <- list()
  if (is.null(plots) || isTRUE(plots)) {
    display_models <- sb_clinician_models(sb_viz_leaves(result))
    for (name in c("performance", if (type == "binary") "roc")) {
      p <- tryCatch(if (name == "performance") {
        sb_plot_performance(
          result,
          if (type == "binary") "roc_auc" else "harrell_c"
        )
      } else {
        sb_plot_roc(result)
      }, error = identity)
      if (inherits(p, "error")) {
        plot_errors[[name]] <- conditionMessage(p)
      } else {
        p <- sb_clinician_clean_plot(sb_compact_plot(p), display_models)
        if (name == "performance") p <- p + ggplot2::labs(title = "Model comparison", x = "Model")
        overview[[name]] <- p
      }
    }
  }
  details <- list()
  if (isTRUE(plots) || (is.null(plots) && length(plan$pipelines) == 1L)) {
    sb_with_selector_reuse(for (name in names(plan$pipelines)) {
      fit <- sb_fit_final(task, plan$pipelines[[name]], capture_diagnostics = TRUE)
      figures <- NULL
      if (fit$status == "success") {
        figures <- tryCatch(sb_plots(fit), error = function(error) {
          list(plots = list(), unavailable = list(algorithm = conditionMessage(error)))
        })
      }
      details[[name]] <- list(fit = fit, figures = figures)
    })
  }
  output <- list(
    comparison = result, pipelines = plan$pipelines, unavailable = plan$unavailable,
    details = details, plots = overview, plot_errors = plot_errors, task_type = type,
    development_task = task
  )
  output$issues <- sb_run_issues(output)
  structure(output, class = c("stormbreaker_run", "list"))
}

#' Run binary models with short names and fixed algorithm defaults
#' @param data Data frame or an existing task. Repeated patients require a task
#'   created with its group column; automatic row IDs do not infer patients.
#' @param outcome Outcome column (binary 0/1, or multiclass labels).
#' @param predictors Feature column names. NULL uses all non-role columns.
#' @param id ID column. NULL generates row IDs for independent rows.
#' @param models Method names from [sb_models()]. NULL requests the whole task
#'   menu. For binary tasks this is one shared combination pool: prediction
#'   models and feature methods are not chosen from separate menus.
#' @param stages Defaults to 1 when models are specified. Otherwise binary
#'   defaults to c(1,2), and survival retains c(1,2,3). Binary 1 fits each
#'   predictive method; 2 enumerates every supported extractor/model pair
#'   within the chosen pool; explicit 3 enumerates two ordered extractors
#'   followed by a model, with repetitions allowed. It is not just presets.
#'   Large pools can be expensive; three-stage expansion is never implicit
#'   for binary tasks. Survival retains its fixed three-stage presets.
#' @param selector Legacy explicit selector specification, overriding stages.
#'   For new binary routes use [bi_pipeline()] with ordered model names.
#' @param pipelines Optional uniquely named list from [bi_pipeline()] or
#'   [surv_pipeline()]. Mix single learners and ordered selector chains in one
#'   comparison with shared folds. Names are used by [sb_complete_run()].
#'   Do not also supply models or selector; stages is ignored for this list.
#' @param folds,seed Outer resampling controls.
#' @param decision_thresholds Optional risk probabilities for decision curves.
#'   Set these when planning the analysis. NULL leaves decision curves unrequested.
#'   For an existing task, its declared settings are used.
#' @param plots NULL (default) draws comparison overviews, but only performs a
#'   final fit and detailed plotting when exactly one pipeline is requested.
#'   TRUE explicitly fits and plots all candidates; FALSE skips all plots and
#'   final fits. Use [sb_complete_run()] for a chosen binary/survival pipeline
#'   after batch comparison. No browser opens; unavailable views stay visible.
#' @return List with CV comparison, requested pipelines, unavailable combinations
#'   and final-fit figures. The issues table offers readable explanations while
#'   preserving original reasons. No winner is automatically selected.
#'   The result retains development data for [sb_complete_run()]; treat saved
#'   run objects as patient-level data, not anonymous reports.
#'   Printing shows a short summary; all existing list fields remain available.
#'   Optional plotting failures are recorded in issues without discarding fits.
#' @examples
#' genes <- c("gene_a", "gene_b")
#' development <- bi_example(40, seed = 42)
#' binary <- bi_run(development, outcome = "outcome", id = "patient_id",
#'   predictors = genes, models = "Logistic", folds = 2, seed = 42,
#'   plots = FALSE)
#' result_table(binary$comparison, metric = "roc_auc")
#'
#' follow_up <- surv_example(40, seed = 42)
#' survival <- surv_run(follow_up, time = "time", event = "event",
#'   id = "patient_id", predictors = genes, time_unit = "month",
#'   models = "Cox", folds = 2, seed = 42, plots = FALSE)
#' result_table(survival$comparison, metric = "harrell_c")
#' @export
bi_run <- function(data, outcome = NULL, predictors = NULL, id = NULL, models = NULL,
                   stages = if (is.null(models)) c(1L, 2L) else 1L,
                   selector = NULL, folds = 5L, seed = 1L, plots = NULL,
                   pipelines = NULL, decision_thresholds = NULL) {
  task <- sb_doctor_data_task(data, "binary", id, predictors, outcome,
    decision_thresholds = decision_thresholds)
  sb_doctor_run(task, "binary", models, stages, selector, folds, seed, plots, pipelines)
}

#' @rdname bi_run
#' @param time,event Survival time and 0/1 event columns.
#' @param time_unit Unit of the supplied follow-up times, default month.
#' @param prediction_horizon Optional survival prediction time in time_unit.
#'   Enables horizon-specific metrics and curves. NULL requests risk ranking only.
#'   For an existing task, its declared horizon is used.
#' @export
surv_run <- function(data, time = NULL, event = NULL, predictors = NULL, id = NULL, models = NULL,
                     stages = if (is.null(models)) c(1L, 2L, 3L) else 1L,
                     selector = NULL, folds = 5L, seed = 1L, plots = NULL,
                     time_unit = "month", pipelines = NULL,
                     prediction_horizon = NULL, decision_thresholds = NULL) {
  task <- sb_doctor_data_task(data, "survival", id, predictors, time = time,
    event = event, time_unit = time_unit, prediction_horizon = prediction_horizon,
    decision_thresholds = decision_thresholds)
  sb_doctor_run(task, "survival", models, stages, selector, folds, seed, plots, pipelines)
}

#' Run existing multiclass learners by short name
#' @inheritParams bi_run
#' @param plots NULL draws detailed figures for one method only; TRUE draws all
#'   requested methods and FALSE skips detailed fits. For selected multiclass
#'   methods after comparison, use [sb_analyze_methods()].
#' @return Existing multiclass CV comparison. Multiclass selector chains and
#'   three-stage presets are not implemented by this entrypoint.
#' @export
multi_run <- function(data, outcome = NULL, predictors = NULL, id = NULL,
                      models = NULL, folds = 5L, seed = 1L, plots = NULL) {
  task <- sb_doctor_data_task(data, "multiclass", id, predictors, outcome)
  menu <- sb_models("multiclass")
  menu <- menu[menu$supported & !menu$method %in% sb_multiclass_standalone_selectors(), ]
  if (is.null(models)) {
    methods <- menu$method
  } else {
    index <- match(tolower(models), tolower(menu$name))
    if (anyNA(index)) stop("Choose supported multiclass learner names from sb_models('multiclass').")
    methods <- menu$method[index]
  }
  list(
    comparison = sb_compare_multiclass(task, methods = methods, folds = folds, seed = seed),
    details = if (isTRUE(plots) || (is.null(plots) && length(methods) == 1L)) {
      sb_analyze_methods(task, methods)
    } else {
      NULL
    },
    methods = methods, task_type = "multiclass"
  )
}
