sb_pipeline_inner_selection_method <- function(pipeline) {
  if (pipeline$learner %in% sb_svm_probability_learner_ids()) {
    if (sb_selector_is_learned(pipeline$selector)) {
      return("svm_fold_local_probability_calibration")
    }
    return(if (identical(pipeline$selector, "fixed_features")) {
      "predeclared_feature_subset"
    } else "none")
  }
  if (!pipeline$learner %in% c("glmnet_cox", "glmnet_logistic")) {
    return(NA_character_)
  }
  if (sb_selector_is_learned(pipeline$selector)) {
    return("glmnet_exclude_function")
  }
  if (identical(pipeline$selector, "fixed_features")) {
    return("predeclared_feature_subset")
  }
  "none"
}

sb_fit_runtime_selector_uncached <- function(selector, x, endpoint) {
  if (identical(selector$id, "ordered_chain")) return(selector$fit(x, endpoint))
  task_type <- endpoint[["task_type"]]
  if (!isTRUE(selector[["needs_outcome"]])) {
    return(selector[["fit"]](x))
  }
  if (identical(task_type, "survival")) {
    if (selector$id %in% c("lasso_nonzero", "enet_cox_nonzero", "ridge_cox_weights")) return(
      selector$fit(x, endpoint$time, endpoint$event, group_ids = endpoint$group_ids))
    return(selector[["fit"]](
      x,
      endpoint[["time"]],
      endpoint[["event"]]
    ))
  }
  if (identical(task_type, "binary")) {
    if (selector$id %in% sb_binary_model_selector_ids()) return(
      selector$fit(x, endpoint$outcome, group_ids = endpoint$group_ids))
    return(selector[["fit"]](x, endpoint[["outcome"]]))
  }
  sb_abort(
    "The runtime endpoint type is not supported.",
    "stormbreaker_bad_runtime_endpoint"
  )
}

sb_make_glmnet_fold_selector_filter <- function(selector, fitted_selector,
                                                patient_ids, fold_id,
                                                task_type, group_ids = NULL) {
  if (utils::packageVersion("glmnet") <
      base::numeric_version("4.1.2")) {
    sb_abort(
      paste0(
        "Fold-local supervised selection requires `glmnet` >= 4.1-2, ",
        "whose function-valued `exclude` filter is refitted inside ",
        "each CV analysis fold."
      ),
      "stormbreaker_glmnet_fold_filter_unavailable"
    )
  }
  patient_ids <- as.character(patient_ids)
  valid_ids <- length(patient_ids) == length(fold_id) &&
    length(patient_ids) > 0L && !anyNA(patient_ids) &&
    all(nzchar(trimws(patient_ids))) && anyDuplicated(patient_ids) == 0L
  if (!valid_ids) {
    sb_abort(
      "Fold-local glmnet selection requires one unique patient ID per row.",
      "stormbreaker_bad_runtime_patient_ids"
    )
  }
  if (!task_type %in% c("survival", "binary")) {
    sb_abort(
      "Fold-local glmnet selection received an unsupported endpoint type.",
      "stormbreaker_bad_runtime_endpoint"
    )
  }

  expected_folds <- sort(unique(fold_id))
  records <- stats::setNames(
    vector("list", length(expected_folds)),
    as.character(expected_folds)
  )
  exclude <- function(x, y, weights = NULL, ...) {
    row_ids <- rownames(x)
    positions <- match(row_ids, patient_ids)
    valid_rows <- !is.null(row_ids) && length(row_ids) == nrow(x) &&
      !anyNA(positions) && anyDuplicated(row_ids) == 0L
    if (!valid_rows) {
      sb_abort(
        paste0(
          "glmnet did not preserve the patient-row identity needed for ",
          "fold-local supervised selection."
        ),
        "stormbreaker_glmnet_fold_filter_protocol"
      )
    }

    x_frame <- as.data.frame(
      x,
      optional = TRUE,
      stringsAsFactors = FALSE
    )
    if (length(positions) == length(patient_ids)) {
      fitted <- fitted_selector
    } else {
      omitted_positions <- setdiff(seq_along(patient_ids), positions)
      omitted_folds <- unique(fold_id[omitted_positions])
      valid_inner_analysis <- length(omitted_folds) == 1L &&
        setequal(positions, which(fold_id != omitted_folds[[1L]]))
      if (!valid_inner_analysis) {
        sb_abort(
          paste0(
            "glmnet supplied rows that do not match one declared inner ",
            "analysis fold."
          ),
          "stormbreaker_glmnet_fold_filter_protocol"
        )
      }
      fold <- as.integer(omitted_folds[[1L]])
      endpoint <- if (identical(task_type, "survival")) {
        list(
          task_type = "survival",
          time = as.numeric(y[, 1L]),
          event = as.integer(y[, 2L])
        )
      } else {
        list(
          task_type = "binary",
          outcome = as.integer(y)
        )
      }
      if (!is.null(group_ids)) endpoint$group_ids <- group_ids[positions]
      fitted <- sb_fit_runtime_selector(selector, x_frame, endpoint)
      records[[as.character(fold)]] <<- list(
        fold = fold,
        analysis_patient_ids = row_ids,
        selected_features = fitted$selected_features
      )
    }
    which(!colnames(x) %in% fitted$selected_features)
  }
  inner_selections <- function() {
    result <- unname(records[as.character(expected_folds)])
    if (any(vapply(result, is.null, logical(1)))) {
      sb_abort(
        paste0(
          "glmnet did not invoke the fold-local selector for every ",
          "declared inner analysis fold."
        ),
        "stormbreaker_glmnet_fold_filter_protocol"
      )
    }
    result
  }
  list(exclude = exclude, inner_selections = inner_selections)
}

sb_prepare_runtime_learner <- function(learner, pipeline, x, endpoint,
                                       patient_ids,
                                       selection_context = NULL,
                                       group_ids = NULL) {
  context <- .sb_selector_reuse_context$active
  if (is.null(context) || !learner$id %in% c("glmnet_logistic", "glmnet_cox",
      sb_svm_probability_learner_ids())) {
    return(sb_prepare_runtime_learner_uncached(learner, pipeline, x, endpoint,
      patient_ids, selection_context, group_ids))
  }
  selection <- NULL
  if (sb_selector_is_learned(pipeline$selector) && !is.null(selection_context)) {
    s <- selection_context$selector
    selection <- list(x = selection_context$x,
      fitted_selector = selection_context$fitted_selector,
      id = s$id, params = s$params, capture_history = s$capture_history,
      upstream_context = s$upstream_context,
      fit_body = body(s$fit), fit_formals = formals(s$fit))
  }
  signature <- list(params = learner$params, x = x, endpoint = endpoint,
    patient_ids = patient_ids, group_ids = group_ids, selection = selection,
    selector_is_learned = sb_selector_is_learned(pipeline$selector),
    rng_kind = RNGkind())
  entries <- context$preparations[[learner$id]]
  for (entry in entries) {
    if (identical(signature, entry$signature)) {
      context$preparation_hits <- context$preparation_hits + 1L
      prepared <- entry$value
      prepared$tuning$inner_selection_method <- sb_pipeline_inner_selection_method(pipeline)
      return(prepared)
    }
  }
  context$preparation_misses <- context$preparation_misses + 1L
  warned <- FALSE
  prepared <- withCallingHandlers(
    sb_prepare_runtime_learner_uncached(learner, pipeline, x, endpoint,
      patient_ids, selection_context, group_ids),
    warning = function(w) warned <<- TRUE)
  if (!warned && length(prepared$tuning$warnings) == 0L) {
    context$preparations[[learner$id]] <- c(entries,
      list(list(signature = signature, value = prepared)))
  }
  prepared
}

sb_prepare_runtime_learner_uncached <- function(learner, pipeline, x, endpoint,
                                                patient_ids,
                                                selection_context = NULL,
                                                group_ids = NULL) {
  if (!identical(length(patient_ids), nrow(x))) {
    sb_abort(
      "Runtime patient IDs must match the fitted analysis rows.",
      "stormbreaker_bad_runtime_patient_ids"
    )
  }
  tuning_x <- x
  tuning_selector <- NULL
  tuning_fitted_selector <- NULL
  if (learner[["id"]] %in% c("glmnet_cox", "glmnet_logistic",
                            sb_svm_probability_learner_ids()) &&
      sb_selector_is_learned(pipeline[["selector"]])) {
    valid_context <- is.list(selection_context) &&
      is.list(selection_context[["selector"]]) &&
      inherits(
        selection_context[["fitted_selector"]],
        "stormbreaker_selector_fit"
      ) && is.data.frame(selection_context[["x"]]) &&
      nrow(selection_context[["x"]]) == nrow(x)
    if (!valid_context) {
      sb_abort(
        paste0(
          "A selected internally resampled learner requires its full outer-analysis ",
          "predictors and fitted selector."
        ),
        if (learner$id %in% sb_svm_probability_learner_ids()) {
          "stormbreaker_bad_svm_selection_context"
        } else "stormbreaker_bad_glmnet_selection_context"
      )
    }
    tuning_x <- selection_context[["x"]]
    tuning_selector <- selection_context[["selector"]]
    tuning_fitted_selector <- selection_context[["fitted_selector"]]
  }
  if (learner$id %in% sb_svm_probability_learner_ids()) {
    tuned <- sb_calibrate_svm_probability(tuning_x, endpoint$outcome,
      learner$id, learner$params, patient_ids, group_ids, tuning_selector)
    tuned$inner_selection_method <- sb_pipeline_inner_selection_method(pipeline)
    return(list(
      learner = sb_svm_calibrated_runtime_learner(learner, tuned$calibration),
      tuning = tuned
    ))
  }
  if (identical(learner[["id"]], "glmnet_cox")) {
    tuned <- sb_select_glmnet_cox_lambda(
      tuning_x,
      endpoint[["time"]],
      endpoint[["event"]],
      inner_folds = learner[["params"]][["inner_folds"]],
      seed = learner[["params"]][["inner_seed"]],
      alpha = learner[["params"]][["alpha"]],
      standardize = learner[["params"]][["standardize"]],
      selector = tuning_selector,
      fitted_selector = tuning_fitted_selector,
      patient_ids = if (is.null(tuning_selector)) NULL else patient_ids,
      selection_rule = learner[["params"]][["selection_rule"]],
      measure = learner[["params"]][["measure"]],
      group_ids = group_ids
    )
    tuned$inner_selection_method <- sb_pipeline_inner_selection_method(
      pipeline
    )
    names(tuned$fold_id) <- as.character(patient_ids)
    return(list(
      learner = sb_glmnet_cox_learner_adapter(
        alpha = tuned$alpha,
        lambda = tuned$lambda,
        standardize = tuned$standardize
      ),
      tuning = tuned
    ))
  }
  if (identical(learner[["id"]], "glmnet_logistic")) {
    tuned <- sb_select_glmnet_logistic_lambda(
      tuning_x,
      endpoint[["outcome"]],
      inner_folds = learner[["params"]][["inner_folds"]],
      seed = learner[["params"]][["inner_seed"]],
      alpha = learner[["params"]][["alpha"]],
      standardize = learner[["params"]][["standardize"]],
      selector = tuning_selector,
      fitted_selector = tuning_fitted_selector,
      patient_ids = if (is.null(tuning_selector)) NULL else patient_ids,
      selection_rule = learner[["params"]][["selection_rule"]],
      measure = learner[["params"]][["measure"]],
      group_ids = group_ids
    )
    tuned$inner_selection_method <- sb_pipeline_inner_selection_method(
      pipeline
    )
    names(tuned$fold_id) <- as.character(patient_ids)
    return(list(
      learner = sb_glmnet_logistic_learner_adapter(
        alpha = tuned$alpha,
        lambda = tuned$lambda,
        standardize = tuned$standardize
      ),
      tuning = tuned
    ))
  }
  list(learner = learner, tuning = NULL)
}

sb_fit_runtime_learner <- function(learner, x, endpoint) {
  task_type <- endpoint[["task_type"]]
  if (identical(task_type, "survival")) {
    return(learner[["fit"]](
      x,
      endpoint[["time"]],
      endpoint[["event"]]
    ))
  }
  if (identical(task_type, "binary")) {
    return(learner[["fit"]](x, endpoint[["outcome"]]))
  }
  sb_abort(
    "The runtime endpoint type is not supported.",
    "stormbreaker_bad_runtime_endpoint"
  )
}
