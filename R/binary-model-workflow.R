sb_binary_unified_menu <- function(menu) {
  selectors <- sb_doctor_selector_menu("binary")
  shared <- c(
    Lasso = "glmnet_logistic", Stepglm = "stepglm_binary",
    glmBoost = "glmboost_binary", GamBoost = "gamboost_binary",
    RF = "random_forest_binary", Enet = "elastic_net_logistic",
    PAMR = "pamr_binary", PLSR = "plsrglm_binary",
    GBM = "gbm_binary", XGBoost = "xgboost_binary",
    ABESS = "abess_binary_gic"
  )
  menu$can_predict <- TRUE
  menu$stage_selector <- unname(selectors[names(shared)[match(menu$model_id, shared)]])
  menu$can_extract_features <- !is.na(menu$stage_selector)
  additional <- setdiff(names(selectors), names(shared))
  registry <- sb_pipeline_registry("binary")
  extra <- menu[rep(1L, length(additional)), , drop = FALSE]
  for (column in names(extra)) extra[[column]][] <- NA
  for (i in seq_along(additional)) {
    id <- unname(selectors[additional[i]])
    row <- registry[registry$selector == id & registry$learner == "logistic", , drop = FALSE]
    extra$model_id[i] <- id
    extra$model_label[i] <- additional[i]
    extra$name[i] <- additional[i]
    extra$task_type[i] <- "binary"
    extra$source[i] <- "stormbreaker_builtin"
    extra$selector[i] <- id
    extra$admission_status[i] <- "admitted"
    extra$tuning_semantics[i] <- "native_feature_method_defaults"
    extra$parameter_semantics[i] <- "Existing native feature extraction; no terminal prediction adapter"
    extra$optional_dependency[i] <- row$optional_dependency[1L]
    if (id %in% sb_forest_selector_ids()) {
      dependency <- sb_forest_selector_dependency(id)
      extra$optional_dependency[i] <- paste(dependency$package, dependency$engine, sep = " + ")
    }
    extra$can_predict[i] <- FALSE
    extra$can_extract_features[i] <- TRUE
    extra$stage_selector[i] <- id
  }
  dependency_columns <- c("r_packages", "python_modules", "missing_r_packages", "dependency_status")
  extra <- sb_add_dependency_status(extra[setdiff(names(extra), dependency_columns)])
  menu <- rbind(menu, extra[names(menu)])
  rownames(menu) <- NULL
  menu
}

sb_binary_method_ids <- function(methods, menu) {
  if (!is.character(methods) || !length(methods) || anyNA(methods) ||
      any(!nzchar(trimws(methods)))) {
    sb_abort("Choose method names from sb_models('binary').", "stormbreaker_bad_model_choice")
  }
  index <- match(tolower(methods), tolower(menu$name))
  missing <- is.na(index)
  index[missing] <- match(sb_canonical_learner_id(methods[missing]), menu$model_id)
  selectors <- sb_doctor_selector_menu("binary")
  aliases <- match(tolower(methods), tolower(names(selectors)))
  missing <- is.na(index) & !is.na(aliases)
  index[missing] <- match(unname(selectors[aliases[missing]]), menu$stage_selector)
  if (anyNA(index)) {
    sb_abort(paste0("Unknown binary method: ", paste(methods[is.na(index)], collapse = ", "),
      ". See sb_models('binary')."), "stormbreaker_bad_model_choice")
  }
  index
}

sb_binary_model_route <- function(methods, menu = sb_models("binary"), bases = NULL) {
  rows <- menu[sb_binary_method_ids(methods, menu), , drop = FALSE]
  last <- nrow(rows)
  if (!rows$can_predict[last]) {
    sb_abort(paste0(rows$name[last], " can extract features but cannot be the final prediction model."),
      "stormbreaker_method_cannot_predict")
  }
  front <- seq_len(last - 1L)
  if (any(!rows$can_extract_features[front])) {
    sb_abort(paste0("No feature-extraction adapter for intermediate method(s): ",
      paste(rows$name[front][!rows$can_extract_features[front]], collapse = ", "), "."),
      "stormbreaker_method_cannot_extract")
  }
  base <- if (is.null(bases)) sb_model_pipelines("binary", rows$model_id[last])[[1L]] else bases[[rows$model_id[last]]]
  if (!length(front)) return(base)
  stages <- lapply(rows$stage_selector[front], sb_doctor_stage, learner = base$learner, base = base)
  if (length(stages) == 1L) return(stages[[1L]])
  sb_new_chain_pipeline(base, stages)
}

sb_binary_model_plan <- function(models, stages, selector = NULL) {
  menu <- sb_models("binary")
  chosen <- if (is.null(models)) seq_len(nrow(menu)) else sb_binary_method_ids(models, menu)
  if (anyDuplicated(chosen)) {
    sb_abort("Choose each binary method once in the combination pool.", "stormbreaker_bad_model_choice")
  }
  if (!is.numeric(stages) || !length(stages) || anyNA(stages) || !all(stages %in% 1:3)) {
    stop("stages must contain 1, 2 or 3.")
  }
  pool <- menu[chosen, , drop = FALSE]
  terminal <- pool[pool$can_predict, , drop = FALSE]
  extractors <- pool[pool$can_extract_features, , drop = FALSE]
  bases <- if (nrow(terminal)) sb_model_pipelines("binary", terminal$model_id) else list()
  pipelines <- unavailable <- list()
  add <- function(methods) {
    label <- paste(methods, collapse = " + ")
    value <- tryCatch(sb_binary_model_route(methods, menu, bases), error = identity)
    if (inherits(value, "error")) {
      unavailable[[label]] <<- conditionMessage(value)
    } else {
      pipelines[[label]] <<- value
    }
  }
  if (!is.null(selector)) {
    for (i in seq_len(nrow(terminal))) {
      name <- terminal$name[i]
      base <- bases[[terminal$model_id[i]]]
      value <- tryCatch(sb_doctor_pipeline_from_base(base, selector, "binary"), error = identity)
      label <- paste(c(selector$methods, name), collapse = " + ")
      if (inherits(value, "error")) unavailable[[label]] <- conditionMessage(value) else pipelines[[label]] <- value
    }
    return(list(pipelines = pipelines, unavailable = unavailable))
  }
  if (1L %in% stages) for (model in terminal$name) add(model)
  if (2L %in% stages) {
    for (first in extractors$name) for (model in terminal$name) add(c(first, model))
  }
  if (3L %in% stages) {
    for (first in extractors$name) for (second in extractors$name) {
      for (model in terminal$name) add(c(first, second, model))
    }
  }
  for (name in pool$name[!pool$can_predict]) {
    if (!nrow(terminal) || !any(stages %in% c(2L, 3L))) {
      unavailable[[name]] <- "Feature-extraction method: choose a final prediction model and stages = 2 or 3."
    }
  }
  list(pipelines = pipelines, unavailable = unavailable)
}
