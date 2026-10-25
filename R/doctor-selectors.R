sb_doctor_selector_menu <- function(task) {
  if (task == "binary") c(Lasso = "lasso_binary_nonzero", Stepglm = "stepglm_binary_aic",
    glmBoost = "glmboost_binary_nonzero", GamBoost = "gamboost_binary_nonzero", RF = "rf_binary_importance",
    Enet = "enet_binary_nonzero", PAMR = "pamr_binary_nonzero", PLSR = "plsrglm_binary_weights", GBM = "gbm_binary_importance", XGBoost = "xgboost_binary_importance", SVM_RFE = "svm_rfe_binary",
    Boruta = "boruta_binary_roughfix", varSelRF = "varselrf_binary_oob_1se",
    VSURF = "vsurf_binary_prediction", ABESS = "abess_binary_gic_nonzero",
    SIS = "sis_binary_nlogn", mRMR = "mrmr_binary_nlogn", JMI = "jmi_binary_nlogn",
    CMIM = "cmim_binary_nlogn", DISR = "disr_binary_nlogn", JMIM = "jmim_binary_nlogn",
    LPC = "lpc_binary_nlogn", Stabsel = "stabsel_binary_lasso") else
    c(Lasso = "lasso_nonzero", Enet = "enet_cox_nonzero", StepCox = "stepcox", CoxBoost = "coxboost_nonzero",
      RF = "rsf_minimal_depth", GBM = "gbm_cox_importance", glmBoost = "glmboost_cox_nonzero",
      GamBoost = "gamboost_cox_nonzero",
      XGBoost = "xgboost_cox_importance", PAMR = "pamr_survival_nonzero",
      SuperPC = "superpc_survival_scores", Ridge = "ridge_cox_weights",
      PLSR = "plsrcox_weights", SVM = "svm_survival_weights", Boruta = "boruta_survival_roughfix",
      ABESS = "abess_cox_gic_nonzero", LPC = "lpc_cox_nlogn", Stabsel = "stabsel_cox_lasso")
}

sb_doctor_selector <- function(methods, task) {
  menu <- sb_doctor_selector_menu(task)
  if (!is.character(methods) || !length(methods) || anyNA(methods)) stop("Supply ordered selector names.")
  methods[tolower(methods) == "elasticnet"] <- "Enet"
  if (task == "survival") methods[tolower(methods) == "mboost"] <- "glmBoost"
  if (task == "survival") methods[tolower(methods) == "rsf"] <- "RF"
  index <- match(tolower(methods), tolower(names(menu)))
  if (anyNA(index)) stop("Unknown selector. Choices: ", paste(names(menu), collapse = ", "))
  structure(list(task_type = task, methods = names(menu)[index], ids = unname(menu[index])),
    class = "stormbreaker_selector_chain_spec")
}

#' Internal compatibility constructor for saved binary workflows
#' @param methods Ordered binary selector names from historical code.
#' @return An internal selector specification. New code uses [bi_pipeline()].
#' @keywords internal
bi_selector <- function(methods) sb_doctor_selector(methods, "binary")

#' Choose ordered survival feature stages by short name
#' @param methods Ordered survival names such as c("Lasso", "StepCox").
#'   Existing native defaults, selection rules and fold-local fitting are kept.
#'   Binary tasks use ordered names directly in [bi_pipeline()], not a separate
#'   binary selector entrypoint.
#' @return A survival-stage specification for [surv_pipeline()].
#' @export
surv_selector <- function(methods) sb_doctor_selector(methods, "survival")

sb_doctor_stage <- function(id, learner, base) {
  args <- list(selector = id, learner = learner)
  if (learner %in% c("glmnet_logistic", "glmnet_cox")) {
    args$inner_folds <- base$learner_params$inner_folds
    args$inner_seed <- base$learner_params$inner_seed
    args$glmnet_alpha <- base$learner_params$alpha
  }
  if (id == "stepcox") args$stepcox_direction <- "both"
  if (id %in% c("lasso_nonzero", "enet_cox_nonzero", "coxboost_nonzero")) {
    args$selector_inner_folds <- 5L
    args$selector_seed <- 1L
  }
  if (id == "rsf_minimal_depth") {
    args$selector_rsf_nodesize <- 5L
    args$selector_seed <- 1L
  }
  do.call(sb_pipeline, args)
}

sb_validate_chain_stages <- function(base, stages) {
  if (!is.list(stages) || !length(stages)) stop("A selector chain requires stages.")
  # Stage admission retains the restrictions on selection inside native CV.
  lapply(stages, function(s) {
    s <- sb_validate_pipeline_contract(s, base$task_type)
    if (!identical(s$learner, base$learner) || s$selector == "ordered_chain") stop("Invalid selector stage.")
    s
  })
}

sb_chain_pipeline_id <- function(base, stages) {
  paste(c("chain", vapply(stages, `[[`, character(1), "id"), base$id), collapse = "__then__")
}

sb_new_chain_pipeline <- function(base, stages) {
  if (!length(stages)) return(base)
  stages <- sb_validate_chain_stages(base, stages)
  result <- base
  result$selector <- "ordered_chain"
  result$selector_params <- list(stages = stages)
  result$chain_base <- base
  result$id <- sb_chain_pipeline_id(base, stages)
  result
}

sb_validate_chain_pipeline <- function(pipeline, task_type) {
  base <- sb_validate_pipeline_contract(pipeline$chain_base, task_type)
  stages <- pipeline$selector_params$stages
  stages <- sb_validate_chain_stages(base, stages)
  fields <- setdiff(sb_pipeline_execution_fields(task_type), c("id", "selector", "selector_params"))
  if (!identical(pipeline$id, sb_chain_pipeline_id(base, stages)) ||
    !all(vapply(fields, function(field) {
      if (field == "learner_params") {
        sb_pipeline_parameters_equal(pipeline[[field]], base[[field]])
      } else {
        identical(pipeline[[field]], base[[field]])
      }
    }, logical(1)))) {
    stop("Selector chain does not match its executable specification.")
  }
  pipeline$chain_base <- base
  pipeline$learner_params <- base$learner_params
  pipeline$selector_params$stages <- stages
  pipeline
}

sb_chain_selector_adapter <- function(params, capture_history = FALSE) {
  adapters <- lapply(params$stages, function(p) {
    if (capture_history) sb_captured_selector_adapter(p$selector, p$selector_params)
    else sb_resolve_selector_adapter(p$selector, p$selector_params)
  })
  list(id = "ordered_chain", params = params, needs_outcome = TRUE,
    fit = function(x, endpoint) {
      if (endpoint$task_type == "survival" && anyDuplicated(endpoint$group_ids) &&
          (any(vapply(params$stages, function(p) p$selector == "coxboost_nonzero", logical(1))) ||
           params$stages[[1L]]$learner %in% sb_analysis_tuned_survival_learner_ids())) {
        stop("This native survival CV engine cannot keep repeated patients together; choose a group-aware learner/selector.")
      }
      original <- names(x)
      source <- x
      states <- vector("list", length(adapters))
      for (i in seq_along(adapters)) {
        fit_stage <- function() {
          if (i == 1L || adapters[[i]]$id != "svm_rfe_binary") {
            return(sb_fit_runtime_selector(adapters[[i]], x, endpoint))
          }
          prefix <- sb_chain_selector_adapter(list(stages = params$stages[seq_len(i - 1L)]))
          prepare <- function(rows) {
            inner <- list(task_type = "binary", outcome = endpoint$outcome[rows],
              group_ids = endpoint$group_ids[rows])
            fitted <- prefix$fit(source[rows, , drop = FALSE], inner)
            list(train = prefix$apply(fitted, source[rows, , drop = FALSE]),
              test = prefix$apply(fitted, source[!rows, , drop = FALSE]))
          }
          controls <- adapters[[i]]$params
          stage <- adapters[[i]]
          # RFE's inner folds need the original input and the complete prefix,
          # not just the features selected on this outer analysis set.
          stage$upstream_context <- list(x = source,
            prefix = lapply(params$stages[seq_len(i - 1L)], function(p) {
              list(selector = p$selector, params = p$selector_params)
            }))
          stage$fit <- function(x, outcome, group_ids = NULL) {
            fitted <- sb_fit_svm_rfe(x, outcome, group_ids,
              controls$inner_folds, controls$seed, prepare_fold = prepare)
            fitted$params <- controls
            fitted
          }
          sb_fit_runtime_selector(stage, x, endpoint)
        }
        states[[i]] <- tryCatch(fit_stage(),
          error = function(e) sb_abort(paste0("Selector stage ", i, " (", adapters[[i]]$id, "): ", conditionMessage(e)),
            "stormbreaker_chain_stage_failed", parent = e))
        x <- adapters[[i]]$apply(states[[i]], x)
        if (!ncol(x)) sb_abort(paste("Selector stage", i, "retained no genes."), "stormbreaker_empty_selection")
      }
      structure(list(id = "ordered_chain", source_features = original,
        selected_features = names(x), params = params, stages = states,
        scores = data.frame(feature = original, selected = original %in% names(x))),
        class = "stormbreaker_selector_fit")
    },
    apply = function(fitted, x) {
      if (!identical(names(x), fitted$source_features)) stop("Predictors differ from the saved selector chain.")
      for (i in seq_along(adapters)) x <- adapters[[i]]$apply(fitted$stages[[i]], x)
      x
    })
}
