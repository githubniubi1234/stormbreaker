.sb_selector_reuse_context <- new.env(parent = emptyenv())
.sb_selector_reuse_context$active <- NULL

#' Reuse identical selector and eligible learner fits in one expression
#'
#' This helper shares successful, warning-free Boruta binary/survival and
#' varSelRF, the high-dimensional omics selectors, ABESS, PAMR, binary model
#' selectors, binary AUC/survival Cox Top-K, survival Lasso/Elastic Net,
#' Ridge/PLS-Cox/survival-SVM weight
#' selection, CoxBoost, RSF minimal-depth, StepCox and survival boosting fits when the
#' complete predictor frame (including row and column identities), endpoint,
#' selector settings (including requested diagnostic history) and RNG kind match.
#' Eligible seeded selectors preserve the caller's random state;
#' deterministic StepCox does not use random draws.
#' Sequential MIME selectors use the same mechanism for their individual
#' stages, allowing identical prefixes to be shared with a single stage.
#' A non-first SVM-RFE stage also matches the original analysis predictors and
#' upstream selector specifications, retaining prefix refitting in its inner folds.
#' Other selectors execute normally.
#' No input or fitted state is persisted.
#'
#' Reuse does not preselect genes on all patients. Different outer-analysis,
#' inner-analysis and final-training inputs receive separate fits. Final fits
#' may also reuse fixed survival-SVM, logistic, Cox, Gaussian naive Bayes and
#' prepared glmnet fits, fixed StepGLM/GLMBoost/PLS-GLM and the classical binary
#' family (including trees, forests, discriminant, kernel and prototype models)
#' and mboost Cox/AFT/ranking or mgcv GAM fits after any selector.
#' Analysis-CV learners reuse their complete search and fit after any legally
#' registered selector; this does not create additional combinations. Models are reused only
#' when ordered predictors, patient IDs, endpoint and learner settings match.
#' Tuned-learner reuse includes the complete native search and selected fit,
#' not a search performed on other data.
#' Matching final inputs are retained within the expression, including when
#' another input for the same learner is fitted between two matching calls.
#' CV and final learner retention share a 256 MiB per-expression payload budget.
#' Stored models are serialized; entry sizes include matching inputs and metadata.
#' Least-recently-used entries are released when needed. An oversized model is
#' returned normally but not retained. Evicted inputs are fitted normally if
#' needed again. This bounds retained learner entries, not total R memory:
#' active fitting, returned results, selectors and preparation have other costs.
#' Each route still constructs its own selection, portable state, horizon
#' baseline, diagnostics and provenance. CV fits of fixed survival-SVM, logistic,
#' Cox, Gaussian naive Bayes, prepared glmnet, fixed StepGLM/GLMBoost/PLS-GLM
#' and the classical binary, mboost Cox/AFT/ranking and mgcv GAM families
#' use a separate in-memory collection of matching
#' analysis inputs; different analysis sets are fitted separately. Other CV learners
#' execute normally; tuning is never transferred to different analysis data.
#' Glmnet and SVM probability preparation also share identical native searches, matching learner
#' controls, patients, endpoints, groups and full upstream selection context.
#' Identical outer-selected columns alone never identify a learned-selector search.
#' Each consumer retains its own inner-selection method description.
#' SVM model matching additionally includes its fitted probability calibration.
#' Native BlackBoost-Cox, RSF and AdaBoost fits also match the caller's initial random state and
#' replay its random-state transition; different starting states are not shared.
#' Mutable mboost model environments are isolated by serialization, so changing
#' one consumer's fitted model cannot change another consumer's model.
#' SurvReg, ranger, flexsurv/spline, timereg, ahaz, Firth Cox, ncvreg and ABESS fits,
#' omics-specialist binary, mclust, mixOmics and HDDA fits also share identical
#' analysis inputs in CV and final fitting. Their native searches, when present,
#' remain confined to those exact analysis data.
#' CatBoost, fixed GBM/CoxBoost/PLS-Cox/SuperPC, MTLR, PAMR survival,
#' Akritas, bagged survival kNN, Bujar, RLT, WSRF, AdaBoost/LogitBoost and
#' rpart/RSF/ctree/cforest/AORSF survival fits are also eligible. Eligibility
#' alone does not guarantee a hit: inputs and required random states must match.
#' SLOPE, randomGLM and cancerclass may share the complete analysis-local fit,
#' including native searches or out-of-fold calibration. This does not add new
#' selector combinations or transfer internal calibration to other patients.
#' BART survival/AFT, PRE, penAFT, CV-tuned CoxBoost/PLS-Cox/SuperPC/GBM-Cox
#' also match their complete analysis inputs. Failed or warning-bearing fits
#' are not retained; every affected route preserves its own failure reporting.
#' XGBoost and LightGBM fits use native model serialization for independent
#' booster handles. LightGBM retains original fit importance values because its
#' text model rounds split gains; predictions and reported importance stay unchanged.
#' Errors are never cached. Ordinary selector fits with warnings are not cached.
#' Final learner fits with warnings are also executed again, including warnings
#' retained as fit metadata from unselected PLS-Cox components.
#' The context is discarded on normal
#' return or error; nested calls have independent contexts.
#'
#' @param expr An expression containing repeated comparisons, CV runs or final
#'   fits. Do not change package implementations during the expression.
#' @param enabled Logical; FALSE executes the expression without reuse, even
#'   inside an enabled outer context.
#' @param diagnostics Logical; TRUE returns a list with `result` and `reuse`
#'   (selector hits, misses and entries, plus separately named
#'   `final_learner_hits`, `final_learner_misses` and `final_learner_entries`);
#'   analogous `cv_learner_*` fields describe CV learner reuse;
#'   `learner_bytes`, `learner_byte_limit`, `learner_evictions` and
#'   `learner_oversized` report the bounded learner collection;
#'   `preparation_hits` and `preparation_misses` describe native tuning reuse;
#'   FALSE returns the expression value unchanged.
#' @return The value of `expr`, or a result/diagnostics list.
#' @examples
#' if (requireNamespace("Boruta", quietly = TRUE) &&
#'     requireNamespace("ranger", quietly = TRUE)) {
#'   d <- sb_example_binary(n = 100, seed = 42)
#'   task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
#'   shared <- sb_with_selector_reuse({
#'     lapply(c("logistic", "gaussian_naive_bayes"), function(learner) {
#'       pipeline <- sb_pipeline("boruta_binary_roughfix", learner)
#'       sb_run_cv(task, pipeline, folds = 3, seed = 42)
#'     })
#'   }, diagnostics = TRUE)
#'   shared$reuse
#' }
#' @export
sb_with_selector_reuse <- function(expr, enabled = TRUE, diagnostics = FALSE) {
  if (!is.logical(enabled) || length(enabled) != 1L || is.na(enabled) ||
      !is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics)) {
    sb_abort("enabled and diagnostics must be single non-missing logical values.",
      "stormbreaker_bad_selector_reuse")
  }
  previous <- .sb_selector_reuse_context$active
  on.exit(.sb_selector_reuse_context$active <- previous, add = TRUE)
  context <- new.env(parent = emptyenv())
  context$entries <- new.env(parent = emptyenv())
  context$hits <- 0L
  context$misses <- 0L
  context$n_entries <- 0L
  context$final_learners <- new.env(parent = emptyenv())
  context$final_learner_hits <- 0L
  context$final_learner_misses <- 0L
  context$cv_learners <- new.env(parent = emptyenv())
  context$cv_learner_hits <- 0L
  context$cv_learner_misses <- 0L
  context$learner_byte_limit <- 256 * 1024^2
  context$learner_bytes <- 0
  context$learner_clock <- 0
  context$learner_evictions <- 0L
  context$learner_oversized <- 0L
  context$preparations <- new.env(parent = emptyenv())
  context$preparation_hits <- 0L
  context$preparation_misses <- 0L
  .sb_selector_reuse_context$active <- if (enabled) context else NULL
  value <- eval(substitute(expr), envir = parent.frame())
  if (diagnostics) {
    return(list(result = value, reuse = list(enabled = enabled,
      hits = context$hits, misses = context$misses, entries = context$n_entries,
      final_learner_hits = context$final_learner_hits,
      final_learner_misses = context$final_learner_misses,
      final_learner_entries = sum(lengths(as.list(context$final_learners))),
      cv_learner_hits = context$cv_learner_hits,
      cv_learner_misses = context$cv_learner_misses,
      cv_learner_entries = sum(lengths(as.list(context$cv_learners))),
      learner_bytes = context$learner_bytes,
      learner_byte_limit = context$learner_byte_limit,
      learner_evictions = context$learner_evictions,
      learner_oversized = context$learner_oversized,
      preparation_hits = context$preparation_hits,
      preparation_misses = context$preparation_misses)))
  }
  value
}

# Internal handoff for bounded Boruta workers. The controller owns persistence;
# these objects contain analysis predictors/endpoints and are not portable models.
# Consumers still go through the ordinary exact-input/signature matcher.
sb_boruta_worker_checkpoint <- function(id) {
  context <- .sb_selector_reuse_context$active
  if (!id %in% c("boruta_binary_roughfix", "boruta_survival_roughfix") ||
      is.null(context)) {
    stop("Boruta worker handoff requires an active reuse context and Boruta selector.")
  }
  list(id = id, entries = context$entries[[id]])
}

sb_restore_boruta_worker_checkpoint <- function(checkpoint) {
  context <- .sb_selector_reuse_context$active
  if (is.null(context) || !checkpoint$id %in%
      c("boruta_binary_roughfix", "boruta_survival_roughfix")) {
    stop("Boruta worker handoff requires an active reuse context and Boruta selector.")
  }
  previous <- length(context$entries[[checkpoint$id]])
  context$entries[[checkpoint$id]] <- checkpoint$entries
  context$n_entries <- context$n_entries - previous + length(checkpoint$entries)
  invisible(NULL)
}

sb_fit_final_runtime_learner <- function(learner, x, endpoint, patient_ids,
                                       selector_id) {
  sb_fit_reusable_runtime_learner(learner, x, endpoint, patient_ids,
    selector_id, phase = "final")
}

sb_fit_reusable_runtime_learner <- function(learner, x, endpoint, patient_ids,
                                           selector_id, phase = "cv") {
  context <- .sb_selector_reuse_context$active
  eligible <- c(sb_survivalsvm_learner_ids(), "logistic", "coxph", "gaussian_naive_bayes",
    "glmnet_logistic", "glmnet_cox", sb_specialized_binary_linear_learner_ids(),
    sb_classical_binary_learner_ids(), sb_mboost_cox_learner_ids(),
    sb_mboost_aft_learner_ids(), sb_mboost_ranking_learner_ids(),
    sb_mgcv_gam_learner_ids(), sb_flexsurv_parametric_learner_ids(),
    sb_flexsurvspline_learner_ids(), sb_timereg_survival_learner_ids(),
    sb_ahaz_survival_learner_ids(), sb_coxphf_learner_ids(),
    sb_ncvreg_learner_ids(), sb_abess_learner_ids(),
    sb_omics_specialist_binary_learner_ids(), sb_mclust_da_learner_ids(),
    sb_mixomics_da_learner_ids(), sb_hdda_binary_learner_ids(),
    sb_ranger_learner_ids(), sb_survreg_learner_ids(),
    sb_xgboost_learner_ids(), sb_lightgbm_binary_learner_ids(),
    "catboost_binary", "catboost_cox", "catboost_aft_survival",
    "gbm_cox", "coxboost", "plsrcox", "superpc",
    sb_mtlr_survival_learner_ids(), sb_pamr_survival_learner_ids(),
    sb_akritas_learner_ids(), sb_bagged_knn_survival_learner_ids(),
    sb_bujar_linear_learner_ids(), sb_rlt_learner_ids(), sb_wsrf_binary_learner_ids(),
    "adaboost_binary", "logitboost_binary", "rpart_survival",
    "rsf_survival", "ctree_survival", "cforest_survival",
    "aorsf_survival", "slope_binary", "random_glm_binary",
    "cancerclass_binary", "bart_survival", "bart_aft_survival",
    "coxboost_cv", "plsrcox_cv", "superpc_cv", "gbm_cox_cv",
    "pre_binary_rule_ensemble", "pre_survival_rule_ensemble",
    "penaft_gehan_elastic_net")
  if (is.null(context) || !learner$id %in% eligible) {
    return(sb_fit_runtime_learner(learner, x, endpoint))
  }
  # These native fits consume caller RNG (including AdaBoost's rpart CV).
  # Matching its starting state and replaying the transition preserves both
  # the model/diagnostics and the random stream of subsequent computations.
  rng_sensitive <- sb_learner_uses_caller_rng(learner$id)
  rng_input <- if (rng_sensitive && exists(".Random.seed", .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", .GlobalEnv, inherits = FALSE)
  } else NULL
  signature <- list(id = learner$id, params = learner$params,
    fit_context = learner$fit_context,
    rng_input = rng_input,
    fit_body = body(learner$fit), fit_formals = formals(learner$fit),
    rng_kind = RNGkind())
  store <- context[[paste0(phase, "_learners")]]
  hits <- paste0(phase, "_learner_hits")
  misses <- paste0(phase, "_learner_misses")
  entries <- store[[learner$id]]
  for (index in seq_along(entries)) {
    entry <- entries[[index]]
    if (identical(signature, entry$signature) &&
      identical(x, entry$x) && identical(endpoint, entry$endpoint) &&
      identical(patient_ids, entry$patient_ids)) {
      context[[hits]] <- context[[hits]] + 1L
      context$learner_clock <- context$learner_clock + 1
      entries[[index]]$last_used <- context$learner_clock
      store[[learner$id]] <- entries
      if (rng_sensitive) assign(".Random.seed", entry$rng_after, .GlobalEnv)
      return(sb_restore_reused_learner_model(unserialize(entry$model), learner$id))
    }
  }
  context[[misses]] <- context[[misses]] + 1L
  warned <- FALSE
  model <- withCallingHandlers(sb_fit_runtime_learner(learner, x, endpoint),
    warning = function(w) warned <<- TRUE)
  if (!warned && length(attr(model, "stormbreaker_fit_warnings", exact = TRUE)) == 0L &&
      (!rng_sensitive || !is.null(rng_input))) {
    entry <- list(signature = signature,
      x = x, endpoint = endpoint, patient_ids = patient_ids,
      model = serialize(sb_store_reused_learner_model(model, learner$id), NULL),
      rng_after = if (rng_sensitive) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL)
    sb_retain_shared_learner(context, phase, learner$id, entry)
  }
  model
}

sb_retain_shared_learner <- function(context, phase, id, entry) {
  context$learner_clock <- context$learner_clock + 1
  entry$last_used <- context$learner_clock
  entry$bytes <- 0
  entry$bytes <- as.numeric(utils::object.size(entry))
  if (entry$bytes > context$learner_byte_limit) {
    context$learner_oversized <- context$learner_oversized + 1L
    return(invisible(NULL))
  }
  while (context$learner_bytes + entry$bytes > context$learner_byte_limit) {
    oldest <- Inf
    for (name in c("cv_learners", "final_learners")) {
      store <- context[[name]]
      for (key in ls(store, all.names = TRUE)) {
        entries <- store[[key]]
        for (i in seq_along(entries)) if (entries[[i]]$last_used < oldest) {
          oldest <- entries[[i]]$last_used
          victim <- list(store = name, key = key, index = i)
        }
      }
    }
    store <- context[[victim$store]]
    entries <- store[[victim$key]]
    context$learner_bytes <- context$learner_bytes - entries[[victim$index]]$bytes
    store[[victim$key]] <- entries[-victim$index]
    context$learner_evictions <- context$learner_evictions + 1L
  }
  store <- context[[paste0(phase, "_learners")]]
  store[[id]] <- c(store[[id]], list(entry))
  context$learner_bytes <- context$learner_bytes + entry$bytes
  invisible(NULL)
}

sb_store_reused_learner_model <- function(model, learner) {
  if (learner %in% sb_lightgbm_binary_learner_ids()) {
    text <- model$booster$save_model_to_string()
    fields <- intersect(c("best_iter", "best_score", "params", "record_evals",
      "data_processor", "raw"), names(model$booster))
    metadata <- stats::setNames(lapply(fields, function(name) model$booster[[name]]), fields)
    # LightGBM's text format rounds split gains, although predictions are exact.
    # Retain the original diagnostic values alongside the independent booster.
    model$reused_native_importance <- list(model_string = text,
      observed = as.data.frame(lightgbm::lgb.importance(model$booster, percentage = FALSE),
        stringsAsFactors = FALSE))
    model$booster <- NULL
    return(list(format = "lightgbm", value = model, text = text, metadata = metadata))
  }
  if (learner %in% sb_xgboost_learner_ids()) {
    raw <- xgboost::xgb.save.raw(model$booster, raw_format = "json")
    # Older XGBoost keeps R-side call/iteration metadata alongside the pointer.
    metadata <- if (is.list(model$booster)) {
      model$booster[setdiff(names(model$booster), "handle")]
    } else NULL
    model$booster <- NULL
    return(list(format = "xgboost", value = model, raw = raw, metadata = metadata))
  }
  if (learner %in% c(sb_mboost_cox_learner_ids(),
      sb_mboost_aft_learner_ids(), sb_mboost_ranking_learner_ids(),
      "gamboost_binary", "blackboost_binary", "catboost_binary",
      "catboost_cox", "catboost_aft_survival")) {
    # mstop<- mutates native environments: every consumer needs an isolated fit.
    return(list(format = "serialized", value = serialize(model, NULL)))
  }
  list(format = "value", value = model)
}

sb_restore_reused_learner_model <- function(stored, learner) {
  if (stored$format == "serialized") return(unserialize(stored$value))
  if (stored$format == "lightgbm") {
    model <- stored$value
    model$booster <- lightgbm::lgb.load(model_str = stored$text)
    for (name in names(stored$metadata)) model$booster[[name]] <- stored$metadata[[name]]
    return(model)
  }
  if (stored$format == "xgboost") {
    model <- stored$value
    model$booster <- sb_load_xgboost_state_model(list(model_raw = stored$raw,
      params = model$params, features = model$features), learner)$booster
    if (is.list(model$booster) && !is.null(stored$metadata)) {
      model$booster[names(stored$metadata)] <- stored$metadata
    }
    return(model)
  }
  stored$value
}

sb_fit_runtime_selector <- function(selector, x, endpoint) {
  context <- .sb_selector_reuse_context$active
  if (is.null(context) ||
      !selector$id %in% c(sb_forest_selector_ids(), "stepcox", "lasso_nonzero",
        "univ_cox_top_k", "univ_binary_auc_top_k",
        sb_high_dimensional_selector_ids(), sb_pamr_selector_ids(),
        sb_binary_model_selector_ids(), sb_abess_selector_ids(),
        "enet_cox_nonzero", sb_survival_weight_selector_ids(),
        "coxboost_nonzero", "rsf_minimal_depth", sb_survival_boosting_selector_ids())) {
    return(sb_fit_runtime_selector_uncached(selector, x, endpoint))
  }
  # Source-reference environments change identity after worker serialization;
  # executable body/formals, rather than debugger metadata, identify the fit.
  fit_definition <- utils::removeSource(selector$fit)
  signature <- list(id = selector$id, adapter_params = selector$params,
    capture_history = selector$capture_history,
    upstream_context = selector$upstream_context,
    method_params = if (selector$id %in% sb_forest_selector_ids()) {
      sb_forest_selector_params(selector$id)
    } else NULL,
    needs_outcome = selector$needs_outcome,
    fit_body = body(fit_definition), fit_formals = formals(fit_definition),
    rng_kind = RNGkind())
  entries <- context$entries[[selector$id]]
  for (entry in entries) {
    if (identical(signature, entry$signature) && identical(x, entry$x) &&
        identical(endpoint, entry$endpoint)) {
      context$hits <- context$hits + 1L
      return(entry$value)
    }
  }
  context$misses <- context$misses + 1L
  warned <- FALSE
  value <- withCallingHandlers(
    sb_fit_runtime_selector_uncached(selector, x, endpoint),
    warning = function(w) warned <<- TRUE
  )
  if (!warned) {
    entries[[length(entries) + 1L]] <- list(signature = signature,
      x = x, endpoint = endpoint, value = value)
    context$entries[[selector$id]] <- entries
    context$n_entries <- context$n_entries + 1L
  }
  value
}
