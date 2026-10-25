test_that("AdaBoost reuse preserves native tree CV and caller random state", {
  skip_if_not_installed("fastAdaboost")
  d <- sb_example_binary(80L, 1208L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "binary", outcome = d$outcome)
  learner <- sb_resolve_pipeline_runtime(sb_pipeline("none", "adaboost_binary"))$learner
  for (phase in c("cv", "final")) {
    native <- withr::with_seed(55L, {
      model <- sb_fit_runtime_learner(learner, x, endpoint)
      list(model = model, rng = .Random.seed)
    })
    out <- sb_with_selector_reuse({
      run <- function() withr::with_seed(55L, {
        model <- sb_fit_reusable_runtime_learner(learner, x, endpoint,
          d$patient_id, "none", phase = phase)
        list(model = model, rng = .Random.seed)
      })
      first <- run(); second <- run()
      expect_equal(first$model, native$model)
      expect_equal(second$model, native$model)
      expect_identical(first$rng, native$rng)
      expect_identical(second$rng, native$rng)
    }, diagnostics = TRUE)
    expect_identical(out$reuse[[paste0(phase, "_learner_hits")]], 1L)
    expect_identical(out$reuse[[paste0(phase, "_learner_misses")]], 1L)
  }
})

check_parametric_learner_reuse <- function(akritas_only = FALSE) {
  dependencies <- if (akritas_only) c("survivalmodels", "distr6") else
    c("flexsurv", "timereg", "ahaz", "coxphf", "ncvreg", "abess")
  for (p in dependencies) {
    skip_if_not_installed(p)
  }
  remaining_binary <- c("catboost_binary", "rlt_reinforcement_binary",
    "adaboost_binary", "logitboost_binary", "wsrf_binary")
  remaining_survival <- c("catboost_cox", "catboost_aft_survival",
    "gbm_cox", "coxboost", "plsrcox", "superpc",
    "mtlr_survival", "pamr_survival",
    "bagged_knn_survival", "bujar_linear_boost",
    "rlt_reinforcement_survival", "rpart_survival", "rsf_survival",
    "ctree_survival", "cforest_survival", "aorsf_survival")
  extra_binary <- c(sb_omics_specialist_binary_learner_ids(), sb_mclust_da_learner_ids(),
    sb_mixomics_da_learner_ids(), sb_hdda_binary_learner_ids())
  ids <- c(remaining_binary, remaining_survival,
    sb_lightgbm_binary_learner_ids(), sb_xgboost_learner_ids(), extra_binary,
    sb_ranger_learner_ids(), sb_survreg_learner_ids(),
    sb_flexsurv_parametric_learner_ids(),
    sb_flexsurvspline_learner_ids(), sb_timereg_survival_learner_ids(),
    sb_ahaz_survival_learner_ids(), sb_coxphf_learner_ids(),
    sb_ncvreg_learner_ids(), sb_abess_learner_ids())
  if (akritas_only) ids <- "akritas_survival"
  features <- c("gene_a", "gene_b", "gene_c")
  calls <- 0L
  native <- sb_fit_runtime_learner
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls <<- calls + 1L
    native(learner, x, endpoint)
  })
  for (id in ids) {
    features <- c("gene_a", "gene_b", "gene_c")
    binary <- id %in% c(remaining_binary, sb_lightgbm_binary_learner_ids(), sb_xgboost_binary_learner_ids(),
      extra_binary, sb_ranger_binary_learner_ids(),
      sb_ncvreg_binary_learner_ids(), sb_abess_binary_learner_ids())
    if (binary) {
      d <- sb_example_binary(180L, 1208L)
      if (id %in% extra_binary) {
        d$gene_d <- withr::with_seed(702L, stats::rnorm(nrow(d)))
        features <- c(features, "gene_d")
      }
      task <- sb_binary_task(d, "patient_id", "outcome", features)
    } else {
      d <- sb_example_survival(180L, 11032L)
      if (id == "aalen_additive_survival") d <- withr::with_seed(9431L, {
        x <- matrix(stats::runif(180L * 3L), ncol = 3L)
        event_time <- stats::rexp(180L, 0.08 + drop(x %*% c(0.08, 0.05, 0.03)))
        censor_time <- stats::rexp(180L, 0.025)
        data.frame(patient_id = paste0("A", seq_len(180L)),
          time = pmax(pmin(event_time, censor_time), 0.01),
          event = as.integer(event_time <= censor_time),
          gene_a = x[, 1L], gene_b = x[, 2L], gene_c = x[, 3L])
      })
      task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
        prediction_horizon = 7)
    }
    pipelines <- list(sb_pipeline("none", id), sb_pipeline("fixed_features", id,
      selector_features = features, feature_set_id = "all"))
    run <- function(enabled, final) withr::with_seed(101L, {
      warnings <- character()
      value <- withCallingHandlers({
        if (final) sb_with_selector_reuse(lapply(pipelines,
          function(p) sb_fit_final(task, p, capture_diagnostics = TRUE)), enabled = enabled)
        else sb_compare_cv(task, pipelines, folds = 3L, seed = 11L, reuse_selectors = enabled)
      }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      })
      list(value = value, warnings = warnings, rng = .Random.seed)
    })
    for (final in c(FALSE, TRUE)) {
      calls <- 0L; plain <- run(FALSE, final); nplain <- calls
      calls <- 0L; shared <- run(TRUE, final)
      comparable <- function(value) {
        if (final && id %in% c("aorsf_survival", "ctree_survival",
            "cforest_survival")) value <- lapply(value, function(fit) {
          # Environment reference encoding can differ; compare the restored
          # native model and, below, exact portable prediction frames.
          fit$learner_state$engine_raw <- unserialize(fit$learner_state$engine_raw)
          fit
        })
        if (id %in% c("catboost_binary", "catboost_cox", "catboost_aft_survival")) {
          clean_cat <- function(x) {
            if (inherits(x, "catboost.Model")) {
              return(list(fields = x[setdiff(names(x), "cpp_obj")],
                parameters = catboost::catboost.get_model_params(x)))
            }
            if (is.list(x)) {
              y <- lapply(seq_along(x), function(i) {
                if (!is.null(names(x)) && names(x)[[i]] == "model_raw" && is.raw(x[[i]])) {
                  clean_cat(unserialize(x[[i]]))
                } else clean_cat(x[[i]])
              })
              attributes(y) <- attributes(x)
              return(y)
            }
            x
          }
          # Native CatBoost model identities are volatile between two fits;
          # compare parameters/audits and exact OOF/portable predictions instead.
          value <- clean_cat(value)
        }
        if (final && id == "rsf_survival") value <- lapply(value, function(fit) {
          clean_rsf <- function(engine) {
            engine[c("ctime.internal", "ctime.external")] <- NULL
            engine$forest$sampsize <- list(body = body(engine$forest$sampsize),
              formals = formals(engine$forest$sampsize))
            engine
          }
          fit$learner_state$engine_raw <- clean_rsf(unserialize(fit$learner_state$engine_raw))
          fit$diagnostic_artifacts$learner$engine <- clean_rsf(fit$diagnostic_artifacts$learner$engine)
          fit
        })
        if (final && id %in% sb_lightgbm_binary_learner_ids()) {
          value <- lapply(value, function(fit) {
            booster <- fit$diagnostic_artifacts$learner$booster
            fields <- intersect(c("best_iter", "best_score", "params", "record_evals",
              "data_processor", "raw"), names(booster))
            fit$diagnostic_artifacts$learner$booster <- list(
              model = booster$save_model_to_string(),
              metadata = setNames(lapply(fields, function(name) booster[[name]]), fields))
            # Internal source-precision audit storage is not a separate model
            # result; the actual public importance tables remain compared below.
            fit$diagnostic_artifacts$learner$reused_native_importance <- NULL
            fit
          })
        }
        if (final && id %in% sb_xgboost_learner_ids()) {
          value <- lapply(value, function(fit) {
            booster <- fit$diagnostic_artifacts$learner$booster
            # Native pointers differ even between two non-shared fits.
            fit$diagnostic_artifacts$learner$booster <- list(
              model = xgboost::xgb.save.raw(booster, raw_format = "json"),
              metadata = if (is.list(booster)) booster[setdiff(names(booster), "handle")] else NULL)
            fit
          })
        }
        value
      }
      expect_equal(comparable(shared$value), comparable(plain$value), info = id)
      expect_identical(shared$warnings, plain$warnings, info = id)
      expect_identical(shared$rng, plain$rng, info = id)
      expect_identical(calls, if (id %in% c("rsf_survival", "adaboost_binary")) nplain else if (final) 1L else 3L, info = id)
      expect_identical(nplain, if (final) 2L else 6L, info = id)
      if (final) {
        expect_identical(vapply(shared$value, `[[`, "", "status"), rep("success", 2L), info = id)
        for (i in 1:2) expect_identical(
          sb_predict_final(unserialize(serialize(shared$value[[i]], NULL)),
            d[c("patient_id", features)], "patient_id"),
          sb_predict_final(plain$value[[i]], d[c("patient_id", features)], "patient_id"), info = id)
      } else {
        expect_identical(unique(shared$value$comparison$status), "success", info = id)
        for (i in 1:2) expect_identical(shared$value$results[[i]]$predictions,
          plain$value$results[[i]]$predictions, info = id)
      }
    }
  }
}

test_that("parametric, penalized and discriminant learners reuse identical analysis fits", {
  check_parametric_learner_reuse()
})

test_that("Akritas reuses identical analysis fits when its dependencies are available", {
  check_parametric_learner_reuse(akritas_only = TRUE)
})

test_that("failed Akritas fits are never reused as successful models", {
  d <- sb_example_survival(48L, 42L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month")
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(...) {
    calls <<- calls + 1L
    stop("synthetic learner failure")
  })
  for (phase in c("cv", "final")) {
    calls <- 0L
    results <- sb_with_selector_reuse(lapply(1:2, function(i) {
      p <- sb_pipeline("none", "akritas_survival")
      if (phase == "final") sb_fit_final(task, p) else
        sb_run_cv(task, p, folds = 3L, seed = 82L)
    }))
    expect_identical(calls, if (phase == "final") 2L else 6L)
    for (result in results) {
      if (phase == "final") expect_identical(result$status, "failed") else {
        expect_identical(result$complete_oof, FALSE)
        expect_true(all(result$runs$status == "failed"))
      }
    }
  }
})
