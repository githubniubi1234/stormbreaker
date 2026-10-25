test_that("analysis-local searches reuse the complete fitted model", {
  for (package in c("SLOPE", "randomGLM", "cancerclass", "BART", "CoxBoost",
    "plsRcox", "superpc", "gbm", "pre", "penAFT")) skip_if_not_installed(package)
  # Load before both arms: an installed-package build-version warning is not
  # a fit warning and is emitted only on the first namespace load.
  if (!"package:gbm" %in% search()) {
    suppressWarnings(suppressPackageStartupMessages(library("gbm", character.only = TRUE)))
    withr::defer(detach("package:gbm", unload = FALSE))
  }
  features <- c("gene_a", "gene_b", "gene_c")
  calls <- 0L
  native <- sb_fit_runtime_learner
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls <<- calls + 1L
    native(learner, x, endpoint)
  })
  binary_ids <- c("random_glm_binary", "cancerclass_binary",
    "pre_binary_rule_ensemble")
  survival_ids <- c("bart_survival", "bart_aft_survival", "coxboost_cv",
    "plsrcox_cv", "superpc_cv", "gbm_cox_cv", "pre_survival_rule_ensemble",
    "penaft_gehan_elastic_net")
  for (id in c(survival_ids, binary_ids)) {
    features <- c("gene_a", "gene_b", "gene_c")
    if (id %in% binary_ids) {
      d <- sb_example_binary(120L, 1208L)
      task <- sb_binary_task(d, "patient_id", "outcome", features)
    } else {
      d <- sb_example_survival(180L, 11032L)
      task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
        prediction_horizon = 7)
    }
    pipelines <- list(sb_pipeline("none", id), sb_pipeline("fixed_features", id,
      selector_features = features, feature_set_id = "all"))
    for (final in c(FALSE, TRUE)) {
      run <- function(enabled) withr::with_seed(101L, {
        warnings <- character()
        value <- withCallingHandlers({
          if (final) sb_with_selector_reuse(lapply(pipelines,
            function(p) sb_fit_final(task, p, capture_diagnostics = TRUE)), enabled = enabled)
          else sb_compare_cv(task, pipelines, folds = 3L, seed = 11L, reuse_selectors = enabled)
        }, warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
        list(value = value, warnings = warnings, rng = .Random.seed)
      })
      calls <- 0L; plain <- run(FALSE); plain_calls <- calls
      calls <- 0L; shared <- run(TRUE)
      expect_equal(shared, plain, info = id)
      expect_identical(plain_calls, if (final) 2L else 6L, info = id)
      # This SuperPC fixture has one genuine native nonconvergence fold.
      # Both routes must retry/report that failed fit, sharing only successes.
      expected_calls <- if (final) 1L else if (id == "superpc_cv") 4L else 3L
      expect_identical(calls, expected_calls, info = id)
      if (final) {
        expect_identical(vapply(shared$value, `[[`, "", "status"), rep("success", 2L), info = id)
        for (i in 1:2) expect_identical(
          sb_predict_final(unserialize(serialize(shared$value[[i]], NULL)),
            d[c("patient_id", features)], "patient_id"),
          sb_predict_final(plain$value[[i]], d[c("patient_id", features)], "patient_id"), info = id)
      } else {
        expect_identical(unique(shared$value$comparison$status),
          if (id == "superpc_cv") "failed" else "success", info = id)
        if (id == "superpc_cv") expect_identical(
          shared$value$results[[1L]]$runs$error_class[[1L]],
          "stormbreaker_superpc_cv_warning")
      }
    }
  }
})
