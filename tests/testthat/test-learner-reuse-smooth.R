test_that("smooth and boosting families preserve CV and final predictions under reuse", {
  for (package in c("mboost", "mgcv", "partykit")) skip_if_not_installed(package)
  ids <- c(sb_mboost_cox_learner_ids(), sb_mboost_aft_learner_ids(),
    sb_mboost_ranking_learner_ids(), sb_mgcv_gam_learner_ids())
  features <- c("gene_a", "gene_b", "gene_c")
  native <- sb_fit_runtime_learner
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls <<- calls + 1L; native(learner, x, endpoint)
  })
  for (id in ids) {
    if (id == "gam_binary_reml") {
      d <- sb_example_binary(120L, 1208L)
      task <- sb_binary_task(d, "patient_id", "outcome", features)
    } else {
      d <- sb_example_survival(120L, 11032L)
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
      # Serializing an isolated mboost environment can change reference encoding,
      # not the model. Compare restored engines rather than their raw bytes.
      comparable <- function(value) {
        if (final) value <- lapply(value, function(fit) {
          if (is.raw(fit$learner_state$engine_raw)) {
            fit$learner_state$engine_raw <- unserialize(fit$learner_state$engine_raw)
          }
          fit
        })
        value
      }
      expect_equal(comparable(shared$value), comparable(plain$value), info = id)
      expect_identical(shared$warnings, plain$warnings, info = id)
      expect_identical(shared$rng, plain$rng, info = id)
      expected_calls <- if (id == "blackboost_cox") nplain else if (final) 1L else 3L
      expect_identical(calls, expected_calls, info = id)
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
})

test_that("native RNG-dependent fits reuse only the same starting random state", {
  skip_if_not_installed("mboost")
  skip_if_not_installed("partykit")
  d <- sb_example_survival(80L, 11032L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "survival", time = d$time, event = d$event)
  learner <- sb_mboost_cox_learner_adapter("blackboost_cox")
  native <- withr::with_seed(55L, {
    model <- sb_fit_runtime_learner(learner, x, endpoint)
    list(model = model, rng = .Random.seed)
  })
  out <- sb_with_selector_reuse({
    run <- function() withr::with_seed(55L, {
      model <- sb_fit_reusable_runtime_learner(learner, x, endpoint, d$patient_id, "none")
      list(model = model, rng = .Random.seed)
    })
    first <- run(); second <- run()
    expect_equal(first$model, native$model)
    expect_equal(second$model, native$model)
    expect_identical(first$rng, native$rng)
    expect_identical(second$rng, native$rng)
  }, diagnostics = TRUE)
  expect_identical(out$reuse$cv_learner_hits, 1L)
  expect_identical(out$reuse$cv_learner_misses, 1L)
})
