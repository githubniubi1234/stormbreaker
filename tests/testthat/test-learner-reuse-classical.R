test_that("classical binary final fits retain native results under reuse", {
  d <- sb_example_binary(120L, 1208L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  native <- sb_fit_runtime_learner
  calls <- 0L
  native_calibration <- sb_calibrate_svm_probability
  calibrations <- 0L
  local_mocked_bindings(sb_calibrate_svm_probability = function(...) {
    calibrations <<- calibrations + 1L
    native_calibration(...)
  })
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls <<- calls + 1L
    native(learner, x, endpoint)
  })
  for (id in sb_classical_binary_learner_ids()) {
    pipelines <- list(sb_pipeline("none", id), sb_pipeline("fixed_features", id,
      selector_features = features, feature_set_id = "all"))
    run <- function(enabled) withr::with_seed(101L, {
      value <- sb_with_selector_reuse(lapply(pipelines, function(p) {
        sb_fit_final(task, p, capture_diagnostics = TRUE)
      }), enabled = enabled, diagnostics = TRUE)
      list(value = value, rng = .Random.seed)
    })
    calls <- 0L; calibrations <- 0L; plain <- run(FALSE)
    plain_calls <- calls
    plain_calibrations <- calibrations
    calls <- 0L; calibrations <- 0L; shared <- run(TRUE)
    comparable <- function(value) {
      if (id %in% c("c50_binary", "c50_rules_binary")) {
        clean_clock <- function(engine) {
          for (field in intersect(c("names", "output"), names(engine))) {
            engine[[field]] <- gsub("(?m)^\\| on [^\\n]+$", "| on <native clock>",
              engine[[field]], perl = TRUE)
            engine[[field]] <- gsub("(?m)^(C5\\.0 \\[Release [^\\]]+\\][ \\t]+)[^\\n]+$",
              "\\1<native clock>", engine[[field]], perl = TRUE)
          }
          engine
        }
        value <- lapply(value, function(fit) {
          fit$learner_state$engine_raw <- clean_clock(unserialize(fit$learner_state$engine_raw))
          fit$diagnostic_artifacts$learner$engine <- clean_clock(fit$diagnostic_artifacts$learner$engine)
          fit
        })
      }
      # A serialization round trip changes environment reference encoding,
      # not the restored model; keep full engine and exact prediction checks.
      if (id %in% c("gamboost_binary", "blackboost_binary",
          "gaussian_process_binary", "neuralnet_multilayer_binary",
          "ctree_binary", "cforest_binary", "loclda_binary")) {
        value <- lapply(value, function(fit) {
          if (is.raw(fit$learner_state$engine_raw)) {
            fit$learner_state$engine_raw <- unserialize(fit$learner_state$engine_raw)
          }
          fit
        })
      }
      value
    }
    expect_equal(comparable(shared$value$result), comparable(plain$value$result), info = id)
    expect_identical(shared$rng, plain$rng, info = id)
    if (all(vapply(plain$value$result, function(f) f$status == "success", logical(1)))) {
      expect_identical(shared$value$reuse$final_learner_hits, 1L, info = id)
      expect_identical(calls, plain_calls - 1L, info = id)
      if (id %in% sb_svm_probability_learner_ids()) {
        expect_identical(plain_calibrations, 2L, info = id)
        expect_identical(calibrations, 1L, info = id)
      }
      for (i in seq_along(pipelines)) expect_identical(
        sb_predict_final(unserialize(serialize(shared$value$result[[i]], NULL)),
          d[c("patient_id", features)], "patient_id"),
        sb_predict_final(plain$value$result[[i]], d[c("patient_id", features)], "patient_id"),
        info = id)
    } else {
      expect_identical(shared$value$reuse$final_learner_hits, 0L, info = id)
      expect_identical(calls, plain_calls, info = id)
    }
  }
})

test_that("classical binary CV reuse preserves held-out predictions and failures", {
  d <- sb_example_binary(120L, 1208L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  native <- sb_fit_runtime_learner
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, x, endpoint) {
    calls <<- calls + 1L
    native(learner, x, endpoint)
  })
  for (id in sb_classical_binary_learner_ids()) {
    pipelines <- list(sb_pipeline("none", id), sb_pipeline("fixed_features", id,
      selector_features = features, feature_set_id = "all"))
    run <- function(enabled) withr::with_seed(101L, {
      warnings <- character()
      value <- withCallingHandlers(sb_compare_cv(task, pipelines, folds = 3L,
        seed = 11L, reuse_selectors = enabled), warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      list(value = value, warnings = warnings, rng = .Random.seed)
    })
    calls <- 0L; plain <- run(FALSE); plain_calls <- calls
    calls <- 0L; shared <- run(TRUE)
    expect_equal(shared$value, plain$value, info = id)
    for (i in seq_along(pipelines)) expect_identical(
      shared$value$results[[i]]$predictions, plain$value$results[[i]]$predictions,
      info = id)
    expect_identical(shared$warnings, plain$warnings, info = id)
    expect_identical(shared$rng, plain$rng, info = id)
    expect_lte(calls, plain_calls)
    if (all(shared$value$comparison$status == "success")) {
      expect_identical(calls, 3L, info = id)
    }
  }
})

test_that("SVM fit reuse includes its inner calibration coefficients", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(100L, 1208L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  endpoint <- list(task_type = "binary", outcome = d$outcome)
  adapter <- sb_classical_binary_learner_adapter("svm_radial_binary")
  first <- sb_svm_calibrated_runtime_learner(adapter, list(A = -1, B = .1, positive_class = "1"))
  second <- sb_svm_calibrated_runtime_learner(adapter, list(A = -2, B = .2, positive_class = "1"))
  out <- sb_with_selector_reuse({
    a <- sb_fit_reusable_runtime_learner(first, x, endpoint, d$patient_id, "none")
    b <- sb_fit_reusable_runtime_learner(second, x, endpoint, d$patient_id, "none")
    expect_false(identical(first$predict(a, x), second$predict(b, x)))
    sb_fit_reusable_runtime_learner(first, x, endpoint, d$patient_id, "none")
  }, diagnostics = TRUE)
  expect_identical(out$reuse$cv_learner_misses, 2L)
  expect_identical(out$reuse$cv_learner_hits, 1L)
})
