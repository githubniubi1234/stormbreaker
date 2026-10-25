make_final_cv_reuse_task <- function() {
  withr::with_seed(8901L, {
    x <- as.data.frame(matrix(stats::rnorm(90L * 12L), nrow = 90L))
    names(x) <- sprintf("gene_%02d", seq_len(12L))
    lp <- x[[1]] - x[[2]] + 0.6 * x[[3]]
    et <- stats::rexp(90L, rate = exp(lp) / 15)
    ct <- stats::rexp(90L, rate = 1 / 25)
    d <- data.frame(patient_id = paste0("patient_", seq_len(90L)),
      time = pmax(pmin(et, ct), 0.05), event = as.integer(et <= ct), x)
    sb_surv_task(d, "patient_id", "time", "event", names(x), "month",
      prediction_horizon = 7)
  })
}

test_that("identical final CV learners share complete fitting and diagnostics", {
  skip_if_not_installed("CoxBoost")
  skip_if_not_installed("plsRcox")
  skip_if_not_installed("penAFT")
  task <- make_final_cv_reuse_task()
  native <- sb_fit_runtime_learner
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  withr::local_seed(38L)
  rng <- .Random.seed
  for (id in c("coxboost_cv", "plsrcox_cv", "penaft_gehan_elastic_net")) {
    args <- if (id == "penaft_gehan_elastic_net") list() else
      list(learner_tuning_folds = 3L, learner_seed = if (id == "coxboost_cv") 8902L else 8903L)
    pp <- list(do.call(sb_pipeline, c(list("none", id), args)),
      do.call(sb_pipeline, c(list("fixed_features", id, selector_features = task$predictors,
        feature_set_id = "all"), args)))
    run <- function(enabled) sb_with_selector_reuse(lapply(pp, function(p) {
      sb_fit_final(task, p, capture_diagnostics = TRUE)
    }), enabled = enabled, diagnostics = TRUE)
    calls <- 0L
    plain <- run(FALSE)
    expect_identical(calls, 2L, info = id)
    calls <- 0L
    shared <- run(TRUE)
    expect_identical(calls, 1L, info = id)
    expect_identical(shared$reuse$final_learner_hits, 1L, info = id)
    for (i in 1:2) {
      a <- plain$result[[i]]; b <- shared$result[[i]]
      expect_identical(b$status, "success", info = id)
      expect_equal(b, a, info = id)
      expect_identical(b$learner_state, a$learner_state, info = id)
      expect_identical(b$tuning, a$tuning, info = id)
      expect_identical(sb_predict_final(unserialize(serialize(b, NULL)),
        task$data[c("patient_id", task$predictors)], "patient_id"),
        sb_predict_final(a, task$data[c("patient_id", task$predictors)], "patient_id"), info = id)
    }
    expect_identical(.Random.seed, rng, info = id)
  }
})

test_that("stored component-search warnings prevent final model reuse", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(...) {
    calls <<- calls + 1L
    structure(list(value = 1L), stormbreaker_fit_warnings = "late-component warning")
  })
  learner <- sb_plsrcox_cv_learner_adapter()
  x <- data.frame(a = 1:6, b = 6:1)
  endpoint <- list(task_type = "survival", time = 1:6, event = rep(0:1, 3))
  z <- sb_with_selector_reuse({
    lapply(1:2, function(i) sb_fit_final_runtime_learner(
      learner, x, endpoint, letters[1:6], "none"))
  }, diagnostics = TRUE)
  expect_identical(calls, 2L)
  expect_identical(z$reuse$final_learner_hits, 0L)
  expect_identical(z$reuse$final_learner_entries, 0L)
  expect_identical(vapply(z$result, function(m) attr(m, "stormbreaker_fit_warnings"), ""),
    rep("late-component warning", 2L))
})

test_that("final CV learner reuse distinguishes seeds and tuning folds", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(...) {
    calls <<- calls + 1L
    list(value = calls)
  })
  x <- data.frame(a = 1:6, b = 6:1)
  e <- list(task_type = "survival", time = 1:6, event = rep(0:1, 3))
  for (l in list(sb_coxboost_cv_learner_adapter(), sb_plsrcox_cv_learner_adapter(),
    sb_rule_aft_learner_adapter("penaft_gehan_elastic_net"))) {
    calls <- 0L
    run <- function(learner = l) sb_fit_final_runtime_learner(learner, x, e, letters[1:6], "none")
    z <- sb_with_selector_reuse({
      first <- run()
      expect_identical(run(), first)
      changed <- l; changed$params$seed <- changed$params$seed + 1L; run(changed)
      changed <- l
      if (l$id == "penaft_gehan_elastic_net") changed$params$internal_folds <- 3L else
        changed$params$inner_folds <- 3L
      run(changed)
    }, diagnostics = TRUE)
    expect_identical(calls, 3L)
    expect_identical(z$reuse$final_learner_hits, 1L)
    expect_identical(z$reuse$final_learner_entries, 3L)
  }
})

test_that("native PLS search warnings remain on every unshared final route", {
  skip_if_not_installed("plsRcox")
  task <- make_final_cv_reuse_task()
  message <- "Loglik converged before variable  3 ; coefficient may be infinite. "
  cv_calls <- 0L
  local_mocked_bindings(
    cv.plsRcox = function(...) {
      cv_calls <<- cv_calls + 1L
      warning(message)
      list(lambda.min5 = 2L, folds = split(seq_len(90L), rep(1:3, each = 30L)),
        cv.error5 = seq(0, 1, length.out = 11L), cv.se5 = rep(0.01, 11L),
        completed.cv5 = matrix(TRUE, 11L, 3L))
    },
    plsRcox = function(Xplan, nt, ...) {
      list(Coeffs = matrix(0.1, nrow = ncol(Xplan), ncol = 1L,
        dimnames = list(names(Xplan), NULL)), computed_nt = nt)
    }, .package = "plsRcox")
  pp <- list(sb_pipeline("none", "plsrcox_cv", learner_tuning_folds = 3L),
    sb_pipeline("fixed_features", "plsrcox_cv", learner_tuning_folds = 3L,
      selector_features = task$predictors, feature_set_id = "all"))
  run <- function(enabled) sb_with_selector_reuse(
    lapply(pp, function(p) sb_fit_final(task, p)), enabled = enabled, diagnostics = TRUE)
  plain <- run(FALSE)
  cv_calls <- 0L
  shared <- run(TRUE)
  expect_identical(cv_calls, 2L)
  expect_identical(shared$reuse$final_learner_entries, 0L)
  expect_identical(shared$reuse$final_learner_hits, 0L)
  expect_equal(shared$result, plain$result)
  for (fit in shared$result) {
    expect_identical(fit$status, "success")
    expect_identical(fit$tuning$warnings, paste0("Learner refit: ", message))
  }
})

test_that("final tuned learners preserve absent RNG and later diagnostic requests", {
  skip_if_not_installed("CoxBoost")
  skip_if_not_installed("plsRcox")
  skip_if_not_installed("penAFT")
  task <- make_final_cv_reuse_task()
  withr::local_seed(39L)
  rm(".Random.seed", envir = .GlobalEnv)
  for (id in c("coxboost_cv", "plsrcox_cv", "penaft_gehan_elastic_net")) {
    args <- if (id == "penaft_gehan_elastic_net") list() else
      list(learner_tuning_folds = 3L, learner_seed = if (id == "coxboost_cv") 8902L else 8903L)
    p <- do.call(sb_pipeline, c(list("none", id), args))
    q <- do.call(sb_pipeline, c(list("fixed_features", id,
      selector_features = task$predictors, feature_set_id = "all"), args))
    expected <- sb_fit_final(task, q, capture_diagnostics = TRUE)
    actual <- sb_with_selector_reuse({
      first <- sb_fit_final(task, p)
      expect_null(first$diagnostic_artifacts)
      sb_fit_final(task, q, capture_diagnostics = TRUE)
    }, diagnostics = TRUE)
    expect_identical(actual$result$status, "success", info = id)
    expect_identical(actual$reuse$final_learner_hits, 1L, info = id)
    expect_equal(actual$result, expected, info = id)
    expect_identical(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), FALSE,
      info = id)
  }
})
