test_that("identical final GBM and RLT routes share fitting without changing outputs", {
  skip_if_not_installed("gbm")
  skip_if_not_installed("RLT")
  d <- sb_example_survival(90L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
    prediction_horizon = 7)
  gbm_fit <- sb_fit_gbm_cox_cv
  rlt_fit <- sb_fit_rlt_learner
  calls <- 0L
  local_mocked_bindings(
    sb_fit_gbm_cox_cv = function(...) { calls <<- calls + 1L; gbm_fit(...) },
    sb_fit_rlt_learner = function(...) { calls <<- calls + 1L; rlt_fit(...) })
  withr::local_seed(35L)
  before <- .Random.seed
  for (id in c("gbm_cox_cv", "rlt_reinforcement_survival")) {
    args <- if (id == "gbm_cox_cv") list(learner_tuning_folds = 3L,
      gbm_cv_max_trees = 30L) else list()
    pipelines <- list(do.call(sb_pipeline, c(list("none", id), args)),
      do.call(sb_pipeline, c(list("fixed_features", id,
        selector_features = features, feature_set_id = "all"), args)))
    run <- function(enabled) sb_with_selector_reuse(
      lapply(pipelines, function(p) sb_fit_final(task, p, capture_diagnostics = TRUE)),
      enabled = enabled, diagnostics = TRUE)
    calls <- 0L
    plain <- run(FALSE)
    expect_identical(calls, 2L, info = id)
    calls <- 0L
    shared <- run(TRUE)
    expect_identical(calls, 1L, info = id)
    expect_identical(shared$reuse$final_learner_hits, 1L, info = id)
    for (i in seq_along(pipelines)) {
      expect_identical(shared$result[[i]]$status, "success", info = id)
      expect_equal(shared$result[[i]], plain$result[[i]], info = id)
      expect_identical(shared$result[[i]]$learner_state,
        plain$result[[i]]$learner_state, info = id)
      expect_identical(sb_predict_final(unserialize(serialize(shared$result[[i]], NULL)),
        d[c("patient_id", features)], "patient_id"),
        sb_predict_final(plain$result[[i]], d[c("patient_id", features)], "patient_id"), info = id)
    }
    expect_identical(.Random.seed, before, info = id)
  }
  expect_null(.sb_selector_reuse_context$active)
})

test_that("final ensemble reuse keeps tuning settings and failures distinct", {
  x <- data.frame(a = 1:6, b = 6:1)
  e <- list(task_type = "survival", time = 1:6, event = rep(0:1, 3))
  ids <- letters[1:6]
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, ...) {
    calls <<- calls + 1L
    if (identical(learner$params$mode, "warning")) warning("native warning")
    if (identical(learner$params$mode, "error")) sb_abort("native failure", "native_failure")
    list(call = calls)
  })
  for (learner in list(sb_gbm_cox_cv_learner_adapter(),
    sb_rlt_learner_adapter("rlt_reinforcement_survival"))) {
    calls <- 0L
    run <- function(l = learner, selector = "none") {
      sb_fit_final_runtime_learner(l, x, e, ids, selector)
    }
    z <- sb_with_selector_reuse({
      first <- run()
      expect_identical(run(selector = "fixed_features"), first)
      altered <- learner; altered$params$seed <- 9L; run(altered)
      altered <- learner
      if (learner$id == "gbm_cox_cv") altered$params$inner_folds <- 3L else
        altered$params$num_trees <- 20L
      run(altered)
      run(selector = "univ_cox_top_k")
      run(selector = "univ_cox_top_k")
      altered <- learner; altered$params$mode <- "warning"
      for (i in 1:2) expect_warning(run(altered), "native warning")
      altered$params$mode <- "error"
      for (i in 1:2) expect_error(run(altered), class = "native_failure")
    }, diagnostics = TRUE)
    expect_identical(calls, 7L)
    expect_identical(z$reuse$final_learner_hits, 3L)
    expect_identical(z$reuse$final_learner_entries, 3L)
  }
})

test_that("GBM CV is not reused across assessment folds or horizon extraction", {
  skip_if_not_installed("gbm")
  d <- sb_example_survival(180L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- function(h) sb_surv_task(d, "patient_id", "time", "event", features,
    "month", prediction_horizon = h)
  p <- sb_pipeline("none", "gbm_cox_cv", learner_tuning_folds = 3L, gbm_cv_max_trees = 30L)
  native <- sb_fit_gbm_cox_cv
  calls <- 0L
  local_mocked_bindings(sb_fit_gbm_cox_cv = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  plain <- lapply(c(3, 9), function(h) sb_fit_final(task(h), p))
  calls <- 0L
  shared <- sb_with_selector_reuse({
    fits <- lapply(c(3, 9), function(h) sb_fit_final(task(h), p))
    cv <- sb_run_cv(task(3), p, folds = 3L, seed = 27L, progress = FALSE)
    list(fits = fits, cv = cv)
  }, diagnostics = TRUE)
  expect_identical(calls, 4L)
  expect_identical(shared$reuse$final_learner_hits, 1L)
  expect_equal(shared$result$fits, plain)
  expect_identical(shared$result$cv$complete_oof, TRUE)
  expect_identical(shared$result$cv$complete_survival_probability_oof, TRUE)
  expect_identical(shared$result$cv,
    sb_run_cv(task(3), p, folds = 3L, seed = 27L, progress = FALSE))
  a <- sb_predict_final(shared$result$fits[[1]], d[c("patient_id", features)], "patient_id")
  b <- sb_predict_final(shared$result$fits[[2]], d[c("patient_id", features)], "patient_id")
  expect_identical(a$predictions$risk, b$predictions$risk)
  expect_identical(identical(a$predictions$survival_probability, b$predictions$survival_probability), FALSE)
})

test_that("cohort comparisons reuse identical CV and final GBM fits separately", {
  skip_if_not_installed("gbm")
  features <- c("gene_a", "gene_b", "gene_c")
  make_task <- function(seed, prefix) {
    d <- sb_example_survival(180L, seed)
    d$patient_id <- paste0(prefix, d$patient_id)
    sb_surv_task(d, "patient_id", "time", "event", features, "month")
  }
  development <- make_task(11032L, "dev")
  external <- make_task(11033L, "external")
  pipelines <- list(sb_pipeline("none", "gbm_cox_cv", learner_tuning_folds = 3L,
      gbm_cv_max_trees = 30L), sb_pipeline("fixed_features", "gbm_cox_cv",
      selector_features = features, feature_set_id = "all", learner_tuning_folds = 3L,
      gbm_cv_max_trees = 30L))
  native <- sb_fit_gbm_cox_cv
  calls <- 0L
  local_mocked_bindings(sb_fit_gbm_cox_cv = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  shared <- sb_compare_cohorts(development, list(External = external),
    pipelines = pipelines, folds = 3L, seed = 27L, progress = FALSE)
  expect_identical(calls, 4L)
  expect_identical(unname(vapply(shared$fits, `[[`, "", "status")), rep("success", 2L))
  calls <- 0L
  local_mocked_bindings(sb_fit_final_runtime_learner = function(learner, x, endpoint, ...) {
    sb_fit_runtime_learner(learner, x, endpoint)
  })
  plain <- sb_compare_cohorts(development, list(External = external),
    pipelines = pipelines, folds = 3L, seed = 27L, progress = FALSE)
  expect_identical(calls, 5L)
  expect_equal(shared$fits, plain$fits)
  expect_identical(shared$predictions, plain$predictions)
  expect_equal(shared$cohorts, plain$cohorts)
})

test_that("seeded final ensembles do not create caller RNG state", {
  skip_if_not_installed("gbm")
  skip_if_not_installed("RLT")
  d <- sb_example_survival(90L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
  withr::local_seed(35L)
  rm(".Random.seed", envir = .GlobalEnv)
  for (p in list(sb_pipeline("none", "gbm_cox_cv", learner_tuning_folds = 3L,
      gbm_cv_max_trees = 30L), sb_pipeline("none", "rlt_reinforcement_survival"))) {
    fits <- sb_with_selector_reuse(list(sb_fit_final(task, p), sb_fit_final(task, p)),
      diagnostics = TRUE)
    expect_identical(fits$result[[1]]$status, "success")
    expect_equal(fits$result[[1]], fits$result[[2]])
    expect_identical(fits$reuse$final_learner_hits, 1L)
    expect_identical(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), FALSE)
  }
})
