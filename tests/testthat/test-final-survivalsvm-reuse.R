test_that("identical none and fixed-feature final SVM fits share only native fitting", {
  skip_if_not_installed("survivalsvm")
  d <- sb_example_survival(72L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
    prediction_horizon = 7)
  pipelines <- list(sb_pipeline("none", "survivalsvm"),
    sb_pipeline("fixed_features", "survivalsvm",
      selector_features = features, feature_set_id = "declared"))
  native <- sb_fit_survivalsvm
  calls <- 0L
  local_mocked_bindings(sb_fit_survivalsvm = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  run <- function(enabled) sb_with_selector_reuse(
    lapply(pipelines, function(p) sb_fit_final(task, p, capture_diagnostics = TRUE)),
    enabled = enabled, diagnostics = TRUE)
  plain <- run(FALSE)
  expect_identical(calls, 2L)
  calls <- 0L
  shared <- run(TRUE)
  expect_identical(calls, 1L)
  expect_identical(shared$reuse$final_learner_hits, 1L)
  for (i in seq_along(pipelines)) {
    a <- plain$result[[i]]; b <- shared$result[[i]]
    expect_identical(b$status, "success")
    expect_identical(b$pipeline, pipelines[[i]])
    expect_equal(b, a)
    expect_identical(sb_predict_final(unserialize(serialize(b, NULL)),
      d[c("patient_id", features)], "patient_id"),
      sb_predict_final(unserialize(serialize(a, NULL)),
        d[c("patient_id", features)], "patient_id"))
  }
  expect_null(.sb_selector_reuse_context$active)
})

test_that("all fixed survival-SVM variants replay exactly without consuming RNG", {
  skip_if_not_installed("survivalsvm")
  d <- sb_example_survival(72L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
    prediction_horizon = 7)
  native <- sb_fit_survivalsvm
  calls <- 0L
  local_mocked_bindings(sb_fit_survivalsvm = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  withr::local_seed(21L)
  before <- .Random.seed
  for (id in sb_survivalsvm_learner_ids()) {
    p <- sb_pipeline("none", id)
    q <- sb_pipeline("fixed_features", id, selector_features = features,
      feature_set_id = "all")
    plain <- sb_fit_final(task, q, capture_diagnostics = TRUE)
    count <- calls
    shared <- sb_with_selector_reuse({
      first <- sb_fit_final(task, p)
      second <- sb_fit_final(task, q, capture_diagnostics = TRUE)
      expect_null(first$diagnostic_artifacts)
      second
    }, diagnostics = TRUE)
    expect_identical(calls - count, 1L, info = id)
    expect_identical(plain$status, "success", info = id)
    expect_identical(shared$result$status, "success", info = id)
    expect_equal(shared$result, plain, info = id)
    expect_identical(sb_predict_final(unserialize(serialize(shared$result, NULL)),
      d[c("patient_id", features)], "patient_id"),
      sb_predict_final(plain, d[c("patient_id", features)], "patient_id"), info = id)
    expect_identical(.Random.seed, before, info = id)
  }
  rm(".Random.seed", envir = .GlobalEnv)
  sb_with_selector_reuse({
    sb_fit_final(task, sb_pipeline("none", "survivalsvm"))
    sb_fit_final(task, sb_pipeline("none", "survivalsvm"))
  })
  expect_identical(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), FALSE)
})

test_that("final learner reuse distinguishes complete data and retains distinct fits", {
  learner <- sb_survivalsvm_learner_adapter()
  x <- data.frame(a = 1:6, b = 6:1)
  endpoint <- list(task_type = "survival", time = 1:6, event = rep(0:1, 3),
    group_ids = NULL)
  ids <- paste0("patient", 1:6)
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(...) {
    calls <<- calls + 1L
    list(value = calls)
  })
  run <- function(l = learner, data = x, e = endpoint, patients = ids,
                  selector = "none") {
    sb_fit_final_runtime_learner(l, data, e, patients, selector)
  }
  out <- sb_with_selector_reuse({
    a <- run()
    expect_identical(run(selector = "fixed_features"), a)
    run(data = x[c("b", "a")])
    run(data = x[6:1, ])
    changed <- x; changed$a[1] <- 99
    run(data = changed)
    changed <- x; rownames(changed) <- ids
    run(data = changed)
    e <- endpoint; e$event[1] <- 1L
    run(e = e)
    e <- endpoint; e$time[1] <- 1.5
    run(e = e)
    e <- endpoint; e$group_ids <- rep(1:3, 2)
    run(e = e)
    run(patients = paste0("other", ids))
    l <- learner; l$params$gamma_mu <- 2
    run(l = l)
    l <- learner; body(l$fit) <- quote(stop("changed implementation"))
    run(l = l)
    run()
  }, diagnostics = TRUE)
  expect_identical(calls, 11L)
  expect_identical(out$reuse$final_learner_hits, 2L)
  expect_identical(out$reuse$final_learner_entries, 11L)
  expect_identical(out$reuse$hits, 0L)
})

test_that("warnings and failures are retried while successful eligible fits reuse", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_learner = function(learner, ...) {
    calls <<- calls + 1L
    if (learner$params$mode == "warning") warning("native warning")
    if (learner$params$mode == "error") stop("native failure")
    list(value = calls)
  })
  l <- sb_survivalsvm_learner_adapter()
  x <- data.frame(a = 1:4); e <- list(task_type = "survival", time = 1:4, event = c(1, 0, 1, 0))
  ids <- letters[1:4]
  run <- function(selector = "none") sb_fit_final_runtime_learner(l, x, e, ids, selector)
  z <- sb_with_selector_reuse({
    l$params$mode <- "warning"
    for (i in 1:2) expect_warning(run(), "native warning")
    l$params$mode <- "error"
    for (i in 1:2) expect_error(run(), "native failure")
    l$params$mode <- "ordinary"
    l$id <- "gbm_cox_cv"
    for (i in 1:2) run("univ_cox_top_k")
    l$id <- "coxboost"
    for (i in 1:2) run()
  }, diagnostics = TRUE)
  expect_identical(calls, 6L)
  expect_identical(z$reuse$final_learner_entries, 2L)
  expect_identical(z$reuse$final_learner_hits, 2L)
})

test_that("final learner scopes are isolated and CV still fits every assessment fold", {
  skip_if_not_installed("survivalsvm")
  d <- sb_example_survival(90L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
  p <- sb_pipeline("none", "survivalsvm")
  native <- sb_fit_survivalsvm
  calls <- 0L
  local_mocked_bindings(sb_fit_survivalsvm = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  out <- sb_with_selector_reuse({
    first <- sb_fit_final(task, p)
    sb_with_selector_reuse(sb_fit_final(task, p), enabled = FALSE)
    sb_with_selector_reuse(sb_fit_final(task, p))
    expect_equal(sb_fit_final(task, p), first)
    cv <- sb_run_cv(task, p, folds = 3L, seed = 21L, progress = FALSE)
    cv
  }, diagnostics = TRUE)
  expect_identical(calls, 6L)
  expect_identical(out$reuse$final_learner_hits, 1L)
  expect_null(.sb_selector_reuse_context$active)
  expect_error(sb_with_selector_reuse({sb_fit_final(task, p); stop("scope exit")}), "scope exit")
  expect_null(.sb_selector_reuse_context$active)
})

test_that("horizon extraction stays route-local after a shared final SVM fit", {
  skip_if_not_installed("survivalsvm")
  d <- sb_example_survival(72L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  make_task <- function(h) sb_surv_task(d, "patient_id", "time", "event", features,
    "month", prediction_horizon = h)
  p <- sb_pipeline("none", "survivalsvm")
  plain <- lapply(c(3, 9), function(h) sb_fit_final(make_task(h), p))
  shared <- sb_with_selector_reuse(
    lapply(c(3, 9), function(h) sb_fit_final(make_task(h), p)), diagnostics = TRUE)
  expect_identical(shared$reuse$final_learner_hits, 1L)
  expect_equal(shared$result, plain)
  a <- sb_predict_final(shared$result[[1]], d[c("patient_id", features)], "patient_id")
  b <- sb_predict_final(shared$result[[2]], d[c("patient_id", features)], "patient_id")
  expect_identical(a$predictions$risk, b$predictions$risk)
  expect_identical(identical(a$predictions$survival_probability, b$predictions$survival_probability), FALSE)
})

test_that("cohort comparisons share matching CV and final fits in separate scopes", {
  skip_if_not_installed("survivalsvm")
  features <- c("gene_a", "gene_b", "gene_c")
  make_task <- function(seed, prefix) {
    d <- sb_example_survival(72L, seed)
    d$patient_id <- paste0(prefix, d$patient_id)
    sb_surv_task(d, "patient_id", "time", "event", features, "month")
  }
  development <- make_task(11032L, "dev")
  external <- make_task(11033L, "external")
  pipelines <- list(sb_pipeline("none", "survivalsvm"),
    sb_pipeline("fixed_features", "survivalsvm",
      selector_features = features, feature_set_id = "all"))
  native <- sb_fit_survivalsvm
  calls <- 0L
  local_mocked_bindings(sb_fit_survivalsvm = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  shared <- sb_compare_cohorts(development, list(External = external),
    pipelines = pipelines, folds = 3L, seed = 21L, progress = FALSE)
  expect_identical(calls, 4L)
  calls <- 0L
  local_mocked_bindings(sb_fit_final_runtime_learner = function(learner, x, endpoint, ...) {
    sb_fit_runtime_learner(learner, x, endpoint)
  })
  plain <- sb_compare_cohorts(development, list(External = external),
    pipelines = pipelines, folds = 3L, seed = 21L, progress = FALSE)
  expect_identical(calls, 5L)
  expect_equal(shared$fits, plain$fits)
  expect_identical(shared$predictions, plain$predictions)
  expect_equal(shared$cohorts, plain$cohorts)
})
