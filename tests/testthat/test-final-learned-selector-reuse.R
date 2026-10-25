test_that("learned and declared selectors share identical final SVM inputs", {
  skip_if_not_installed("survivalsvm")
  d <- sb_example_survival(72L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
    prediction_horizon = 7)
  p <- sb_pipeline("univ_cox_top_k", "survivalsvm", selector_n_features = 2L)
  first <- sb_fit_final(task, p)
  expect_identical(first$status, "success")
  q <- sb_pipeline("fixed_features", "survivalsvm",
    selector_features = first$required_features, feature_set_id = "same_selected")
  native <- sb_fit_survivalsvm
  calls <- 0L
  local_mocked_bindings(sb_fit_survivalsvm = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  run <- function(enabled) sb_with_selector_reuse(
    lapply(list(p, q), function(pipeline) sb_fit_final(task, pipeline,
      capture_diagnostics = TRUE)), enabled = enabled, diagnostics = TRUE)
  withr::local_seed(42L)
  before <- .Random.seed
  plain <- run(FALSE)
  expect_identical(calls, 2L)
  calls <- 0L
  shared <- run(TRUE)
  expect_identical(calls, 1L)
  expect_identical(shared$reuse$final_learner_hits, 1L)
  expect_identical(.Random.seed, before)
  for (i in 1:2) {
    expect_identical(shared$result[[i]]$status, "success")
    expect_equal(shared$result[[i]], plain$result[[i]])
    expect_identical(shared$result[[i]]$learner_state, plain$result[[i]]$learner_state)
    expect_identical(sb_predict_final(unserialize(serialize(shared$result[[i]], NULL)),
      d[c("patient_id", features)], "patient_id"),
      sb_predict_final(plain$result[[i]], d[c("patient_id", features)], "patient_id"))
  }
})
