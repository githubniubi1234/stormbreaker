test_that("equivalent parameter representations remain executable", {
  p <- sb_pipeline(learner = "nnet_binary")
  canonical <- p
  p$learner_params$maxit <- as.numeric(p$learner_params$maxit)
  p$learner_params <- rev(p$learner_params)
  expect_identical(sb_validate_pipeline_contract(p, "binary"), canonical)
  chain <- sb_new_chain_pipeline(p, list(p))
  chain$learner_params <- rev(chain$learner_params)
  expect_identical(sb_validate_pipeline_contract(chain, "binary"),
    sb_new_chain_pipeline(canonical, list(canonical)))
  p$learner_params$maxit <- 200
  canonical$learner_params$maxit <- 200L
  expect_identical(sb_validate_pipeline_contract(p, "binary"), canonical)
})

test_that("equivalent parameters reach the native forest with usable types", {
  skip_if_not_installed("ranger")
  canonical <- sb_pipeline(learner = "ranger_binary", ranger_num_trees = 5L)
  p <- canonical
  p$learner_params$num_trees <- 5
  p$learner_params <- rev(p$learner_params)
  validated <- sb_validate_pipeline_contract(p, "binary")
  x <- data.frame(a = sin(1:30), b = cos(1:30))
  y <- rep(0:1, 15L)
  fitted <- sb_fit_ranger_binary(x, y, validated$learner_params)
  reference <- sb_fit_ranger_binary(x, y, canonical$learner_params)
  expect_identical(sb_predict_ranger_binary(fitted, x),
    sb_predict_ranger_binary(reference, x))
  chain <- sb_validate_pipeline_contract(sb_new_chain_pipeline(p, list(p)), "binary")
  expect_identical(chain$learner_params, canonical$learner_params)
  expect_identical(chain$selector_params$stages[[1L]]$learner_params,
    canonical$learner_params)
})

test_that("parameter equivalence does not hide scientific changes", {
  expect_identical(sb_pipeline_parameters_equal(
    list(a = 1L, b = list(x = 2L, y = 3L)),
    list(b = list(y = 3, x = 2), a = 1)), TRUE)
  expect_identical(sb_pipeline_parameters_equal(1, 1 + 1e-12), FALSE)
  expect_identical(sb_pipeline_parameters_equal(c(3L, 2L), c(2, 3)), FALSE)
  expect_identical(sb_pipeline_parameters_equal(list(a = 1L), list(b = 1)), FALSE)
  expect_identical(sb_pipeline_parameters_equal(list(1L, 2L), list(2, 1)), FALSE)
  expect_identical(sb_pipeline_parameters_equal(c(x = 1L), c(y = 1)), FALSE)
  expect_identical(sb_pipeline_parameters_equal(TRUE, 1L), FALSE)
  p <- sb_pipeline(learner = "nnet_binary")
  p$learner_params$decay <- p$learner_params$decay + 1e-12
  expect_s3_class(tryCatch(sb_validate_pipeline_contract(p, "binary"),
    error = identity), "stormbreaker_unresolved_pipeline")
})

test_that("chain annotations do not change the executable specification", {
  base <- sb_pipeline(learner = "logistic")
  stage <- sb_pipeline(learner = "logistic")
  chain <- sb_new_chain_pipeline(base, list(stage))
  chain$note <- "Clinician annotation"
  chain$chain_base$note <- "Development cohort"
  chain$selector_params$stages[[1L]]$note <- "First stage"
  expect_identical(sb_validate_pipeline_contract(chain, "binary"), chain)

  chain$positive_class <- 0L
  expect_identical(conditionMessage(tryCatch(
    sb_validate_pipeline_contract(chain, "binary"), error = identity)),
    "Selector chain does not match its executable specification.")
})

test_that("chain validation does not reconstruct an already built chain", {
  base <- sb_pipeline(learner = "logistic")
  chain <- sb_new_chain_pipeline(base, list(base))
  local_mocked_bindings(sb_new_chain_pipeline = function(...) {
    stop("constructor should not run during validation")
  })
  expect_identical(sb_validate_pipeline_contract(chain, "binary"), chain)
})

test_that("invalid chain stages and missing pipeline class still fail", {
  base <- sb_pipeline(learner = "logistic")
  chain <- sb_new_chain_pipeline(base, list(base))
  expect_s3_class(tryCatch(sb_validate_pipeline_contract(unclass(chain), "binary"),
    error = identity), "stormbreaker_bad_pipeline")
  empty <- chain
  empty$selector_params$stages <- list()
  expect_identical(conditionMessage(tryCatch(
    sb_validate_pipeline_contract(empty, "binary"), error = identity)),
    "A selector chain requires stages.")
  chain$selector_params$stages[[1L]] <- sb_pipeline()
  expect_s3_class(tryCatch(sb_validate_pipeline_contract(chain, "binary"),
    error = identity), "stormbreaker_pipeline_task_mismatch")
})

test_that("ordinary C50 warnings remain warnings and preserve the fit", {
  skip_if_not_installed("C50")
  native <- structure(list(size = 1L), class = "C5.0")
  local_mocked_bindings(C5.0 = function(...) {
    warning("informational native warning")
    native
  }, .package = "C50")
  seen <- character()
  result <- withCallingHandlers(sb_fit_c50_binary(
    data.frame(x = 1:4), factor(c(0, 0, 1, 1)), "c50_binary",
    sb_classical_binary_params("c50_binary")), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(result, native)
  expect_identical(seen, "informational native warning")
})

test_that("robust logistic still rejects unconverged native fits", {
  skip_if_not_installed("robustbase")
  local_mocked_bindings(glmrob = function(...) {
    list(converged = FALSE, coefficients = c(0, 1), iter = 1L)
  }, .package = "robustbase")
  result <- tryCatch(sb_fit_robust_logistic_mqle(
    data.frame(x = 1:4), factor(c(0, 0, 1, 1)),
    "robust_logistic_mqle_binary_fixed",
    sb_classical_binary_params("robust_logistic_mqle_binary_fixed")),
    error = identity)
  expect_s3_class(result, "stormbreaker_incomplete_robust_logistic_fit")
})

test_that("nnet warnings do not replace its convergence check", {
  skip_if_not_installed("nnet")
  local_mocked_bindings(nnet = function(...) {
    warning("native diagnostic")
    list(convergence = 1L)
  }, .package = "nnet")
  seen <- character()
  result <- tryCatch(withCallingHandlers(sb_fit_classical_binary(
    data.frame(x = 1:8), rep(c(0L, 1L), 4L), "nnet_binary",
    sb_classical_binary_params("nnet_binary")), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  expect_identical(seen, "native diagnostic")
  expect_s3_class(result, "stormbreaker_nnet_nonconvergence")
})

test_that("BART warnings do not bypass the posterior-tree check", {
  skip_if_not_installed("BART")
  local_mocked_bindings(pbart = function(...) {
    warning("native diagnostic")
    list()
  }, .package = "BART")
  seen <- character()
  result <- tryCatch(withCallingHandlers(sb_fit_bart_binary(
    matrix(1:8, ncol = 1L), rep(c(0L, 1L), 4L), "bart_binary",
    sb_classical_binary_params("bart_binary")), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  expect_identical(seen, "native diagnostic")
  expect_s3_class(result, "stormbreaker_incomplete_bart_binary_fit")
})
