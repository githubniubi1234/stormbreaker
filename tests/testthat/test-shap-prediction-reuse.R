test_that("SHAP validates the saved fit once without changing predictions", {
  d <- sb_example_binary(60, 42)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, bi_pipeline("Logistic"))
  validate <- sb_validate_final_fit_for_prediction
  calls <- 0L
  local_mocked_bindings(sb_validate_final_fit_for_prediction = function(fit) {
    calls <<- calls + 1L
    validate(fit)
  })

  result <- sb_model_shap(fit, d[1:2, ], "patient_id", d,
    permutations = 2, background_n = 3, seed = 12)

  expect_identical(calls, 1L)
  expected <- sb_predict_final(fit, d[1:2, ], "patient_id")$predictions$probability
  expect_identical(result$scores$score, expected)
  expect_equal(unname(result$reference_value + rowSums(result$contributions)), expected)
  expect_identical(result, model_shap(fit, d[1:2, ], "patient_id", d,
    permutations = 2, background_n = 3, seed = 12))
})

test_that("tree SHAP does not validate a second time through the generic wrapper", {
  skip_if_not_installed("rpart")
  d <- sb_example_binary(60, 42)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "decision_tree_binary"))
  expect_identical(fit$status, "success")
  validate <- sb_validate_final_fit_for_prediction
  calls <- 0L
  local_mocked_bindings(sb_validate_final_fit_for_prediction = function(fit) {
    calls <<- calls + 1L
    validate(fit)
  })

  result <- sb_tree_shap(fit, d[1:2, ], "patient_id", d,
    permutations = 2, background_n = 3)

  expect_identical(calls, 1L)
  expect_s3_class(result, "stormbreaker_tree_shap")
  expect_equal(unname(result$reference_value + rowSums(result$contributions)), result$scores$score)
})

test_that("the SHAP inner predictor still rejects invalid probability values", {
  d <- sb_example_binary(60, 42)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, bi_pipeline("Logistic"))
  local_mocked_bindings(sb_predict_final_state = function(state, x) {
    list(probability = rep(1.1, nrow(x)))
  })

  error <- tryCatch(sb_model_shap(fit, d[1:2, ], "patient_id", d,
    permutations = 2, background_n = 3), error = identity)

  expect_s3_class(error, "stormbreaker_invalid_final_probability")
})

test_that("survival SHAP explains the saved risk and preserves patient order", {
  d <- sb_example_survival(60, 42)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b"), "month")
  fit <- sb_fit_final(task, surv_pipeline("Cox"))
  query <- d[c(4, 1, 3), ]

  result <- sb_model_shap(fit, query, "patient_id", d,
    permutations = 2, background_n = 3)

  expect_identical(result$scores$patient_id, query$patient_id)
  expect_identical(result$scores$score,
    sb_predict_final(fit, query, "patient_id")$predictions$risk)
  expect_equal(unname(result$reference_value + rowSums(result$contributions)), result$scores$score)
  expect_identical(result$score_scale, "saved_risk_score")
})

test_that("survival SHAP still rejects non-finite inner predictions", {
  d <- sb_example_survival(60, 42)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b"), "month")
  fit <- sb_fit_final(task, surv_pipeline("Cox"))
  local_mocked_bindings(sb_predict_final_state = function(state, x) {
    list(risk = rep(Inf, nrow(x)))
  })

  error <- tryCatch(sb_model_shap(fit, d[1:2, ], "patient_id", d,
    permutations = 2, background_n = 3), error = identity)

  expect_s3_class(error, "stormbreaker_invalid_final_risk")
})
