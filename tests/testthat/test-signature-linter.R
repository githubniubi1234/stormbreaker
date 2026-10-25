test_that("the signature linter admits complete raw and constructed formulas", {
  raw_report <- sb_lint_signature(list(
    signature_id = "published_binary",
    task_type = "binary",
    coefficients = c(gene_a = 0.7, gene_b = -0.4),
    intercept = -0.2,
    source = "DOI supplied by the reviewer"
  ))
  survival <- sb_linear_signature(
    "published_survival",
    "survival",
    c(gene_a = 0.5),
    prediction_horizon = 12,
    time_unit = "month",
    baseline_survival = 0.8
  )
  object_report <- sb_lint_signature(survival)
  risk_report <- sb_lint_signature(list(
    signature_id = "published_risk_score",
    task_type = "survival",
    coefficients = c(gene_a = 0.5)
  ))

  expect_identical(
    names(raw_report),
    c(
      "status", "signature_id", "task_type", "n_features", "output",
      "source", "issue_code", "message"
    )
  )
  expect_identical(raw_report$status, "admitted")
  expect_identical(raw_report$signature_id, "published_binary")
  expect_identical(raw_report$task_type, "binary")
  expect_identical(raw_report$n_features, 2L)
  expect_identical(raw_report$output, "probability")
  expect_identical(raw_report$source, "DOI supplied by the reviewer")
  expect_true(is.na(raw_report$issue_code))
  expect_true(is.na(raw_report$message))

  expect_identical(object_report$status, "admitted")
  expect_identical(object_report$signature_id, "published_survival")
  expect_identical(object_report$n_features, 1L)
  expect_identical(object_report$output, "risk+survival_probability")
  expect_identical(risk_report$status, "admitted")
  expect_identical(risk_report$output, "risk")
})

test_that("the signature linter retains incomplete literature formulas", {
  missing_coefficients <- sb_lint_signature(list(
    signature_id = "missing_coefficients",
    task_type = "survival"
  ))
  missing_binary_intercept <- sb_lint_signature(list(
    signature_id = "missing_binary_intercept",
    task_type = "binary",
    coefficients = c(gene_a = 0.7)
  ))
  partial_horizon <- sb_lint_signature(list(
    signature_id = "partial_horizon",
    task_type = "survival",
    coefficients = c(gene_a = 0.5),
    prediction_horizon = 12,
    time_unit = "month"
  ))
  unknown_argument <- sb_lint_signature(list(
    signature_id = "unknown_argument",
    task_type = "survival",
    coefficients = c(gene_a = 0.5),
    published_cutoff = 0.3
  ))
  bad_time_unit <- sb_lint_signature(list(
    signature_id = "bad_time_unit",
    task_type = "survival",
    coefficients = c(gene_a = 0.5),
    prediction_horizon = 12,
    time_unit = "week",
    baseline_survival = 0.8
  ))
  unnamed <- sb_lint_signature(list(
    "unnamed_specification",
    "survival",
    c(gene_a = 0.5)
  ))
  bad_container <- sb_lint_signature("not a signature specification")

  expect_identical(missing_coefficients$status, "rejected")
  expect_identical(
    missing_coefficients$issue_code,
    "stormbreaker_missing_signature_argument"
  )
  expect_identical(missing_coefficients$signature_id, "missing_coefficients")
  expect_true(is.na(missing_coefficients$n_features))

  expect_identical(missing_binary_intercept$status, "rejected")
  expect_identical(
    missing_binary_intercept$issue_code,
    "stormbreaker_missing_signature_intercept"
  )
  expect_identical(missing_binary_intercept$n_features, 1L)
  expect_identical(missing_binary_intercept$task_type, "binary")

  expect_identical(partial_horizon$status, "rejected")
  expect_identical(
    partial_horizon$issue_code,
    "stormbreaker_bad_signature_horizon"
  )
  expect_identical(unknown_argument$status, "rejected")
  expect_identical(
    unknown_argument$issue_code,
    "stormbreaker_unknown_signature_argument"
  )
  expect_identical(bad_time_unit$status, "rejected")
  expect_identical(bad_time_unit$issue_code, "stormbreaker_bad_time_unit")
  expect_identical(unnamed$status, "rejected")
  expect_identical(unnamed$issue_code, "stormbreaker_bad_signature_spec")
  expect_identical(bad_container$status, "rejected")
  expect_identical(
    bad_container$issue_code,
    "stormbreaker_bad_signature_spec"
  )
})

test_that("the signature linter does not disguise unexpected internal errors", {
  local_mocked_bindings(
    sb_signature_from_lint_spec = function(spec) {
      sb_abort(
        "injected internal signature failure",
        "stormbreaker_internal_signature_failure"
      )
    }
  )

  error <- tryCatch(
    sb_lint_signature(list(
      signature_id = "otherwise_complete",
      task_type = "binary",
      coefficients = c(gene_a = 0.7),
      intercept = -0.2
    )),
    error = function(error) error
  )

  expect_s3_class(error, "stormbreaker_internal_signature_failure")
  expect_identical(
    conditionMessage(error),
    "injected internal signature failure"
  )
})
