test_that("a signature inventory retains every literature record in order", {
  records <- list(
    paper_a = list(
      signature_id = "complete_binary",
      task_type = "binary",
      coefficients = c(gene_a = 0.7),
      intercept = -0.2,
      source = "doi:paper-a"
    ),
    paper_b = list(
      signature_id = "binary_without_intercept",
      task_type = "binary",
      coefficients = c(gene_a = 0.5),
      source = "doi:paper-b"
    ),
    paper_c = list(
      task_type = "survival",
      coefficients = c(gene_a = 0.4),
      source = "doi:paper-c"
    ),
    paper_d = sb_linear_signature(
      "complete_survival",
      "survival",
      c(gene_a = 0.3)
    ),
    paper_e = "unparseable record"
  )

  observed <- sb_signature_inventory(records)

  expect_identical(
    names(observed),
    c(
      "record_id", "status", "signature_id", "task_type", "n_features",
      "output", "source", "issue_code", "message"
    )
  )
  expect_identical(observed$record_id, names(records))
  expect_identical(
    observed$status,
    c("admitted", "rejected", "rejected", "admitted", "rejected")
  )
  expect_identical(observed$source[1:3], c(
    "doi:paper-a", "doi:paper-b", "doi:paper-c"
  ))
  expect_identical(
    observed$issue_code[2L],
    "stormbreaker_missing_signature_intercept"
  )
  expect_identical(
    observed$issue_code[3L],
    "stormbreaker_missing_signature_argument"
  )
  expect_true(is.na(observed$signature_id[3L]))
  expect_identical(observed$output[4L], "risk")
  expect_identical(
    observed$issue_code[5L],
    "stormbreaker_bad_signature_spec"
  )
  expect_identical(nrow(observed), length(records))
})

test_that("a signature inventory requires unique non-empty record IDs", {
  valid_spec <- list(
    signature_id = "complete",
    task_type = "survival",
    coefficients = c(gene_a = 0.4)
  )
  empty <- capture_sb_error(sb_signature_inventory(list()))
  unnamed <- capture_sb_error(sb_signature_inventory(list(valid_spec)))
  duplicated <- capture_sb_error(sb_signature_inventory(structure(
    list(valid_spec, valid_spec),
    names = c("paper", "paper")
  )))
  blank <- capture_sb_error(sb_signature_inventory(structure(
    list(valid_spec),
    names = " "
  )))
  bad_container <- capture_sb_error(sb_signature_inventory("not a list"))
  bad_frame <- capture_sb_error(sb_signature_inventory(data.frame(record = 1)))
  single_spec <- capture_sb_error(sb_signature_inventory(valid_spec))
  single_object <- capture_sb_error(sb_signature_inventory(
    sb_linear_signature(
      "single_object",
      "survival",
      c(gene_a = 0.4)
    )
  ))
  legal_core_named_record <- sb_signature_inventory(list(
    task_type = valid_spec
  ))

  expect_s3_class(empty, "stormbreaker_bad_signature_inventory")
  expect_s3_class(unnamed, "stormbreaker_bad_signature_inventory")
  expect_s3_class(duplicated, "stormbreaker_bad_signature_inventory")
  expect_s3_class(blank, "stormbreaker_bad_signature_inventory")
  expect_s3_class(bad_container, "stormbreaker_bad_signature_inventory")
  expect_s3_class(bad_frame, "stormbreaker_bad_signature_inventory")
  expect_s3_class(single_spec, "stormbreaker_bad_signature_inventory")
  expect_s3_class(single_object, "stormbreaker_bad_signature_inventory")
  expect_identical(legal_core_named_record$record_id, "task_type")
  expect_identical(legal_core_named_record$status, "admitted")
})

test_that("a signature inventory never disguises an internal linter error", {
  local_mocked_bindings(
    sb_lint_signature = function(spec) {
      sb_abort(
        "injected inventory failure",
        "stormbreaker_internal_inventory_failure"
      )
    }
  )

  error <- tryCatch(
    sb_signature_inventory(list(record = list(
      signature_id = "complete",
      task_type = "survival",
      coefficients = c(gene_a = 0.4)
    ))),
    error = function(error) error
  )

  expect_s3_class(error, "stormbreaker_internal_inventory_failure")
  expect_identical(conditionMessage(error), "injected inventory failure")
})
