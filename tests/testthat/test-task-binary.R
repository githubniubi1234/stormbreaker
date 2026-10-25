test_that("group roles are explicit without weakening unique observation IDs", {
  d <- sb_example_binary(20L, 3755L)
  d$group <- rep(1:10, each = 2)
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a", group_col = "group")
  expect_identical(task$data, d)
  expect_identical(task$n, 20L)
  expect_identical(task$group_col, "group")
  for (column in c("patient_id", "outcome", "gene_a")) {
    error <- capture_sb_error(sb_binary_task(
      d, "patient_id", "outcome", "gene_a", group_col = column
    ))
    expect_s3_class(error, "stormbreaker_overlapping_group_roles")
  }
  missing <- capture_sb_error(sb_binary_task(
    d, "patient_id", "outcome", "gene_a", group_col = "absent"
  ))
  expect_s3_class(missing, "stormbreaker_missing_group_column")
  for (bad in list(rep(NA_character_, 20), rep(" ", 20),
                  rep(Inf, 20), as.list(1:20))) {
    changed <- d
    changed$group <- bad
    error <- capture_sb_error(sb_binary_task(
      changed, "patient_id", "outcome", "gene_a", group_col = "group"
    ))
    expect_s3_class(error, "stormbreaker_bad_group_ids")
  }
  d$patient_id[2] <- d$patient_id[1]
  error <- capture_sb_error(sb_binary_task(
    d, "patient_id", "outcome", "gene_a", group_col = "group"
  ))
  expect_s3_class(error, "stormbreaker_duplicate_id")
})

test_that("a valid binary task records explicit class semantics", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b")
  )

  expect_s3_class(task, "stormbreaker_binary_task")
  expect_identical(task$n, nrow(d))
  expect_identical(task$n_positive, sum(d$outcome == 1L))
  expect_identical(task$n_negative, sum(d$outcome == 0L))
  expect_identical(task$positive_class, 1L)
  expect_identical(task$negative_class, 0L)
  expect_identical(task$predictors, c("gene_a", "gene_b"))
  expect_identical(task$independent_unit, "patient")
})

test_that("logical binary outcomes preserve explicit zero-one meaning", {
  d <- make_sb_binary_fixture()
  d$outcome <- as.logical(d$outcome)
  task <- sb_binary_task(d, "patient_id", "outcome", "gene_a")

  expect_identical(task$n_positive, sum(d$outcome))
  expect_identical(task$positive_class, 1L)
})

test_that("binary outcomes require finite zero and one with both classes", {
  d <- make_sb_binary_fixture()
  d$outcome[1] <- 2L
  value_error <- capture_sb_error(
    sb_binary_task(d, "patient_id", "outcome", "gene_a")
  )
  expect_s3_class(value_error, "stormbreaker_bad_binary_outcome")

  d <- make_sb_binary_fixture()
  d$outcome <- 1L
  negative_error <- capture_sb_error(
    sb_binary_task(d, "patient_id", "outcome", "gene_a")
  )
  expect_s3_class(negative_error, "stormbreaker_no_negative_class")

  d$outcome <- 0L
  positive_error <- capture_sb_error(
    sb_binary_task(d, "patient_id", "outcome", "gene_a")
  )
  expect_s3_class(positive_error, "stormbreaker_no_positive_class")
})

test_that("binary task roles and numeric predictor values are explicit", {
  d <- make_sb_binary_fixture()
  role_error <- capture_sb_error(
    sb_binary_task(d, "patient_id", "outcome", "outcome")
  )
  expect_s3_class(role_error, "stormbreaker_reserved_predictor")

  d$group <- factor(rep(c("A", "B"), length.out = nrow(d)))
  type_error <- capture_sb_error(
    sb_binary_task(d, "patient_id", "outcome", "group")
  )
  expect_s3_class(type_error, "stormbreaker_bad_predictor_type")

  d$gene_a[1] <- NA_real_
  value_error <- capture_sb_error(
    sb_binary_task(d, "patient_id", "outcome", "gene_a")
  )
  expect_s3_class(value_error, "stormbreaker_bad_predictor_values")
})

test_that("binary tasks reserve the model intercept name", {
  d <- make_sb_binary_fixture()
  d[["(Intercept)"]] <- d$gene_a

  error <- capture_sb_error(sb_binary_task(
    d, "patient_id", "outcome", c("(Intercept)", "gene_b")
  ))

  expect_s3_class(error, "stormbreaker_reserved_intercept_name")
  expect_match(conditionMessage(error), "reserved for the model intercept")
})

test_that("binary patient identifiers remain unique and complete", {
  d <- make_sb_binary_fixture()
  d$patient_id[2] <- d$patient_id[1]
  duplicate_error <- capture_sb_error(
    sb_binary_task(d, "patient_id", "outcome", "gene_a")
  )
  expect_s3_class(duplicate_error, "stormbreaker_duplicate_id")

  d <- make_sb_binary_fixture()
  d$patient_id[1] <- ""
  blank_error <- capture_sb_error(
    sb_binary_task(d, "patient_id", "outcome", "gene_a")
  )
  expect_s3_class(blank_error, "stormbreaker_bad_id")
})

test_that("binary tasks preserve study-declared decision thresholds", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b"),
    decision_thresholds = c(low = 0.10, high = 0.35)
  )

  expect_identical(task$decision_thresholds, c(0.10, 0.35))
  expect_identical(
    stormbreaker:::sb_validate_binary_task_contract(task)$decision_thresholds,
    c(0.10, 0.35)
  )
})

test_that("binary tasks omit decision thresholds by default", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b")
  )

  expect_false("decision_thresholds" %in% names(task))
})

test_that("binary decision thresholds must be unique increasing probabilities", {
  d <- make_sb_binary_fixture()
  bad <- list(
    numeric(),
    c(0, 0.2),
    c(0.2, 1),
    c(NA_real_, 0.2),
    c(Inf, 0.2),
    c(0.2, 0.2),
    c(0.4, 0.2),
    "0.2"
  )

  for (value in bad) {
    error <- capture_sb_error(sb_binary_task(
      d,
      "patient_id",
      "outcome",
      c("gene_a", "gene_b"),
      decision_thresholds = value
    ))
    expect_s3_class(error, "stormbreaker_bad_decision_thresholds")
  }
})

test_that("binary task reconstruction rejects invalid threshold mutation", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b"),
    decision_thresholds = c(0.10, 0.30)
  )
  task$decision_thresholds <- c(0.10, 1)

  error <- capture_sb_error(
    stormbreaker:::sb_validate_binary_task_contract(task)
  )
  expect_s3_class(error, "stormbreaker_bad_decision_thresholds")
})

test_that("default binary tasks still enter CV and final fitting", {
  d <- make_sb_binary_fixture()
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b")
  )
  pipeline <- sb_pipeline(learner = "logistic")

  cv <- sb_run_cv(task, pipeline, folds = 3, seed = 2026)
  fit <- sb_fit_final(task, pipeline)

  expect_identical(cv$status, "success")
  expect_identical(fit$status, "success")
  expect_false("decision_thresholds" %in% names(cv$task))
  expect_false("decision_thresholds" %in% names(fit$prediction))
})
