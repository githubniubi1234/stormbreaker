test_that("a valid survival task records endpoint semantics", {
  d <- make_sb_fixture()
  task <- sb_surv_task(
    d,
    id_col = "patient_id",
    time_col = "time",
    event_col = "event",
    predictors = c("gene_a", "gene_b"),
    time_unit = "month"
  )

  expect_s3_class(task, "stormbreaker_surv_task")
  expect_identical(task$n, nrow(d))
  expect_identical(task$n_events, sum(d$event))
  expect_identical(task$predictors, c("gene_a", "gene_b"))
  expect_identical(task$time_unit, "month")
  expect_null(task$prediction_horizon)
  expect_named(
    task,
    c(
      "data", "id_col", "time_col", "event_col", "predictors",
      "time_unit", "n", "n_events", "n_censored", "independent_unit"
    ),
    ignore.order = FALSE
  )
})

test_that("a survival task records one predeclared prediction horizon", {
  d <- make_sb_fixture()
  task <- sb_surv_task(
    d,
    id_col = "patient_id",
    time_col = "time",
    event_col = "event",
    predictors = c("gene_a", "gene_b"),
    time_unit = "month",
    prediction_horizon = 24L
  )

  expect_identical(task$prediction_horizon, 24)
  expect_identical(
    names(task),
    c(
      "data", "id_col", "time_col", "event_col", "predictors",
      "time_unit", "prediction_horizon", "n", "n_events", "n_censored",
      "independent_unit"
    )
  )
  printed <- paste(capture.output(print(task)), collapse = "\n")
  expect_match(printed, "Prediction horizon: 24 \\[ month \\]")
})

test_that("a survival task preserves declared event-risk thresholds", {
  d <- make_sb_fixture()
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 24,
    decision_thresholds = c(low = 0.10, high = 0.35)
  )

  expect_identical(task$decision_thresholds, c(0.10, 0.35))
  expect_identical(
    sb_validate_task_contract(task)$decision_thresholds,
    c(0.10, 0.35)
  )
  printed <- paste(capture.output(print(task)), collapse = "\n")
  expect_match(printed, "Decision thresholds: 2")
})

test_that("survival decision thresholds require a declared horizon", {
  d <- make_sb_fixture()
  error <- capture_sb_error(sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    "gene_a",
    "month",
    decision_thresholds = c(0.10, 0.35)
  ))

  expect_s3_class(
    error,
    "stormbreaker_decision_thresholds_require_prediction_horizon"
  )
})

test_that("survival decision thresholds use the shared strict validator", {
  d <- make_sb_fixture()
  invalid <- list(numeric(), c(0, 0.2), c(0.2, 1), c(0.3, 0.2), c(0.2, 0.2))

  for (value in invalid) {
    error <- capture_sb_error(sb_surv_task(
      d,
      "patient_id",
      "time",
      "event",
      "gene_a",
      "month",
      prediction_horizon = 24,
      decision_thresholds = value
    ))
    expect_s3_class(error, "stormbreaker_bad_decision_thresholds")
  }
})

test_that("survival task reconstruction rejects threshold mutation", {
  d <- make_sb_fixture()
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    "gene_a",
    "month",
    prediction_horizon = 24,
    decision_thresholds = c(0.10, 0.35)
  )
  task$decision_thresholds <- c(0.10, 1)

  error <- capture_sb_error(sb_validate_task_contract(task))
  expect_s3_class(error, "stormbreaker_bad_decision_thresholds")

  no_horizon <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  no_horizon$decision_thresholds <- 0.20
  horizon_error <- capture_sb_error(sb_validate_task_contract(no_horizon))
  expect_s3_class(
    horizon_error,
    "stormbreaker_decision_thresholds_require_prediction_horizon"
  )
})

test_that("a similarly named field cannot activate a prediction horizon", {
  d <- make_sb_fixture()
  task <- sb_surv_task(
    d,
    id_col = "patient_id",
    time_col = "time",
    event_col = "event",
    predictors = c("gene_a", "gene_b"),
    time_unit = "month"
  )
  task$prediction_horizons <- 5

  validated <- sb_validate_task_contract(task)

  expect_identical("prediction_horizon" %in% names(validated), FALSE)
})

test_that("prediction horizon must be one finite positive number", {
  d <- make_sb_fixture()
  invalid <- list(
    0,
    -1,
    NA_real_,
    Inf,
    "24",
    c(12, 24)
  )

  for (value in invalid) {
    error <- capture_sb_error(
      sb_surv_task(
        d, "patient_id", "time", "event", "gene_a", "month",
        prediction_horizon = value
      )
    )
    expect_s3_class(error, "stormbreaker_bad_prediction_horizon")
  }
})

test_that("duplicate patients and invalid endpoints are rejected", {
  d <- make_sb_fixture()
  d$patient_id[2] <- d$patient_id[1]
  duplicate_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "gene_a", "month"
    )
  )
  expect_s3_class(duplicate_error, "stormbreaker_duplicate_id")

  d <- make_sb_fixture()
  d$event[1] <- 2
  event_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "gene_a", "month"
    )
  )
  expect_s3_class(event_error, "stormbreaker_bad_event")

  d <- make_sb_fixture()
  d$time[1] <- -1
  time_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "gene_a", "month"
    )
  )
  expect_s3_class(time_error, "stormbreaker_bad_time")
})

test_that("predictor roles, values and time units are explicit", {
  d <- make_sb_fixture()
  role_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "event", "month"
    )
  )
  expect_s3_class(role_error, "stormbreaker_reserved_predictor")

  d$gene_a[1] <- NA_real_
  missing_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "gene_a", "month"
    )
  )
  expect_s3_class(missing_error, "stormbreaker_bad_predictor_values")

  d <- make_sb_fixture()
  unit_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "gene_a", "months"
    )
  )
  expect_s3_class(unit_error, "stormbreaker_bad_time_unit")
})

test_that("all-event, all-censored, and nonnumeric predictors are rejected", {
  d <- make_sb_fixture()
  d$event <- 1L
  censoring_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "gene_a", "month"
    )
  )
  expect_s3_class(censoring_error, "stormbreaker_no_censoring")

  d <- make_sb_fixture()
  d$event <- 0L
  event_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "gene_a", "month"
    )
  )
  expect_s3_class(event_error, "stormbreaker_no_events")

  d <- make_sb_fixture()
  d$gene_text <- as.character(d$gene_a)
  predictor_error <- capture_sb_error(
    sb_surv_task(
      d, "patient_id", "time", "event", "gene_text", "month"
    )
  )
  expect_s3_class(predictor_error, "stormbreaker_bad_predictor_type")
})
test_that("survival task validation and subsetting retain patient grouping", {
  d <- sb_example_survival(40L, 3756L)
  d$group <- rep(1:20, each = 2)
  task <- sb_surv_task(d, "patient_id", "time", "event", "gene_a", "month",
    prediction_horizon = 6, decision_thresholds = c(0.1, 0.5), group_col = "group")
  expect_identical(sb_validate_task_contract(task), task)
  subset <- sb_subset_task(task, 1:30)
  expect_identical(subset$group_col, "group")
  expect_identical(subset$data$group, d$group[1:30])
  expect_identical(subset$prediction_horizon, 6)
  for (column in c("patient_id", "time", "event", "gene_a")) {
    error <- capture_sb_error(sb_surv_task(
      d, "patient_id", "time", "event", "gene_a", "month", group_col = column
    ))
    expect_s3_class(error, "stormbreaker_overlapping_group_roles")
  }
})
