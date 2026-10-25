test_that("binary endpoints retain only patient ID and outcome", {
  expect_identical(
    names(formals(sb_binary_endpoint)),
    c("data", "id_col", "outcome_col")
  )
  data <- data.frame(
    patient_id = c("e1", "e2", "e3"),
    outcome = c(1L, 1L, 0L),
    gene_a = c(10, 20, 30),
    note = c("a", "b", "c"),
    stringsAsFactors = FALSE
  )

  endpoint <- sb_binary_endpoint(data, "patient_id", "outcome")

  expect_s3_class(endpoint, "stormbreaker_binary_endpoint")
  expect_identical(names(endpoint$data), c("patient_id", "outcome"))
  expect_identical(endpoint$data, data[c("patient_id", "outcome")])
  expect_identical(endpoint$n, 3L)
  expect_identical(endpoint$n_positive, 2L)
  expect_identical(endpoint$n_negative, 1L)
  expect_identical(endpoint$positive_class, 1L)
  expect_identical(endpoint$negative_class, 0L)
  expect_identical(endpoint$independent_unit, "patient")
})

test_that("binary endpoints record one patient and one observed class", {
  one_positive <- data.frame(patient_id = "e1", outcome = 1L)
  one_negative_class <- data.frame(
    patient_id = c("e1", "e2"),
    outcome = c(0L, 0L)
  )

  positive <- sb_binary_endpoint(
    one_positive,
    "patient_id",
    "outcome"
  )
  negative <- sb_binary_endpoint(
    one_negative_class,
    "patient_id",
    "outcome"
  )

  expect_identical(positive$n, 1L)
  expect_identical(positive$n_positive, 1L)
  expect_identical(positive$n_negative, 0L)
  expect_identical(negative$n_positive, 0L)
  expect_identical(negative$n_negative, 2L)
})

test_that("logical binary endpoints retain explicit zero-one semantics", {
  data <- data.frame(
    patient_id = c("e1", "e2"),
    outcome = c(TRUE, FALSE)
  )

  endpoint <- sb_binary_endpoint(data, "patient_id", "outcome")

  expect_identical(endpoint$data$outcome, c(TRUE, FALSE))
  expect_identical(endpoint$n_positive, 1L)
  expect_identical(endpoint$n_negative, 1L)
})

test_that("survival endpoints retain only ID time and event semantics", {
  expect_identical(
    names(formals(sb_surv_endpoint)),
    c(
      "data", "id_col", "time_col", "event_col", "time_unit",
      "prediction_horizon"
    )
  )
  data <- data.frame(
    patient_id = c("e1", "e2", "e3"),
    follow_up = c(2, 5, 7),
    died = c(1L, 0L, 1L),
    gene_a = c(4, 5, 6),
    stringsAsFactors = FALSE
  )

  endpoint <- sb_surv_endpoint(
    data,
    "patient_id",
    "follow_up",
    "died",
    "month",
    prediction_horizon = 6
  )

  expect_s3_class(endpoint, "stormbreaker_surv_endpoint")
  expect_identical(
    names(endpoint$data),
    c("patient_id", "follow_up", "died")
  )
  expect_identical(
    endpoint$data,
    data[c("patient_id", "follow_up", "died")]
  )
  expect_identical(endpoint$n, 3L)
  expect_identical(endpoint$n_events, 2L)
  expect_identical(endpoint$n_censored, 1L)
  expect_identical(endpoint$time_unit, "month")
  expect_identical(endpoint$prediction_horizon, 6)
  expect_identical(endpoint$independent_unit, "patient")
})

test_that("survival endpoints allow one row all events and all censoring", {
  one <- data.frame(patient_id = "e1", time = 2, event = 1L)
  all_events <- data.frame(
    patient_id = paste0("e", 1:3),
    time = c(1, 2, 3),
    event = 1L
  )
  all_censored <- data.frame(
    patient_id = paste0("c", 1:3),
    time = c(4, 5, 6),
    event = 0L
  )

  one_endpoint <- sb_surv_endpoint(
    one, "patient_id", "time", "event", "day"
  )
  event_endpoint <- sb_surv_endpoint(
    all_events, "patient_id", "time", "event", "month"
  )
  censored_endpoint <- sb_surv_endpoint(
    all_censored, "patient_id", "time", "event", "year"
  )

  expect_identical(one_endpoint$n, 1L)
  expect_identical(event_endpoint$n_events, 3L)
  expect_identical(event_endpoint$n_censored, 0L)
  expect_identical(censored_endpoint$n_events, 0L)
  expect_identical(censored_endpoint$n_censored, 3L)
})

test_that("endpoint constructors require one nonambiguous endpoint row", {
  empty_binary <- data.frame(
    patient_id = character(),
    outcome = integer()
  )
  empty_error <- capture_sb_error(
    sb_binary_endpoint(empty_binary, "patient_id", "outcome")
  )
  expect_s3_class(empty_error, "stormbreaker_empty_endpoint")

  binary <- data.frame(patient_id = "e1", outcome = 1L)
  overlap_error <- capture_sb_error(
    sb_binary_endpoint(binary, "patient_id", "patient_id")
  )
  expect_s3_class(
    overlap_error,
    "stormbreaker_overlapping_endpoint_roles"
  )

  duplicate <- data.frame(
    patient_id = c("e1", "e2"),
    outcome = c(0L, 1L),
    second = c(1L, 0L),
    check.names = FALSE
  )
  names(duplicate)[3L] <- "outcome"
  duplicate_error <- capture_sb_error(
    sb_binary_endpoint(duplicate, "patient_id", "outcome")
  )
  expect_s3_class(
    duplicate_error,
    "stormbreaker_duplicate_endpoint_column"
  )

  survival <- data.frame(patient_id = "e1", time = 1, event = 1L)
  survival_role_error <- capture_sb_error(
    sb_surv_endpoint(
      survival, "patient_id", "time", "time", "month"
    )
  )
  expect_s3_class(
    survival_role_error,
    "stormbreaker_overlapping_endpoint_roles"
  )
})

test_that("endpoint patient identifiers remain complete and unique", {
  duplicate <- data.frame(
    patient_id = c("e1", "e1"),
    outcome = c(0L, 1L)
  )
  duplicate_error <- capture_sb_error(
    sb_binary_endpoint(duplicate, "patient_id", "outcome")
  )
  expect_s3_class(duplicate_error, "stormbreaker_duplicate_id")

  blank <- data.frame(patient_id = " ", outcome = 0L)
  blank_error <- capture_sb_error(
    sb_binary_endpoint(blank, "patient_id", "outcome")
  )
  expect_s3_class(blank_error, "stormbreaker_bad_id")
})

test_that("endpoint values and time declarations are explicit", {
  binary <- data.frame(patient_id = "e1", outcome = 2L)
  binary_error <- capture_sb_error(
    sb_binary_endpoint(binary, "patient_id", "outcome")
  )
  expect_s3_class(binary_error, "stormbreaker_bad_binary_outcome")

  survival <- data.frame(patient_id = "e1", time = 1, event = 1L)
  bad_time <- survival
  bad_time$time <- -1
  time_error <- capture_sb_error(
    sb_surv_endpoint(
      bad_time, "patient_id", "time", "event", "month"
    )
  )
  expect_s3_class(time_error, "stormbreaker_bad_time")

  bad_event <- survival
  bad_event$event <- 2L
  event_error <- capture_sb_error(
    sb_surv_endpoint(
      bad_event, "patient_id", "time", "event", "month"
    )
  )
  expect_s3_class(event_error, "stormbreaker_bad_event")

  unit_error <- capture_sb_error(
    sb_surv_endpoint(
      survival, "patient_id", "time", "event", "months"
    )
  )
  expect_s3_class(unit_error, "stormbreaker_bad_time_unit")

  for (horizon in list(0, -1, NA_real_, Inf, "1", c(1, 2))) {
    horizon_error <- capture_sb_error(
      sb_surv_endpoint(
        survival,
        "patient_id",
        "time",
        "event",
        "month",
        prediction_horizon = horizon
      )
    )
    expect_s3_class(
      horizon_error,
      "stormbreaker_bad_prediction_horizon"
    )
  }
})

test_that("endpoint contract validation detects changed endpoint metadata", {
  binary <- sb_binary_endpoint(
    data.frame(patient_id = c("e1", "e2"), outcome = c(0L, 1L)),
    "patient_id",
    "outcome"
  )
  binary$n_positive <- 2L
  binary_error <- capture_sb_error(
    sb_validate_binary_endpoint_contract(binary)
  )
  expect_s3_class(binary_error, "stormbreaker_mutated_endpoint")

  survival <- sb_surv_endpoint(
    data.frame(
      patient_id = c("e1", "e2"),
      time = c(1, 2),
      event = c(1L, 0L)
    ),
    "patient_id",
    "time",
    "event",
    "month"
  )
  survival$n_events <- 2L
  survival_error <- capture_sb_error(
    sb_validate_surv_endpoint_contract(survival)
  )
  expect_s3_class(survival_error, "stormbreaker_mutated_endpoint")
})

test_that("endpoint print methods describe outcome-only objects", {
  binary <- sb_binary_endpoint(
    data.frame(patient_id = "e1", outcome = 1L),
    "patient_id",
    "outcome"
  )
  survival <- sb_surv_endpoint(
    data.frame(patient_id = "e1", time = 1, event = 1L),
    "patient_id",
    "time",
    "event",
    "month",
    prediction_horizon = 2
  )

  binary_output <- paste(capture.output(print(binary)), collapse = "\n")
  survival_output <- paste(
    capture.output(print(survival)),
    collapse = "\n"
  )

  expect_match(binary_output, "outcome-only evaluation endpoint")
  expect_match(binary_output, "Patients: 1")
  expect_match(survival_output, "outcome-only evaluation endpoint")
  expect_match(survival_output, "Prediction horizon: 2 \\[ month \\]")
  expect_no_match(binary_output, "fitted")
  expect_no_match(survival_output, "validated")
})
