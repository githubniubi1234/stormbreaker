make_stepcox_fixture <- function(n = 120L, seed = 2401L) {
  set.seed(seed)
  signal_a <- stats::rnorm(n)
  signal_b <- stats::rnorm(n)
  noise <- stats::rnorm(n)
  event_time <- stats::rexp(
    n,
    rate = exp(1.5 * signal_a - 1.2 * signal_b)
  )
  censor_time <- stats::rexp(n, rate = 0.12)
  data.frame(
    patient_id = sprintf("stepcox_%03d", seq_len(n)),
    time = pmin(event_time, censor_time),
    event = as.integer(event_time <= censor_time),
    signal_a = signal_a,
    signal_b = signal_b,
    noise = noise,
    stringsAsFactors = FALSE
  )
}

test_that("StepCox directions select from analysis outcomes and preserve names", {
  data <- make_stepcox_fixture()
  x <- data[c("signal_a", "signal_b", "noise")]
  names(x) <- c("signal A", "signal-B", "noise")

  for (direction in c("forward", "backward", "both")) {
    adapter <- sb_stepcox_selector_adapter(direction)
    fitted <- adapter$fit(x, data$time, data$event)
    selected <- adapter$apply(fitted, x)

    expect_s3_class(fitted, "stormbreaker_selector_fit")
    expect_identical(fitted$id, "stepcox")
    expect_identical(fitted$direction, direction)
    expect_identical(fitted$source_features, names(x))
    expect_setequal(
      intersect(fitted$selected_features, c("signal A", "signal-B")),
      c("signal A", "signal-B")
    )
    expect_identical(names(selected), fitted$selected_features)
    expect_identical(selected, x[fitted$selected_features])
    expect_identical(adapter$params, list(direction = direction))
    expect_identical(adapter$needs_outcome, TRUE)
    expect_identical(fitted$scores$feature, names(x))
    expect_identical(
      fitted$scores$selected,
      fitted$scores$feature %in% fitted$selected_features
    )
    expect_true(all(is.finite(
      fitted$scores$coefficient[fitted$scores$selected]
    )))
    expect_true(is.finite(fitted$final_aic))
    expect_s3_class(fitted$selection_path, "data.frame")
    expect_false(any(grepl(
      ".stormbreaker_stepcox_",
      fitted$selection_path$Step,
      fixed = TRUE
    )))
  }

  fitted <- sb_stepcox_selector_adapter("forward")$fit(
    x,
    data$time,
    data$event
  )
  schema_error <- capture_sb_error(sb_apply_stepcox(
    fitted,
    x[c("noise", "signal A", "signal-B")]
  ))
  expect_s3_class(schema_error, "stormbreaker_selector_schema_mismatch")
})

test_that("StepCox pipeline parameters are explicit and task-specific", {
  pipelines <- lapply(c("forward", "backward", "both"), function(direction) {
    sb_pipeline(
      selector = "stepcox",
      learner = "coxph",
      stepcox_direction = direction
    )
  })

  expect_identical(
    vapply(pipelines, function(pipeline) pipeline$id, character(1)),
    paste0("stepcox__coxph__direction", c("forward", "backward", "both"))
  )
  expect_identical(
    lapply(pipelines, function(pipeline) pipeline$selector_params),
    lapply(
      c("forward", "backward", "both"),
      function(direction) list(direction = direction)
    )
  )

  missing <- capture_sb_error(sb_pipeline(
    selector = "stepcox",
    learner = "coxph"
  ))
  invalid <- capture_sb_error(sb_pipeline(
    selector = "stepcox",
    learner = "coxph",
    stepcox_direction = "sideways"
  ))
  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "coxph",
    stepcox_direction = "forward"
  ))
  incompatible <- capture_sb_error(sb_pipeline(
    selector = "stepcox",
    learner = "logistic",
    stepcox_direction = "forward"
  ))

  expect_s3_class(missing, "stormbreaker_missing_selector_parameter")
  expect_s3_class(invalid, "stormbreaker_bad_stepcox_direction")
  expect_s3_class(irrelevant, "stormbreaker_irrelevant_selector_parameter")
  expect_s3_class(incompatible, "stormbreaker_incompatible_binary_selector")
})

test_that("StepCox canonicalizes named and classed direction values", {
  data <- make_stepcox_fixture(n = 120L, seed = 2405L)
  x <- data[c("signal_a", "signal_b", "noise")]
  plain_direction <- "forward"
  named_direction <- c(requested = "forward")
  classed_direction <- structure(
    "forward",
    class = c("stormbreaker_test_direction", "character")
  )

  plain_pipeline <- sb_pipeline(
    selector = "stepcox",
    learner = "coxph",
    stepcox_direction = plain_direction
  )
  named_pipeline <- sb_pipeline(
    selector = "stepcox",
    learner = "coxph",
    stepcox_direction = named_direction
  )
  classed_pipeline <- sb_pipeline(
    selector = "stepcox",
    learner = "coxph",
    stepcox_direction = classed_direction
  )
  plain_fit <- sb_fit_stepcox(
    x,
    data$time,
    data$event,
    plain_direction
  )
  named_fit <- sb_fit_stepcox(
    x,
    data$time,
    data$event,
    named_direction
  )
  classed_fit <- sb_fit_stepcox(
    x,
    data$time,
    data$event,
    classed_direction
  )

  expect_identical(named_pipeline, plain_pipeline)
  expect_identical(classed_pipeline, plain_pipeline)
  expect_identical(named_fit, plain_fit)
  expect_identical(classed_fit, plain_fit)
  expect_identical(plain_fit$start_model, "null")
})

test_that("StepCox runs fold-locally and supports final Cox prediction", {
  data <- make_stepcox_fixture(n = 150L, seed = 2411L)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    c("signal_a", "signal_b", "noise"),
    "month"
  )
  pipeline <- sb_pipeline(
    selector = "stepcox",
    learner = "coxph",
    stepcox_direction = "forward"
  )

  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 2412L)
  fit <- sb_fit_final(task, pipeline)
  new_data <- data[1:12, c(
    "patient_id",
    fit$required_features
  )]
  prediction <- sb_predict_final(fit, new_data, "patient_id")

  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_identical(result$selection_stability$selector, "stepcox")
  expect_identical(result$selection_stability$status, "success")
  expect_true(all(vapply(result$selections, function(selection) {
    selection$status == "success" &&
      all(c("signal_a", "signal_b") %in% selection$selected_features)
  }, logical(1))))
  expect_identical(fit$status, "success")
  expect_identical(fit$selection$state$direction, "forward")
  expect_setequal(
    intersect(fit$required_features, c("signal_a", "signal_b")),
    c("signal_a", "signal_b")
  )
  expect_identical(prediction$status, "success")
  expect_identical(
    prediction$predictions$patient_id,
    new_data$patient_id
  )
  expect_true(all(is.finite(prediction$predictions$risk)))
})

test_that("outer assessment outcomes cannot change StepCox selection", {
  data <- make_stepcox_fixture(n = 120L, seed = 2415L)
  split <- list(
    fold = 1L,
    analysis_rows = 1:90,
    assessment_rows = 91:120,
    analysis_ids = data$patient_id[1:90],
    assessment_ids = data$patient_id[91:120]
  )
  pipeline <- sb_pipeline(
    selector = "stepcox",
    learner = "coxph",
    stepcox_direction = "forward"
  )
  make_task <- function(frame) {
    sb_surv_task(
      frame,
      "patient_id",
      "time",
      "event",
      c("signal_a", "signal_b", "noise"),
      "month"
    )
  }

  original <- sb_run_one_fold(
    make_task(data),
    pipeline,
    sb_resolve_pipeline_runtime(pipeline),
    split,
    seed = 2416L
  )
  changed <- data
  changed$time[91:120] <- changed$time[91:120] + 1000
  changed$event[91:120] <- 1L - changed$event[91:120]
  perturbed <- sb_run_one_fold(
    make_task(changed),
    pipeline,
    sb_resolve_pipeline_runtime(pipeline),
    split,
    seed = 2416L
  )

  expect_identical(
    original$selection$selected_features,
    perturbed$selection$selected_features
  )
  expect_identical(original$selection$scores, perturbed$selection$scores)
  expect_equal(
    original$predictions$risk,
    perturbed$predictions$risk,
    tolerance = 0
  )
})

test_that("StepCox is repeated inside glmnet inner analysis folds", {
  skip_if_not_installed("glmnet")

  data <- make_stepcox_fixture(n = 150L, seed = 2421L)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    c("signal_a", "signal_b", "noise"),
    "month"
  )
  pipeline <- sb_pipeline(
    selector = "stepcox",
    learner = "glmnet_cox",
    stepcox_direction = "backward",
    glmnet_alpha = 0,
    inner_folds = 3L,
    inner_seed = 2422L
  )

  fit <- sb_fit_final(task, pipeline)

  expect_identical(fit$status, "success")
  expect_identical(
    fit$tuning$inner_selection_method,
    "glmnet_exclude_function"
  )
  expect_length(fit$tuning$inner_selections, 3L)
  expect_true(all(vapply(fit$tuning$inner_selections, function(selection) {
    length(selection$analysis_patient_ids) < nrow(data) &&
      all(c("signal_a", "signal_b") %in% selection$selected_features)
  }, logical(1))))
})
