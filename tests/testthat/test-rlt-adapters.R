make_rlt_binary_fixture <- function(n = 100L, p = 6L, seed = 43101L) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    eta <- x[, 1L] - 0.7 * x[, min(2L, p)]
    outcome <- stats::rbinom(n, 1L, stats::plogis(eta))
    feature_names <- c(
      "gene one",
      if (p >= 2L) "2gene",
      if (p > 2L) paste0("gene-", seq_len(p - 2L))
    )
    data.frame(
      patient_id = sprintf("RB%03d", seq_len(n)),
      outcome = outcome,
      x,
      check.names = FALSE
    ) |>
      stats::setNames(c("patient_id", "outcome", feature_names))
  })
}

make_rlt_survival_fixture <- function(n = 120L, p = 6L, seed = 43201L) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    eta <- 0.7 * x[, 1L] - 0.5 * x[, min(2L, p)]
    event_time <- stats::rexp(n, rate = exp(eta))
    censor_time <- stats::rexp(n, rate = 0.45)
    time <- pmin(event_time, censor_time) + 0.01
    feature_names <- c(
      "gene one",
      if (p >= 2L) "2gene",
      if (p > 2L) paste0("gene-", seq_len(p - 2L))
    )
    data.frame(
      patient_id = sprintf("RS%03d", seq_len(n)),
      time = time,
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    ) |>
      stats::setNames(c("patient_id", "time", "event", feature_names))
  })
}

rlt_native_fit <- function(x, endpoint, params) {
  matrix_x <- as.matrix(x)
  storage.mode(matrix_x) <- "double"
  arguments <- list(
    x = matrix_x,
    model = params$model,
    print.summary = 0L,
    use.cores = params$cores,
    ntrees = params$num_trees,
    nmin = params$native_nmin,
    reinforcement = params$reinforcement,
    combsplit = params$combination_split_size,
    embed.n.th = params$embedded_stop_size,
    importance = params$importance,
    track.obs = params$track_observations
  )
  if (identical(params$model, "classification")) {
    arguments$y <- factor(
      endpoint$outcome,
      levels = c(0L, 1L),
      labels = c("control", "case")
    )
  } else {
    arguments$y <- endpoint$time
    arguments$censor <- endpoint$event
  }
  withr::with_seed(params$seed, do.call(RLT::RLT, arguments))
}

test_that("two RLT reinforcement learners expose their legal routes", {
  expect_identical(
    sb_rlt_learner_ids(),
    c(
      "rlt_reinforcement_binary",
      "rlt_reinforcement_survival"
    )
  )
  for (learner in sb_rlt_learner_ids()) {
    params <- sb_rlt_learner_params(learner)
    expect_identical(params$num_trees, 100L)
    expect_identical(params$reinforcement, TRUE)
    expect_identical(params$combination_split_size, 1L)
    expect_identical(params$native_nmin, 5L)
    expect_identical(params$embedded_stop_size, 1L)
    expect_identical(params$cores, 1L)
    expect_identical(params$seed, 1L)
  }

  binary_selectors <- expected_fixed_learner_selectors("binary")
  survival_selectors <- expected_fixed_learner_selectors("survival")
  registry <- sb_pipeline_registry()
  rows <- registry[registry$learner %in% sb_rlt_learner_ids(), , drop = FALSE]
  expect_identical(
    nrow(rows),
    length(binary_selectors) + length(survival_selectors)
  )
  expect_true(all(rows$optional_dependency == "RLT" |
    grepl("RLT", rows$optional_dependency, fixed = TRUE)))
  expect_true(all(rows$tuning_semantics == "none"))
  expect_setequal(
    rows$selector[rows$task_type == "binary"],
    binary_selectors
  )
  expect_setequal(
    rows$selector[rows$task_type == "survival"],
    survival_selectors
  )

  for (index in seq_len(nrow(rows))) {
    spec <- list(
      selector = rows$selector[[index]],
      learner = rows$learner[[index]]
    )
    if (rows$selector[[index]] %in% c(
      "univ_binary_auc_top_k", "univ_cox_top_k"
    )) {
      spec$selector_n_features <- 3L
    }
    if (identical(rows$selector[[index]], "stepcox")) {
      spec$stepcox_direction <- "both"
    }
    if (identical(rows$selector[[index]], "fixed_features")) {
      spec$selector_features <- c("gene one", "2gene", "gene-1")
      spec$feature_set_id <- "rlt_test"
    }
    if (spec$selector %in% c("lasso_nonzero", "enet_cox_nonzero", "coxboost_nonzero")) {
      spec$selector_inner_folds <- 3L
      spec$selector_seed <- 1L
    }
    if (spec$selector == "rsf_minimal_depth") {
      spec$selector_rsf_nodesize <- 5L
      spec$selector_seed <- 1L
    }
    expect_identical(sb_lint_pipeline(spec)$status, "admitted")
  }
})

test_that("binary RLT matches native probabilities and portable replay", {
  skip_if_not_installed("RLT", minimum_version = "3.2.6")

  data <- make_rlt_binary_fixture()
  features <- setdiff(names(data), c("patient_id", "outcome"))
  x <- data[features]
  assessment <- x[c(2L, 19L, 61L), , drop = FALSE]
  adapter <- sb_rlt_learner_adapter("rlt_reinforcement_binary")

  set.seed(43102L)
  caller_rng <- .Random.seed
  model <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, caller_rng)
  native <- rlt_native_fit(
    x,
    list(outcome = data$outcome),
    adapter$params
  )
  native_probability <- as.matrix(
    stats::predict(native, assessment)$ProbPrediction
  )[, 2L]
  probability <- adapter$predict(model, assessment)
  expect_identical(probability, unname(native_probability))
  expect_true(all(probability >= 0 & probability <= 1))
  expect_identical(model$engine$ylevels, c("control", "case"))
  expect_false(any(c(
    "AllError", "VarImp", "ObsTrack", "n", "subjectnames",
    "subject.weight", "x", "y", "outcome"
  ) %in% names(model$engine)))

  state <- unserialize(serialize(adapter$extract_state(model), NULL))
  compact <- unserialize(state$engine_raw)
  replay <- adapter$predict_state(state, assessment)
  expect_identical(replay$probability, probability)
  expect_false(any(c(
    "AllError", "VarImp", "ObsTrack", "n", "subjectnames",
    "subject.weight", "x", "y", "outcome"
  ) %in% names(compact)))
  expect_error(
    adapter$predict(model, assessment[rev(names(assessment))]),
    class = "stormbreaker_learner_schema_mismatch"
  )
})

test_that("survival RLT reconstructs native hazards and direct survival", {
  skip_if_not_installed("RLT", minimum_version = "3.2.6")

  data <- make_rlt_survival_fixture()
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[c(3L, 27L, 92L), , drop = FALSE]
  adapter <- sb_rlt_learner_adapter("rlt_reinforcement_survival")

  set.seed(43202L)
  caller_rng <- .Random.seed
  model <- adapter$fit(x, data$time, data$event)
  expect_identical(.Random.seed, caller_rng)
  native <- rlt_native_fit(
    x,
    list(time = data$time, event = data$event),
    adapter$params
  )
  native_prediction <- stats::predict(native, assessment)
  hazard <- as.matrix(native_prediction$SurvPred)
  survival <- matrix(NA_real_, nrow(hazard), ncol(hazard))
  for (row in seq_len(nrow(hazard))) {
    survival[row, ] <- cumprod(1 - hazard[row, ])
  }
  interval <- diff(c(0, native_prediction$timepoints))
  previous <- cbind(1, survival[, -ncol(survival), drop = FALSE])
  native_risk <- -as.vector(previous %*% interval)
  curves <- sb_predict_rlt_survival_curves(model, assessment)
  expect_identical(curves$hazard, hazard)
  expect_identical(curves$survival, survival)
  expect_identical(curves$risk, unname(native_risk))
  expect_true(all(apply(
    curves$survival,
    1L,
    function(values) all(diff(values) <= 0)
  )))

  early_horizon <- min(curves$time) / 2
  early <- adapter$predict_survival(
    model, x, data$time, data$event, assessment, early_horizon
  )
  later_horizon <- stats::median(curves$time)
  later <- adapter$predict_survival(
    model, x, data$time, data$event, assessment, later_horizon
  )
  expect_identical(early$probability, rep(1, nrow(assessment)))
  expect_true(all(later$probability >= 0 & later$probability <= 1))
  expect_identical(later$source, sb_rlt_direct_survival_source())
  expect_identical(later$ties, "not_applicable")

  state <- unserialize(serialize(adapter$extract_state(model, list(
    x = x,
    time = data$time,
    event = data$event,
    prediction_horizon = later_horizon,
    time_unit = "month"
  )), NULL))
  replay <- adapter$predict_state(state, assessment)
  expect_identical(replay$risk, curves$risk)
  expect_identical(replay$survival_probability, later$probability)
  expect_identical(
    state$risk_formula,
    "negative_step_RMST_over_native_event_grid"
  )
  compact <- unserialize(state$engine_raw)
  expect_false(any(c(
    "AllError", "VarImp", "ObsTrack", "n", "subjectnames",
    "subject.weight", "x", "y", "outcome", "time", "event"
  ) %in% names(compact)))
})

test_that("RLT native-domain failures remain typed and do not enter C code", {
  skip_if_not_installed("RLT", minimum_version = "3.2.6")

  data <- make_rlt_survival_fixture(n = 30L, p = 2L, seed = 43301L)
  x <- data[setdiff(names(data), c("patient_id", "time", "event"))]
  adapter <- sb_rlt_learner_adapter("rlt_reinforcement_survival")
  event <- integer(nrow(data))
  event[[1L]] <- 1L
  expect_error(
    adapter$fit(x, data$time, event),
    class = "stormbreaker_rlt_insufficient_event_times"
  )
  event[[2L]] <- 1L
  tied_time <- data$time
  tied_time[1:2] <- tied_time[[1L]]
  expect_error(
    adapter$fit(x, tied_time, event),
    class = "stormbreaker_rlt_insufficient_event_times"
  )
  zero_time <- data$time
  zero_time[[1L]] <- 0
  expect_error(
    adapter$fit(x, zero_time, data$event),
    class = "stormbreaker_bad_rlt_time"
  )
  expect_error(
    sb_rlt_stop_on_warning(
      {
        warning("Subject weights length must be n, reset to equal weights")
        TRUE
      },
      "rlt_reinforcement_binary",
      "fitting"
    ),
    class = "stormbreaker_rlt_warning"
  )
})

test_that("RLT dependency failure is explicit", {
  local_mocked_bindings(sb_rlt_available = function() FALSE)
  expect_error(
    sb_require_rlt_dependency(),
    class = "stormbreaker_missing_rlt"
  )
})

test_that("RLT executes through CV, final fit, and outcome-blind prediction", {
  skip_if_not_installed("RLT", minimum_version = "3.2.6")

  binary <- make_rlt_binary_fixture(n = 90L, seed = 43401L)
  binary_features <- setdiff(names(binary), c("patient_id", "outcome"))
  binary_task <- sb_binary_task(
    binary,
    "patient_id",
    "outcome",
    binary_features
  )
  binary_pipeline <- sb_pipeline(
    learner = "rlt_reinforcement_binary"
  )
  binary_cv <- sb_run_cv(
    binary_task,
    binary_pipeline,
    folds = 3L,
    seed = 43402L
  )
  binary_fit <- sb_fit_final(binary_task, binary_pipeline)
  binary_prediction <- sb_predict_final(
    unserialize(serialize(binary_fit, NULL)),
    binary[1:6, c("patient_id", rev(binary_features)), drop = FALSE],
    "patient_id"
  )
  rlt_version <- as.character(utils::packageVersion("RLT"))
  expect_identical(binary_cv$status, "success")
  expect_true(binary_cv$complete_oof)
  expect_identical(binary_cv$provenance$rlt_version, rlt_version)
  expect_identical(binary_fit$status, "success")
  expect_identical(binary_fit$provenance$rlt_version, rlt_version)
  expect_identical(binary_prediction$status, "success")

  survival <- make_rlt_survival_fixture(n = 105L, seed = 43403L)
  survival_features <- setdiff(
    names(survival),
    c("patient_id", "time", "event")
  )
  horizon <- as.numeric(stats::quantile(survival$time, 0.45))
  survival_task <- sb_surv_task(
    survival,
    "patient_id",
    "time",
    "event",
    survival_features,
    "month",
    prediction_horizon = horizon
  )
  survival_pipeline <- sb_pipeline(
    learner = "rlt_reinforcement_survival"
  )
  survival_cv <- sb_run_cv(
    survival_task,
    survival_pipeline,
    folds = 3L,
    seed = 43404L
  )
  survival_fit <- sb_fit_final(survival_task, survival_pipeline)
  survival_prediction <- sb_predict_final(
    unserialize(serialize(survival_fit, NULL)),
    survival[1:6, c("patient_id", rev(survival_features)), drop = FALSE],
    "patient_id"
  )
  expect_identical(survival_cv$status, "success")
  expect_true(survival_cv$complete_oof)
  expect_true(survival_cv$complete_survival_probability_oof)
  expect_identical(survival_cv$provenance$rlt_version, rlt_version)
  expect_identical(survival_fit$status, "success")
  expect_identical(survival_fit$provenance$rlt_version, rlt_version)
  expect_identical(survival_prediction$status, "success")
  expect_true(all(is.finite(survival_prediction$predictions$risk)))
  expect_true(all(
    survival_prediction$predictions$survival_probability >= 0 &
      survival_prediction$predictions$survival_probability <= 1
  ))
})

test_that("RLT ordinary warnings preserve native fit and predictions", {
  skip_if_not_installed("RLT")
  d <- make_rlt_binary_fixture(n = 40L, p = 3L)
  x <- d[setdiff(names(d), c("patient_id", "outcome"))]
  id <- "rlt_reinforcement_binary"
  params <- sb_rlt_learner_params(id)
  params$num_trees <- 5L
  adapter <- sb_rlt_learner_adapter(id, params)
  reference <- adapter$fit(x, d$outcome)
  expected <- adapter$predict(reference, x)
  native_fit <- RLT::RLT
  native_predict <- get("predict.RLT", asNamespace("RLT"))
  local_mocked_bindings(RLT = function(...) {
    warning("informational fit warning")
    native_fit(...)
  }, .package = "RLT")
  local_mocked_bindings(predict.RLT = function(...) {
    warning("informational prediction warning")
    native_predict(...)
  }, .package = "RLT")
  seen <- character()
  observed <- withCallingHandlers({
    fitted <- adapter$fit(x, d$outcome)
    adapter$predict(fitted, x)
  }, warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(observed, expected)
  expect_identical(all(c("informational fit warning", "informational prediction warning") %in% seen), TRUE)
})

test_that("RLT retains weight reset and feature mismatch protection", {
  messages <- c(
    "Subject weights length must be n, reset to equal weights",
    "Subject weights cannot be 0 or negative, reset to equal weights",
    "Variable weights length must by p, reset to equal weights",
    "Variable weights cannot be negative, reset to equal weights",
    "test data variables names does not match training data...")
  for (message in messages) {
    error <- tryCatch(sb_rlt_stop_on_warning({ warning(message); TRUE },
      "rlt_reinforcement_binary", "fitting"), error = identity)
    expect_s3_class(error, "stormbreaker_rlt_warning")
  }
})
