make_ahaz_bujar_fixture <- function(
  n = 90L, p = 6L, seed = 601L, tied = FALSE,
  heavy_censoring = FALSE, signal = c(0.8, -0.6, 0.35)
) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    beta <- rep(0, p)
    beta[seq_len(min(length(signal), p))] <- signal[
      seq_len(min(length(signal), p))
    ]
    risk <- drop(x %*% beta)
    event_time <- stats::rexp(n, exp(risk) / 10)
    censor_time <- stats::rexp(
      n,
      if (heavy_censoring) 0.16 else 0.04
    )
    time <- pmin(event_time, censor_time)
    if (tied) {
      time <- pmax(0.25, round(time * 2) / 2)
    } else {
      time <- time + seq_len(n) * 1e-10
    }
    names <- c(
      "gene one", "2gene", "if", "gene-4",
      paste0("gene_", seq_len(max(0L, p - 4L)))
    )[seq_len(p)]
    data.frame(
      patient_id = paste0("PRIVATE_PATIENT_", seq_len(n)),
      time = time,
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    ) |>
      stats::setNames(c("patient_id", "time", "event", names))
  })
}

ahaz_bujar_features <- function(data) {
  setdiff(names(data), c("patient_id", "time", "event"))
}

test_that("four additive-hazard and Buckley-James learners add 32 routes", {
  ids <- c(sb_ahaz_survival_learner_ids(), sb_bujar_linear_learner_ids())
  expect_identical(ids, c(
    "ahaz_lin_ying",
    "ahaz_lasso_bic",
    "ahaz_sscad_bic",
    "bujar_linear_boost"
  ))
  rows <- sb_pipeline_registry()
  rows <- rows[rows$learner %in% ids, , drop = FALSE]
  expect_identical(nrow(rows), 4L * length(expected_fixed_learner_selectors("survival")))
  for (learner in ids) {
    expect_setequal(
      rows$selector[rows$learner == learner],
      expected_fixed_learner_selectors("survival")
    )
  }
  expect_true(all(
    rows$tuning_semantics[
      rows$learner %in% c("ahaz_lasso_bic", "ahaz_sscad_bic")
    ] == "analysis_fold_bic"
  ))
  expect_true(all(
    rows$tuning_semantics[
      rows$learner %in% c(
        "ahaz_lin_ying",
        "bujar_linear_boost"
      )
    ] == "none"
  ))
  expect_route_dependencies(rows, c(ahaz_lin_ying = "ahaz",
    ahaz_lasso_bic = "ahaz", ahaz_sscad_bic = "ahaz", bujar_linear_boost = "bujar"))
})

test_that("ahaz learners reproduce native coefficients and BIC choices", {
  skip_if_not_installed("ahaz", minimum_version = "1.15.1")
  data <- make_ahaz_bujar_fixture()
  features <- ahaz_bujar_features(data)
  x <- data[features]
  assessment <- x[c(2L, 31L, 76L), , drop = FALSE]
  engine_x <- sb_ahaz_engine_x(x)
  surv <- survival::Surv(data$time, data$event)

  for (learner in sb_ahaz_survival_learner_ids()) {
    adapter <- sb_ahaz_survival_learner_adapter(learner)
    model <- adapter$fit(x, data$time, data$event)
    native <- sb_fit_ahaz_engine(surv, engine_x, learner, adapter$params)
    native_coefficient <- unname(as.numeric(stats::coef(native)))
    expect_equal(model$coefficient, native_coefficient, tolerance = 0)
    expect_equal(
      adapter$predict(model, assessment),
      unname(as.vector(as.matrix(assessment) %*% native_coefficient)),
      tolerance = 0
    )
    if (identical(learner, "ahaz_lin_ying")) {
      expect_true(is.na(model$audit$selected_lambda))
    } else {
      expect_equal(
        model$audit$selected_lambda,
        as.numeric(native$lambda.min),
        tolerance = 0
      )
      expect_identical(
        model$audit$path_length,
        as.integer(length(native$lambda))
      )
    }
  }
})

test_that("ahaz horizon output is analysis-fold Cox calibrated and portable", {
  skip_if_not_installed("ahaz", minimum_version = "1.15.1")
  data <- make_ahaz_bujar_fixture(seed = 602L)
  features <- ahaz_bujar_features(data)
  x <- data[features]
  rownames(x) <- data$patient_id
  assessment <- x[c(3L, 40L, 81L), , drop = FALSE]
  horizon <- as.numeric(stats::median(data$time[data$event == 1L]))

  for (learner in sb_ahaz_survival_learner_ids()) {
    adapter <- sb_ahaz_survival_learner_adapter(learner)
    model <- adapter$fit(x, data$time, data$event)
    probability <- adapter$predict_survival(
      model,
      x,
      data$time,
      data$event,
      assessment,
      horizon
    )
    expect_identical(
      probability$source,
      paste0(learner, "_analysis_cox_breslow_calibration")
    )
    expect_identical(probability$ties, "breslow")
    expect_gte(probability$calibration_coefficient, 0)
    expect_true(all(probability$probability >= 0 &
      probability$probability <= 1))

    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = horizon,
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, assessment)
    expect_equal(
      replay$risk,
      adapter$predict(model, assessment),
      tolerance = 0
    )
    expect_equal(
      replay$survival_probability,
      probability$probability,
      tolerance = 0
    )
    expect_false(any(c(
      "engine", "x", "X", "y", "outcome", "time", "event"
    ) %in% names(state)))
    expect_length(
      grepRaw(
        "PRIVATE_PATIENT_",
        serialize(state, NULL),
        fixed = TRUE,
        all = TRUE
      ),
      0L
    )
  }
})

test_that("ahaz resolves event ties while retaining other native limitations", {
  skip_if_not_installed("ahaz", minimum_version = "1.15.1")
  tied <- make_ahaz_bujar_fixture(tied = TRUE, seed = 603L)
  tied_x <- tied[ahaz_bujar_features(tied)]
  original <- tied
  for (learner in sb_ahaz_survival_learner_ids()) {
    adapter <- sb_ahaz_survival_learner_adapter(learner)
    model <- adapter$fit(tied_x, tied$time, tied$event)
    fit_time <- sb_ahaz_fit_times(tied$time, tied$event)
    native <- sb_fit_ahaz_engine(
      survival::Surv(fit_time$time, tied$event),
      sb_ahaz_engine_x(tied_x), learner, adapter$params
    )
    expect_equal(model$coefficient, unname(as.numeric(stats::coef(native))))
    expect_equal(all(is.finite(adapter$predict(model, tied_x))), TRUE)
    expect_identical(tied, original)
    expect_identical(model$audit$time_ties, fit_time$audit)
    risk <- adapter$predict(model, tied_x)
    horizon <- stats::median(tied$time)
    expected <- sb_fit_ahaz_calibration(tied$time, tied$event, risk, horizon)
    probability <- sb_predict_ahaz_survival_probability(
      model, tied_x, tied$time, tied$event, tied_x, horizon
    )
    expect_equal(probability$calibration_coefficient, expected$coefficient)
    expect_equal(probability$cumulative_baseline_hazard,
      expected$cumulative_baseline_hazard)
  }

  high_dimensional <- make_ahaz_bujar_fixture(
    n = 45L,
    p = 80L,
    seed = 704L
  )
  high_x <- high_dimensional[ahaz_bujar_features(high_dimensional)]
  expect_error(
    sb_ahaz_survival_learner_adapter("ahaz_lin_ying")$fit(
      high_x,
      high_dimensional$time,
      high_dimensional$event
    ),
    class = "stormbreaker_invalid_ahaz_fit"
  )

  expect_error(
    sb_fit_ahaz_calibration(
      seq_len(20),
      rep(c(1L, 0L), 10L),
      rep(0, 20L),
      10
    ),
    class = "stormbreaker_ahaz_calibration_failed"
  )
})

test_that("Bujar linear boosting is deterministic, portable, and side-effect free", {
  skip_if_not_installed("bujar", minimum_version = "0.2-11")
  data <- make_ahaz_bujar_fixture(seed = 604L)
  features <- ahaz_bujar_features(data)
  x <- data[features]
  rownames(x) <- data$patient_id
  assessment <- x[c(4L, 38L, 79L), , drop = FALSE]
  adapter <- sb_bujar_linear_learner_adapter()

  set.seed(9604L)
  caller_rng <- .Random.seed
  withr::local_dir(withr::local_tempdir())
  model <- adapter$fit(x, data$time, data$event)
  expect_identical(.Random.seed, caller_rng)
  expect_false(file.exists("file"))
  link <- unname(
    model$intercept + as.vector(as.matrix(assessment) %*% model$coefficient)
  )
  expect_equal(adapter$predict(model, assessment), -link, tolerance = 0)

  horizon <- as.numeric(stats::median(data$time[data$event == 1L]))
  probability <- adapter$predict_survival(
    model,
    x,
    data$time,
    data$event,
    assessment,
    horizon
  )
  expect_identical(
    probability$source,
    "bujar_linear_boost_analysis_residual_kaplan_meier"
  )
  expect_identical(probability$ties, "not_applicable")
  expect_true(all(probability$probability >= 0 &
    probability$probability <= 1))

  state <- unserialize(serialize(adapter$extract_state(model, list(
    x = x,
    time = data$time,
    event = data$event,
    prediction_horizon = horizon,
    time_unit = "month"
  )), NULL))
  replay <- adapter$predict_state(state, assessment)
  expect_equal(replay$risk, -link, tolerance = 0)
  expect_equal(
    replay$survival_probability,
    probability$probability,
    tolerance = 0
  )
  expect_false(any(c(
    "engine", "x", "X", "y", "outcome", "time", "event"
  ) %in% names(state)))
  expect_length(
    grepRaw(
      "PRIVATE_PATIENT_",
      serialize(state, NULL),
      fixed = TRUE,
      all = TRUE
    ),
    0L
  )
})

test_that("Bujar nonconvergence is a typed failure without a fallback", {
  skip_if_not_installed("bujar", minimum_version = "0.2-11")
  data <- make_ahaz_bujar_fixture(seed = 603L)
  x <- data[ahaz_bujar_features(data)]
  expect_error(
    sb_bujar_linear_learner_adapter()$fit(x, data$time, data$event),
    class = "stormbreaker_bujar_nonconvergence"
  )
})

test_that("new optional dependency failures are explicit", {
  local_mocked_bindings(sb_ahaz_available = function() FALSE)
  expect_error(sb_require_ahaz_dependency(), class = "stormbreaker_missing_ahaz")
  local_mocked_bindings(sb_bujar_available = function() FALSE)
  expect_error(
    sb_require_bujar_dependency(),
    class = "stormbreaker_missing_bujar"
  )
})

test_that("ahaz and Bujar run through final outcome-blind prediction", {
  skip_if_not_installed("ahaz", minimum_version = "1.15.1")
  skip_if_not_installed("bujar", minimum_version = "0.2-11")
  data <- make_ahaz_bujar_fixture(n = 100L, p = 4L, seed = 609L)
  features <- ahaz_bujar_features(data)
  horizon <- as.numeric(stats::median(data$time[data$event == 1L]))
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    features,
    "month",
    prediction_horizon = horizon
  )
  new_data <- data[1:5, c(
    "patient_id", "time", "event", rev(features)
  )]
  altered <- new_data
  altered$time <- rev(altered$time) + 1000
  altered$event <- 1L - altered$event

  for (learner in c(
    sb_ahaz_survival_learner_ids(),
    sb_bujar_linear_learner_ids()
  )) {
    pipeline <- sb_pipeline(learner = learner)
    fit <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success")
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      new_data,
      "patient_id"
    )
    altered_prediction <- sb_predict_final(fit, altered, "patient_id")
    expect_identical(prediction$status, "success")
    expect_equal(
      prediction$predictions$risk,
      altered_prediction$predictions$risk,
      tolerance = 0
    )
    expect_equal(
      prediction$predictions$survival_probability,
      altered_prediction$predictions$survival_probability,
      tolerance = 0
    )
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(
      prediction$predictions$survival_probability >= 0 &
        prediction$predictions$survival_probability <= 1
    ))
  }
})

test_that("ahaz BIC and Bujar execute inside outer analysis folds", {
  skip_if_not_installed("ahaz", minimum_version = "1.15.1")
  skip_if_not_installed("bujar", minimum_version = "0.2-11")
  data <- make_ahaz_bujar_fixture(n = 100L, p = 4L, seed = 609L)
  features <- ahaz_bujar_features(data)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    features,
    "month"
  )

  ahaz <- sb_run_cv(
    task,
    sb_pipeline(learner = "ahaz_lasso_bic"),
    folds = 3L,
    seed = 17L
  )
  expect_identical(ahaz$status, "success")
  expect_true(ahaz$complete_oof)
  expect_identical(
    ahaz$provenance$ahaz_version,
    as.character(utils::packageVersion("ahaz"))
  )
  expect_true(all(vapply(
    ahaz$coefficients,
    function(fold) is.finite(fold$lambda) && fold$lambda > 0,
    logical(1)
  )))

  bujar <- sb_run_cv(
    task,
    sb_pipeline(learner = "bujar_linear_boost"),
    folds = 3L,
    seed = 17L
  )
  expect_identical(bujar$status, "success")
  expect_true(bujar$complete_oof)
  expect_identical(
    bujar$provenance$bujar_version,
    as.character(utils::packageVersion("bujar"))
  )
})

test_that("ahaz tied-event CV retains original assessment endpoints", {
  skip_if_not_installed("ahaz", minimum_version = "1.15.1")
  data <- make_ahaz_bujar_fixture(n = 120L, p = 4L, seed = 609L, tied = TRUE)
  original <- data
  task <- sb_surv_task(data, "patient_id", "time", "event",
    ahaz_bujar_features(data), "month")
  for (learner in sb_ahaz_survival_learner_ids()) {
    result <- sb_run_cv(task, sb_pipeline(learner = learner), folds = 3L, seed = 17L)
    expect_identical(result$status, "success")
    expect_identical(result$complete_oof, TRUE)
    expect_identical(data, original)
    rows <- match(result$predictions$patient_id, data$patient_id)
    expect_equal(result$predictions$time, data$time[rows])
    expect_equal(result$predictions$event, data$event[rows])
  }
})

test_that("ahaz CV audit reconstructs calibrated horizon probabilities", {
  skip_if_not_installed("ahaz", minimum_version = "1.15.1")
  data <- make_ahaz_bujar_fixture(
    n = 160L,
    p = 6L,
    seed = 645L,
    signal = c(1.4, -1.1, 0.8)
  )
  features <- ahaz_bujar_features(data)
  horizon <- as.numeric(stats::median(data$time[data$event == 1L]))
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    features,
    "month",
    prediction_horizon = horizon
  )
  result <- sb_run_cv(
    task,
    sb_pipeline(learner = "ahaz_lasso_bic"),
    folds = 3L,
    seed = 17L
  )

  expect_identical(result$status, "success")
  expect_identical(result$complete_survival_probability_oof, TRUE)
  for (fold in seq_along(result$folds)) {
    rows <- result$predictions$fold == fold
    audit <- result$baseline_hazards[[fold]]
    coefficient <- audit$calibration_coefficient
    expect_length(coefficient, 1L)
    if (length(coefficient) != 1L) {
      next
    }
    reconstructed <- exp(
      -audit$cumulative_baseline_hazard *
        exp(coefficient * result$predictions$risk[rows])
    )
    expect_equal(
      result$predictions$survival_probability[rows],
      reconstructed,
      tolerance = 1e-14
    )
  }
})

for (warning_stage in c("fit", "prediction")) {
  test_that(paste("Bujar ordinary", warning_stage, "warnings preserve predictions"), {
    skip_if_not_installed("bujar")
    d <- make_ahaz_bujar_fixture(n = 40L, p = 3L)
    x <- d[ahaz_bujar_features(d)]
    adapter <- sb_bujar_linear_learner_adapter()
    reference <- adapter$fit(x, d$time, d$event)
    expected <- adapter$predict(reference, x)
    message <- paste("informational", warning_stage, "warning")
    if (warning_stage == "fit") {
      native_fit <- bujar::bujar
      local_mocked_bindings(bujar = function(...) {
        warning(message)
        native_fit(...)
      }, .package = "bujar")
    } else {
      native_predict <- stats::predict
      local_mocked_bindings(predict = function(...) {
        warning(message)
        native_predict(...)
      }, .package = "stats")
    }
    seen <- character()
    observed <- withCallingHandlers({
      fitted <- adapter$fit(x, d$time, d$event)
      adapter$predict(fitted, x)
    }, warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    expect_identical(observed, expected)
    expect_identical(seen, message)
  })
}

test_that("Bujar warnings do not bypass convergence or coefficient checks", {
  skip_if_not_installed("bujar")
  d <- make_ahaz_bujar_fixture(n = 40L, p = 3L)
  x <- d[ahaz_bujar_features(d)]
  adapter <- sb_bujar_linear_learner_adapter()
  native_fit <- bujar::bujar
  defect <- "contype"
  local_mocked_bindings(bujar = function(...) {
    engine <- native_fit(...)
    warning("informational fit warning")
    if (defect == "contype") engine$contype <- 1
    if (defect == "nonconv") engine$nonconv <- TRUE
    if (defect == "coefficient") engine$coef.bj[1] <- NA_real_
    engine
  }, .package = "bujar")
  for (defect in c("contype", "nonconv", "coefficient")) {
    error <- tryCatch(
      suppressWarnings(adapter$fit(x, d$time, d$event)),
      error = identity
    )
    expected <- if (defect == "coefficient") {
      "stormbreaker_invalid_bujar_fit"
    } else {
      "stormbreaker_bujar_nonconvergence"
    }
    expect_s3_class(error, expected)
  }
})
