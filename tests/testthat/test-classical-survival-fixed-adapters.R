test_that("classical survival pipeline settings are fixed and explicit", {
  coxboost <- sb_pipeline(learner = "coxboost")
  plsrcox <- sb_pipeline(learner = "plsrcox")
  superpc <- sb_pipeline(learner = "superpc")
  survivalsvm <- sb_pipeline(learner = "survivalsvm")

  expect_identical(
    coxboost$id,
    "none__coxboost__steps100__penaltyevent9"
  )
  expect_identical(plsrcox$id, "none__plsrcox__components2")
  expect_identical(
    superpc$id,
    "none__superpc__threshold1__components1"
  )
  expect_identical(survivalsvm$id, "none__survivalsvm__gamma1")
  expect_identical(coxboost$risk_formula, "x_beta")
  expect_identical(plsrcox$risk_reference, "zero")
  expect_identical(superpc$risk_direction, "higher_is_higher_hazard")
  expect_identical(survivalsvm$learner_params$kernel, "lin_kernel")

  expect_error(
    sb_pipeline(learner = "coxboost", coxboost_penalty = 0),
    class = "stormbreaker_bad_coxboost_penalty"
  )
  expect_error(
    sb_pipeline(learner = "superpc", superpc_threshold = -1),
    class = "stormbreaker_bad_superpc_threshold"
  )
  expect_error(
    sb_pipeline(learner = "superpc", superpc_components = 4L),
    class = "stormbreaker_bad_superpc_components"
  )
  expect_error(
    sb_pipeline(learner = "survivalsvm", survivalsvm_gamma_mu = 0),
    class = "stormbreaker_bad_survivalsvm_gamma_mu"
  )
  expect_error(
    sb_pipeline(learner = "coxph", plsrcox_components = 1L),
    class = "stormbreaker_irrelevant_classical_survival_parameter"
  )
})

test_that("classical survival adapters preserve their native linear scores", {
  skip_if_not_installed("CoxBoost", minimum_version = "1.5.1")
  skip_if_not_installed("plsRcox", minimum_version = "1.8.0")
  skip_if_not_installed("superpc", minimum_version = "1.12")
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  d <- make_sb_fixture(n = 72L, seed = 7811L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(3L, 18L, 41L), , drop = FALSE]

  adapters <- list(
    coxboost = sb_coxboost_learner_adapter(stepno = 10L),
    plsrcox = sb_plsrcox_learner_adapter(components = 2L),
    superpc = sb_superpc_learner_adapter(
      threshold = 0,
      components = 1L
    ),
    survivalsvm = sb_survivalsvm_learner_adapter(gamma_mu = 1)
  )
  models <- lapply(adapters, function(adapter) {
    adapter$fit(x, d$time, d$event)
  })

  coxboost_native <- as.numeric(stats::predict(
    models$coxboost$engine,
    newdata = as.matrix(new_x),
    at.step = 10L,
    type = "lp"
  ))
  plsrcox_native <- as.numeric(stats::predict(
    models$plsrcox$engine,
    newdata = new_x,
    comps = 2L,
    type = "lp",
    verbose = FALSE
  ))
  superpc_native <- as.numeric(
    models$superpc$engine$prediction$v.pred.1df
  )
  svm_center <- models$survivalsvm$audit$center
  svm_scale <- models$survivalsvm$audit$scale
  scaled_new_x <- as.data.frame(
    sweep(sweep(as.matrix(new_x), 2L, svm_center, "-"),
          2L, svm_scale, "/"),
    optional = TRUE
  )
  names(scaled_new_x) <- names(new_x)
  survivalsvm_native <- -as.numeric(stats::predict(
    models$survivalsvm$engine,
    newdata = scaled_new_x
  )$predicted) * models$survivalsvm$audit$calibration_coefficient

  centered <- function(value) value - value[[1L]]
  expect_equal(
    centered(adapters$coxboost$predict(
      models$coxboost,
      new_x
    )),
    centered(coxboost_native),
    tolerance = 1e-10
  )
  expect_equal(
    centered(adapters$plsrcox$predict(
      models$plsrcox,
      new_x
    )),
    centered(plsrcox_native),
    tolerance = 1e-10
  )
  expect_equal(
    centered(adapters$superpc$predict(models$superpc, x)),
    centered(superpc_native),
    tolerance = 1e-10
  )
  expect_equal(
    centered(adapters$survivalsvm$predict(
      models$survivalsvm,
      new_x
    )),
    centered(survivalsvm_native),
    tolerance = 1e-8
  )

  for (learner in names(adapters)) {
    adapter <- adapters[[learner]]
    model <- models[[learner]]
    risk <- adapter$predict(model, new_x)
    probability <- adapter$predict_survival(
      model,
      x,
      d$time,
      d$event,
      new_x,
      prediction_horizon = 7
    )
    state <- adapter$extract_state(model, list(
      x = x,
      time = d$time,
      event = d$event,
      prediction_horizon = 7,
      time_unit = "month"
    ))
    restored <- unserialize(serialize(state, NULL))
    state_prediction <- adapter$predict_state(restored, new_x)

    expect_identical(state$learner, learner)
    expect_false("engine" %in% names(state))
    expect_identical(state$audit, model$audit)
    expect_equal(state_prediction$risk, risk, tolerance = 0)
    expect_equal(
      state_prediction$survival_probability,
      probability$probability,
      tolerance = 0
    )
    expect_true(all(probability$probability >= 0 &
      probability$probability <= 1))
    expect_identical(probability$ties, "breslow")
  }
})

test_that("classical survival learners run CV and portable final prediction", {
  skip_if_not_installed("CoxBoost", minimum_version = "1.5.1")
  skip_if_not_installed("plsRcox", minimum_version = "1.8.0")
  skip_if_not_installed("superpc", minimum_version = "1.12")
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  d <- make_sb_fixture(n = 72L, seed = 7821L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )
  pipelines <- list(
    sb_pipeline(learner = "coxboost", coxboost_stepno = 10L),
    sb_pipeline(learner = "plsrcox", plsrcox_components = 1L),
    sb_pipeline(
      learner = "superpc",
      superpc_threshold = 0,
      superpc_components = 1L
    ),
    sb_pipeline(learner = "survivalsvm", survivalsvm_gamma_mu = 1)
  )
  version_fields <- c(
    "coxboost_version", "plsrcox_version", "superpc_version",
    "survivalsvm_version"
  )

  for (index in seq_along(pipelines)) {
    pipeline <- pipelines[[index]]
    result <- sb_run_cv(task, pipeline, folds = 3L, seed = 7822L)
    fit <- sb_fit_final(task, pipeline)
    new_data <- d[c(5L, 29L, 61L), c("patient_id", predictors)]
    prediction <- sb_predict_final(fit, new_data, "patient_id")
    direct <- sb_predict_cox_learner_state(
      fit$learner_state,
      new_data[fit$required_features]
    )
    explanation <- sb_additive_shap(
      fit,
      new_data,
      "patient_id",
      d[1:24, fit$required_features, drop = FALSE]
    )

    expect_identical(result$status, "success")
    expect_true(result$complete_oof)
    expect_true(result$complete_survival_probability_oof)
    expect_length(result$coefficients, 3L)
    expect_true(all(vapply(
      result$coefficients,
      function(value) identical(value$status, "success"),
      logical(1)
    )))
    expect_true(version_fields[[index]] %in% names(result$provenance))
    expect_identical(fit$status, "success")
    expect_true(version_fields[[index]] %in% names(fit$provenance))
    expect_identical(prediction$status, "success")
    expect_equal(prediction$predictions$risk, direct$risk, tolerance = 0)
    expect_equal(
      prediction$predictions$survival_probability,
      direct$survival_probability,
      tolerance = 0
    )
    expect_identical(explanation$learner, pipeline$learner)
    expect_equal(
      explanation$scores$score,
      prediction$predictions$risk,
      tolerance = 1e-14
    )
  }
})

test_that("classical survival fit warnings are explicit failures", {
  skip_if_not_installed("CoxBoost", minimum_version = "1.5.1")
  skip_if_not_installed("plsRcox", minimum_version = "1.8.0")
  skip_if_not_installed("superpc", minimum_version = "1.12")
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  d <- make_sb_fixture(n = 36L, seed = 7829L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  injected_warning <- function(...) {
    warning("INJECTED_CLASSICAL_SURVIVAL_WARNING")
  }

  coxboost_error <- capture_sb_error(testthat::with_mocked_bindings(
    CoxBoost = injected_warning,
    sb_coxboost_learner_adapter(5L)$fit(x, d$time, d$event),
    .package = "CoxBoost"
  ))
  expect_s3_class(coxboost_error, "stormbreaker_coxboost_fit_warning")
  expect_match(
    conditionMessage(coxboost_error),
    "INJECTED_CLASSICAL_SURVIVAL_WARNING",
    fixed = TRUE
  )
  expect_s3_class(coxboost_error$parent, "warning")

  plsrcox_error <- capture_sb_error(testthat::with_mocked_bindings(
    plsRcox = injected_warning,
    sb_plsrcox_learner_adapter(1L)$fit(x, d$time, d$event),
    .package = "plsRcox"
  ))
  expect_s3_class(plsrcox_error, "stormbreaker_plsrcox_fit_warning")
  expect_match(
    conditionMessage(plsrcox_error),
    "INJECTED_CLASSICAL_SURVIVAL_WARNING",
    fixed = TRUE
  )
  expect_s3_class(plsrcox_error$parent, "warning")

  superpc_train_error <- capture_sb_error(testthat::with_mocked_bindings(
    `superpc.train` = injected_warning,
    sb_superpc_learner_adapter(0, 1L)$fit(x, d$time, d$event),
    .package = "superpc"
  ))
  expect_s3_class(
    superpc_train_error,
    "stormbreaker_superpc_train_warning"
  )
  expect_s3_class(superpc_train_error$parent, "warning")

  superpc_predict_error <- capture_sb_error(testthat::with_mocked_bindings(
    `superpc.predict` = injected_warning,
    sb_superpc_learner_adapter(0, 1L)$fit(x, d$time, d$event),
    .package = "superpc"
  ))
  expect_s3_class(
    superpc_predict_error,
    "stormbreaker_superpc_predict_warning"
  )
  expect_s3_class(superpc_predict_error$parent, "warning")

  survivalsvm_error <- capture_sb_error(testthat::with_mocked_bindings(
    survivalsvm = injected_warning,
    sb_survivalsvm_learner_adapter(1)$fit(x, d$time, d$event),
    .package = "survivalsvm"
  ))
  expect_s3_class(
    survivalsvm_error,
    "stormbreaker_survivalsvm_fit_warning"
  )
  expect_match(
    conditionMessage(survivalsvm_error),
    "INJECTED_CLASSICAL_SURVIVAL_WARNING",
    fixed = TRUE
  )
  expect_s3_class(survivalsvm_error$parent, "warning")

  calibration_error <- capture_sb_error(sb_fit_survivalsvm_calibration(
    time = 1:12,
    event = rep(1L, 12L),
    raw_risk = -(1:12),
    ties = "breslow"
  ))
  expect_s3_class(
    calibration_error,
    "stormbreaker_survivalsvm_calibration_warning"
  )
  expect_match(
    conditionMessage(calibration_error),
    "did not converge",
    fixed = TRUE
  )
  expect_s3_class(calibration_error$parent, "warning")
})

test_that("single-feature SuperPC preserves the native one-component score", {
  skip_if_not_installed("superpc")
  d <- sb_example_survival(n = 200L, seed = 20260905L)
  noise <- withr::with_seed(20260906L, matrix(rnorm(200L * 17L), nrow = 200L))
  for (j in seq_len(ncol(noise))) d[[paste0("noise_", j)]] <- noise[, j]
  d <- d[1:160, ]
  x <- d[setdiff(names(d), c("patient_id", "time", "event"))]
  adapter <- sb_superpc_learner_adapter(1, 1L)
  model <- adapter$fit(x, d$time, d$event)

  native_data <- list(x = t(as.matrix(x)), y = d$time, censoring.status = d$event)
  trained <- superpc::superpc.train(native_data, type = "survival", s0.perc = 0.5)
  selected <- abs(trained$feature.scores) >= 1
  expect_equal(sum(selected), 1L)
  # Oracle only: duplicate the same signal so native PCA retains matrix shape.
  # Its normalized rank-one score is unchanged; production adds no predictors.
  native_data$x <- rbind(native_data$x, native_data$x[selected, , drop = FALSE])
  trained$feature.scores <- c(trained$feature.scores, trained$feature.scores[selected])
  native <- superpc::superpc.predict(
    trained, native_data, native_data, threshold = 1,
    n.components = 1L, prediction.type = "continuous"
  )
  risk <- adapter$predict(model, x)
  expected <- as.numeric(native$v.pred.1df)
  expect_equal(risk - risk[1L], expected - expected[1L], tolerance = 1e-10)
  expect_identical(model$audit$selected_features, "gene_a")
  expect_identical(model$params$threshold, 1)
  expect_identical(model$params$components, 1L)
})

test_that("fixed SuperPC reports its native one-input training boundary", {
  skip_if_not_installed("superpc", minimum_version = "1.12")
  d <- make_sb_fixture(n = 72L, seed = 7830L)

  expect_error(
    sb_superpc_learner_adapter(1, 1L)$fit(
      d["gene_a"], d$time, d$event
    ),
    class = "stormbreaker_superpc_too_few_features"
  )
})

test_that("SuperPC native prediction errors explain the selected feature count", {
  skip_if_not_installed("superpc")
  d <- make_sb_fixture(n = 36L, seed = 7829L)
  local_mocked_bindings(
    `superpc.train` = function(...) {
      list(feature.scores = c(1.5, 0.5, 0.2))
    },
    `superpc.predict` = function(...) stop("INJECTED_ARRAY_ERROR"),
    .package = "superpc"
  )

  error <- capture_sb_error(
    sb_superpc_learner_adapter(1, 2L)$fit(
      d[c("gene_a", "gene_b", "gene_c")], d$time, d$event
    )
  )
  expect_s3_class(error, "stormbreaker_superpc_predict_error")
  expect_match(conditionMessage(error), "selecting 1 of 3 predictors", fixed = TRUE)
  expect_match(conditionMessage(error), "threshold=1, components=2", fixed = TRUE)
  expect_s3_class(error$parent, "simpleError")
  expect_identical(conditionMessage(error$parent), "INJECTED_ARRAY_ERROR")
})

test_that("missing classical survival dependencies fail without fallback", {
  d <- make_sb_fixture(n = 30L, seed = 7831L)
  x <- d[c("gene_a", "gene_b", "gene_c")]

  local_mocked_bindings(sb_coxboost_available = function() FALSE)
  expect_error(
    sb_coxboost_learner_adapter(5L)$fit(x, d$time, d$event),
    class = "stormbreaker_missing_coxboost"
  )

  local_mocked_bindings(sb_plsrcox_available = function() FALSE)
  expect_error(
    sb_plsrcox_learner_adapter(1L)$fit(x, d$time, d$event),
    class = "stormbreaker_missing_plsrcox"
  )

  local_mocked_bindings(sb_superpc_available = function() FALSE)
  expect_error(
    sb_superpc_learner_adapter(0, 1L)$fit(x, d$time, d$event),
    class = "stormbreaker_missing_superpc"
  )

  local_mocked_bindings(sb_survivalsvm_available = function() FALSE)
  expect_error(
    sb_survivalsvm_learner_adapter(1)$fit(x, d$time, d$event),
    class = "stormbreaker_missing_survivalsvm"
  )
})
