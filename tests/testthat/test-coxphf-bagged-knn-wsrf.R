make_native_expansion_binary_fixture <- function(
  n = 90L, p = 5L, seed = 44101L
) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    eta <- 0.8 * x[, 1L] - 0.6 * x[, min(2L, p)]
    outcome <- stats::rbinom(n, 1L, stats::plogis(eta))
    data.frame(
      patient_id = sprintf("NB%03d", seq_len(n)),
      outcome = outcome,
      x,
      check.names = FALSE
    ) |>
      stats::setNames(c(
        "patient_id", "outcome", "gene one",
        if (p >= 2L) "2gene",
        if (p > 2L) paste0("gene-", seq_len(p - 2L))
      ))
  })
}

make_native_expansion_survival_fixture <- function(
  n = 100L, p = 5L, seed = 44201L
) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    eta <- 0.7 * x[, 1L] - 0.5 * x[, min(2L, p)]
    event_time <- stats::rexp(n, rate = exp(eta))
    censor_time <- stats::rexp(n, rate = 0.45)
    data.frame(
      patient_id = sprintf("NS%03d", seq_len(n)),
      time = pmin(event_time, censor_time) + 0.02,
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    ) |>
      stats::setNames(c(
        "patient_id", "time", "event", "gene one",
        if (p >= 2L) "2gene",
        if (p > 2L) paste0("gene-", seq_len(p - 2L))
      ))
  })
}

test_that("three native learners expose their declared selector routes", {
  ids <- c(
    sb_wsrf_binary_learner_ids(),
    sb_coxphf_learner_ids(),
    sb_bagged_knn_survival_learner_ids()
  )
  binary_selectors <- expected_fixed_learner_selectors("binary")
  survival_selectors <- expected_fixed_learner_selectors("survival")
  registry <- sb_pipeline_registry()
  rows <- registry[registry$learner %in% ids, , drop = FALSE]
  expect_identical(
    nrow(rows),
    length(binary_selectors) + 2L * length(survival_selectors)
  )
  expect_setequal(
    rows$selector[rows$learner == sb_wsrf_binary_learner_ids()],
    binary_selectors
  )
  for (learner in c(
    sb_coxphf_learner_ids(),
    sb_bagged_knn_survival_learner_ids()
  )) {
    expect_setequal(
      rows$selector[rows$learner == learner],
      survival_selectors
    )
  }
  expect_route_dependencies(rows, c(wsrf_binary = "wsrf",
    coxphf_firth = "coxphf", bagged_knn_survival = "bnnSurvival"))

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
      spec$feature_set_id <- "native_expansion_test"
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

test_that("Firth Cox matches native coefficients and portable x-beta replay", {
  skip_if_not_installed("coxphf", minimum_version = "1.13.4")
  data <- make_native_expansion_survival_fixture()
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[c(2L, 19L, 73L), , drop = FALSE]
  adapter <- sb_coxphf_learner_adapter("coxphf_firth")

  model <- adapter$fit(x, data$time, data$event)
  prepared <- sb_coxphf_formula_data(list(
    x = x,
    time = data$time,
    event = data$event
  ))
  native <- coxphf::coxphf(
    prepared$formula,
    prepared$data,
    pl = FALSE,
    maxit = 100L,
    maxhs = 5L,
    epsilon = 1e-6,
    gconv = 1e-4,
    maxstep = 0.5,
    firth = TRUE,
    penalty = 0.5
  )
  expect_identical(model$coefficient, unname(native$coefficients))
  risk <- adapter$predict(model, assessment)
  expect_identical(
    risk,
    unname(as.vector(as.matrix(assessment) %*% native$coefficients))
  )

  horizon <- as.numeric(stats::median(data$time))
  probability <- adapter$predict_survival(
    model, x, data$time, data$event, assessment, horizon
  )
  state <- adapter$extract_state(model, list(
    x = x,
    time = data$time,
    event = data$event,
    prediction_horizon = horizon,
    time_unit = "month"
  ))
  replay <- adapter$predict_state(
    unserialize(serialize(state, NULL)),
    assessment
  )
  expect_identical(replay$risk, risk)
  expect_identical(replay$survival_probability, probability$probability)
  expect_identical(probability$ties, "breslow")
  expect_false(any(c("engine", "time", "event", "outcome") %in% names(state)))
  expect_error(
    adapter$predict(model, assessment[rev(names(assessment))]),
    class = "stormbreaker_learner_schema_mismatch"
  )
})

test_that("Firth Cox handles monotone likelihood and exposes warnings", {
  skip_if_not_installed("coxphf", minimum_version = "1.13.4")
  n <- 30L
  x <- data.frame(marker = rev(seq_len(n)))
  adapter <- sb_coxphf_learner_adapter("coxphf_firth")
  model <- adapter$fit(x, seq_len(n) + 0.01, rep(1L, n))
  expect_true(all(is.finite(model$coefficient)))
  expect_true(model$iterations < adapter$params$max_iterations)
  expect_error(
    sb_coxphf_stop_on_warning(
      {
        warning("native convergence warning")
        TRUE
      },
      "fitting"
    ),
    class = "stormbreaker_coxphf_warning"
  )
})

test_that("bagged kNN survival matches native curves and retains required data", {
  skip_if_not_installed("bnnSurvival", minimum_version = "0.1.5")
  data <- make_native_expansion_survival_fixture(seed = 44301L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  rownames(x) <- paste0("PRIVATE_PATIENT_", seq_len(nrow(x)))
  assessment <- x[c(4L, 28L, 77L), , drop = FALSE]
  adapter <- sb_bagged_knn_survival_learner_adapter(
    "bagged_knn_survival"
  )
  old_cores <- getOption("mc.cores")
  set.seed(44302L)
  caller_rng <- .Random.seed
  model <- adapter$fit(x, data$time, data$event)
  expect_identical(.Random.seed, caller_rng)
  expect_identical(getOption("mc.cores"), old_cores)

  prepared <- sb_bagged_knn_formula_data(list(
    x = x,
    time = data$time,
    event = data$event
  ))
  native <- sb_with_bagged_knn_single_core(withr::with_seed(
    1L,
    bnnSurvival::bnnSurvival(
      prepared$formula,
      prepared$analysis,
      k = max(1L, floor(nrow(x) / 10L)),
      num_base_learners = 50L,
      num_features_per_base_learner = length(prepared$internal_features),
      metric = "mahalanobis",
      weighting_function = sb_bagged_knn_constant_weight,
      replace = TRUE,
      sample_fraction = 1
    )
  ))
  predict_generic <- get(
    "predict",
    envir = asNamespace("bnnSurvival"),
    inherits = FALSE
  )
  native_result <- sb_with_bagged_knn_single_core(
    predict_generic(native, sb_bagged_knn_internal_data(assessment)$data)
  )
  curves <- sb_predict_bagged_knn_survival_curves(model, assessment)
  expect_identical(
    curves$survival,
    as.matrix(bnnSurvival::predictions(native_result))
  )
  expect_true(all(apply(curves$survival, 1L, function(values) {
    all(diff(values) <= 1e-10)
  })))

  horizon <- as.numeric(stats::median(data$time))
  probability <- adapter$predict_survival(
    model, x, data$time, data$event, assessment, horizon
  )
  state <- adapter$extract_state(model, list(
    x = x,
    time = data$time,
    event = data$event,
    prediction_horizon = horizon,
    time_unit = "month"
  ))
  replay <- adapter$predict_state(
    unserialize(serialize(state, NULL)),
    assessment
  )
  expect_identical(replay$risk, curves$risk)
  expect_identical(replay$survival_probability, probability$probability)
  expect_identical(
    state$training_data_retention,
    "required_by_bagged_knn_prediction_state"
  )
  restored_engine <- unserialize(state$engine_raw)
  training_data <- methods::slot(restored_engine, "train_data")
  expect_gt(nrow(training_data), 0L)
  expect_false(any(startsWith(
    rownames(training_data),
    "PRIVATE_PATIENT_"
  )))
  expect_length(
    grepRaw("PRIVATE_PATIENT_", state$engine_raw, fixed = TRUE, all = TRUE),
    0L
  )
  expect_identical(
    environment(methods::slot(restored_engine, "weighting_function")),
    environment(sb_bagged_knn_constant_weight)
  )

  one_feature <- x[1L]
  one_model <- adapter$fit(one_feature, data$time, data$event)
  expect_true(one_model$uses_zero_variance_bridge)
  expect_true(all(is.finite(adapter$predict(
    one_model,
    one_feature[1:2, , drop = FALSE]
  ))))
})

test_that("weighted-subspace forest matches native votes and compact replay", {
  skip_if_not_installed("wsrf", minimum_version = "1.7.32")
  data <- make_native_expansion_binary_fixture(seed = 44401L)
  features <- setdiff(names(data), c("patient_id", "outcome"))
  x <- data[features]
  assessment <- x[c(3L, 31L, 82L), , drop = FALSE]
  adapter <- sb_wsrf_binary_learner_adapter("wsrf_binary")

  set.seed(44402L)
  caller_rng <- .Random.seed
  model <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, caller_rng)
  prepared <- sb_wsrf_internal_data(x)
  native <- withr::with_seed(
    1L,
    wsrf::wsrf(
      prepared$data,
      factor(
        data$outcome,
        levels = c(0L, 1L),
        labels = c("control", "case")
      ),
      mtry = floor(log2(ncol(x)) + 1),
      ntree = 500L,
      weights = TRUE,
      parallel = FALSE,
      importance = FALSE,
      nodesize = 2L
    )
  )
  native_probability <- stats::predict(
    native,
    sb_wsrf_internal_data(assessment)$data,
    type = "prob"
  )$prob[, "case"]
  probability <- adapter$predict(model, assessment)
  expect_identical(probability, unname(native_probability))
  expect_null(model$engine$targetData)
  expect_true(all(lengths(model$engine$OOBSets) == 0L))
  expect_true(all(lengths(model$engine$OOBPredictLabels) == 0L))

  state <- adapter$extract_state(model)
  replay <- adapter$predict_state(
    unserialize(serialize(state, NULL)),
    assessment
  )
  expect_identical(replay$probability, probability)
  expect_false(state$training_outcome_retained)
  expect_false(state$oob_patient_indices_retained)
  expect_error(
    adapter$predict(model, assessment[rev(names(assessment))]),
    class = "stormbreaker_learner_schema_mismatch"
  )

  one_feature <- x[1L]
  one_model <- adapter$fit(one_feature, data$outcome)
  expect_true(one_model$uses_zero_variance_bridge)
  expect_true(all(adapter$predict(
    one_model,
    one_feature[1:2, , drop = FALSE]
  ) >= 0))
})

test_that("new native dependency failures remain explicit", {
  local_mocked_bindings(sb_coxphf_available = function() FALSE)
  expect_error(
    sb_require_coxphf_dependency(),
    class = "stormbreaker_missing_coxphf"
  )
  local_mocked_bindings(sb_bagged_knn_survival_available = function() FALSE)
  expect_error(
    sb_require_bagged_knn_survival_dependency(),
    class = "stormbreaker_missing_bnn_survival"
  )
  local_mocked_bindings(sb_wsrf_available = function() FALSE)
  expect_error(sb_require_wsrf_dependency(), class = "stormbreaker_missing_wsrf")
})

test_that("three native learners execute CV and outcome-blind final prediction", {
  skip_if_not_installed("coxphf", minimum_version = "1.13.4")
  skip_if_not_installed("bnnSurvival", minimum_version = "0.1.5")
  skip_if_not_installed("wsrf", minimum_version = "1.7.32")

  binary <- make_native_expansion_binary_fixture(n = 84L, seed = 44501L)
  binary_features <- setdiff(names(binary), c("patient_id", "outcome"))
  binary_task <- sb_binary_task(
    binary,
    "patient_id",
    "outcome",
    binary_features
  )
  binary_pipeline <- sb_pipeline(learner = "wsrf_binary")
  binary_cv <- sb_run_cv(binary_task, binary_pipeline, folds = 3L, seed = 1L)
  binary_fit <- sb_fit_final(binary_task, binary_pipeline)
  binary_prediction <- sb_predict_final(
    unserialize(serialize(binary_fit, NULL)),
    binary[1:5, c("patient_id", rev(binary_features)), drop = FALSE],
    "patient_id"
  )
  expect_identical(binary_cv$status, "success")
  expect_true(binary_cv$complete_oof)
  expect_identical(
    binary_cv$provenance$wsrf_version,
    as.character(utils::packageVersion("wsrf"))
  )
  expect_identical(binary_fit$status, "success")
  expect_identical(binary_prediction$status, "success")

  survival <- make_native_expansion_survival_fixture(
    n = 90L,
    seed = 44502L
  )
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
  for (learner in c(
    "coxphf_firth",
    "bagged_knn_survival"
  )) {
    pipeline <- sb_pipeline(learner = learner)
    cv <- sb_run_cv(survival_task, pipeline, folds = 3L, seed = 2L)
    fit <- sb_fit_final(survival_task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      survival[1:5, c("patient_id", rev(survival_features)), drop = FALSE],
      "patient_id"
    )
    expect_identical(cv$status, "success")
    expect_true(cv$complete_oof)
    expect_true(cv$complete_survival_probability_oof)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(
      prediction$predictions$survival_probability >= 0 &
        prediction$predictions$survival_probability <= 1
    ))
  }
})

test_that("bagged kNN ordinary warnings preserve native fit and predictions", {
  skip_if_not_installed("bnnSurvival")
  d <- make_native_expansion_survival_fixture(n = 40L, p = 3L)
  x <- d[setdiff(names(d), c("patient_id", "time", "event"))]
  id <- "bagged_knn_survival"
  params <- sb_bagged_knn_survival_params(id)
  params$num_base_learners <- 5L
  adapter <- sb_bagged_knn_survival_learner_adapter(id, params)
  reference <- adapter$fit(x, d$time, d$event)
  expected <- adapter$predict(reference, x)
  native_fit <- bnnSurvival::bnnSurvival
  native_predict <- sb_predict_bagged_knn_engine
  local_mocked_bindings(bnnSurvival = function(...) {
    warning("informational fit warning")
    native_fit(...)
  }, .package = "bnnSurvival")
  local_mocked_bindings(sb_predict_bagged_knn_engine = function(...) { warning("informational prediction warning"); native_predict(...) })
  seen <- character()
  observed <- withCallingHandlers({
    fitted <- adapter$fit(x, d$time, d$event)
    adapter$predict(fitted, x)
  }, warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(observed, expected)
  expect_identical(all(c("informational fit warning", "informational prediction warning") %in% seen), TRUE)
})

test_that("WSRF ordinary warnings preserve native fit and predictions", {
  skip_if_not_installed("wsrf")
  d <- make_native_expansion_binary_fixture(n = 40L, p = 3L)
  x <- d[setdiff(names(d), c("patient_id", "outcome"))]
  id <- "wsrf_binary"
  params <- sb_wsrf_binary_params(id)
  params$num_trees <- 5L
  adapter <- sb_wsrf_binary_learner_adapter(id, params)
  reference <- adapter$fit(x, d$outcome)
  expected <- adapter$predict(reference, x)
  native_fit <- wsrf::wsrf
  native_predict <- stats::predict
  local_mocked_bindings(wsrf = function(...) {
    warning("informational fit warning")
    native_fit(...)
  }, .package = "wsrf")
  local_mocked_bindings(predict = function(...) { warning("informational prediction warning"); native_predict(...) }, .package = "stats")
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
