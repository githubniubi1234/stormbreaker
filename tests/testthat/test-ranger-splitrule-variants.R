make_ranger_splitrule_binary_fixture <- function(n = 84L, seed = 8841L) {
  withr::with_seed(seed, {
    outcome <- sample(rep(c(0L, 1L), length.out = n))
    data.frame(
      patient_id = sprintf("RSB%03d", seq_len(n)),
      outcome = outcome,
      "gene one" = stats::rnorm(n, 1.1 * outcome),
      "2gene" = stats::rnorm(n, -0.8 * outcome),
      "gene-three" = stats::rnorm(n),
      check.names = FALSE
    )
  })
}

make_ranger_splitrule_survival_fixture <- function(n = 84L, seed = 8842L) {
  withr::with_seed(seed, {
    x1 <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    x3 <- stats::rnorm(n)
    event_time <- stats::rexp(
      n,
      rate = 0.08 * exp(x1 - 0.7 * x2)
    ) + 0.1
    censor_time <- stats::rexp(n, rate = 0.05) + 0.1
    data.frame(
      patient_id = sprintf("RSS%03d", seq_len(n)),
      time = pmin(event_time, censor_time),
      event = as.integer(event_time <= censor_time),
      "gene one" = x1,
      "2gene" = x2,
      "gene-three" = x3,
      check.names = FALSE
    )
  })
}

test_that("ranger split-rule learners declare five distinct specifications", {
  binary <- c(
    "ranger_binary",
    "ranger_extratrees_binary",
    "ranger_hellinger_binary"
  )
  survival <- c(
    "ranger_survival",
    "ranger_extratrees_survival",
    "ranger_cindex_survival",
    "ranger_maxstat_survival"
  )

  expect_identical(sb_ranger_binary_learner_ids(), binary)
  expect_identical(sb_ranger_survival_learner_ids(), survival)
  expect_identical(sb_ranger_learner_ids(), c(binary, survival))

  binary_splitrules <- c("gini", "extratrees", "hellinger")
  for (i in seq_along(binary)) {
    params <- sb_ranger_binary_params(learner = binary[[i]])
    expect_identical(params$splitrule, binary_splitrules[[i]])
    expect_identical(params$replace, !identical(binary_splitrules[[i]], "extratrees"))
    expect_identical(params$probability, TRUE)
  }
  expect_identical(
    sb_ranger_binary_params(
      learner = "ranger_extratrees_binary"
    )$num_random_splits,
    1L
  )

  survival_splitrules <- c("logrank", "extratrees", "C", "maxstat")
  for (i in seq_along(survival)) {
    params <- sb_ranger_survival_params(learner = survival[[i]])
    expect_identical(params$splitrule, survival_splitrules[[i]])
    expect_identical(
      params$replace,
      !identical(survival_splitrules[[i]], "extratrees")
    )
  }
  extra <- sb_ranger_survival_params(
    learner = "ranger_extratrees_survival"
  )
  maxstat <- sb_ranger_survival_params(
    learner = "ranger_maxstat_survival"
  )
  expect_identical(extra$num_random_splits, 1L)
  expect_identical(maxstat$alpha, 0.5)
  expect_identical(maxstat$minprop, 0.1)

  binary_selectors <- expected_fixed_learner_selectors("binary")
  survival_selectors <- expected_fixed_learner_selectors("survival")
  registry <- sb_pipeline_registry()
  rows <- registry[
    registry$learner %in% c(binary[-1L], survival[-1L]), ,
    drop = FALSE
  ]
  expect_identical(
    nrow(rows),
    length(binary[-1L]) * length(binary_selectors) +
      length(survival[-1L]) * length(survival_selectors)
  )
  expect_identical(
    sort(unique(rows$selector[rows$task_type == "binary"])),
    sort(binary_selectors)
  )
  expect_identical(
    sort(unique(rows$selector[rows$task_type == "survival"])),
    sort(survival_selectors)
  )
  expect_route_dependencies(rows, stats::setNames(
    rep("ranger", length(c(binary[-1L], survival[-1L]))), c(binary[-1L], survival[-1L])))
})

test_that("new ranger split rules support one non-syntactic feature", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  binary_data <- make_ranger_splitrule_binary_fixture(n = 48L, seed = 8849L)
  binary_x <- binary_data["gene one"]
  for (learner in sb_ranger_binary_learner_ids()[-1L]) {
    adapter <- sb_ranger_binary_learner_adapter(
      num_trees = 12L,
      learner = learner
    )
    model <- adapter$fit(binary_x, binary_data$outcome)
    prediction <- adapter$predict(model, binary_x[1:4, , drop = FALSE])

    expect_length(prediction, 4L)
    expect_true(all(is.finite(prediction)))
  }

  survival_data <- make_ranger_splitrule_survival_fixture(
    n = 48L,
    seed = 8850L
  )
  survival_x <- survival_data["gene one"]
  for (learner in sb_ranger_survival_learner_ids()[-1L]) {
    adapter <- sb_ranger_survival_learner_adapter(
      num_trees = 12L,
      learner = learner
    )
    model <- adapter$fit(
      survival_x,
      survival_data$time,
      survival_data$event
    )
    prediction <- adapter$predict(
      model,
      survival_x[1:4, , drop = FALSE]
    )

    expect_length(prediction, 4L)
    expect_true(all(is.finite(prediction)))
  }
})

test_that("binary ranger split rules replay native probabilities", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  data <- make_ranger_splitrule_binary_fixture()
  x <- data[c("gene one", "2gene", "gene-three")]
  new_x <- x[c(3L, 27L, 68L), , drop = FALSE]
  withr::local_seed(8843L)

  for (learner in sb_ranger_binary_learner_ids()[-1L]) {
    adapter <- sb_ranger_binary_learner_adapter(
      num_trees = 24L,
      mtry = 1L,
      min_node_size = 10L,
      seed = 1L,
      learner = learner
    )
    rng_before <- .Random.seed
    model <- adapter$fit(x, data$outcome)
    expect_identical(.Random.seed, rng_before)
    repeated <- adapter$fit(x, data$outcome)
    expect_identical(
      serialize(model$forest, NULL),
      serialize(repeated$forest, NULL)
    )

    native <- stats::predict(
      model$forest,
      data = new_x,
      num.trees = 24L,
      type = "response",
      seed = 0L,
      num.threads = 1L,
      verbose = FALSE
    )$predictions[, "1"]
    observed <- adapter$predict(model, new_x)
    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    replay <- adapter$predict_state(state, new_x)

    expect_identical(model$learner, learner)
    expect_equal(observed, unname(native), tolerance = 0)
    expect_equal(replay$probability, observed, tolerance = 0)
    expect_identical(state$learner, learner)
    expect_identical(state$params$splitrule, adapter$params$splitrule)
    expect_type(state$model_raw, "raw")
    expect_false("forest" %in% names(state))
  }
})

test_that("survival ranger split rules replay native mortality and survival", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  data <- make_ranger_splitrule_survival_fixture()
  x <- data[c("gene one", "2gene", "gene-three")]
  new_x <- x[c(3L, 27L, 68L), , drop = FALSE]
  horizon <- as.numeric(stats::quantile(data$time, 0.4))
  withr::local_seed(8844L)

  for (learner in sb_ranger_survival_learner_ids()[-1L]) {
    adapter <- sb_ranger_survival_learner_adapter(
      num_trees = 24L,
      mtry = 1L,
      min_node_size = 3L,
      seed = 1L,
      learner = learner
    )
    rng_before <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, rng_before)
    repeated <- adapter$fit(x, data$time, data$event)
    expect_identical(
      serialize(model$forest, NULL),
      serialize(repeated$forest, NULL)
    )

    native <- stats::predict(
      model$forest,
      data = new_x,
      num.trees = 24L,
      type = "response",
      seed = 0L,
      num.threads = 1L,
      verbose = FALSE
    )
    risk <- adapter$predict(model, new_x)
    probability <- adapter$predict_survival(
      model,
      x,
      data$time,
      data$event,
      new_x,
      horizon
    )
    eligible <- which(native$unique.death.times <= horizon)
    native_probability <- if (length(eligible) == 0L) {
      rep(1, nrow(new_x))
    } else {
      unname(native$survival[, max(eligible)])
    }
    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = horizon,
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, new_x)

    expect_identical(model$learner, learner)
    expect_equal(risk, unname(rowSums(native$chf)), tolerance = 0)
    expect_equal(probability$probability, native_probability, tolerance = 0)
    expect_equal(replay$risk, risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      probability$probability,
      tolerance = 0
    )
    expect_identical(state$learner, learner)
    expect_identical(state$params$splitrule, adapter$params$splitrule)
    expect_null(state$baseline_hazard)
    expect_type(state$model_raw, "raw")
  }
})

test_that("new ranger split rules execute through CV and final prediction", {
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  binary_data <- make_ranger_splitrule_binary_fixture(n = 72L, seed = 8845L)
  binary_features <- c("gene one", "2gene", "gene-three")
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    binary_features
  )
  for (learner in sb_ranger_binary_learner_ids()[-1L]) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 2L,
      ranger_num_trees = 18L
    )
    comparison <- sb_run_cv(binary_task, pipeline, folds = 3L, seed = 8846L)
    fit <- sb_fit_final(binary_task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      binary_data[1:5, c("patient_id", rev(binary_features))],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_identical(fit$learner_state$learner, learner)
  }

  survival_data <- make_ranger_splitrule_survival_fixture(
    n = 72L,
    seed = 8847L
  )
  survival_features <- c("gene one", "2gene", "gene-three")
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    survival_features,
    "month",
    prediction_horizon = as.numeric(stats::quantile(survival_data$time, 0.4))
  )
  for (learner in sb_ranger_survival_learner_ids()[-1L]) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 2L,
      ranger_num_trees = 18L
    )
    comparison <- sb_run_cv(
      survival_task,
      pipeline,
      folds = 3L,
      seed = 8848L
    )
    fit <- sb_fit_final(survival_task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      survival_data[1:5, c("patient_id", rev(survival_features))],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(
      comparison$complete_survival_probability_oof,
      TRUE
    )
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_identical(fit$learner_state$learner, learner)
  }
})
