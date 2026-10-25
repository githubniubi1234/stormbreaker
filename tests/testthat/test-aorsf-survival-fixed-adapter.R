test_that("aorsf survival declares one fixed native forest", {
  expect_true(
    "aorsf_survival" %in% sb_additional_survival_learner_ids()
  )
  expect_identical(
    sb_additional_survival_dependency("aorsf_survival"),
    list(package = "aorsf", packages = "aorsf",
         available = requireNamespace("aorsf", quietly = TRUE))
  )
  expect_identical(
    sb_additional_survival_params("aorsf_survival"),
    list(
      num_trees = 500L,
      num_split = 5L,
      num_retry = 3L,
      mtry_rule = "ceiling_sqrt_p",
      sample_with_replacement = TRUE,
      sample_fraction = 0.632,
      control_method = "glm",
      control_scale_x = TRUE,
      control_ties = "efron",
      control_max_iter = 1L,
      seed = 1L,
      threads = 1L,
      importance = "none",
      attach_data = FALSE
    )
  )
})

test_that("aorsf survival risk and horizon survival equal native predictions", {
  skip_if_not_installed("aorsf")
  d <- make_sb_fixture(n = 60L, seed = 8431L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:7, , drop = FALSE]
  params <- sb_additional_survival_params("aorsf_survival")
  params$num_trees <- 40L
  adapter <- sb_additional_survival_learner_adapter(
    "aorsf_survival",
    params
  )
  model <- adapter$fit(x, d$time, d$event)
  risk <- adapter$predict(model, assessment)
  survival <- adapter$predict_survival(
    model,
    x,
    d$time,
    d$event,
    assessment,
    prediction_horizon = 7
  )
  native_risk <- stats::predict(
    model$engine,
    new_data = assessment,
    pred_type = "mort",
    n_thread = 1L
  )
  native_survival <- stats::predict(
    model$engine,
    new_data = assessment,
    pred_type = "surv",
    pred_horizon = 7,
    n_thread = 1L
  )
  state <- adapter$extract_state(model, list(
    x = x,
    time = d$time,
    event = d$event,
    prediction_horizon = 7,
    time_unit = "month"
  ))
  replay <- adapter$predict_state(
    unserialize(serialize(state, NULL)),
    assessment
  )

  expect_s3_class(model$engine, "ObliqueForestSurvival")
  expect_null(model$engine$data)
  expect_equal(risk, unname(as.numeric(native_risk)), tolerance = 1e-14)
  expect_equal(
    survival$probability,
    unname(as.numeric(native_survival)),
    tolerance = 1e-14
  )
  expect_equal(replay$risk, risk, tolerance = 0)
  expect_equal(
    replay$survival_probability,
    survival$probability,
    tolerance = 0
  )
  expect_identical(state$risk_reference, "aorsf_native")
  expect_identical(state$risk_formula, "native_mortality_over_event_times")
  expect_identical(
    state$survival_probability$source,
    "aorsf_oblique_survival_forest_direct"
  )
  expect_false("baseline_hazard" %in% names(Filter(Negate(is.null), state)))
})

test_that("aorsf keeps raw feature names while using unique native aliases", {
  skip_if_not_installed("aorsf")
  d <- make_sb_fixture(n = 60L, seed = 8431L)
  x <- d[c("gene_c", "gene_a", "gene_b")]
  names(x) <- c("NKX3-1", "NKX3.1", "if")
  native_x <- x
  names(native_x) <- make.names(names(x), unique = TRUE)
  params <- sb_additional_survival_params("aorsf_survival")
  params$num_trees <- 40L
  adapter <- sb_additional_survival_learner_adapter("aorsf_survival", params)
  model <- adapter$fit(x, d$time, d$event)
  renamed <- adapter$fit(native_x, d$time, d$event)
  assessment <- x[1:7, , drop = FALSE]
  native_assessment <- native_x[1:7, , drop = FALSE]
  state <- adapter$extract_state(model, list(
    x = x, time = d$time, event = d$event,
    prediction_horizon = 7, time_unit = "month"
  ))
  replay <- adapter$predict_state(unserialize(serialize(state, NULL)), assessment)
  native_risk <- unname(as.numeric(stats::predict(
    renamed$engine, new_data = native_assessment, pred_type = "mort", n_thread = 1L
  )))
  native_survival <- unname(as.numeric(stats::predict(
    renamed$engine, new_data = native_assessment,
    pred_type = "surv", pred_horizon = 7, n_thread = 1L
  )))

  expect_identical(model$features, names(x))
  expect_identical(state$features, names(x))
  expect_identical(model$engine$get_names_x(), names(native_x))
  expect_identical(model$engine$tree_seeds, renamed$engine$tree_seeds)
  expect_identical(adapter$predict(model, assessment), native_risk)
  expect_identical(replay$risk, native_risk)
  expect_identical(replay$survival_probability, native_survival)
  expect_identical(adapter$predict_survival(model, x, d$time, d$event,
    assessment, prediction_horizon = 7)$probability, native_survival)
})

test_that("aorsf survival executes through CV and final prediction", {
  skip_if_not_installed("aorsf")
  d <- make_sb_fixture(n = 60L, seed = 8432L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  original_params <- sb_additional_survival_params
  local_mocked_bindings(
    sb_additional_survival_params = function(learner) {
      params <- original_params(learner)
      if (identical(learner, "aorsf_survival")) {
        params$num_trees <- 25L
      }
      params
    }
  )
  pipeline <- sb_pipeline(learner = "aorsf_survival")
  cv <- sb_run_cv(task, pipeline, folds = 3L, seed = 8433L)
  fit <- sb_fit_final(task, pipeline)
  prediction <- sb_predict_final(
    fit,
    d[1:6, c("patient_id", "gene_c", "gene_a", "gene_b")],
    "patient_id"
  )

  expect_identical(cv$status, "success")
  expect_identical(cv$complete_oof, TRUE)
  expect_identical(cv$complete_survival_probability_oof, TRUE)
  expect_identical(
    cv$provenance$aorsf_version,
    as.character(utils::packageVersion("aorsf"))
  )
  expect_identical(fit$status, "success")
  expect_identical(
    fit$provenance$aorsf_version,
    as.character(utils::packageVersion("aorsf"))
  )
  expect_identical(fit$learner_state$learner, "aorsf_survival")
  expect_identical(prediction$status, "success")
  expect_true(all(is.finite(prediction$predictions$risk)))
  expect_true(all(prediction$predictions$survival_probability >= 0 &
    prediction$predictions$survival_probability <= 1))
})
