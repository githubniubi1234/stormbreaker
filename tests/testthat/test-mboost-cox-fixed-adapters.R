test_that("mboost Cox learners keep three distinct native algorithms", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  expect_identical(
    sb_mboost_cox_learner_ids(),
    c(
      "gamboost_cox",
      "glmboost_cox",
      "blackboost_cox"
    )
  )

  d <- make_sb_fixture(n = 60L, seed = 7961L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  expected_class <- c(
    gamboost_cox = "gamboost",
    glmboost_cox = "glmboost",
    blackboost_cox = "blackboost"
  )

  for (learner in sb_mboost_cox_learner_ids()) {
    params <- sb_mboost_cox_params(learner)
    adapter <- sb_mboost_cox_learner_adapter(learner, params)
    model <- adapter$fit(x, d$time, d$event)

    expect_s3_class(model$engine, expected_class[[learner]])
    expect_identical(as.integer(mboost::mstop(model$engine)), 100L)
    expect_identical(model$params$mstop, 100L)
    expect_identical(model$params$nu, 0.1)
    expect_identical(model$params$family, "CoxPH")
  }
})

test_that("mboost Cox adapters return native risk, survival, and replayable state", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  d <- make_sb_fixture(n = 60L, seed = 7962L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(3L, 22L, 51L), , drop = FALSE]

  for (learner in sb_mboost_cox_learner_ids()) {
    adapter <- sb_mboost_cox_learner_adapter(learner)
    model <- adapter$fit(x, d$time, d$event)
    risk <- adapter$predict(model, new_x)
    native_risk <- as.numeric(stats::predict(
      model$engine,
      newdata = new_x,
      type = "link"
    ))
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
    replay <- adapter$predict_state(
      unserialize(serialize(state, NULL)),
      new_x
    )

    expect_equal(risk, native_risk, tolerance = 1e-14)
    expect_true(all(is.finite(risk)))
    expect_true(all(probability$probability >= 0 &
      probability$probability <= 1))
    expect_equal(replay$risk, risk, tolerance = 1e-14)
    expect_equal(
      replay$survival_probability,
      probability$probability,
      tolerance = 1e-14
    )
    expect_identical(
      state$training_data_retention,
      "required_by_mboost_prediction_state"
    )
  }
})

test_that("mboost Cox routes complete CV and final outcome-blind prediction", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")

  d <- make_sb_fixture(n = 60L, seed = 7963L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = 7
  )
  new_data <- d[
    c(5L, 27L, 54L),
    c("patient_id", "gene_c", "gene_a", "gene_b")
  ]

  for (learner in sb_mboost_cox_learner_ids()) {
    pipeline <- sb_pipeline(learner = learner)
    result <- sb_run_cv(task, pipeline, folds = 3L, seed = 7964L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(fit, new_data, "patient_id")

    expect_identical(result$status, "success")
    expect_true(result$complete_oof)
    expect_true(result$complete_survival_probability_oof)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(prediction$predictions$survival_probability >= 0 &
      prediction$predictions$survival_probability <= 1))
  }
})
