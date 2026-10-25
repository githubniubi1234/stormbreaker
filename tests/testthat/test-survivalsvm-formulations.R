make_survivalsvm_formulation_fixture <- function(n = 72L, seed = 9303L) {
  withr::with_seed(seed, {
    gene_a <- stats::rnorm(n)
    gene_b <- stats::rnorm(n)
    gene_c <- stats::rnorm(n)
    linear_predictor <- 0.7 * gene_a - 0.5 * gene_b + 0.25 * gene_c
    event_time <- stats::rexp(n, rate = 0.08 * exp(linear_predictor))
    censor_time <- stats::rexp(n, rate = 0.035)
    data.frame(
      patient_id = sprintf("SVM%03d", seq_len(n)),
      time = pmax(pmin(event_time, censor_time), 0.01),
      event = as.integer(event_time <= censor_time),
      `gene A` = gene_a,
      `gene-B` = gene_b,
      gene_c = gene_c,
      check.names = FALSE
    )
  })
}

test_that("survival-SVM formulations are distinct fixed native models", {
  learners <- c(
    "survivalsvm_vanbelle1",
    "survivalsvm_vanbelle2",
    "survivalsvm_hybrid"
  )
  expect_setequal(
    sb_survivalsvm_linear_learner_ids(),
    c("survivalsvm", learners)
  )

  expected_types <- c("vanbelle1", "vanbelle2", "hybrid")
  for (index in seq_along(learners)) {
    params <- sb_survivalsvm_params(
      gamma_mu = 0.5,
      learner = learners[[index]]
    )
    expect_identical(params$gamma_mu, 0.5)
    expect_identical(params$type, expected_types[[index]])
    expect_identical(params$diff_meth, "makediff3")
    expect_identical(params$kernel, "lin_kernel")
    expect_identical(params$opt_meth, "quadprog")
    expect_identical(params$standardize, TRUE)
    expect_identical(params$calibration_ties, "breslow")
  }

  regression <- sb_survivalsvm_params(gamma_mu = 0.5)
  expect_identical(regression$type, "regression")
  expect_false("diff_meth" %in% names(regression))

  registry <- sb_pipeline_registry()
  expected_selectors <- expected_fixed_learner_selectors("survival")
  for (index in seq_along(learners)) {
    learner <- learners[[index]]
    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_identical(unique(rows$task_type), "survival")
    expect_setequal(rows$selector, expected_selectors)
    expect_true(all(grepl(
      "survivalsvm", rows$optional_dependency,
      fixed = TRUE
    )))
    expect_identical(unique(rows$tuning_semantics), "none")

    pipeline <- sb_pipeline(
      learner = learner,
      survivalsvm_gamma_mu = 0.5
    )
    expect_identical(pipeline$learner, learner)
    expect_identical(pipeline$learner_params$type, expected_types[[index]])
    expect_identical(pipeline$learner_params$diff_meth, "makediff3")
    expect_identical(
      pipeline$id,
      paste0("none__", learner, "__gamma0.5")
    )
  }
})

test_that("new survival-SVM formulations preserve native linear scores", {
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  data <- make_survivalsvm_formulation_fixture()
  predictors <- c("gene A", "gene-B", "gene_c")
  x <- data[predictors]
  assessment <- x[c(3L, 18L, 41L, 67L), , drop = FALSE]
  learners <- c(
    "survivalsvm_vanbelle1",
    "survivalsvm_vanbelle2",
    "survivalsvm_hybrid"
  )
  native_types <- c("vanbelle1", "vanbelle2", "hybrid")
  risks <- list()
  withr::local_seed(9304L)

  for (index in seq_along(learners)) {
    learner <- learners[[index]]
    adapter <- sb_survivalsvm_learner_adapter(
      gamma_mu = 1,
      learner = learner
    )
    rng_before <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, rng_before)
    repeated <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, rng_before)
    expect_identical(model$coefficient, repeated$coefficient)

    scaled_assessment <- as.data.frame(
      sweep(
        sweep(as.matrix(assessment), 2L, model$audit$center, "-"),
        2L,
        model$audit$scale,
        "/"
      ),
      optional = TRUE,
      stringsAsFactors = FALSE
    )
    names(scaled_assessment) <- predictors
    native <- -as.numeric(stats::predict(
      model$engine,
      newdata = scaled_assessment
    )$predicted) * model$audit$calibration_coefficient
    observed <- adapter$predict(model, assessment)
    centered <- function(value) value - value[[1L]]

    expect_equal(centered(observed), centered(native), tolerance = 1e-10)
    expect_identical(model$learner, learner)
    expect_identical(model$params$type, native_types[[index]])
    expect_identical(model$audit$diff_meth, "makediff3")

    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = 7,
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, assessment)
    expected_survival <- adapter$predict_survival(
      model,
      x,
      data$time,
      data$event,
      assessment,
      prediction_horizon = 7
    )$probability

    expect_identical(state$learner, learner)
    expect_false("engine" %in% names(state))
    expect_equal(replay$risk, observed, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      expected_survival,
      tolerance = 0
    )
    risks[[learner]] <- observed
  }

  for (pair in utils::combn(names(risks), 2L, simplify = FALSE)) {
    expect_gt(max(abs(risks[[pair[[1L]]]] - risks[[pair[[2L]]]])), 1e-6)
  }
})

test_that("new survival-SVM formulations execute fold-local CV and final fit", {
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  data <- make_survivalsvm_formulation_fixture(seed = 9305L)
  predictors <- c("gene A", "gene-B", "gene_c")
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )
  learners <- c(
    "survivalsvm_vanbelle1",
    "survivalsvm_vanbelle2",
    "survivalsvm_hybrid"
  )

  for (learner in learners) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 2L,
      survivalsvm_gamma_mu = 1
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 9306L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      data[1:6, c("patient_id", rev(predictors)), drop = FALSE],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_identical(comparison$complete_survival_probability_oof, TRUE)
    expect_identical(sum(comparison$runs$status == "success"), 3L)
    expect_identical(
      comparison$provenance$survivalsvm_version,
      as.character(utils::packageVersion("survivalsvm"))
    )
    expect_identical(fit$status, "success")
    expect_identical(fit$learner_state$learner, learner)
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$risk, 6L)
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(prediction$predictions$survival_probability >= 0 &
      prediction$predictions$survival_probability <= 1))
  }
})
