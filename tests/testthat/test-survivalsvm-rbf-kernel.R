make_survivalsvm_rbf_fixture <- function(n = 72L, seed = 11032L) {
  withr::with_seed(seed, {
    gene_a <- stats::rnorm(n)
    gene_b <- stats::rnorm(n)
    gene_c <- stats::rnorm(n)
    gene_d <- stats::rnorm(n)
    nonlinear_predictor <- 0.7 * gene_a - 0.5 * gene_b +
      0.35 * gene_c^2 - 0.25 * gene_a * gene_d
    event_time <- stats::rexp(
      n,
      rate = 0.07 * exp(nonlinear_predictor)
    )
    censor_time <- stats::rexp(n, rate = 0.03)
    data.frame(
      patient_id = sprintf("RBF%03d", seq_len(n)),
      time = pmax(pmin(event_time, censor_time), 0.01),
      event = as.integer(event_time <= censor_time),
      `gene A` = gene_a,
      `gene-B` = gene_b,
      gene_c = gene_c,
      gene_d = gene_d,
      check.names = FALSE
    )
  })
}

test_that("kernel survival-SVM one-input boundaries are exact", {
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")
  data <- sb_example_survival(n = 200L, seed = 20260905L)[1:160, ]
  x <- data["gene_a"]
  unsupported <- c(
    sb_survivalsvm_additive_learner_ids(),
    setdiff(
      sb_survivalsvm_rbf_learner_ids(),
      "survivalsvm_rbf"
    )
  )

  for (learner in unsupported) {
    expect_error(
      sb_survivalsvm_learner_adapter(1, learner)$fit(
        x, data$time, data$event
      ),
      class = "stormbreaker_survivalsvm_too_few_features",
      info = learner
    )
  }

  expect_s3_class(
    sb_survivalsvm_learner_adapter(
      1, "survivalsvm_rbf"
    )$fit(x, data$time, data$event),
    "stormbreaker_survivalsvm_rbf"
  )
  expect_s3_class(
    sb_survivalsvm_learner_adapter(
      1, "survivalsvm_vanbelle1"
    )$fit(x, data$time, data$event),
    "stormbreaker_survivalsvm_vanbelle1"
  )
})

test_that("RBF survival-SVM formulations are registered native learners", {
  learners <- c(
    "survivalsvm_rbf",
    "survivalsvm_vanbelle1_rbf",
    "survivalsvm_vanbelle2_rbf",
    "survivalsvm_hybrid_rbf"
  )
  expect_setequal(sb_survivalsvm_rbf_learner_ids(), learners)
  expect_setequal(
    sb_survivalsvm_kernel_learner_ids(),
    c(sb_survivalsvm_additive_learner_ids(), learners)
  )
  expect_setequal(
    sb_survivalsvm_learner_ids(),
    c(
      sb_survivalsvm_linear_learner_ids(),
      sb_survivalsvm_kernel_learner_ids()
    )
  )
  expect_length(intersect(learners, sb_linear_survival_learner_ids()), 0L)

  expected_types <- c("regression", "vanbelle1", "vanbelle2", "hybrid")
  expected_selectors <- expected_fixed_learner_selectors("survival")
  registry <- sb_pipeline_registry()
  for (index in seq_along(learners)) {
    learner <- learners[[index]]
    params <- sb_survivalsvm_params(0.5, learner)
    expect_identical(params$type, expected_types[[index]])
    expect_identical(params$kernel, "rbf_kernel")
    expect_null(params$kernel_pars)
    expect_identical(
      params$kernel_parameter_rule,
      "analysis_median_squared_distance"
    )
    expect_identical(params$gamma_mu, 0.5)
    expect_identical(params$standardize, TRUE)
    expect_identical(params$opt_meth, "quadprog")
    if (identical(params$type, "regression")) {
      expect_false("diff_meth" %in% names(params))
    } else {
      expect_identical(params$diff_meth, "makediff3")
    }

    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_setequal(rows$selector, expected_selectors)
    expect_true(all(grepl(
      "survivalsvm", rows$optional_dependency,
      fixed = TRUE
    )))
    expect_identical(unique(rows$tuning_semantics), "none")
    expect_true(all(grepl(
      "analysis-median-squared-distance", rows$parameter_semantics,
      fixed = TRUE
    )))

    pipeline <- sb_pipeline(
      learner = learner,
      survivalsvm_gamma_mu = 0.5
    )
    expect_identical(pipeline$learner_params$kernel, "rbf_kernel")
    expect_identical(
      pipeline$risk_type,
      "calibrated_survival_svm_kernel_score"
    )
    expect_identical(
      pipeline$risk_reference,
      "survivalsvm_rbf_kernel_native"
    )
    expect_identical(
      pipeline$risk_formula,
      "negative_native_kernel_score_times_analysis_cox_coefficient"
    )
  }
})

test_that("RBF survival-SVM kernel matrix follows the native formula", {
  basis <- matrix(c(0, 0, 1, -1), nrow = 2L, byrow = TRUE)
  assessment <- matrix(c(0, 1, 2, -1), nrow = 2L, byrow = TRUE)
  observed <- sb_survivalsvm_rbf_kernel_matrix(
    basis,
    assessment,
    kernel_pars = c(0.5, 0.5)
  )
  expected <- matrix(
    c(exp(-1), exp(-5), exp(-5), exp(-1)),
    nrow = 2L,
    byrow = TRUE
  )

  expect_equal(observed, expected, tolerance = 1e-15)
})

test_that("RBF survival-SVM states exactly replay native kernel scores", {
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  data <- make_survivalsvm_rbf_fixture()
  predictors <- c("gene A", "gene-B", "gene_c", "gene_d")
  x <- data[predictors]
  rownames(x) <- data$patient_id
  assessment <- x[c(4L, 19L, 42L, 68L), , drop = FALSE]
  learners <- sb_survivalsvm_rbf_learner_ids()
  risks <- list()
  withr::local_seed(11033L)

  for (learner in learners) {
    adapter <- sb_survivalsvm_learner_adapter(
      gamma_mu = 1,
      learner = learner
    )
    rng_before <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, rng_before)

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
    expect_equal(observed, native, tolerance = 1e-10)
    expect_equal(
      model$audit$kernel_pars,
      stats::median(as.numeric(stats::dist(scale(x)))^2),
      tolerance = 0
    )

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

    expect_s3_class(state, "stormbreaker_survivalsvm_rbf_learner_state")
    expect_false("engine" %in% names(state))
    expect_null(rownames(state$basis))
    expect_identical(
      state$training_data_retention,
      "required_scaled_kernel_basis_without_outcomes"
    )
    expect_identical(dim(state$basis), c(length(state$row_weight), 4L))
    expect_equal(state$kernel_pars, model$audit$kernel_pars, tolerance = 0)
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

test_that("RBF survival-SVM formulations execute fold-local CV and final fit", {
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  data <- make_survivalsvm_rbf_fixture(seed = 11034L)
  predictors <- c("gene A", "gene-B", "gene_c", "gene_d")
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    predictors,
    "month",
    prediction_horizon = 7
  )

  for (learner in sb_survivalsvm_rbf_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 3L,
      survivalsvm_gamma_mu = 1
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 11035L)
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
    expect_identical(fit$status, "success")
    expect_identical(fit$learner_state$learner, learner)
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$risk, 6L)
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(prediction$predictions$survival_probability >= 0 &
      prediction$predictions$survival_probability <= 1))
  }
})
