test_that("additive kernel query batching preserves feature and sample ordering", {
  basis <- matrix(c(-2, 0, 3, 1, 4, 7), nrow = 3L)
  colnames(basis) <- c("gene A", "gene-B")
  x <- matrix(c(-4, 2, 6, 0, 5, 9), nrow = 3L)
  colnames(x) <- colnames(basis)
  for (columns in list(1L, 1:2)) {
    b <- basis[, columns, drop = FALSE]
    queries <- x[, columns, drop = FALSE]
    result <- sb_survivalsvm_additive_kernel_matrix(b, queries)
    individual <- do.call(cbind, lapply(seq_len(nrow(queries)), function(i) {
      sb_survivalsvm_additive_kernel_matrix(b, queries[i, , drop = FALSE])
    }))
    expect_identical(result, individual)
    expect_identical(
      sb_survivalsvm_additive_kernel_matrix(b, queries[3:1, , drop = FALSE]),
      result[, 3:1, drop = FALSE]
    )
  }
})

make_survivalsvm_additive_fixture <- function(n = 72L, seed = 10032L) {
  withr::with_seed(seed, {
    gene_a <- stats::rnorm(n)
    gene_b <- stats::rnorm(n)
    gene_c <- stats::rnorm(n)
    gene_d <- stats::rnorm(n)
    linear_predictor <- 0.8 * gene_a - 0.6 * gene_b +
      0.35 * gene_c + 0.2 * gene_d
    event_time <- stats::rexp(n, rate = 0.08 * exp(linear_predictor))
    censor_time <- stats::rexp(n, rate = 0.03)
    data.frame(
      patient_id = sprintf("ASVM%03d", seq_len(n)),
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

test_that("additive survival-SVM formulations are registered native learners", {
  learners <- c(
    "survivalsvm_additive",
    "survivalsvm_vanbelle1_additive",
    "survivalsvm_vanbelle2_additive",
    "survivalsvm_hybrid_additive"
  )
  expect_setequal(sb_survivalsvm_additive_learner_ids(), learners)
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
    expect_identical(params$kernel, "add_kernel")
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

    pipeline <- sb_pipeline(
      learner = learner,
      survivalsvm_gamma_mu = 0.5
    )
    expect_identical(pipeline$learner_params$kernel, "add_kernel")
    expect_identical(
      pipeline$risk_type,
      "calibrated_survival_svm_kernel_score"
    )
    expect_identical(
      pipeline$risk_reference,
      "survivalsvm_additive_kernel_native"
    )
    expect_identical(
      pipeline$risk_formula,
      "negative_native_kernel_score_times_analysis_cox_coefficient"
    )
  }
})

test_that("additive survival-SVM states exactly replay native kernel scores", {
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  data <- make_survivalsvm_additive_fixture()
  predictors <- c("gene A", "gene-B", "gene_c", "gene_d")
  x <- data[predictors]
  rownames(x) <- data$patient_id
  assessment <- x[c(3L, 18L, 41L, 67L), , drop = FALSE]
  learners <- sb_survivalsvm_additive_learner_ids()
  risks <- list()
  withr::local_seed(10033L)

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

    expect_s3_class(state, "stormbreaker_survivalsvm_additive_learner_state")
    expect_false("engine" %in% names(state))
    expect_null(rownames(state$basis))
    expect_identical(
      state$training_data_retention,
      "required_scaled_kernel_basis_without_outcomes"
    )
    expect_identical(dim(state$basis), c(length(state$row_weight), 4L))
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

test_that("additive survival-SVM formulations execute fold-local CV and final fit", {
  skip_if_not_installed("survivalsvm", minimum_version = "0.0.6")

  data <- make_survivalsvm_additive_fixture(seed = 10034L)
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

  for (learner in sb_survivalsvm_additive_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 3L,
      survivalsvm_gamma_mu = 1
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 10035L)
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
