make_new_method_binary_fixture <- function(n = 72L, p = 8L, seed = 1701L) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n)
    colnames(x) <- paste0("gene_", seq_len(p))
    outcome <- as.integer(
      x[, 1L] - 0.7 * x[, 2L] + 0.3 * x[, 3L] +
        stats::rnorm(n, sd = 0.6) > 0
    )
    data.frame(
      patient_id = sprintf("NM%03d", seq_len(n)),
      outcome = outcome,
      x,
      check.names = FALSE
    )
  })
}

make_new_method_survival_fixture <- function(n = 72L, p = 5L,
                                             seed = 1703L) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n)
    colnames(x) <- paste0("gene_", seq_len(p))
    event_time <- stats::rexp(n, rate = exp(0.45 * x[, 1L] - 0.3 * x[, 2L]))
    censor_time <- stats::rexp(n, rate = 0.35)
    data.frame(
      patient_id = sprintf("MS%03d", seq_len(n)),
      time = pmin(event_time, censor_time) + 0.05,
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    )
  })
}

test_that("K-TSP keeps the declared class-one orientation", {
  skip_if_not_installed("switchBox", minimum_version = "1.44.0")

  data <- make_new_method_binary_fixture()
  x <- data[setdiff(names(data), c("patient_id", "outcome"))]
  adapter <- sb_omics_specialist_binary_learner_adapter("ktsp_binary")
  model <- adapter$fit(x, data$outcome)
  native <- switchBox::SWAP.Train.KTSP(
    inputMat = t(data.matrix(x)),
    phenoGroup = factor(data$outcome, levels = c(0L, 1L)),
    krange = seq_len(4L),
    FilterFunc = switchBox::SWAP.Filter.Wilcoxon,
    handleTies = FALSE,
    disjoint = TRUE,
    featureNo = 8L,
    verbose = FALSE
  )
  assessment <- x[1:7, , drop = FALSE]
  native_statistic <- switchBox::SWAP.KTSP.Statistics(
    t(data.matrix(assessment)),
    native
  )$statistics
  expected <- unname(native_statistic / nrow(native$TSPs) + 0.5)

  expect_identical(model$engine$labels, c("0", "1"))
  expect_equal(model$engine$pairs, unname(native$TSPs), tolerance = 0)
  expect_equal(model$engine$pair_scores, unname(as.numeric(native$score)), tolerance = 0)
  expect_equal(adapter$predict(model, assessment), expected, tolerance = 0)
  expect_identical(
    adapter$predict(model, assessment) > 0.5,
    unname(switchBox::SWAP.KTSP.Classify(
      t(data.matrix(assessment)),
      native
    ) == "1")
  )
})

test_that("MTLR accepts its native interval-by-coefficient weight matrix", {
  skip_if_not_installed("MTLR", minimum_version = "0.2.2")

  data <- make_new_method_survival_fixture()
  x <- data[setdiff(names(data), c("patient_id", "time", "event"))]
  adapter <- sb_mtlr_survival_learner_adapter()
  model <- adapter$fit(x, data$time, data$event)

  expect_identical(ncol(model$engine$weight_matrix), ncol(x) + 1L)
  expect_identical(
    nrow(model$engine$weight_matrix),
    length(model$engine$time_points)
  )
})

test_that("new methods add their generated legal routes", {
  registry <- sb_pipeline_registry()
  binary_selector_count <- length(expected_fixed_learner_selectors("binary"))
  survival_selector_count <- length(expected_fixed_learner_selectors("survival"))
  expected_task <- c(
    stats::setNames(
      rep("binary", length(sb_omics_specialist_binary_learner_ids())),
      sb_omics_specialist_binary_learner_ids()
    ),
    mtlr_survival = "survival"
  )
  expected_routes <- c(
    stats::setNames(
      rep(
        binary_selector_count,
        length(sb_omics_specialist_binary_learner_ids())
      ),
      sb_omics_specialist_binary_learner_ids()
    ),
    mtlr_survival = survival_selector_count
  )

  for (learner in names(expected_task)) {
    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_identical(nrow(rows), unname(expected_routes[[learner]]))
    expect_identical(unique(rows$task_type), unname(expected_task[[learner]]))
    expect_identical(anyDuplicated(rows$selector), 0L)
    expect_true(all(rows$admission_status == "admitted"))
    expect_s3_class(sb_pipeline(learner = learner), "stormbreaker_pipeline")
  }
})

test_that("bias-reduced logistic remains finite under complete separation", {
  skip_if_not_installed("brglm2", minimum_version = "1.1.0")

  n <- 40L
  separating_feature <- seq(-1, 1, length.out = n)
  data <- data.frame(
    patient_id = sprintf("NM%03d", seq_len(n)),
    outcome = as.integer(separating_feature > 0),
    gene_1 = separating_feature,
    gene_2 = sin(seq_len(n)),
    gene_3 = cos(seq_len(n) / 3),
    gene_4 = rep(c(-1, 1), length.out = n)
  )
  x <- data[setdiff(names(data), c("patient_id", "outcome"))]
  adapter <- sb_omics_specialist_binary_learner_adapter(
    "brglm2_logistic"
  )
  model <- adapter$fit(x, data$outcome)
  probability <- adapter$predict(model, x[1:8, , drop = FALSE])
  state <- unserialize(serialize(adapter$extract_state(model), NULL))

  expect_true(all(is.finite(model$engine$coefficient)))
  expect_true(all(probability > 0 & probability < 1))
  expect_equal(
    adapter$predict_state(state, x[1:8, , drop = FALSE])$probability,
    probability,
    tolerance = 0
  )
  expect_identical(state$training_data_retention, "none")
})

test_that("Bayesian Cauchy logistic replays its native posterior mode", {
  skip_if_not_installed("arm", minimum_version = "1.15.3")

  data <- make_new_method_binary_fixture(n = 48L, p = 70L, seed = 1712L)
  x <- data[setdiff(names(data), c("patient_id", "outcome"))]
  adapter <- sb_omics_specialist_binary_learner_adapter(
    "bayesglm_logistic_cauchy"
  )
  params <- adapter$params
  set.seed(2719)
  rng_before <- .Random.seed
  model <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, rng_before)

  design <- cbind("(Intercept)" = 1, data.matrix(x))
  native <- suppressWarnings(arm::bayesglm.fit(
    x = design,
    y = data$outcome,
    family = stats::binomial(link = params$link),
    control = list(maxit = params$maximum_iterations),
    intercept = TRUE,
    prior.mean = params$prior_mean,
    prior.scale = params$prior_scale,
    prior.df = params$prior_df,
    prior.mean.for.intercept = params$intercept_prior_mean,
    prior.scale.for.intercept = params$intercept_prior_scale,
    prior.df.for.intercept = params$intercept_prior_df,
    scaled = params$scaled,
    Warning = FALSE
  ))
  assessment <- x[1:7, , drop = FALSE]
  expected <- stats::plogis(
    unname(native$coefficients[[1L]]) +
      as.vector(data.matrix(assessment) %*%
        unname(native$coefficients[-1L]))
  )
  state <- unserialize(serialize(adapter$extract_state(model), NULL))

  expect_equal(
    unname(model$engine$coefficient),
    unname(native$coefficients[-1L]),
    tolerance = 0
  )
  expect_equal(adapter$predict(model, assessment), expected, tolerance = 0)
  expect_equal(
    adapter$predict_state(state, assessment)$probability,
    expected,
    tolerance = 0
  )
  expect_identical(state$training_data_retention, "none")
})

test_that("hinge boosting stores a portable native margin", {
  skip_if_not_installed("bst", minimum_version = "0.3-24")

  data <- make_new_method_binary_fixture(seed = 1713L)
  x <- data[setdiff(names(data), c("patient_id", "outcome"))]
  adapter <- sb_omics_specialist_binary_learner_adapter(
    "bst_hinge_binary"
  )
  params <- adapter$params
  model <- adapter$fit(x, data$outcome)
  labels <- ifelse(data$outcome == 1L, 1, -1)
  native <- bst::bst(
    data.matrix(x),
    labels,
    cost = params$cost,
    family = params$family,
    learner = params$base_learner,
    ctrl = bst::bst_control(
      mstop = params$iterations,
      nu = params$learning_rate,
      center = params$center,
      intercept = params$fit_intercept,
      trace = FALSE
    )
  )
  assessment <- x[1:7, , drop = FALSE]
  native_margin <- as.numeric(stats::predict(
    native,
    newdata = data.matrix(assessment),
    type = "response"
  ))
  portable_margin <- model$engine$intercept +
    as.vector(data.matrix(assessment) %*% model$engine$coefficient)
  state <- unserialize(serialize(adapter$extract_state(model), NULL))

  expect_equal(portable_margin, native_margin, tolerance = 1e-12)
  expect_equal(
    adapter$predict(model, assessment),
    stats::plogis(2 * native_margin),
    tolerance = 1e-12
  )
  expect_equal(
    adapter$predict_state(state, assessment)$probability,
    adapter$predict(model, assessment),
    tolerance = 0
  )
  expect_identical(state$training_data_retention, "none")
})

test_that("HDA preserves its native one-dimensional posterior", {
  skip_if_not_installed("hda", minimum_version = "0.2-14")

  data <- make_new_method_binary_fixture(n = 96L, p = 8L, seed = 1714L)
  x <- data[setdiff(names(data), c("patient_id", "outcome"))]
  adapter <- sb_omics_specialist_binary_learner_adapter("hda_binary")
  params <- adapter$params
  model <- adapter$fit(x, data$outcome)
  native <- NULL
  invisible(utils::capture.output(
    native <- hda::hda(
      data.matrix(x),
      grouping = factor(data$outcome, levels = c(0L, 1L)),
      newdim = params$reduced_dimension,
      crule = TRUE,
      reg.lamb = params$regularization_lambda,
      reg.gamm = params$regularization_gamma,
      initial.loadings = NULL,
      noutit = params$outer_iterations,
      ninit = params$inner_iterations,
      verbose = FALSE
    ),
    type = "output"
  ))
  assessment <- x[1:7, , drop = FALSE]
  native_prediction <- suppressWarnings(stats::predict(
    native,
    data.matrix(assessment),
    task = "c"
  ))
  expected <- unname(native_prediction$posteriors[, 2L])
  state <- unserialize(serialize(adapter$extract_state(model), NULL))

  expect_equal(adapter$predict(model, assessment), expected, tolerance = 0)
  expect_equal(
    adapter$predict_state(state, assessment)$probability,
    expected,
    tolerance = 0
  )
  expect_setequal(
    names(model$engine),
    c("hda.loadings", "reduced.dimension", "naivebayes")
  )
  expect_identical(model$engine$naivebayes$levels, c("0", "1"))
  expect_identical(state$training_data_retention, "none")
})

test_that("DWD supports p greater than n and replays its native margin", {
  skip_if_not_installed("DWDLargeR", minimum_version = "0.2.0")

  data <- make_new_method_binary_fixture(n = 28L, p = 90L, seed = 1706L)
  x <- data[setdiff(names(data), c("patient_id", "outcome"))]
  x[] <- sweep(
    data.matrix(x),
    2L,
    10^seq(-3, 3, length.out = ncol(x)),
    "*"
  )
  adapter <- sb_omics_specialist_binary_learner_adapter("dwd_binary")
  set.seed(2718)
  rng_before <- .Random.seed
  model <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, rng_before)
  assessment <- x[1:6, , drop = FALSE]
  native_margin <- model$engine$intercept +
    as.vector(data.matrix(assessment) %*% model$engine$coefficient)
  probability <- adapter$predict(model, assessment)

  expect_equal(probability, stats::plogis(2 * native_margin), tolerance = 0)
  expect_true(is.finite(model$engine$penalty))
  training_probability <- adapter$predict(model, x)
  expect_equal(
    mean((training_probability > 0.5) != data$outcome) * 100,
    model$engine$training_error,
    tolerance = 1e-12
  )
  expect_equal(
    adapter$predict_state(
      unserialize(serialize(adapter$extract_state(model), NULL)),
      assessment
    )$probability,
    probability,
    tolerance = 0
  )
})

test_that("DWD maps removed all-zero features back to the original schema", {
  skip_if_not_installed("DWDLargeR", minimum_version = "0.2.0")

  data <- make_new_method_binary_fixture(n = 28L, p = 45L, seed = 1710L)
  x <- data[setdiff(names(data), c("patient_id", "outcome"))]
  x$gene_10 <- 0
  x$gene_20 <- 0
  adapter <- sb_omics_specialist_binary_learner_adapter("dwd_binary")
  model <- adapter$fit(x, data$outcome)

  expect_identical(
    unname(model$engine$coefficient[c("gene_10", "gene_20")]),
    c(0, 0)
  )
  expect_length(adapter$predict(model, x[1:5, , drop = FALSE]), 5L)
})

test_that("MTLR risk, horizon probability, and lean state match native output", {
  skip_if_not_installed("MTLR", minimum_version = "0.2.2")

  data <- make_new_method_survival_fixture(seed = 1707L)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[1:7, , drop = FALSE]
  horizon <- unname(stats::quantile(data$time, 0.35))
  adapter <- sb_mtlr_survival_learner_adapter()
  model <- adapter$fit(x, data$time, data$event)
  internal <- assessment
  names(internal) <- model$internal_features
  native_risk <- -as.numeric(stats::predict(
    model$engine,
    newdata = internal,
    type = "mean_time"
  ))
  native_curve <- stats::predict(
    model$engine,
    newdata = internal,
    type = "prob_times",
    times = horizon
  )
  native_probability <- unname(as.numeric(
    native_curve[nrow(native_curve), -1L]
  ))
  probability <- adapter$predict_survival(
    model,
    x,
    data$time,
    data$event,
    assessment,
    horizon
  )
  state <- unserialize(serialize(adapter$extract_state(model, list(
    x = x,
    time = data$time,
    event = data$event,
    prediction_horizon = horizon,
    time_unit = "month"
  )), NULL))
  replay <- adapter$predict_state(state, assessment)

  expect_equal(adapter$predict(model, assessment), native_risk, tolerance = 0)
  expect_equal(probability$probability, native_probability, tolerance = 0)
  expect_equal(replay$risk, native_risk, tolerance = 0)
  expect_equal(
    replay$survival_probability,
    native_probability,
    tolerance = 0
  )
  expect_false(any(c("x", "y", "response", "Call") %in%
    names(state$engine)))
  expect_identical(state$training_data_retention, "none")
})

test_that("all four methods execute through CV and final state prediction", {
  binary_data <- make_new_method_binary_fixture(seed = 1708L)
  binary_features <- setdiff(
    names(binary_data),
    c("patient_id", "outcome")
  )
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    binary_features
  )
  dependency <- c(
    brglm2_logistic = "brglm2",
    bayesglm_logistic_cauchy = "arm",
    bst_hinge_binary = "bst",
    hda_binary = "hda",
    dwd_binary = "DWDLargeR",
    ktsp_binary = "switchBox"
  )
  for (learner in names(dependency)) {
    skip_if_not_installed(dependency[[learner]])
    pipeline <- sb_pipeline(learner = learner)
    comparison <- sb_run_binary_cv(
      binary_task,
      pipeline,
      folds = 3L,
      seed = 1709L
    )
    fit <- sb_fit_final(binary_task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      binary_data[1:6, c("patient_id", rev(binary_features)), drop = FALSE],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_true(comparison$complete_oof)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
  }

  skip_if_not_installed("MTLR", minimum_version = "0.2.2")
  survival_data <- make_new_method_survival_fixture(seed = 1710L)
  survival_features <- setdiff(
    names(survival_data),
    c("patient_id", "time", "event")
  )
  horizon <- unname(stats::quantile(survival_data$time, 0.3))
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    survival_features,
    "month",
    prediction_horizon = horizon
  )
  pipeline <- sb_pipeline(learner = "mtlr_survival")
  comparison <- sb_run_cv(survival_task, pipeline, folds = 3L, seed = 1711L)
  fit <- sb_fit_final(survival_task, pipeline)
  prediction <- sb_predict_final(
    unserialize(serialize(fit, NULL)),
    survival_data[
      1:6,
      c("patient_id", rev(survival_features)),
      drop = FALSE
    ],
    "patient_id"
  )

  expect_identical(comparison$status, "success")
  expect_true(comparison$complete_oof)
  expect_true(comparison$complete_survival_probability_oof)
  expect_identical(fit$status, "success")
  expect_identical(prediction$status, "success")
})
