test_that("additional survival specifications are explicit", {
  survreg <- sb_additional_survival_params("survreg_weibull")
  rpart <- sb_additional_survival_params("rpart_survival")
  rsf <- sb_additional_survival_params("rsf_survival")
  ctree <- sb_additional_survival_params("ctree_survival")
  cforest <- sb_additional_survival_params("cforest_survival")
  bart <- sb_additional_survival_params("bart_survival")

  expect_identical(survreg, list(
    distribution = "weibull",
    scale = "estimated"
  ))
  expect_identical(rpart$method, "exp")
  expect_identical(rpart[c("cp", "minsplit", "maxdepth", "xval")], list(
    cp = 0.01,
    minsplit = 20L,
    maxdepth = 30L,
    xval = 0L
  ))
  expect_identical(rsf$num_trees, 1000L)
  expect_identical(rsf$mtry_rule, "randomForestSRC_default")
  expect_identical(rsf$min_node_size, 5L)
  expect_identical(rsf$splitrule, "logrank")
  expect_identical(rsf$seed, 1L)
  expect_identical(ctree$mincriterion, 0.95)
  expect_identical(ctree$minsplit, 20L)
  expect_identical(ctree$minbucket, 7L)
  expect_identical(cforest$num_trees, 500L)
  expect_identical(cforest$mtry_rule, "ceiling_sqrt_p")
  expect_identical(cforest$replace, TRUE)
  expect_identical(cforest$sample_fraction, 1)
  expect_identical(cforest$seed, 1L)
  expect_identical(bart, list(
    type = "pbart",
    num_trees = 50L,
    num_cutpoints = 100L,
    time_grid_size = 50L,
    posterior_draws = 250L,
    burn_in = 100L,
    keep_every = 1L,
    prior_k = 2,
    tree_power = 2,
    tree_base = 0.95,
    sparse = FALSE,
    seed = 1L,
    cores = 1L
  ))
})

test_that("Weibull SurvReg nonconvergence is a typed learner failure", {
  value <- seq(-1, 1, length.out = 18L)
  data <- data.frame(
    patient_id = sprintf("patient_%02d", seq_along(value)),
    time = exp(value),
    event = as.integer(!seq_along(value) %in% c(6L, 12L, 18L)),
    gene_a = value
  )
  adapter <- sb_additional_survival_learner_adapter(
    "survreg_weibull"
  )

  error <- capture_sb_error(
    adapter$fit(data["gene_a"], data$time, data$event)
  )
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    "gene_a",
    "month",
    prediction_horizon = 1
  )
  pipeline <- sb_pipeline(learner = "survreg_weibull")
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 7920L)
  fit <- sb_fit_final(task, pipeline)

  expect_s3_class(error, "stormbreaker_survreg_nonconvergence")
  expect_match(
    conditionMessage(error),
    "Ran out of iterations and did not converge",
    fixed = TRUE
  )
  expect_s3_class(error$parent, "warning")
  expect_identical(result$status, "failed")
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_survreg_nonconvergence"
  )
  expect_identical(fit$status, "failed")
  expect_identical(fit$error_stage, "learner_fit")
  expect_identical(fit$error_class, "stormbreaker_survreg_nonconvergence")
})

test_that("survival BART uses native posterior curves and declared risk", {
  skip_if_not_installed("BART", minimum_version = "2.9.0")

  d <- make_sb_fixture(n = 54L, seed = 7922L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(4L, 21L, 47L), , drop = FALSE]
  params <- sb_additional_survival_params("bart_survival")
  params$time_grid_size <- 20L
  params$posterior_draws <- 40L
  params$burn_in <- 20L
  adapter <- sb_additional_survival_learner_adapter(
    "bart_survival",
    params
  )
  model <- adapter$fit(x, d$time, d$event)
  curves <- sb_bart_survival_curves(model, new_x)

  grid <- model$engine$times
  augmented <- cbind(
    t = rep(grid, times = nrow(new_x)),
    as.matrix(new_x)[
      rep(seq_len(nrow(new_x)), each = length(grid)),
      ,
      drop = FALSE
    ]
  )
  if (length(model$engine$rm.const) > 0L) {
    augmented <- augmented[, model$engine$rm.const, drop = FALSE]
  }
  native <- NULL
  method <- get(
    "predict.survbart",
    envir = asNamespace("BART"),
    inherits = FALSE
  )
  invisible(utils::capture.output(
    native <- method(
      model$engine,
      augmented,
      mc.cores = 1L,
      openmp = FALSE
    ),
    type = "output"
  ))
  native_survival <- matrix(
    native$surv.test.mean,
    nrow = nrow(new_x),
    ncol = length(grid),
    byrow = TRUE
  )
  interval <- diff(c(0, grid))
  native_risk <- -as.vector(cbind(
    1,
    native_survival[, -ncol(native_survival), drop = FALSE]
  ) %*% interval)
  horizon <- stats::median(grid)
  horizon_column <- max(which(grid <= horizon))
  probability <- adapter$predict_survival(
    model,
    x,
    d$time,
    d$event,
    new_x,
    horizon
  )

  expect_s3_class(model$engine, "survbart")
  expect_equal(curves$survival, native_survival, tolerance = 1e-14)
  expect_equal(adapter$predict(model, new_x), native_risk, tolerance = 1e-14)
  expect_equal(
    probability$probability,
    native_survival[, horizon_column],
    tolerance = 1e-14
  )
  expect_identical(
    probability$source,
    "BART_discrete_time_posterior_survival_direct"
  )
})

test_that("additional survival adapters return risk, horizon survival, and state", {
  skip_if_not_installed("rpart", minimum_version = "4.1")
  skip_if_not_installed("randomForestSRC", minimum_version = "3.3.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")
  skip_if_not_installed("coin", minimum_version = "1.1")
  skip_if_not_installed("BART", minimum_version = "2.9.0")

  d <- make_sb_fixture(n = 60L, seed = 7921L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(4L, 23L, 49L), , drop = FALSE]
  learners <- sb_additional_survival_learner_ids()

  for (learner in learners) {
    params <- sb_additional_survival_params(learner)
    if (identical(learner, "rsf_survival")) {
      params$num_trees <- 20L
    }
    if (identical(learner, "cforest_survival")) {
      params$num_trees <- 20L
    }
    if (identical(learner, "bart_survival")) {
      params$time_grid_size <- 20L
      params$posterior_draws <- 40L
      params$burn_in <- 20L
    }
    if (identical(learner, "aorsf_survival")) {
      params$num_trees <- 20L
    }
    adapter <- sb_additional_survival_learner_adapter(learner, params)
    model <- adapter$fit(x, d$time, d$event)
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
    replay <- adapter$predict_state(
      unserialize(serialize(state, NULL)),
      new_x
    )

    expect_identical(state$learner, learner)
    expect_length(risk, nrow(new_x))
    expect_true(all(is.finite(risk)))
    expect_true(all(probability$probability >= 0 &
      probability$probability <= 1))
    expect_equal(replay$risk, risk, tolerance = 1e-14)
    expect_equal(
      replay$survival_probability,
      probability$probability,
      tolerance = 1e-14
    )
  }
})

test_that("additional survival states retain only their required engine form", {
  skip_if_not_installed("rpart", minimum_version = "4.1")
  skip_if_not_installed("partykit", minimum_version = "1.2")
  skip_if_not_installed("coin", minimum_version = "1.1")

  d <- make_sb_fixture(n = 48L, seed = 7931L)
  x <- d[c("gene_a", "gene_b", "gene_c")]

  survreg_adapter <- sb_additional_survival_learner_adapter(
    "survreg_weibull"
  )
  survreg_model <- survreg_adapter$fit(x, d$time, d$event)
  survreg_state <- survreg_adapter$extract_state(survreg_model)
  expect_false("engine_raw" %in% names(survreg_state))
  expect_true(all(is.finite(c(
    survreg_state$intercept,
    survreg_state$coefficient,
    survreg_state$scale
  ))))

  rpart_adapter <- sb_additional_survival_learner_adapter(
    "rpart_survival"
  )
  rpart_model <- rpart_adapter$fit(x, d$time, d$event)
  rpart_state <- rpart_adapter$extract_state(rpart_model)
  rpart_engine <- unserialize(rpart_state$engine_raw)
  expect_identical(
    intersect(c("model", "x", "y", "where"), names(rpart_engine)),
    character()
  )
  expect_identical(environment(rpart_engine$terms), baseenv())

  rsf_params <- sb_additional_survival_params("rsf_survival")
  rsf_params$num_trees <- 10L
  rsf_adapter <- sb_additional_survival_learner_adapter(
    "rsf_survival",
    rsf_params
  )
  rsf_state <- rsf_adapter$extract_state(rsf_adapter$fit(
    x,
    d$time,
    d$event
  ))
  expect_identical(
    rsf_state$training_data_retention,
    "event_time_metadata_required_by_randomForestSRC_forest"
  )
  rsf_engine <- unserialize(rsf_state$engine_raw)
  expect_identical(intersect(c("xvar", "yvar"), names(rsf_engine)), character())

  for (learner in c("ctree_survival", "cforest_survival")) {
    params <- sb_additional_survival_params(learner)
    if (identical(learner, "cforest_survival")) {
      params$num_trees <- 10L
    }
    adapter <- sb_additional_survival_learner_adapter(learner, params)
    state <- adapter$extract_state(adapter$fit(x, d$time, d$event))
    expect_identical(
      state$training_data_retention,
      "required_by_partykit_conditional_distribution_state"
    )
  }

  if (requireNamespace("BART", quietly = TRUE)) {
    params <- sb_additional_survival_params("bart_survival")
    params$time_grid_size <- 20L
    params$posterior_draws <- 30L
    params$burn_in <- 15L
    adapter <- sb_additional_survival_learner_adapter(
      "bart_survival",
      params
    )
    state <- adapter$extract_state(adapter$fit(x, d$time, d$event))
    engine <- unserialize(state$engine_raw)
    expect_true(all(names(engine) %in% c(
      "treedraws", "rm.const", "binaryOffset", "offset", "type", "times", "K"
    )))
    expect_false(any(c("yhat.train", "tx.train", "prob.train") %in%
      names(engine)))
  }
})

test_that("RSF fit and prediction restore the caller core option", {
  skip_if_not_installed("randomForestSRC", minimum_version = "3.3.0")

  d <- make_sb_fixture(n = 42L, seed = 7932L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  params <- sb_additional_survival_params("rsf_survival")
  params$num_trees <- 10L
  adapter <- sb_additional_survival_learner_adapter(
    "rsf_survival",
    params
  )
  withr::local_options(rf.cores = 7L, mc.cores = 6L)
  original_rfsrc <- randomForestSRC::rfsrc
  original_predict <- stats::predict
  fit_options <- integer()
  prediction_options <- integer()

  model <- testthat::with_mocked_bindings(
    rfsrc = function(...) {
      fit_options <<- c(
        rf.cores = getOption("rf.cores"),
        mc.cores = getOption("mc.cores")
      )
      original_rfsrc(...)
    },
    adapter$fit(x, d$time, d$event),
    .package = "randomForestSRC"
  )
  expect_identical(fit_options, c(rf.cores = 1L, mc.cores = 1L))
  expect_identical(
    c(rf.cores = getOption("rf.cores"), mc.cores = getOption("mc.cores")),
    c(rf.cores = 7L, mc.cores = 6L)
  )

  testthat::with_mocked_bindings(
    predict = function(...) {
      prediction_options <<- c(
        rf.cores = getOption("rf.cores"),
        mc.cores = getOption("mc.cores")
      )
      original_predict(...)
    },
    adapter$predict(model, x[1:3, , drop = FALSE]),
    .package = "stats"
  )
  expect_identical(
    prediction_options,
    c(rf.cores = 1L, mc.cores = 1L)
  )
  expect_identical(
    c(rf.cores = getOption("rf.cores"), mc.cores = getOption("mc.cores")),
    c(rf.cores = 7L, mc.cores = 6L)
  )
})

test_that("additional survival routes run CV and final prediction", {
  skip_if_not_installed("rpart", minimum_version = "4.1")
  skip_if_not_installed("randomForestSRC", minimum_version = "3.3.0")
  skip_if_not_installed("partykit", minimum_version = "1.2")
  skip_if_not_installed("coin", minimum_version = "1.1")
  skip_if_not_installed("BART", minimum_version = "2.9.0")

  d <- make_sb_fixture(n = 60L, seed = 7941L)
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
      if (identical(learner, "rsf_survival")) {
        params$num_trees <- 20L
      }
      if (identical(learner, "cforest_survival")) {
        params$num_trees <- 20L
      }
      if (identical(learner, "bart_survival")) {
        params$time_grid_size <- 20L
        params$posterior_draws <- 40L
        params$burn_in <- 20L
      }
      if (identical(learner, "aorsf_survival")) {
        params$num_trees <- 20L
      }
      params
    }
  )

  for (learner in sb_additional_survival_learner_ids()) {
    pipeline <- sb_pipeline(learner = learner)
    result <- sb_run_cv(task, pipeline, folds = 3L, seed = 7942L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      d[c(6L, 31L, 52L), c("patient_id", "gene_c", "gene_a", "gene_b")],
      "patient_id"
    )

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
