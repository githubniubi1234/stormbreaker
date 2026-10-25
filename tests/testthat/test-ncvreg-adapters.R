make_ncvreg_binary_fixture <- function(
  n = 140L, p = 6L, seed = 42101L, correlated = FALSE
) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    if (correlated && p > 1L) {
      rho <- 0.8
      sigma <- outer(seq_len(p), seq_len(p), function(i, j) {
        rho^abs(i - j)
      })
      x <- x %*% chol(sigma)
    }
    beta <- numeric(p)
    beta[seq_len(min(p, 4L))] <- c(1.0, -0.8, 0.6, -0.4)[
      seq_len(min(p, 4L))
    ]
    outcome <- stats::rbinom(n, 1L, stats::plogis(drop(x %*% beta)))
    names <- if (p == 1L) {
      "gene one"
    } else {
      c("gene one", "2gene", paste0("gene-", seq_len(p - 2L)))
    }
    data.frame(
      patient_id = sprintf("NB%03d", seq_len(n)),
      outcome = outcome,
      x,
      check.names = FALSE
    ) |>
      stats::setNames(c("patient_id", "outcome", names))
  })
}

make_ncvreg_survival_fixture <- function(
  n = 160L, p = 6L, seed = 42201L, tied = FALSE,
  censor_rate = 0.04, correlated = FALSE
) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    if (correlated && p > 1L) {
      rho <- 0.8
      sigma <- outer(seq_len(p), seq_len(p), function(i, j) {
        rho^abs(i - j)
      })
      x <- x %*% chol(sigma)
    }
    beta <- numeric(p)
    beta[seq_len(min(p, 4L))] <- c(0.9, -0.7, 0.5, -0.3)[
      seq_len(min(p, 4L))
    ]
    event_time <- stats::rexp(n, rate = 0.06 * exp(drop(x %*% beta)))
    censor_time <- stats::rexp(n, rate = censor_rate)
    time <- pmin(event_time, censor_time)
    if (tied) {
      time <- pmax(round(time, 1), 0.1)
    }
    names <- if (p == 1L) {
      "gene one"
    } else {
      c("gene one", "2gene", paste0("gene-", seq_len(p - 2L)))
    }
    data.frame(
      patient_id = sprintf("NS%03d", seq_len(n)),
      time = time,
      event = as.integer(event_time <= censor_time),
      x,
      check.names = FALSE
    ) |>
      stats::setNames(c("patient_id", "time", "event", names))
  })
}

ncvreg_native_fit <- function(x, endpoint, params) {
  matrix_x <- as.matrix(x)
  storage.mode(matrix_x) <- "double"
  lambda_min <- if (nrow(matrix_x) > ncol(matrix_x)) {
    params$lambda_min_ratio_n_gt_p
  } else {
    params$lambda_min_ratio_n_le_p
  }
  dfmax <- ncol(matrix_x) + if (identical(params$family, "binomial")) 1L else 0L
  if (identical(params$family, "binomial")) {
    return(ncvreg::ncvreg(
      X = matrix_x,
      y = as.numeric(endpoint$outcome),
      family = "binomial",
      penalty = params$penalty,
      gamma = params$gamma,
      alpha = params$alpha,
      lambda.min = lambda_min,
      nlambda = params$nlambda,
      eps = params$eps,
      max.iter = params$max_iter,
      convex = params$convex_diagnostics,
      dfmax = dfmax,
      penalty.factor = rep(1, ncol(matrix_x)),
      warn = TRUE,
      returnX = FALSE
    ))
  }
  ncvreg::ncvsurv(
    X = matrix_x,
    y = cbind(endpoint$time, endpoint$event),
    penalty = params$penalty,
    gamma = params$gamma,
    alpha = params$alpha,
    lambda.min = lambda_min,
    nlambda = params$nlambda,
    eps = params$eps,
    max.iter = params$max_iter,
    convex = params$convex_diagnostics,
    dfmax = dfmax,
    penalty.factor = rep(1, ncol(matrix_x)),
    warn = TRUE,
    returnX = FALSE
  )
}

test_that("four MCP and SCAD BIC learners are admitted without grids", {
  expected <- c(
    "ncvreg_mcp_binary_bic", "ncvreg_scad_binary_bic",
    "ncvreg_mcp_cox_bic", "ncvreg_scad_cox_bic"
  )
  expect_identical(sb_ncvreg_learner_ids(), expected)

  for (learner in expected) {
    params <- sb_ncvreg_learner_params(learner)
    expect_identical(params$alpha, 1)
    expect_identical(params$nlambda, 100L)
    expect_identical(params$selection, "analysis_fold_bic")
    expect_identical(params$warning_policy, "fail_fit")
    expect_identical(params$return_x, FALSE)
    if (grepl("_mcp_", learner, fixed = TRUE)) {
      expect_identical(params$penalty, "MCP")
      expect_identical(params$gamma, 3)
    } else {
      expect_identical(params$penalty, "SCAD")
      expect_identical(params$gamma, 3.7)
    }
  }

  binary_selectors <- expected_fixed_learner_selectors("binary")
  survival_selectors <- expected_fixed_learner_selectors("survival")
  rows <- sb_pipeline_registry()
  rows <- rows[rows$learner %in% expected, , drop = FALSE]
  expect_identical(
    nrow(rows),
    2L * length(binary_selectors) + 2L * length(survival_selectors)
  )
  expect_true(all(rows$tuning_semantics == "analysis_fold_bic"))
  expect_true(all(grepl("ncvreg", rows$optional_dependency, fixed = TRUE)))
  expect_setequal(
    unique(rows$selector[rows$task_type == "binary"]),
    binary_selectors
  )
  expect_setequal(
    unique(rows$selector[rows$task_type == "survival"]),
    survival_selectors
  )
})

test_that("binary MCP and SCAD reproduce the native BIC-selected model", {
  skip_if_not_installed("ncvreg", minimum_version = "3.16.0")

  data <- make_ncvreg_binary_fixture(correlated = TRUE)
  features <- setdiff(names(data), c("patient_id", "outcome"))
  x <- data[features]
  assessment <- x[c(2L, 31L, 88L), , drop = FALSE]

  for (learner in sb_ncvreg_binary_learner_ids()) {
    adapter <- sb_ncvreg_learner_adapter(learner)
    set.seed(42102L)
    caller_rng <- .Random.seed
    model <- adapter$fit(x, data$outcome)
    expect_identical(.Random.seed, caller_rng)

    native <- ncvreg_native_fit(
      x,
      list(outcome = data$outcome),
      adapter$params
    )
    bic <- stats::BIC(native)
    selected <- which.min(bic)
    native_coefficient <- drop(stats::coef(
      native,
      lambda = native$lambda[[selected]]
    ))
    native_probability <- drop(stats::predict(
      native,
      X = as.matrix(assessment),
      type = "response",
      lambda = native$lambda[[selected]]
    ))

    expect_equal(model$intercept, native_coefficient[[1L]], tolerance = 0)
    expect_equal(
      unname(model$coefficient),
      unname(native_coefficient[-1L]),
      tolerance = 0
    )
    expect_equal(
      adapter$predict(model, assessment),
      unname(native_probability),
      tolerance = 1e-14
    )
    expect_identical(model$effective$selected_index, as.integer(selected))
    expect_equal(
      model$effective$selected_lambda,
      native$lambda[[selected]],
      tolerance = 0
    )
    expect_false(any(c(
      "engine", "x", "X", "y", "outcome", "time", "event", "fail"
    ) %in% names(model)))

    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    replay <- adapter$predict_state(state, assessment)
    expect_equal(
      replay$probability,
      unname(native_probability),
      tolerance = 1e-14
    )
    expect_false(any(c(
      "engine", "x", "X", "y", "outcome", "time", "event", "fail"
    ) %in% names(state)))
  }
})

test_that("survival MCP and SCAD reproduce native risk and portable survival", {
  skip_if_not_installed("ncvreg", minimum_version = "3.16.0")

  data <- make_ncvreg_survival_fixture(tied = TRUE, correlated = TRUE)
  features <- setdiff(names(data), c("patient_id", "time", "event"))
  x <- data[features]
  assessment <- x[c(3L, 41L, 109L), , drop = FALSE]
  horizons <- as.numeric(stats::quantile(data$time, c(0.35, 0.65)))

  for (learner in sb_ncvreg_survival_learner_ids()) {
    adapter <- sb_ncvreg_learner_adapter(learner)
    set.seed(42202L)
    caller_rng <- .Random.seed
    model <- adapter$fit(x, data$time, data$event)
    expect_identical(.Random.seed, caller_rng)

    native <- ncvreg_native_fit(
      x,
      list(time = data$time, event = data$event),
      adapter$params
    )
    selected <- which.min(stats::BIC(native))
    native_risk <- drop(stats::predict(
      native,
      X = as.matrix(assessment),
      type = "link",
      lambda = native$lambda[[selected]]
    ))
    risk <- adapter$predict(model, assessment)
    early <- adapter$predict_survival(
      model, x, data$time, data$event, assessment, horizons[[1L]]
    )
    late <- adapter$predict_survival(
      model, x, data$time, data$event, assessment, horizons[[2L]]
    )

    expect_equal(risk, unname(native_risk), tolerance = 1e-14)
    expect_true(all(early$probability >= late$probability))
    expect_true(all(early$probability >= 0 & early$probability <= 1))
    expect_identical(early$ties, "breslow")
    expect_identical(
      early$source,
      paste0(learner, "_breslow_native_link")
    )

    state <- unserialize(serialize(adapter$extract_state(model, list(
      x = x,
      time = data$time,
      event = data$event,
      prediction_horizon = horizons[[1L]],
      time_unit = "month"
    )), NULL))
    replay <- adapter$predict_state(state, assessment)
    expect_equal(replay$risk, risk, tolerance = 0)
    expect_equal(
      replay$survival_probability,
      early$probability,
      tolerance = 0
    )
    expect_false(any(c(
      "engine", "x", "X", "y", "outcome", "time", "event", "fail"
    ) %in% names(state)))
  }
})

test_that("one-feature and heavy-censoring fits stay finite", {
  skip_if_not_installed("ncvreg", minimum_version = "3.16.0")

  binary <- make_ncvreg_binary_fixture(p = 1L, seed = 42301L)
  survival <- make_ncvreg_survival_fixture(
    p = 1L,
    seed = 42302L,
    censor_rate = 0.22
  )
  binary_x <- binary[setdiff(names(binary), c("patient_id", "outcome"))]
  survival_x <- survival[setdiff(
    names(survival),
    c("patient_id", "time", "event")
  )]

  for (learner in sb_ncvreg_binary_learner_ids()) {
    adapter <- sb_ncvreg_learner_adapter(learner)
    model <- adapter$fit(binary_x, binary$outcome)
    probability <- adapter$predict(model, binary_x[1:5, , drop = FALSE])
    expect_true(all(is.finite(probability)))
    expect_true(all(probability >= 0 & probability <= 1))
  }
  for (learner in sb_ncvreg_survival_learner_ids()) {
    adapter <- sb_ncvreg_learner_adapter(learner)
    model <- adapter$fit(survival_x, survival$time, survival$event)
    risk <- adapter$predict(model, survival_x[1:5, , drop = FALSE])
    expect_true(all(is.finite(risk)))
  }
})

test_that("native path warnings are preserved as typed fit failures", {
  skip_if_not_installed("ncvreg", minimum_version = "3.16.0")

  cases <- list(
    binary = make_ncvreg_binary_fixture(n = 80L, p = 120L, seed = 41001L),
    survival = make_ncvreg_survival_fixture(
      n = 80L,
      p = 120L,
      seed = 42001L
    )
  )
  for (task in names(cases)) {
    data <- cases[[task]]
    endpoints <- if (identical(task, "binary")) {
      c("patient_id", "outcome")
    } else {
      c("patient_id", "time", "event")
    }
    x <- data[setdiff(names(data), endpoints)]
    learners <- if (identical(task, "binary")) {
      sb_ncvreg_binary_learner_ids()
    } else {
      sb_ncvreg_survival_learner_ids()
    }
    endpoint <- if (identical(task, "binary")) {
      list(outcome = data$outcome)
    } else {
      list(time = data$time, event = data$event)
    }
    for (learner in learners) {
      params <- sb_ncvreg_learner_params(learner)
      native_warnings <- character()
      native <- withCallingHandlers(
        ncvreg_native_fit(x, endpoint, params),
        warning = function(warning) {
          native_warnings <<- c(native_warnings, conditionMessage(warning))
          invokeRestart("muffleWarning")
        }
      )
      expect_s3_class(native, "ncvreg")

      adapter <- sb_ncvreg_learner_adapter(learner)
      fit_call <- if (identical(task, "binary")) {
        quote(adapter$fit(x, data$outcome))
      } else {
        quote(adapter$fit(x, data$time, data$event))
      }
      if (length(native_warnings)) {
        expect_error(
          eval(fit_call),
          class = "stormbreaker_ncvreg_fit_warning"
        )
      } else {
        expect_s3_class(eval(fit_call), "stormbreaker_ncvreg_model")
      }
    }
  }
})

test_that("MCP and SCAD execute through CV and final prediction", {
  skip_if_not_installed("ncvreg", minimum_version = "3.16.0")

  binary <- make_ncvreg_binary_fixture(seed = 42401L)
  binary_features <- setdiff(names(binary), c("patient_id", "outcome"))
  binary_task <- sb_binary_task(
    binary,
    "patient_id",
    "outcome",
    binary_features
  )
  survival <- make_ncvreg_survival_fixture(seed = 42402L)
  survival_features <- setdiff(
    names(survival),
    c("patient_id", "time", "event")
  )
  survival_task <- sb_surv_task(
    survival,
    "patient_id",
    "time",
    "event",
    survival_features,
    "month",
    prediction_horizon = as.numeric(stats::quantile(survival$time, 0.5))
  )

  for (learner in sb_ncvreg_binary_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 4L
    )
    comparison <- sb_run_cv(binary_task, pipeline, folds = 3L, seed = 42403L)
    fit <- sb_fit_final(binary_task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      binary[1:6, c("patient_id", rev(binary_features)), drop = FALSE],
      "patient_id"
    )
    expect_identical(comparison$status, "success")
    expect_true(comparison$complete_oof)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_true(all(prediction$predictions$probability >= 0 &
      prediction$predictions$probability <= 1))
  }

  for (learner in sb_ncvreg_survival_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_cox_top_k",
      learner = learner,
      selector_n_features = 4L
    )
    comparison <- sb_run_cv(
      survival_task,
      pipeline,
      folds = 3L,
      seed = 42404L
    )
    fit <- sb_fit_final(survival_task, pipeline)
    prediction <- sb_predict_final(
      unserialize(serialize(fit, NULL)),
      survival[1:6, c("patient_id", rev(survival_features)), drop = FALSE],
      "patient_id"
    )
    expect_identical(comparison$status, "success")
    expect_true(comparison$complete_oof)
    expect_true(comparison$complete_survival_probability_oof)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_true(all(is.finite(prediction$predictions$risk)))
    expect_true(all(prediction$predictions$survival_probability >= 0 &
      prediction$predictions$survival_probability <= 1))
  }
})
