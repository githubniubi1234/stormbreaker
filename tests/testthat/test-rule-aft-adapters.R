make_rule_aft_fixture <- function(n = 72L, p = 12L, seed = 9101L) {
  withr::with_seed(seed, {
    expression <- matrix(stats::rnorm(n * p), nrow = n)
    colnames(expression) <- paste0("gene_", seq_len(p))
    nonlinear <- 1.3 * expression[, 1L] * expression[, 2L] +
      0.8 * expression[, 3L] - 0.5 * expression[, 4L]^2
    outcome <- as.integer(nonlinear > stats::median(nonlinear))
    log_time <- 2 + 0.6 * expression[, 1L] -
      0.5 * expression[, 2L] +
      0.7 * (expression[, 3L] > 0) + stats::rnorm(n, sd = 0.35)
    data.frame(
      patient_id = paste0("RULE_AFT_PRIVATE_", seq_len(n)),
      time = exp(log_time),
      event = rep(c(1L, 1L, 0L, 1L, 0L), length.out = n),
      outcome = outcome,
      as.data.frame(expression, check.names = FALSE),
      check.names = FALSE
    )
  })
}

rule_aft_features <- function(data) {
  setdiff(names(data), c("patient_id", "time", "event", "outcome"))
}

rule_aft_token_params <- function(learner) {
  params <- sb_rule_aft_params(learner)
  if (identical(learner, "bart_aft_survival")) {
    params$num_trees <- 15L
    params$posterior_draws <- 40L
    params$burn_in <- 20L
    params$time_grid_size <- 20L
  } else if (learner %in% sb_pre_rule_ensemble_learner_ids()) {
    params$num_trees <- 80L
    params$internal_folds <- 3L
  } else {
    params$num_lambda <- 8L
    params$internal_folds <- 3L
  }
  params
}

test_that("four rule/AFT learners expose only meaningful routes", {
  ids <- sb_rule_aft_learner_ids()
  expect_identical(ids, c(
    "bart_aft_survival",
    "pre_binary_rule_ensemble",
    "pre_survival_rule_ensemble",
    "penaft_gehan_elastic_net"
  ))
  registry <- sb_pipeline_registry()
  rows <- registry[registry$learner %in% ids, , drop = FALSE]
  expect_setequal(
    rows$selector[rows$learner == "bart_aft_survival"],
    expected_fixed_learner_selectors("survival")
  )
  for (learner in setdiff(ids, "bart_aft_survival")) {
    expect_setequal(
      rows$selector[rows$learner == learner],
      c("none", "fixed_features")
    )
  }
  pre_rows <- rows$learner %in% sb_pre_rule_ensemble_learner_ids()
  expect_identical(
    unique(rows$tuning_semantics[pre_rows]),
    "analysis_fold_native_coefficient_cv_after_rule_induction"
  )
  expect_true(all(grepl(
    "complete analysis fold",
    rows$parameter_semantics[pre_rows],
    fixed = TRUE
  )))
  expect_identical(
    unique(rows$optional_dependency[
      rows$learner == "pre_survival_rule_ensemble"
    ]),
    "pre + mboost"
  )
  expect_error(
    sb_pipeline(
      selector = "univ_cox_top_k",
      learner = "pre_survival_rule_ensemble",
      selector_n_features = 5L
    ),
    class = "stormbreaker_unadmitted_pipeline"
  )
})

test_that("native rule/AFT fits are high-dimensional, portable, and local-RNG", {
  skip_if_not_installed("BART", minimum_version = "2.9.0")
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("pre", minimum_version = "1.1.0")
  skip_if_not_installed("penAFT", minimum_version = "0.3.2")
  data <- make_rule_aft_fixture(n = 54L, p = 70L, seed = 9102L)
  features <- rule_aft_features(data)
  x <- data[features]
  rownames(x) <- data$patient_id
  assessment <- x[c(2L, 19L, 47L), , drop = FALSE]
  horizon <- as.numeric(stats::median(data$time))

  for (learner in sb_rule_aft_learner_ids()) {
    adapter <- sb_rule_aft_learner_adapter(
      learner,
      rule_aft_token_params(learner)
    )
    set.seed(9200L)
    caller_rng <- .Random.seed
    model <- if (identical(learner, "pre_binary_rule_ensemble")) {
      adapter$fit(x, data$outcome)
    } else {
      adapter$fit(x, data$time, data$event)
    }
    expect_identical(.Random.seed, caller_rng)
    prediction <- adapter$predict(model, assessment)
    expect_length(prediction, nrow(assessment))
    expect_equal(sum(is.finite(prediction)), nrow(assessment))
    if (identical(learner, "pre_binary_rule_ensemble")) {
      native <- stats::predict(
        model$engine,
        newdata = assessment,
        type = "response",
        penalty.par.val = model$params$penalty_rule
      )
      expect_equal(prediction, unname(as.numeric(native)), tolerance = 0)
    }
    if (identical(learner, "penaft_gehan_elastic_net")) {
      manual <- -unname(
        model$intercept +
          as.vector(as.matrix(assessment) %*% model$coefficient)
      )
      expect_equal(prediction, manual, tolerance = 0)
    }

    context <- if (learner %in% sb_rule_aft_survival_learner_ids()) {
      list(
        x = x,
        time = data$time,
        event = data$event,
        prediction_horizon = horizon,
        time_unit = "month"
      )
    } else {
      list()
    }
    state <- unserialize(serialize(
      adapter$extract_state(model, context),
      NULL
    ))
    replay <- adapter$predict_state(state, assessment)
    replay_value <- if (identical(learner, "pre_binary_rule_ensemble")) {
      replay$probability
    } else {
      replay$risk
    }
    expect_equal(replay_value, prediction, tolerance = 1e-10)
    if (learner %in% sb_rule_aft_survival_learner_ids()) {
      expect_equal(
        sum(replay$survival_probability >= 0 &
          replay$survival_probability <= 1),
        nrow(assessment)
      )
    }
    expect_length(
      grepRaw(
        "RULE_AFT_PRIVATE_",
        serialize(state, NULL),
        fixed = TRUE,
        all = TRUE
      ),
      0L
    )
  }
})

test_that("AFT-BART and penAFT horizon formulas match their fitted states", {
  skip_if_not_installed("BART", minimum_version = "2.9.0")
  skip_if_not_installed("penAFT", minimum_version = "0.3.2")
  data <- make_rule_aft_fixture(n = 60L, p = 18L, seed = 9103L)
  features <- rule_aft_features(data)
  x <- data[features]
  assessment <- x[c(3L, 28L, 55L), , drop = FALSE]
  horizon <- as.numeric(stats::median(data$time))

  bart <- sb_rule_aft_learner_adapter(
    "bart_aft_survival",
    rule_aft_token_params("bart_aft_survival")
  )
  bart_model <- bart$fit(x, data$time, data$event)
  draws <- sb_bart_aft_log_time_draws(bart_model, assessment)
  manual_bart <- colMeans(stats::pnorm(
    (log(horizon) - draws) / bart_model$engine$sigma,
    lower.tail = FALSE
  ))
  bart_probability <- bart$predict_survival(
    bart_model,
    x,
    data$time,
    data$event,
    assessment,
    horizon
  )
  expect_equal(
    bart_probability$probability,
    unname(manual_bart),
    tolerance = 0
  )

  penaft <- sb_rule_aft_learner_adapter(
    "penaft_gehan_elastic_net",
    rule_aft_token_params("penaft_gehan_elastic_net")
  )
  penaft_model <- penaft$fit(x, data$time, data$event)
  link <- sb_penaft_log_time_link(penaft_model, x)
  curve <- sb_fit_mboost_gehan_residual_survival(
    data$time,
    data$event,
    link
  )
  manual_penaft <- sb_predict_mboost_gehan_residual_survival(
    curve,
    sb_penaft_log_time_link(penaft_model, assessment),
    horizon
  )
  penaft_probability <- penaft$predict_survival(
    penaft_model,
    x,
    data$time,
    data$event,
    assessment,
    horizon
  )
  expect_equal(
    penaft_probability$probability,
    manual_penaft,
    tolerance = 0
  )
})

test_that("rule/AFT learners run through CV and outcome-blind final prediction", {
  skip_on_cran()
  skip_if_not_installed("BART", minimum_version = "2.9.0")
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("pre", minimum_version = "1.1.0")
  skip_if_not_installed("penAFT", minimum_version = "0.3.2")
  data <- make_rule_aft_fixture(n = 75L, p = 6L, seed = 9104L)
  features <- rule_aft_features(data)
  horizon <- as.numeric(stats::median(data$time))

  for (learner in sb_rule_aft_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "fixed_features",
      learner = learner,
      selector_features = features,
      feature_set_id = "rule_aft_test"
    )
    task <- if (identical(learner, "pre_binary_rule_ensemble")) {
      sb_binary_task(
        data,
        "patient_id",
        "outcome",
        features
      )
    } else {
      sb_surv_task(
        data,
        "patient_id",
        "time",
        "event",
        features,
        "month",
        prediction_horizon = horizon
      )
    }
    cv <- sb_run_cv(task, pipeline, folds = 3L, seed = 9301L)
    expect_identical(unique(cv$runs$status), "success")
    expect_identical(unique(cv$predictions$status), "success")
    if (identical(learner, "pre_survival_rule_ensemble")) {
      expect_true(all(c("pre_version", "mboost_version") %in%
        names(cv$provenance)))
      expect_false(anyNA(unlist(cv$provenance[
        c("pre_version", "mboost_version")
      ])))
    }

    fit <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success")
    if (identical(learner, "pre_survival_rule_ensemble")) {
      expect_true(all(c("pre_version", "mboost_version") %in%
        names(fit$provenance)))
      expect_false(anyNA(unlist(fit$provenance[
        c("pre_version", "mboost_version")
      ])))
    }
    new_data <- data[1:5, , drop = FALSE]
    altered <- new_data
    altered$outcome <- 1L - altered$outcome
    altered$time <- altered$time + 1000
    altered$event <- 1L - altered$event
    prediction <- sb_predict_final(fit, new_data, "patient_id")
    changed <- sb_predict_final(fit, altered, "patient_id")
    expect_identical(prediction$status, "success")
    value <- if (identical(learner, "pre_binary_rule_ensemble")) {
      "probability"
    } else {
      "risk"
    }
    expect_equal(
      prediction$predictions[[value]],
      changed$predictions[[value]],
      tolerance = 0
    )
    if (learner %in% sb_rule_aft_survival_learner_ids()) {
      expect_equal(
        prediction$predictions$survival_probability,
        changed$predictions$survival_probability,
        tolerance = 0
      )
    }
  }
})

test_that("rule/AFT dependencies include the survival PRE engine", {
  dependency <- sb_rule_aft_dependency("pre_survival_rule_ensemble")
  expect_identical(dependency$packages, c("pre", "mboost"))

  local_mocked_bindings(
    sb_rule_aft_dependency = function(learner) {
      list(
        package = "pre",
        packages = c("pre", "mboost"),
        installed = c(pre = TRUE, mboost = FALSE),
        available = FALSE
      )
    }
  )
  expect_error(
    sb_require_rule_aft_dependency("pre_survival_rule_ensemble"),
    "mboost",
    class = "stormbreaker_missing_rule_aft_dependency"
  )
})
