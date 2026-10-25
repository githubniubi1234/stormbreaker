test_that("saved BART and cforest states evaluate a native curve once", {
  skip_if_not_installed("BART")
  skip_if_not_installed("partykit")
  skip_if_not_installed("coin")
  d <- sb_example_survival(60L, 992L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  query <- x[c(3L, 12L, 27L, 59L), , drop = FALSE]
  bart <- sb_bart_survival_curves
  party <- sb_partykit_survival_distributions
  calls <- 0L
  local_mocked_bindings(
    sb_bart_survival_curves = function(...) { calls <<- calls + 1L; bart(...) },
    sb_partykit_survival_distributions = function(...) { calls <<- calls + 1L; party(...) })
  withr::local_seed(93L)
  for (id in c("bart_survival", "cforest_survival")) {
    a <- sb_additional_survival_learner_adapter(id)
    model <- a$fit(x, d$time, d$event)
    for (h in c(0.001, as.numeric(median(d$time)), as.numeric(max(d$time)))) {
      state <- unserialize(serialize(a$extract_state(model, list(x = x, time = d$time,
        event = d$event, prediction_horizon = h, time_unit = "month")), NULL))
      expected <- list(risk = a$predict(model, query), prediction_horizon = h,
        time_unit = "month", survival_probability = a$predict_survival(
          model, x, d$time, d$event, query, h)$probability)
      calls <- 0L; rng <- .Random.seed
      actual <- a$predict_state(state, query)
      expect_identical(actual, expected, info = id)
      expect_identical(calls, 1L, info = id)
      expect_identical(.Random.seed, rng, info = id)
    }
  }
})

test_that("curve-capable adapters expose exact CV details only for pure prediction", {
  skip_if_not_installed("BART")
  skip_if_not_installed("partykit")
  skip_if_not_installed("coin")
  d <- sb_example_survival(60L, 992L)
  x <- d[c("gene_a", "gene_b")]
  for (id in c("bart_survival", "cforest_survival")) {
    a <- sb_additional_survival_learner_adapter(id)
    expect_identical(is.function(a$predict_details), TRUE, info = id)
    if (!is.function(a$predict_details)) next
    model <- a$fit(x, d$time, d$event)
    details <- a$predict_details(model, x[1:3, , drop = FALSE])
    expect_identical(details$risk, a$predict(model, x[1:3, , drop = FALSE]))
    for (h in c(0.001, as.numeric(median(d$time)))) {
      expected <- a$predict_survival(model, x, d$time, d$event, x[1:3, , drop = FALSE], h)
      actual <- a$predict_survival(model, x, d$time, d$event, x[1:3, , drop = FALSE], h,
        prediction_details = details)
      expect_identical(actual, expected, info = id)
    }
  }
  for (id in c("rsf_survival", "ctree_survival", "aorsf_survival", "rpart_survival",
    "survreg_weibull")) {
    expect_null(sb_additional_survival_learner_adapter(id)$predict_details)
  }
})

test_that("shared native curves preserve CV and horizon-only failure boundaries", {
  skip_if_not_installed("BART")
  skip_if_not_installed("partykit")
  skip_if_not_installed("coin")
  d <- sb_example_survival(72L, 992L)
  d$time[1L] <- max(d$time) + 1
  factory <- sb_additional_survival_learner_adapter
  bart <- sb_bart_survival_curves
  party <- sb_partykit_survival_distributions
  calls <- 0L
  local_mocked_bindings(
    sb_bart_survival_curves = function(...) { calls <<- calls + 1L; bart(...) },
    sb_partykit_survival_distributions = function(...) { calls <<- calls + 1L; party(...) })
  withr::local_seed(94L)
  before <- .Random.seed
  for (id in c("bart_survival", "cforest_survival")) {
    for (h in c(as.numeric(median(d$time)), as.numeric(max(d$time)))) {
      task <- sb_surv_task(d, "patient_id", "time", "event",
        c("gene_a", "gene_b", "gene_c"), "month", prediction_horizon = h)
      p <- sb_pipeline("univ_cox_top_k", id, selector_n_features = 2L)
      calls <- 0L
      actual <- sb_run_cv(task, p, folds = 3L, seed = 82L)
      expect_identical(calls, 6L, info = id)
      calls <- 0L
      expected <- with_mocked_bindings(sb_run_cv(task, p, folds = 3L, seed = 82L),
        sb_additional_survival_learner_adapter = function(...) {
          a <- factory(...); a$predict_details <- NULL; a
        })
      expect_identical(calls, if (h < max(d$time)) 9L else 8L, info = id)
      expect_identical(actual, expected, info = id)
      expect_identical(actual$complete_oof, TRUE, info = id)
      expect_identical(actual$complete_survival_probability_oof, h < max(d$time), info = id)
      expect_identical(.Random.seed, before, info = id)
    }
  }
})

test_that("public final predictions replay shared distributions without outcomes", {
  skip_if_not_installed("BART")
  skip_if_not_installed("partykit")
  skip_if_not_installed("coin")
  d <- sb_example_survival(60L, 992L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
    prediction_horizon = as.numeric(median(d$time)))
  for (id in c("bart_survival", "ctree_survival", "cforest_survival")) {
    fit <- sb_fit_final(task, sb_pipeline("fixed_features", id,
      selector_features = features[1:2], feature_set_id = "two"))
    expect_identical(fit$status, "success", info = id)
    blind <- d[c("patient_id", rev(features))]
    actual <- sb_predict_final(unserialize(serialize(fit, NULL)), blind, "patient_id")
    expect_identical(actual$status, "success", info = id)
    expect_identical(actual, sb_predict_final(fit, blind, "patient_id"), info = id)
    expect_identical(actual$predictions$patient_id, d$patient_id, info = id)
  }
})

test_that("a warned native distribution is not reused during state prediction", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("coin")
  d <- sb_example_survival(60L, 992L)
  x <- d[c("gene_a", "gene_b")]
  a <- sb_additional_survival_learner_adapter("cforest_survival")
  model <- a$fit(x, d$time, d$event)
  state <- a$extract_state(model, list(x = x, time = d$time, event = d$event,
    prediction_horizon = 5, time_unit = "month"))
  expected <- a$predict_state(state, x[1:3, , drop = FALSE])
  original <- sb_partykit_survival_distributions
  calls <- 0L
  local_mocked_bindings(sb_partykit_survival_distributions = function(...) {
    calls <<- calls + 1L
    warning("native prediction warning", call. = FALSE)
    original(...)
  })
  warnings <- character()
  actual <- withCallingHandlers(a$predict_state(state, x[1:3, , drop = FALSE]),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  expect_identical(actual, expected)
  expect_identical(calls, 2L)
  expect_identical(warnings, rep("native prediction warning", 2L))
})

test_that("risk-only state calls once and native errors remain failures", {
  skip_if_not_installed("partykit")
  skip_if_not_installed("coin")
  d <- sb_example_survival(60L, 992L)
  x <- d[c("gene_a", "gene_b")]
  a <- sb_additional_survival_learner_adapter("ctree_survival")
  model <- a$fit(x, d$time, d$event)
  state <- a$extract_state(model)
  original <- sb_partykit_survival_distributions
  calls <- 0L
  local_mocked_bindings(sb_partykit_survival_distributions = function(...) {
    calls <<- calls + 1L; original(...)
  })
  expected <- a$predict(model, x)
  calls <- 0L
  withr::local_seed(95L)
  rm(".Random.seed", envir = .GlobalEnv)
  expect_identical(a$predict_state(state, x), list(risk = expected))
  expect_identical(calls, 1L)
  expect_identical(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), FALSE)
  local_mocked_bindings(sb_partykit_survival_distributions = function(...) {
    sb_abort("native distribution failed", "stormbreaker_test_prediction_failure")
  })
  error <- capture_sb_error(a$predict_state(state, x))
  expect_s3_class(error, "stormbreaker_test_prediction_failure")
})
