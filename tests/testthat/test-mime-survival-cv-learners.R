make_mime_cv_learner_fixture <- function(n = 90L, p = 12L,
                                         seed = 8901L) {
  set.seed(seed)
  x <- as.data.frame(matrix(stats::rnorm(n * p), nrow = n, ncol = p))
  names(x) <- sprintf("gene_%02d", seq_len(p))
  linear_predictor <- x[[1L]] - x[[2L]] + 0.6 * x[[3L]]
  event_time <- stats::rexp(n, rate = exp(linear_predictor) / 15)
  censor_time <- stats::rexp(n, rate = 1 / 25)
  list(
    x = x,
    time = pmax(pmin(event_time, censor_time), 0.05),
    event = as.integer(event_time <= censor_time)
  )
}

test_that("reference CV survival learners fit through their native engines", {
  skip_if_not_installed("CoxBoost")
  skip_if_not_installed("plsRcox")
  skip_if_not_installed("superpc")
  skip_if_not_installed("gbm")

  data <- make_mime_cv_learner_fixture()
  adapters <- list(
    coxboost_cv = sb_coxboost_cv_learner_adapter(
      inner_folds = 3L,
      seed = 8902L,
      maximum_steps = 80L
    ),
    plsrcox_cv = sb_plsrcox_cv_learner_adapter(
      inner_folds = 3L,
      seed = 8903L,
      maximum_components = 10L
    ),
    superpc_cv = sb_superpc_cv_learner_adapter(
      inner_folds = 3L,
      seed = 8904L
    ),
    gbm_cox_cv = sb_gbm_cox_cv_learner_adapter(
      inner_folds = 3L,
      seed = 8905L,
      maximum_trees = 60L
    )
  )

  models <- lapply(adapters, function(adapter) {
    model <- adapter$fit(data$x, data$time, data$event)
    observed <- adapter$predict(model, data$x[1:5, , drop = FALSE])
    state <- adapter$extract_state(model, list())
    restored <- adapter$predict_state(state, data$x[1:5, , drop = FALSE])
    expect_identical(adapter$id, state$learner)
    expect_equal(restored$risk, observed, tolerance = 0)
    expect_true(all(is.finite(observed)))
    model
  })

  expect_true(is.finite(models$coxboost_cv$audit$effective_penalty))
  expect_gte(models$coxboost_cv$audit$optimal_step, 1L)
  expect_lte(models$coxboost_cv$audit$optimal_step, 80L)
  expect_setequal(models$coxboost_cv$audit$fold_id, 1:3)
  expect_gte(models$plsrcox_cv$audit$selected_components, 1L)
  expect_lte(models$plsrcox_cv$audit$selected_components, 10L)
  expect_length(models$plsrcox_cv$audit$folds, 3L)
  expect_true(is.finite(models$superpc_cv$audit$selected_threshold))
  expect_identical(models$superpc_cv$audit$requested_folds, 3L)
  expect_gte(models$gbm_cox_cv$n_trees, 1L)
  expect_lte(models$gbm_cox_cv$n_trees, 60L)
  expect_length(models$gbm_cox_cv$tuning$cv_error, 60L)
})

test_that("reference CV learner controls reconstruct pipeline contracts", {
  for (learner in sb_analysis_tuned_survival_learner_ids()) {
    pipeline <- sb_pipeline(
      learner = learner,
      learner_tuning_folds = 3L,
      learner_seed = 8911L,
      gbm_cv_max_trees = if (identical(learner, "gbm_cox_cv")) {
        60L
      } else {
        NULL
      }
    )
    expect_identical(sb_lint_pipeline(pipeline)$status, "admitted")
    expect_identical(
      sb_validate_pipeline_contract(pipeline, "survival"),
      pipeline
    )
  }
})

test_that("PLS-Cox CV uses the sole feasible component for one predictor", {
  skip_if_not_installed("plsRcox")
  data <- make_mime_cv_learner_fixture(p = 3L)
  data$x <- data$x[1L]
  adapter <- sb_plsrcox_cv_learner_adapter(
    inner_folds = 3L, seed = 8912L, maximum_components = 10L
  )
  model <- adapter$fit(data$x, data$time, data$event)
  observed <- adapter$predict(model, data$x[1:5, , drop = FALSE])
  native <- plsRcox::plsRcox(
    data$x,
    time = data$time,
    event = data$event,
    nt = 1L,
    scaleX = TRUE,
    scaleY = FALSE,
    sparseStop = FALSE,
    allres = FALSE,
    verbose = FALSE
  )
  expected <- as.numeric(stats::predict(
    native, newdata = data$x[1:5, , drop = FALSE], type = "lp"
  ))
  expect_identical(model$audit$selected_components, 1L)
  expect_identical(
    model$audit$component_selection,
    "single_predictor_has_one_feasible_component"
  )
  expect_length(model$audit$folds, 0L)
  expect_true(all(is.finite(observed)))
  expect_equal(
    observed - mean(observed), expected - mean(expected), tolerance = 1e-12
  )
})

test_that("SuperPC CV reports its two-predictor tuning boundary", {
  skip_if_not_installed("superpc")
  data <- make_mime_cv_learner_fixture(p = 3L)
  data$x <- data$x[1L]
  error <- capture_sb_error(
    sb_superpc_cv_learner_adapter(inner_folds = 3L, seed = 8913L)$fit(
      data$x, data$time, data$event
    )
  )
  expect_s3_class(error, "stormbreaker_superpc_cv_too_few_features")
  expect_match(
    conditionMessage(error),
    "requires at least 2 predictors but received 1",
    fixed = TRUE
  )
})

test_that("PLS-Cox CV retains warnings beyond its selected component count", {
  skip_if_not_installed("plsRcox")
  data <- make_mime_cv_learner_fixture(n = 30L, p = 3L)
  message <- "Loglik converged before variable  3 ; coefficient may be infinite. "
  cv_calls <- 0L
  final_components <- NULL
  native_cv <- function(...) {
    cv_calls <<- cv_calls + 1L
    warning(message)
    warning(message)
    list(lambda.min5 = 2L, folds = split(seq_len(30L), rep(1:3, each = 10L)),
         cv.error5 = c(0.5, 0.6, 0.7, 0.65), cv.se5 = rep(0.01, 4L),
         completed.cv5 = matrix(TRUE, nrow = 4L, ncol = 3L))
  }
  native_fit <- function(Xplan, nt, ...) {
    final_components <<- nt
    list(Coeffs = matrix(0.1, nrow = ncol(Xplan), ncol = 1L,
                         dimnames = list(names(Xplan), NULL)), computed_nt = nt)
  }
  model <- testthat::with_mocked_bindings(
    sb_plsrcox_cv_learner_adapter(inner_folds = 3L, maximum_components = 3L)$fit(
      data$x, data$time, data$event
    ),
    cv.plsRcox = native_cv, plsRcox = native_fit, .package = "plsRcox"
  )
  expect_identical(cv_calls, 1L)
  expect_identical(final_components, 2L)
  expect_identical(model$audit$selected_components, 2L)
  expect_identical(model$audit$component_search_warnings, rep(message, 2L))
  expect_identical(attr(model, "stormbreaker_fit_warnings"), rep(message, 2L))
})

test_that("PLS-Cox CV still rejects selected and unclassified warnings", {
  skip_if_not_installed("plsRcox")
  data <- make_mime_cv_learner_fixture(n = 30L, p = 3L)
  messages <- c(
    "Loglik converged before variable  2 ; coefficient may be infinite. ",
    "Loglik converged before variable  2,3 ; coefficient may be infinite. ",
    "Loglik converged before variable  0 ; coefficient may be infinite. ",
    "Loglik converged before variable  4 ; coefficient may be infinite. ",
    "Ran out of iterations and did not converge",
    "INJECTED_UNKNOWN_PLSRCOX_WARNING"
  )
  for (message in messages) {
    error <- capture_sb_error(testthat::with_mocked_bindings(
      sb_plsrcox_cv_learner_adapter(inner_folds = 3L, maximum_components = 3L)$fit(
        data$x, data$time, data$event
      ),
      cv.plsRcox = function(...) {
        warning(message)
        list(lambda.min5 = 2L)
      },
      .package = "plsRcox"
    ))
    expect_s3_class(error, "stormbreaker_plsrcox_cv_warning")
    expect_match(conditionMessage(error), message, fixed = TRUE)
    expect_s3_class(error$parent, "warning")
  }
})

test_that("PLS-Cox CV does not accept incomplete selected folds or refit warnings", {
  skip_if_not_installed("plsRcox")
  data <- make_mime_cv_learner_fixture(n = 30L, p = 3L)
  message <- "Loglik converged before variable  3 ; coefficient may be infinite. "
  for (complete in c(FALSE, TRUE)) {
    error <- capture_sb_error(testthat::with_mocked_bindings(
      sb_plsrcox_cv_learner_adapter(inner_folds = 3L, maximum_components = 3L)$fit(
        data$x, data$time, data$event
      ),
      cv.plsRcox = function(...) {
        warning(message)
        completed <- matrix(TRUE, nrow = 4L, ncol = 3L)
        completed[3L, 1L] <- complete
        list(lambda.min5 = 2L, completed.cv5 = completed)
      },
      plsRcox = function(...) warning(message),
      .package = "plsRcox"
    ))
    expect_s3_class(error, if (complete) {
      "stormbreaker_plsrcox_cv_warning"
    } else {
      "stormbreaker_invalid_plsrcox_cv"
    })
    expect_match(conditionMessage(error), if (complete) {
      "final fit raised a warning"
    } else {
      "did not complete every fold"
    }, fixed = TRUE)
  }
})

test_that("PLS-Cox CV fits the original 20-feature case without changing selection", {
  skip_if_not_installed("plsRcox", minimum_version = "1.8.0")
  data <- sb_example_survival(n = 200L, seed = 20260905L)
  noise <- withr::with_seed(20260906L, matrix(rnorm(200L * 17L), nrow = 200L))
  for (j in seq_len(ncol(noise))) {
    data[[sprintf("gene_%02d", j + 3L)]] <- noise[, j]
  }
  features <- grep("^gene_", names(data), value = TRUE)
  training <- data[1:160, ]
  task <- sb_surv_task(training, "patient_id", "time", "event", features,
                       "month", prediction_horizon = 6)
  fit <- sb_fit_final(task, sb_pipeline(learner = "plsrcox_cv"))
  expect_identical(fit$status, "success")
  expect_identical(fit$learner_state$audit$selected_components, 4L)
  native_warnings <- character()
  native_cv <- withCallingHandlers(withr::with_seed(1L, plsRcox::cv.plsRcox(
    list(x = training[features], time = training$time, status = training$event),
    nfold = 5L, nt = 10L, plot.it = FALSE, folddetails = TRUE, verbose = FALSE,
    control = survival::coxph.control(eps = 1e-11, iter.max = 50L,
      toler.inf = survival::coxph.control()$toler.inf)
  )), warning = function(w) {
    native_warnings <<- c(native_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(as.integer(native_cv$lambda.min5), 4L)
  expect_identical(fit$learner_state$audit$component_search_warnings, native_warnings)
  expect_identical(fit$tuning$warnings, paste0("Learner refit: ", native_warnings))
  expect_match(fit$tuning$warnings, "coefficient may be infinite", fixed = TRUE)
  assessment <- data[161:200, c("patient_id", features)]
  prediction <- sb_predict_final(fit, assessment, "patient_id")
  expect_identical(prediction$status, "success")
  expect_identical(prediction$predictions$patient_id, assessment$patient_id)
  expect_identical(all(is.finite(prediction$predictions$risk)), TRUE)
  expect_identical(all(prediction$predictions$survival_probability >= 0 &
                         prediction$predictions$survival_probability <= 1), TRUE)
  native <- plsRcox::plsRcox(training[features], time = training$time,
                            event = training$event, nt = 4L, scaleX = TRUE,
                            scaleY = FALSE, sparseStop = FALSE,
                            allres = FALSE, verbose = FALSE)
  native_risk <- as.numeric(stats::predict(native,
    newdata = assessment[features], type = "lp"
  ))
  observed <- prediction$predictions$risk
  expect_equal(observed - mean(observed), native_risk - mean(native_risk),
               tolerance = 1e-10)
})
