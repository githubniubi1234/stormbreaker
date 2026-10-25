test_that("inner CV selects lambda reproducibly from its own patients", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 80, seed = 20260813)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  first <- sb_select_glmnet_cox_lambda(
    x, d$time, d$event, inner_folds = 5L, seed = 1701L
  )
  repeated <- sb_select_glmnet_cox_lambda(
    x, d$time, d$event, inner_folds = 5L, seed = 1701L
  )

  expect_identical(first$fold_id, repeated$fold_id)
  expect_equal(first$lambda, repeated$lambda, tolerance = 0)
  expect_equal(
    first$lambda,
    first$curve$lambda[which.min(first$curve$mean_deviance)],
    tolerance = 0
  )
  expect_identical(first$selection_rule, "lambda.min")

  fold_has_both_outcomes <- vapply(seq_len(5L), function(fold) {
    setequal(unique(d$event[first$fold_id == fold]), c(0L, 1L))
  }, logical(1))
  expect_identical(fold_has_both_outcomes, rep(TRUE, 5L))
})

test_that("Cox inner tuning passes the declared standardization setting", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 80, seed = 1)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  tuned <- sb_select_glmnet_cox_lambda(
    x,
    d$time,
    d$event,
    inner_folds = 3L,
    seed = 1701L,
    standardize = FALSE,
    selection_rule = "lambda.1se",
    measure = "C"
  )
  manual <- glmnet::cv.glmnet(
    x = as.matrix(x),
    y = survival::Surv(d$time, d$event),
    family = "cox",
    cox.ties = "breslow",
    alpha = 1,
    foldid = tuned$fold_id,
    type.measure = "C",
    standardize = FALSE
  )

  expect_identical(tuned$standardize, FALSE)
  expect_identical(tuned$selection_rule, "lambda.1se")
  expect_identical(tuned$measure, "C")
  expect_equal(tuned$lambda, unname(manual$lambda.1se), tolerance = 0)
  expect_gt(tuned$lambda, unname(manual$lambda.min))
  expect_equal(tuned$curve$lambda, manual$lambda, tolerance = 0)
  expect_equal(tuned$curve$mean_measure, manual$cvm, tolerance = 0)
})

test_that("outer assessment changes cannot alter inner lambda selection", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 90, seed = 20260813)
  analysis_rows <- 1:70
  assessment_rows <- 71:90
  predictors <- c("gene_a", "gene_b", "gene_c")
  select_from_analysis <- function(data) {
    sb_select_glmnet_cox_lambda(
      data[analysis_rows, predictors, drop = FALSE],
      data$time[analysis_rows],
      data$event[analysis_rows],
      inner_folds = 5L,
      seed = 1701L
    )
  }
  original <- select_from_analysis(d)

  changed <- d
  changed$time[assessment_rows] <- changed$time[assessment_rows] + 10000
  changed$event[assessment_rows] <- 1L - changed$event[assessment_rows]
  changed$gene_a[assessment_rows] <- seq_along(assessment_rows) * 10000
  perturbed <- select_from_analysis(changed)

  expect_equal(original$lambda, perturbed$lambda, tolerance = 0)
  expect_identical(original$fold_id, perturbed$fold_id)
  expect_identical(original$curve, perturbed$curve)
})

test_that("inner CV reports missing glmnet and infeasible folds", {
  d <- make_sb_fixture(n = 20, seed = 20260813)
  x <- d[c("gene_a", "gene_b")]
  local_mocked_bindings(sb_glmnet_available = function() FALSE)

  missing_error <- capture_sb_error(sb_select_glmnet_cox_lambda(
    x, d$time, d$event, inner_folds = 5L, seed = 1701L
  ))
  expect_s3_class(missing_error, "stormbreaker_missing_glmnet")

  local_mocked_bindings(sb_glmnet_available = function() TRUE)
  small <- make_sb_fixture(n = 8, seed = 20260813)
  fold_error <- capture_sb_error(sb_select_glmnet_cox_lambda(
    small[c("gene_a", "gene_b")],
    small$time,
    small$event,
    inner_folds = 5L,
    seed = 1701L
  ))
  expect_s3_class(fold_error, "stormbreaker_infeasible_inner_folds")
})

test_that("inner CV records a glmnet fitting error explicitly", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 20, seed = 20260813)
  x <- data.frame(constant_a = 1, constant_b = 2)[rep(1L, 20L), ]

  error <- capture_sb_error(sb_select_glmnet_cox_lambda(
    x, d$time, d$event, inner_folds = 5L, seed = 1701L
  ))

  expect_s3_class(error, "stormbreaker_glmnet_tuning_failure")
  expect_s3_class(error$parent, "simpleError")
  expect_match(conditionMessage(error), "zero variance")
})

test_that("inner CV preserves its dedicated warning failure class", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 20, seed = 20260813)
  original_cv_glmnet <- glmnet::cv.glmnet
  local_mocked_bindings(
    cv.glmnet = function(...) {
      warning("synthetic cv warning from glmnet", call. = FALSE)
      original_cv_glmnet(...)
    },
    .package = "glmnet"
  )

  error <- capture_sb_error(sb_select_glmnet_cox_lambda(
    d[c("gene_a", "gene_b")],
    d$time,
    d$event,
    inner_folds = 5L,
    seed = 1701L
  ))

  expect_s3_class(error, "stormbreaker_glmnet_tuning_warning")
  expect_s3_class(error$parent, "simpleWarning")
  expect_match(conditionMessage(error), "synthetic cv warning from glmnet")
})

test_that("a validated convergence-tail warning remains visible", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60, seed = 20260813)
  tail_warning <- paste0(
    "from glmnet C++ code (error code -27); Convergence for 27th ",
    "lambda value not reached after maxit=100000 iterations; ",
    "solutions for larger lambdas returned"
  )
  original_cv_glmnet <- glmnet::cv.glmnet
  original_glmnet <- glmnet::glmnet
  observed_cv_ties <- character()
  observed_path_ties <- character()
  local_mocked_bindings(
    cv.glmnet = function(..., cox.ties = NULL) {
      observed_cv_ties <<- c(observed_cv_ties, cox.ties)
      warning(tail_warning, call. = FALSE)
      original_cv_glmnet(..., cox.ties = cox.ties)
    },
    glmnet = function(..., cox.ties = NULL) {
      observed_path_ties <<- c(observed_path_ties, cox.ties)
      original_glmnet(..., cox.ties = cox.ties)
    },
    .package = "glmnet"
  )

  tuned <- sb_select_glmnet_cox_lambda(
    d[c("gene_a", "gene_b", "gene_c")],
    d$time,
    d$event,
    inner_folds = 5L,
    seed = 1701L
  )
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = tuned$alpha,
    lambda = tuned$lambda,
    standardize = tuned$standardize
  )
  learner$fit(
    d[c("gene_a", "gene_b", "gene_c")],
    d$time,
    d$event
  )

  expect_identical(tuned$warnings, tail_warning)
  expect_identical(observed_cv_ties, "breslow")
  expect_gte(length(observed_path_ties), 6L)
  expect_identical(unique(observed_path_ties), "breslow")
  expect_identical(nrow(tuned$warning_validation), 5L)
  expect_identical(tuned$warning_validation$fold, seq_len(5L))
  expect_true(all(tuned$lambda >= tuned$warning_validation$min_lambda))
  expect_true(all(tuned$warning_validation$assessment_risk_finite))
})

test_that("Cox refits accept only usable, recorded convergence tails", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 80, seed = 20260813)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  original_glmnet <- glmnet::glmnet
  native <- original_glmnet(
    as.matrix(x), survival::Surv(d$time, d$event),
    family = "cox", cox.ties = "breslow"
  )
  tail_warning <- paste0(
    "from glmnet C++ code (error code -88); Convergence for 88th ",
    "lambda value not reached after maxit=100000 iterations; ",
    "solutions for larger lambdas returned"
  )
  local_mocked_bindings(
    glmnet = function(...) {
      warning(tail_warning, call. = FALSE)
      original_glmnet(...)
    },
    .package = "glmnet"
  )
  learner <- sb_glmnet_cox_learner_adapter(1, native$lambda[[5L]])
  fitted <- learner$fit(x, d$time, d$event)
  expect_identical(attr(fitted, "stormbreaker_fit_warnings"), tail_warning)
  expect_equal(
    learner$predict(fitted, x),
    unname(drop(predict(
      native, as.matrix(x), s = native$lambda[[5L]], type = "link"
    ))),
    tolerance = 0
  )

  outside <- sb_glmnet_cox_learner_adapter(1, min(native$lambda) / 2)
  error <- capture_sb_error(outside$fit(x, d$time, d$event))
  expect_s3_class(error, "stormbreaker_glmnet_lambda_outside_path")

  local_mocked_bindings(
    glmnet = function(...) warning("unrecognized fit warning", call. = FALSE),
    .package = "glmnet"
  )
  error <- capture_sb_error(learner$fit(x, d$time, d$event))
  expect_s3_class(error, "stormbreaker_glmnet_fit_warning")
})

test_that("Cox refit warnings reach CV and final-fit results", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 90, seed = 20260813)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_cox", inner_folds = 3L, inner_seed = 1701L
  )
  original_glmnet <- glmnet::glmnet
  tail_warning <- paste0(
    "from glmnet C++ code (error code -88); Convergence for 88th ",
    "lambda value not reached after maxit=100000 iterations; ",
    "solutions for larger lambdas returned"
  )
  local_mocked_bindings(
    glmnet = function(...) {
      warning(tail_warning, call. = FALSE)
      original_glmnet(...)
    },
    .package = "glmnet"
  )
  cv <- sb_run_cv(task, pipeline, folds = 3L, seed = 1701L)
  final <- sb_fit_final(task, pipeline)
  expect_identical(cv$status, "success")
  expect_identical(final$status, "success")
  expect_identical(
    vapply(
      cv$tuning,
      function(z) any(grepl("Learner refit:", z$warnings)),
      logical(1)
    ),
    rep(TRUE, 3L)
  )
  expect_match(tail(final$tuning$warnings, 1L), "Learner refit:")
})

test_that("an inner-selected lambda refits and predicts held-out patients", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 90, seed = 20260813)
  analysis_rows <- 1:70
  assessment_rows <- 71:90
  predictors <- c("gene_a", "gene_b", "gene_c")
  tuned <- sb_select_glmnet_cox_lambda(
    d[analysis_rows, predictors, drop = FALSE],
    d$time[analysis_rows],
    d$event[analysis_rows],
    inner_folds = 5L,
    seed = 1701L
  )
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = tuned$alpha,
    lambda = tuned$lambda,
    standardize = tuned$standardize
  )
  model <- learner$fit(
    d[analysis_rows, predictors, drop = FALSE],
    d$time[analysis_rows],
    d$event[analysis_rows]
  )
  risk <- learner$predict(
    model,
    d[assessment_rows, predictors, drop = FALSE]
  )

  expect_length(risk, length(assessment_rows))
  expect_identical(all(is.finite(risk)), TRUE)
  expect_gt(stats::sd(risk), 0)
})

test_that("glmnet Cox preserves native predictions at rounded path endpoints", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 80, seed = 20260813)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  native <- glmnet::glmnet(as.matrix(x), survival::Surv(d$time, d$event),
    family = "cox", cox.ties = "breslow")
  endpoints <- range(native$lambda)
  near <- endpoints * c(1 - 2 * .Machine$double.eps,
    1 + 2 * .Machine$double.eps)
  for (lambda in near) {
    learner <- sb_glmnet_cox_learner_adapter(1, lambda)
    model <- learner$fit(x, d$time, d$event)
    expect_identical(model$lambda, native$lambda)
    expect_identical(attr(model, "stormbreaker_lambda"), lambda)
    expect_equal(learner$predict(model, x),
      as.numeric(predict(native, as.matrix(x), s = lambda, type = "link")),
      tolerance = 0)
  }
  for (lambda in endpoints * c(1 - 1e-6, 1 + 1e-6)) {
    learner <- sb_glmnet_cox_learner_adapter(1, lambda)
    error <- capture_sb_error(learner$fit(x, d$time, d$event))
    expect_s3_class(error, "stormbreaker_glmnet_lambda_outside_path")
  }
})

test_that("Cox convergence validation distinguishes endpoint roundoff from a missing solution", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 80, seed = 20260813)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  folds <- rep(1:4, length.out = nrow(d))
  minima <- vapply(1:4, function(f) {
    ix <- folds != f
    model <- glmnet::glmnet(as.matrix(x[ix, ]),
      survival::Surv(d$time[ix], d$event[ix]), family = "cox", cox.ties = "breslow")
    min(model$lambda)
  }, numeric(1))
  lambda <- max(minima) * (1 - 2 * .Machine$double.eps)
  validation <- sb_validate_glmnet_selected_lambda(x, d$time, d$event,
    folds, lambda, alpha = 1, standardize = TRUE)
  expect_equal(nrow(validation), 4L)
  expect_true(all(validation$assessment_risk_finite))
  error <- capture_sb_error(sb_validate_glmnet_selected_lambda(
    x, d$time, d$event, folds, max(minima) * (1 - 1e-6),
    alpha = 1, standardize = TRUE))
  expect_s3_class(error, "stormbreaker_glmnet_selected_lambda_not_converged")
})
