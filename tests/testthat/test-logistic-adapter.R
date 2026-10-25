test_that("the logistic adapter matches an independent glm fit", {
  d <- make_sb_binary_fixture(n = 60)
  analysis_rows <- 1:48
  assessment_rows <- 49:60
  predictors <- c("gene_a", "gene_b", "gene_c")
  learner <- sb_logistic_learner_adapter()

  model <- learner$fit(
    d[analysis_rows, predictors, drop = FALSE],
    d$outcome[analysis_rows]
  )
  observed <- learner$predict(
    model,
    d[assessment_rows, predictors, drop = FALSE]
  )
  manual <- stats::glm(
    outcome ~ gene_a + gene_b + gene_c,
    data = d[analysis_rows, , drop = FALSE],
    family = stats::binomial(link = "logit")
  )
  expected <- unname(stats::predict(
    manual,
    newdata = d[assessment_rows, , drop = FALSE],
    type = "response"
  ))

  expect_identical(learner$id, "logistic")
  expect_equal(observed, expected, tolerance = 1e-12)
  expect_identical(all(observed >= 0 & observed <= 1), TRUE)

  coefficients <- learner$extract_coefficients(model)
  expect_identical(
    coefficients$term,
    c("(Intercept)", predictors)
  )
  expect_equal(
    coefficients$coefficient,
    unname(stats::coef(manual)),
    tolerance = 1e-12
  )
})

test_that("the logistic adapter locks its numeric predictor schema", {
  d <- make_sb_binary_fixture(n = 40)
  learner <- sb_logistic_learner_adapter()
  model <- learner$fit(d[c("gene_a", "gene_b")], d$outcome)

  error <- capture_sb_error(
    learner$predict(model, d[c("gene_b", "gene_a")])
  )

  expect_s3_class(error, "stormbreaker_learner_schema_mismatch")
})

test_that("rank-deficient logistic fits fail without dropping terms", {
  d <- make_sb_binary_fixture(n = 40)
  x <- data.frame(first = d$gene_a, duplicate = d$gene_a)
  learner <- sb_logistic_learner_adapter()

  error <- capture_sb_error(learner$fit(x, d$outcome))

  expect_s3_class(error, "stormbreaker_invalid_logistic_fit")
})

test_that("the logistic adapter rejects an ambiguous intercept feature", {
  d <- make_sb_binary_fixture(n = 40)
  x <- d[c("gene_a", "gene_b")]
  names(x)[1L] <- "(Intercept)"
  learner <- sb_logistic_learner_adapter()

  error <- capture_sb_error(learner$fit(x, d$outcome))

  expect_s3_class(error, "stormbreaker_reserved_intercept_name")
})

test_that("saturated probabilities alone do not invalidate an overlapping fit", {
  d <- data.frame(
    gene = c(rep(-1, 20), rep(1, 20), -100, 100),
    outcome = c(rep(0, 18), rep(1, 2), rep(0, 2), rep(1, 18), 0, 1)
  )
  native <- suppressWarnings(stats::glm(
    outcome ~ gene, d, family = stats::binomial()
  ))
  model <- sb_fit_logistic(d["gene"], d$outcome)
  expect_equal(stats::coef(model), stats::coef(native), tolerance = 1e-12)
  expect_equal(
    sb_predict_logistic(model, d["gene"]),
    unname(stats::predict(native, d, type = "response")),
    tolerance = 1e-12
  )
  expect_match(
    attr(model, "stormbreaker_fit_warnings"),
    "fitted probabilities numerically 0 or 1 occurred", fixed = TRUE
  )
})

test_that("a probability warning does not admit a separated logistic fit", {
  d <- data.frame(gene = rep(c(-1, 1), each = 20), outcome = rep(0:1, each = 20))
  native <- suppressWarnings(stats::glm(
    outcome ~ gene, d, family = stats::binomial(), x = TRUE, y = TRUE
  ))
  native$converged <- TRUE
  local_mocked_bindings(
    glm = function(...) {
      warning("glm.fit: fitted probabilities numerically 0 or 1 occurred")
      native
    },
    .package = "stats"
  )
  error <- capture_sb_error(sb_fit_logistic(d["gene"], d$outcome))
  expect_s3_class(error, "stormbreaker_logistic_separation")
})

test_that("logistic still rejects unrelated warnings and inconclusive separation", {
  d <- make_sb_binary_fixture(n = 40)
  original_glm <- stats::glm
  message <- "unexpected fitting warning"
  local_mocked_bindings(
    glm = function(...) {
      warning(message)
      original_glm(...)
    },
    .package = "stats"
  )
  error <- capture_sb_error(sb_fit_logistic(d["gene_a"], d$outcome))
  expect_s3_class(error, "stormbreaker_logistic_warning")

  message <- "glm.fit: fitted probabilities numerically 0 or 1 occurred"
  local_mocked_bindings(
    sb_binary_moderate_calibration_separation = function(...) {
      list(status = "failed", detected = NA, error_message = "solver failed")
    }
  )
  error <- capture_sb_error(sb_fit_logistic(d["gene_a"], d$outcome))
  expect_s3_class(error, "stormbreaker_logistic_separation_check_failed")
})

test_that("logistic fit warnings remain in CV and final-fit results", {
  d <- make_sb_binary_fixture(n = 90)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  original_fit <- sb_fit_logistic
  message <- "glm.fit: fitted probabilities numerically 0 or 1 occurred"
  local_mocked_bindings(
    sb_fit_logistic = function(x, outcome) {
      model <- original_fit(x, outcome)
      attr(model, "stormbreaker_fit_warnings") <- message
      model
    }
  )
  cv <- sb_run_cv(task, sb_pipeline(learner = "logistic"), folds = 3L, seed = 1L)
  final <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  expect_identical(cv$status, "success")
  expect_identical(final$status, "success")
  expect_identical(
    vapply(
      cv$tuning, function(x) any(grepl("Learner refit:", x$warnings)), logical(1)
    ),
    rep(TRUE, 3L)
  )
  expect_match(tail(final$tuning$warnings, 1L), "Learner refit:")
})
