test_that("Gaussian Naive Bayes matches independent population calculations", {
  x <- data.frame(
    gene_a = c(0, 1, 2, 2, 4, 6),
    gene_b = c(1, 3, 5, 0, 1, 2),
    check.names = FALSE
  )
  outcome <- c(0L, 0L, 0L, 1L, 1L, 1L)
  new_x <- data.frame(
    gene_a = c(1.5, 5),
    gene_b = c(2, 1),
    check.names = FALSE
  )
  var_smoothing <- 0.1
  epsilon <- var_smoothing * max(vapply(x, population_variance, numeric(1)))
  means <- rbind(
    `0` = vapply(x[outcome == 0L, , drop = FALSE], mean, numeric(1)),
    `1` = vapply(x[outcome == 1L, , drop = FALSE], mean, numeric(1))
  )
  variances <- rbind(
    `0` = vapply(
      x[outcome == 0L, , drop = FALSE], population_variance, numeric(1)
    ) + epsilon,
    `1` = vapply(
      x[outcome == 1L, , drop = FALSE], population_variance, numeric(1)
    ) + epsilon
  )
  priors <- c(`0` = 0.5, `1` = 0.5)
  adapter <- sb_gaussian_naive_bayes_learner_adapter(var_smoothing)
  model <- adapter$fit(x, outcome)
  observed <- adapter$predict(model, new_x)
  one_at_a_time <- vapply(seq_len(nrow(new_x)), function(row) {
    adapter$predict(model, new_x[row, , drop = FALSE])
  }, numeric(1))
  reversed <- adapter$predict(model, new_x[2:1, , drop = FALSE])[2:1]
  expected <- manual_gaussian_nb_probability(
    as.matrix(new_x), priors, means, variances
  )
  audit <- adapter$extract_parameters(model)

  expect_s3_class(model, "stormbreaker_gaussian_naive_bayes")
  expect_equal(model$epsilon, epsilon, tolerance = 0)
  expect_equal(model$means, means, tolerance = 0)
  expect_equal(model$variances, variances, tolerance = 0)
  expect_equal(observed, expected, tolerance = 1e-15)
  expect_equal(observed, one_at_a_time, tolerance = 0)
  expect_equal(observed, reversed, tolerance = 0)
  expect_identical(all(observed >= 0 & observed <= 1), TRUE)
  expect_identical(
    audit$priors,
    data.frame(
      class = c(0L, 1L),
      n = c(3L, 3L),
      prior = c(0.5, 0.5),
      stringsAsFactors = FALSE
    )
  )
  expect_identical(audit$var_smoothing, var_smoothing)
  expect_equal(audit$epsilon, epsilon, tolerance = 0)
  expect_identical(
    audit$values$feature,
    rep(names(x), times = 2L)
  )
  expect_identical(audit$values$class, rep(c(0L, 1L), each = ncol(x)))
  expect_equal(
    audit$values$mean,
    unname(c(means["0", ], means["1", ])),
    tolerance = 0
  )
  expect_equal(
    audit$values$variance,
    unname(c(variances["0", ], variances["1", ])),
    tolerance = 0
  )
  expect_equal(
    audit$values$standard_deviation,
    sqrt(audit$values$variance),
    tolerance = 0
  )
})

test_that("Gaussian Naive Bayes locks predictor schema and smoothing", {
  x <- data.frame(
    gene_a = c(0, 1, 2, 3),
    gene_b = c(1, 3, 2, 4)
  )
  outcome <- c(0L, 0L, 1L, 1L)
  adapter <- sb_gaussian_naive_bayes_learner_adapter(1e-9)
  model <- adapter$fit(x, outcome)

  schema_error <- capture_sb_error(
    adapter$predict(model, x[c("gene_b", "gene_a")])
  )
  expect_s3_class(schema_error, "stormbreaker_learner_schema_mismatch")

  invalid_values <- list(-1, NA_real_, Inf, "1e-9", c(0, 1e-9))
  for (value in invalid_values) {
    error <- capture_sb_error(
      sb_gaussian_naive_bayes_learner_adapter(value)
    )
    expect_s3_class(
      error,
      "stormbreaker_bad_naive_bayes_var_smoothing"
    )
  }
})

test_that("Gaussian Naive Bayes predicts a one-patient assessment fold", {
  x <- data.frame(
    gene_a = c(0, 1, 2, 3, 4, 5),
    gene_b = c(1, 3, 5, 0, 2, 4)
  )
  outcome <- rep(c(0L, 1L), each = 3L)
  new_x <- x[4L, , drop = FALSE]
  adapter <- sb_gaussian_naive_bayes_learner_adapter(1e-9)
  model <- adapter$fit(x, outcome)

  probability <- adapter$predict(model, new_x)

  expect_length(probability, 1L)
  expect_true(is.finite(probability))
  expect_true(probability >= 0 && probability <= 1)
})

test_that("Gaussian Naive Bayes exposes invalid fitted variances", {
  all_constant <- data.frame(
    gene_a = rep(1, 6),
    gene_b = rep(2, 6)
  )
  outcome <- rep(c(0L, 1L), each = 3L)
  adapter <- sb_gaussian_naive_bayes_learner_adapter(1e-9)

  constant_error <- capture_sb_error(adapter$fit(all_constant, outcome))
  expect_s3_class(
    constant_error,
    "stormbreaker_invalid_naive_bayes_variance"
  )

  within_class_constant <- data.frame(
    gene_a = c(rep(0, 3), rep(1, 3)),
    gene_b = c(rep(2, 3), rep(3, 3))
  )
  no_smoothing <- sb_gaussian_naive_bayes_learner_adapter(0)
  zero_error <- capture_sb_error(
    no_smoothing$fit(within_class_constant, outcome)
  )
  expect_s3_class(
    zero_error,
    "stormbreaker_invalid_naive_bayes_variance"
  )
})
