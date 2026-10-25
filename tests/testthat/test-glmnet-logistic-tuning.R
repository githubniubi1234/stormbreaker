test_that("binary inner CV selects lambda reproducibly", {
  skip_if_not_installed("glmnet")

  d <- make_sb_high_dim_binary_fixture(n = 60, p = 30)
  predictors <- grep("^gene_", names(d), value = TRUE)
  first <- sb_select_glmnet_logistic_lambda(
    d[predictors], d$outcome,
    inner_folds = 5L, seed = 1701L, alpha = 0.5
  )
  repeated <- sb_select_glmnet_logistic_lambda(
    d[predictors], d$outcome,
    inner_folds = 5L, seed = 1701L, alpha = 0.5
  )

  expect_identical(first$fold_id, repeated$fold_id)
  expect_equal(first$lambda, repeated$lambda, tolerance = 0)
  expect_identical(first$curve, repeated$curve)
  expect_equal(
    first$lambda,
    first$curve$lambda[which.min(first$curve$mean_deviance)],
    tolerance = 0
  )
  expect_identical(first$selection_rule, "lambda.min")
  expect_identical(first$measure, "deviance")
  fold_has_both_classes <- vapply(seq_len(5L), function(fold) {
    setequal(unique(d$outcome[first$fold_id == fold]), c(0L, 1L))
  }, logical(1))
  expect_identical(fold_has_both_classes, rep(TRUE, 5L))
})

test_that("binary inner tuning matches direct cv.glmnet", {
  skip_if_not_installed("glmnet")

  d <- make_sb_high_dim_binary_fixture(n = 60, p = 20)
  predictors <- grep("^gene_", names(d), value = TRUE)
  tuned <- sb_select_glmnet_logistic_lambda(
    d[predictors], d$outcome,
    inner_folds = 5L, seed = 101L, alpha = 0.3,
    selection_rule = "lambda.1se", measure = "auc"
  )
  manual <- glmnet::cv.glmnet(
    x = as.matrix(d[predictors]),
    y = d$outcome,
    family = "binomial",
    alpha = 0.3,
    foldid = tuned$fold_id,
    type.measure = "auc",
    standardize = TRUE
  )

  expect_identical(tuned$selection_rule, "lambda.1se")
  expect_identical(tuned$measure, "auc")
  expect_equal(tuned$lambda, unname(manual$lambda.1se), tolerance = 0)
  expect_gt(tuned$lambda, unname(manual$lambda.min))
  expect_equal(tuned$curve$lambda, manual$lambda, tolerance = 0)
  expect_equal(tuned$curve$mean_measure, manual$cvm, tolerance = 0)
  expect_equal(tuned$curve$standard_error, manual$cvsd, tolerance = 0)
  expect_equal(tuned$curve$nonzero, unname(manual$nzero), tolerance = 0)
})

test_that("outer assessment data cannot alter binary inner tuning", {
  skip_if_not_installed("glmnet")

  d <- make_sb_high_dim_binary_fixture(n = 70, p = 20)
  analysis_rows <- 1:55
  assessment_rows <- 56:70
  predictors <- grep("^gene_", names(d), value = TRUE)
  select_from_analysis <- function(data) {
    sb_select_glmnet_logistic_lambda(
      data[analysis_rows, predictors, drop = FALSE],
      data$outcome[analysis_rows],
      inner_folds = 5L,
      seed = 101L,
      alpha = 1
    )
  }
  original <- select_from_analysis(d)
  changed <- d
  changed$outcome[assessment_rows] <- 1L - changed$outcome[assessment_rows]
  changed[assessment_rows, predictors] <-
    changed[assessment_rows, predictors] * 10000
  perturbed <- select_from_analysis(changed)

  expect_identical(original, perturbed)
})

test_that("binary inner tuning records dependency, folds, and warnings", {
  d <- make_sb_binary_fixture(n = 20)
  local_mocked_bindings(sb_glmnet_available = function() FALSE)
  missing_error <- capture_sb_error(sb_select_glmnet_logistic_lambda(
    d[c("gene_a", "gene_b")], d$outcome,
    inner_folds = 5L, seed = 1L
  ))
  expect_s3_class(missing_error, "stormbreaker_missing_glmnet")

  local_mocked_bindings(sb_glmnet_available = function() TRUE)
  small_outcome <- c(rep(1L, 2), rep(0L, 8))
  fold_error <- capture_sb_error(sb_select_glmnet_logistic_lambda(
    data.frame(a = 1:10, b = 11:20), small_outcome,
    inner_folds = 3L, seed = 1L
  ))
  expect_s3_class(fold_error, "stormbreaker_infeasible_inner_folds")

  skip_if_not_installed("glmnet")
  original_cv_glmnet <- glmnet::cv.glmnet
  local_mocked_bindings(
    cv.glmnet = function(...) {
      warning("synthetic binary tuning warning", call. = FALSE)
      original_cv_glmnet(...)
    },
    .package = "glmnet"
  )
  warning_error <- capture_sb_error(sb_select_glmnet_logistic_lambda(
    d[c("gene_a", "gene_b")], d$outcome,
    inner_folds = 5L, seed = 1L
  ))
  expect_s3_class(
    warning_error,
    "stormbreaker_glmnet_logistic_tuning_warning"
  )
  expect_match(conditionMessage(warning_error), "synthetic binary tuning")
})

test_that("binary inner tuning exposes glmnet fitting errors", {
  skip_if_not_installed("glmnet")

  d <- make_sb_binary_fixture(n = 20)
  x <- data.frame(constant_a = rep(1, 20), constant_b = rep(2, 20))
  error <- capture_sb_error(sb_select_glmnet_logistic_lambda(
    x, d$outcome, inner_folds = 5L, seed = 1L
  ))

  expect_s3_class(error, "stormbreaker_glmnet_logistic_tuning_failure")
  expect_s3_class(error$parent, "simpleError")
})
