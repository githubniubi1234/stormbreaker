test_that("locLDA unnamed posterior follows native grouping order", {
  skip_if_not_installed("klaR")
  posterior <- matrix(c(0.8, 0.3, 0.2, 0.7), nrow = 2L)
  local_mocked_bindings(
    predict.loclda = function(object, newdata) list(posterior = posterior),
    .package = "klaR"
  )
  engine <- list(grouping = factor(c(0, 1, 0), levels = c(0, 1)))
  expect_equal(sb_predict_loclda_binary(engine, matrix(1, 2, 1)),
               c(0.2, 0.7))
  engine$grouping <- factor(c(0, 1, 0), levels = c(1, 0))
  expect_equal(sb_predict_loclda_binary(engine, matrix(1, 2, 1)),
               c(0.8, 0.3))
})

test_that("locLDA named columns and invalid native values are preserved", {
  skip_if_not_installed("klaR")
  posterior <- matrix(c(0.2, NaN, 0.8, NaN), nrow = 2L,
                      dimnames = list(NULL, c("1", "0")))
  local_mocked_bindings(
    predict.loclda = function(object, newdata) list(posterior = posterior),
    .package = "klaR"
  )
  engine <- list(grouping = factor(c(0, 1)))
  expect_equal(sb_predict_loclda_binary(engine, matrix(1, 2, 1)),
               c(0.2, NaN))
  colnames(posterior) <- NULL
  expect_equal(sb_predict_loclda_binary(engine, matrix(1, 2, 1)),
               c(0.8, NaN))
})

test_that("locLDA native finite unnamed probabilities remain usable", {
  skip_if_not_installed("klaR")
  x <- withr::with_seed(1907L, matrix(rnorm(160L), ncol = 2L))
  grouping <- factor(rep(c(0, 1), 40L), levels = c(0, 1))
  engine <- klaR::loclda(x, grouping, k = 30L, weighted.apriori = TRUE)
  assessment <- matrix(c(1000, 1001, 1000, 1001), ncol = 2L)
  native <- get("predict.loclda", asNamespace("klaR"))(
    engine, assessment
  )$posterior
  expect_null(colnames(native))
  expect_equal(unname(native[, 2L]), c(0.5, 0.5))
  expect_equal(unname(sb_predict_loclda_binary(engine, assessment)),
               unname(native[, 2L]))
})

test_that("locLDA does not guess unidentified posterior classes", {
  skip_if_not_installed("klaR")
  local_mocked_bindings(
    predict.loclda = function(object, newdata) {
      list(posterior = matrix(0.5, nrow = 2L, ncol = 2L))
    },
    .package = "klaR"
  )
  expect_error(
    sb_predict_loclda_binary(list(grouping = factor("0")),
                                   matrix(1, 2, 1)),
    class = "stormbreaker_invalid_classical_binary_prediction"
  )
})
