manual_binary_auc <- function(outcome, score) {
  positive <- score[outcome == 1L]
  negative <- score[outcome == 0L]
  credits <- unlist(lapply(positive, function(pos) {
    (pos > negative) + 0.5 * (pos == negative)
  }))
  mean(credits)
}

test_that("binary AUC Top-K matches independently enumerated pairs", {
  outcome <- c(0L, 0L, 1L, 1L)
  x <- data.frame(
    higher = c(0, 1, 2, 3),
    lower = c(3, 2, 1, 0),
    tied = c(0, 1, 0, 1),
    constant = rep(4, 4),
    check.names = FALSE
  )
  fitted <- sb_fit_univ_binary_auc_top_k(x, outcome, n_features = 2L)
  scores <- fitted$scores

  expect_s3_class(fitted, "stormbreaker_selector_fit")
  expect_identical(fitted$id, "univ_binary_auc_top_k")
  expect_identical(fitted$source_features, names(x))
  expect_identical(fitted$selected_features, c("higher", "lower"))
  expect_identical(scores$feature, names(x))
  expect_equal(
    scores$auc[1:3],
    unname(vapply(
      x[1:3], manual_binary_auc, numeric(1), outcome = outcome
    )),
    tolerance = 0
  )
  expect_equal(scores$auc[1:3], c(1, 0, 0.5), tolerance = 0)
  expect_equal(scores$discrimination[1:3], c(1, 1, 0.5), tolerance = 0)
  expect_identical(
    scores$direction[1:3],
    c("higher_in_positive", "lower_in_positive", "tied")
  )
  expect_identical(scores$comparable_pairs[1:3], rep(4, 3L))
  expect_identical(scores$rank[1:3], 1:3)
  expect_identical(scores$selected, c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(scores$status[1:3], rep("success", 3L))
  expect_identical(scores$status[4], "failed")
  expect_identical(
    scores$error_class[4],
    "stormbreaker_selector_constant_feature"
  )
})

test_that("binary AUC Top-K applies only its fitted columns in fitted order", {
  x <- data.frame(
    reverse = c(4, 3, 2, 1, 0, -1),
    weak = c(0, 1, 1, 0, 1, 0),
    direct = c(-1, 0, 1, 2, 3, 4),
    check.names = FALSE
  )
  outcome <- c(0L, 0L, 0L, 1L, 1L, 1L)
  adapter <- sb_univ_binary_auc_top_k_selector_adapter(2L)
  fitted <- adapter$fit(x, outcome)
  applied <- adapter$apply(fitted, x)

  expect_identical(adapter$id, "univ_binary_auc_top_k")
  expect_identical(adapter$params, list(n_features = 2L))
  expect_identical(adapter$needs_outcome, TRUE)
  expect_identical(names(applied), fitted$selected_features)
  expect_identical(applied, x[fitted$selected_features])

  schema_error <- capture_sb_error(
    adapter$apply(fitted, x[c("weak", "reverse", "direct")])
  )
  expect_s3_class(schema_error, "stormbreaker_selector_schema_mismatch")
})

test_that("complementary AUC directions break ties by original order", {
  outcome <- c(1L, 0L, 0L, 0L)
  first <- c(3, 1, 2, 4)
  x <- data.frame(
    first = first,
    second = -first,
    check.names = FALSE
  )

  fitted <- sb_fit_univ_binary_auc_top_k(x, outcome, n_features = 1L)

  expect_equal(fitted$scores$auc, c(2 / 3, 1 / 3), tolerance = 1e-15)
  expect_identical(
    fitted$scores$discrimination[[1L]],
    fitted$scores$discrimination[[2L]]
  )
  expect_identical(fitted$scores$rank, 1:2)
  expect_identical(fitted$selected_features, "first")
})

test_that("binary AUC Top-K fails when fewer than K scores are valid", {
  x <- data.frame(
    informative = c(0, 1, 2, 3, 4, 5),
    constant_a = rep(1, 6),
    constant_b = rep(2, 6)
  )
  outcome <- c(0L, 0L, 0L, 1L, 1L, 1L)

  error <- capture_sb_error(
    sb_fit_univ_binary_auc_top_k(x, outcome, n_features = 2L)
  )

  expect_s3_class(error, "stormbreaker_insufficient_valid_features")
  expect_match(conditionMessage(error), "only 1 features")
})

test_that("binary AUC Top-K resolves through the admitted runtime", {
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 2L
  )
  runtime <- sb_resolve_pipeline_runtime(pipeline)

  expect_identical(runtime$id, pipeline$id)
  expect_identical(runtime$selector$id, pipeline$selector)
  expect_identical(runtime$selector$params, pipeline$selector_params)
  expect_type(runtime$selector$fit, "closure")
  expect_type(runtime$selector$apply, "closure")
  expect_identical(runtime$learner$id, "logistic")
  expect_identical(sb_validate_pipeline_runtime(runtime, pipeline), runtime)
})
