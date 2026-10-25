make_selection_audit <- function(
    fold,
    selected_features,
    status = "success",
    error_class = NA_character_,
    error_message = NA_character_) {
  list(
    fold = as.integer(fold),
    selector = "test_selector",
    status = status,
    selected_features = selected_features,
    scores = NULL,
    error_class = error_class,
    error_message = error_message
  )
}

test_that("selection frequency and Jaccard match a hand calculation", {
  selections <- list(
    make_selection_audit(1L, c("a", "b")),
    make_selection_audit(2L, c("a", "c")),
    make_selection_audit(3L, c("a", "b"))
  )

  audit <- sb_selection_stability(
    "test_selector",
    selections,
    c("d", "b", "a", "c")
  )

  expect_named(
    audit,
    c(
      "selector", "status", "scope", "aggregation",
      "frequency_denominator", "n_candidate_features", "n_folds",
      "n_successful_folds", "n_pairs", "n_evaluable_pairs",
      "mean_pairwise_jaccard", "feature_frequency", "per_fold",
      "pairwise", "error_message"
    ),
    ignore.order = FALSE
  )
  expect_identical(audit$status, "success")
  expect_identical(audit$scope, "outer_analysis_fold_selected_feature_sets")
  expect_identical(
    audit$aggregation,
    "mean_of_unique_fold_pair_jaccard"
  )
  expect_identical(audit$frequency_denominator, "successful_selection_folds")
  expect_identical(audit$n_candidate_features, 4L)
  expect_identical(audit$n_folds, 3L)
  expect_identical(audit$n_successful_folds, 3L)
  expect_identical(audit$n_pairs, 3L)
  expect_identical(audit$n_evaluable_pairs, 3L)
  expect_identical(audit$feature_frequency$feature, c("d", "b", "a", "c"))
  expect_identical(audit$feature_frequency$selected_count, c(0L, 2L, 3L, 1L))
  expect_equal(
    audit$feature_frequency$selection_frequency,
    c(0, 2 / 3, 1, 1 / 3),
    tolerance = 0
  )
  expect_identical(audit$pairwise$fold_i, c(1L, 1L, 2L))
  expect_identical(audit$pairwise$fold_j, c(2L, 3L, 3L))
  expect_identical(audit$pairwise$intersection, c(1L, 2L, 1L))
  expect_identical(audit$pairwise$union, c(3L, 2L, 3L))
  expect_equal(audit$pairwise$jaccard, c(1 / 3, 1, 1 / 3), tolerance = 0)
  expect_equal(audit$mean_pairwise_jaccard, 5 / 9, tolerance = 0)
  expect_true(all(audit$pairwise$status == "success"))
  expect_true(is.na(audit$error_message))
})

test_that("fold and selected-feature order cannot alter stability", {
  original <- list(
    make_selection_audit(1L, c("a", "b")),
    make_selection_audit(2L, c("a", "c")),
    make_selection_audit(3L, c("a", "b"))
  )
  reordered <- list(
    make_selection_audit(3L, c("b", "a")),
    make_selection_audit(1L, c("b", "a")),
    make_selection_audit(2L, c("c", "a"))
  )

  expected <- sb_selection_stability(
    "test_selector",
    original,
    c("d", "b", "a", "c")
  )
  observed <- sb_selection_stability(
    "test_selector",
    reordered,
    c("d", "b", "a", "c")
  )

  expect_identical(observed, expected)
})

test_that("no learned selector is explicitly not applicable", {
  selections <- list(
    make_selection_audit(1L, c("a", "b")),
    make_selection_audit(2L, c("a", "b"))
  )

  audit <- sb_selection_stability("none", selections, c("a", "b"))

  expect_identical(audit$status, "not_applicable")
  expect_identical(audit$n_folds, 2L)
  expect_identical(audit$n_successful_folds, 0L)
  expect_identical(audit$n_pairs, 1L)
  expect_identical(audit$n_evaluable_pairs, 0L)
  expect_true(is.na(audit$mean_pairwise_jaccard))
  expect_identical(audit$feature_frequency$selected_count, c(0L, 0L))
  expect_true(all(is.na(audit$feature_frequency$selection_frequency)))
  expect_true(all(audit$per_fold$status == "not_applicable"))
  expect_true(all(audit$pairwise$status == "not_applicable"))
})

test_that("failed selector folds retain partial evidence without a headline", {
  selections <- list(
    make_selection_audit(3L, c("a", "c")),
    make_selection_audit(
      2L,
      character(),
      status = "failed",
      error_class = "stormbreaker_selection_failure",
      error_message = "selector failed"
    ),
    make_selection_audit(1L, c("a", "b"))
  )

  audit <- sb_selection_stability(
    "test_selector",
    selections,
    c("a", "b", "c", "d")
  )

  expect_identical(audit$status, "incomplete")
  expect_identical(audit$n_successful_folds, 2L)
  expect_identical(audit$n_pairs, 3L)
  expect_identical(audit$n_evaluable_pairs, 1L)
  expect_true(is.na(audit$mean_pairwise_jaccard))
  expect_identical(audit$feature_frequency$selected_count, c(2L, 1L, 1L, 0L))
  expect_equal(
    audit$feature_frequency$selection_frequency,
    c(1, 1 / 2, 1 / 2, 0),
    tolerance = 0
  )
  expect_identical(audit$per_fold$fold, c(1L, 2L, 3L))
  expect_identical(audit$per_fold$status, c("success", "failed", "success"))
  expect_identical(
    audit$per_fold$error_class,
    c(NA_character_, "stormbreaker_selection_failure", NA_character_)
  )
  expect_identical(audit$per_fold$error_message[[2L]], "selector failed")
  expect_true(is.na(audit$per_fold$n_selected[[2L]]))
  expect_identical(audit$pairwise$status, c("incomplete", "success", "incomplete"))
  expect_equal(audit$pairwise$jaccard, c(NA, 1 / 3, NA), tolerance = 0)
  expect_match(audit$error_message, "not every outer-fold selection succeeded")
})
