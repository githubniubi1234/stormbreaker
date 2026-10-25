reverse_fixture <- function(probability, outcome = c(0L, 0L, 1L, 1L)) {
  ids <- paste0("p", seq_along(outcome))
  prediction <- structure(list(
    status = "success", task_type = "binary", pipeline_id = "test",
    learner = "logistic", id_col = "id",
    predictions = data.frame(
      patient_id = ids, pipeline_id = "test", learner = "logistic",
      positive_class = 1L, probability = probability, status = "success"
    )
  ), class = "stormbreaker_final_prediction")
  endpoint <- sb_binary_endpoint(
    data.frame(id = rev(ids), outcome = rev(outcome)), "id", "outcome"
  )
  list(prediction = prediction, endpoint = endpoint)
}

test_that("reverse flips only low AUC and retains aligned originals", {
  x <- reverse_fixture(c(.9, .8, .2, .1))
  before <- x$prediction
  original <- sb_evaluate_external(x$prediction, x$endpoint, reverse = FALSE)
  off <- sb_evaluate_external(x$prediction, x$endpoint, reverse = FALSE)
  on <- sb_evaluate_external(x$prediction, x$endpoint, reverse = TRUE)
  expect_identical(sb_evaluate_external(x$prediction, x$endpoint), on)
  expect_identical(off, original)
  expect_identical(x$prediction, before)
  expect_identical(on$original_metrics, original$metrics)
  expect_identical(on$original_scored_predictions, original$scored_predictions)
  expect_equal(on$direction_correction$auc_original, 0)
  expect_equal(on$direction_correction$auc_processed, 1)
  expect_identical(on$direction_correction$reversed, TRUE)
  expect_equal(on$scored_predictions$probability,
               1 - original$scored_predictions$probability)
  expect_equal(on$metrics$brier_score$estimate,
    mean((on$scored_predictions$outcome - on$scored_predictions$probability)^2))
  expect_identical(on$evaluation_boundary$outcome_use,
                  "metrics_and_post_hoc_direction_selection")
})

test_that("good, tied and unavailable AUC do not reverse", {
  for (p in list(c(.1, .2, .8, .9), rep(.5, 4), c(.2, .8, .2, .8))) {
    x <- reverse_fixture(p)
    result <- sb_evaluate_external(x$prediction, x$endpoint, reverse = TRUE)
    expect_identical(result$direction_correction$reversed, FALSE)
    expect_equal(result$scored_predictions$probability,
                 result$original_scored_predictions$probability)
  }
  x <- reverse_fixture(c(.1, .2, .3, .4), rep(1L, 4))
  result <- sb_evaluate_external(x$prediction, x$endpoint, reverse = TRUE)
  expect_identical(result$direction_correction$status, "not_evaluable")
  expect_identical(result$direction_correction$reversed, FALSE)
})

test_that("reverse validates its flag and preserves prediction failures", {
  x <- reverse_fixture(c(.1, .2, .8, .9))
  for (flag in list(NA, 1, "TRUE", logical(), c(TRUE, FALSE))) {
    error <- tryCatch(
      sb_evaluate_external(x$prediction, x$endpoint, reverse = flag),
      error = identity
    )
    expect_s3_class(error, "stormbreaker_bad_reverse")
  }
  result <- sb_evaluate_external(list(), x$endpoint, reverse = TRUE)
  expect_identical(result$status, "failed")
})

test_that("tied scores retain half credit after reversal", {
  x <- reverse_fixture(c(.8, .8, .8, .2))
  result <- sb_evaluate_external(x$prediction, x$endpoint, reverse = TRUE)
  expect_equal(result$direction_correction$auc_original, .25)
  expect_equal(result$direction_correction$auc_processed, .75)
})

test_that("omitted reverse preserves survival evaluation", {
  x <- reverse_fixture(c(.8, .7, .3, .2))
  p <- x$prediction
  p$task_type <- "survival"
  p$predictions$risk <- c(.8, .7, .3, .2)
  e <- sb_surv_endpoint(data.frame(id = paste0("p", 1:4),
    time = 1:4, event = c(1, 0, 1, 0)), "id", "time", "event", "day")
  expect_identical(sb_evaluate_external(p, e),
                   sb_evaluate_external(p, e, reverse = FALSE))
})
