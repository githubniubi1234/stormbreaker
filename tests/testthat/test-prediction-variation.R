variation_comparison <- function() {
  make_task <- function(prefix, seed) {
    d <- sb_example_binary(60, seed)
    d$patient_id <- paste0(prefix, d$patient_id)
    sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  }
  sb_compare_cohorts(make_task("dev", 1),
    list(Validation = make_task("val", 2), Test = make_task("test", 3)),
    roles = c(Validation = "validation", Test = "test"),
    models = c("logistic", "gaussian_naive_bayes"), folds = 3,
    view = FALSE, progress = FALSE)
}

test_that("variation distinguishes exact ties, small differences and singleton batches", {
  rows <- data.frame(status = "success", risk = c(1, 1))
  expect_identical(sb_prediction_variation(rows, "survival")$status, "constant")
  rows$risk[2] <- 1 + 1e-12
  expect_identical(sb_prediction_variation(rows, "survival")$status, "variable")
  expect_identical(sb_prediction_variation(rows[1, ], "survival")$status,
    "insufficient_patients")
  rows$risk[2] <- NA_real_
  expect_identical(sb_prediction_variation(rows, "survival")$status, "incomplete")
})

test_that("different constants in different CV folds do not count as discrimination", {
  rows <- data.frame(status = "success", probability = rep(c(.2, .4, .6), each = 3),
    fold = rep(1:3, each = 3))
  diagnostic <- sb_prediction_variation(rows, "binary")
  expect_identical(diagnostic$status, "constant")
  expect_identical(diagnostic$constant_folds, 1:3)
  rows$probability[1] <- .3
  expect_identical(sb_prediction_variation(rows, "binary")$status, "partly_constant")
})

test_that("constant validation predictions cannot beat a worse-than-chance model", {
  compared <- variation_comparison()
  ids <- unname(vapply(compared$model_pipelines, function(x) x$id, ""))
  compared$cohorts$Validation$comparison$roc_auc <- c(.5, .4)
  compared$cohorts$Validation$evaluations[[ids[1]]]$scored_predictions$probability <- .5
  selected <- sb_select_model(compared)
  expect_identical(selected$pipeline$id, ids[2])
  ranking <- selected$model_selection$ranking
  expect_identical(ranking$eligible[match(ids[1], ranking$pipeline_id)], FALSE)
  expect_match(ranking$exclusion_reason[match(ids[1], ranking$pipeline_id)],
    "constant_validation_predictions:Validation", fixed = TRUE)
  expect_equal(compared$cohorts$Validation$comparison$roc_auc[1], .5)
})

test_that("constant development folds are ineligible without reading test predictions", {
  compared <- variation_comparison()
  ids <- unname(vapply(compared$model_pipelines, function(x) x$id, ""))
  compared$cohorts$Validation$comparison$roc_auc <- c(.9, .6)
  rows <- compared$development$results[[ids[1]]]$predictions
  rows$probability <- rows$fold / 10
  compared$development$results[[ids[1]]]$predictions <- rows
  selected <- sb_select_model(compared)
  expect_identical(selected$pipeline$id, ids[2])
  compared$cohorts$Test$evaluations[[ids[2]]]$scored_predictions$probability <- .5
  expect_identical(sb_select_model(compared)$pipeline$id, ids[2])
  compared$cohorts$Validation$evaluations[[ids[2]]]$scored_predictions$probability <- .5
  error <- tryCatch(sb_select_model(compared), error = identity)
  expect_s3_class(error, "stormbreaker_no_eligible_model")
  expect_match(conditionMessage(error), "constant", fixed = TRUE)
})

test_that("constant final predictions stay available and explicitly labelled", {
  compared <- variation_comparison()
  fit <- compared$fits[[1]]
  d <- sb_example_binary(10, 4)
  d$gene_a <- 0
  d$gene_b <- 0
  prediction <- sb_predict_final(fit, d, "patient_id")
  expect_identical(prediction$status, "success")
  expect_identical(prediction$prediction_variation$status, "constant")
  endpoint <- sb_binary_endpoint(d, "patient_id", "outcome")
  external <- sb_compare_external(list(prediction), endpoint)
  expect_identical(external$comparison$prediction_variation, "constant")
  expect_identical(external$comparison$roc_auc_status, "success")
  expect_equal(external$comparison$roc_auc, .5)
})
