merge_reverse_fixture <- function() {
  ids <- paste0("p", 1:8)
  data <- data.frame(id = ids, gene = seq_along(ids), y = rep(0:1, each = 4))
  prediction <- structure(list(
    status = "success", task_type = "binary", pipeline_id = "merge-test",
    learner = "logistic", id_col = "id",
    predictions = data.frame(patient_id = ids, pipeline_id = "merge-test",
      learner = "logistic", positive_class = 1L,
      probability = c(.9, .8, .7, .6, .4, .3, .2, .1), status = "success")
  ), class = "stormbreaker_final_prediction")
  list(data = data, prediction = prediction,
    endpoint = sb_binary_endpoint(data[8:1, ], "id", "y"))
}

test_that("merged batch comparison retains raw direction and variation diagnostics", {
  x <- merge_reverse_fixture()
  point <- sb_evaluate_external(x$prediction, x$endpoint)
  batch <- sb_compare_external(list(x$prediction), x$endpoint)
  expect_identical(formals(sb_evaluate_external)$reverse, TRUE)
  expect_equal(point$metrics$roc_auc$estimate, 1)
  expect_equal(batch$comparison$roc_auc, 0)
  expect_identical(batch$comparison$prediction_variation,
    sb_prediction_variation(x$prediction$predictions, "binary")$status)
  expect_null(batch$evaluations[[1]]$direction_correction)
  x$prediction$predictions$probability[] <- .5
  tied <- sb_compare_external(list(x$prediction), x$endpoint)
  expect_identical(tied$comparison$prediction_variation, "constant")
  expect_equal(tied$comparison$roc_auc, .5)
})

test_that("completion and export keep raw and processed values explicitly paired", {
  x <- merge_reverse_fixture()
  run <- list(
    task_type = "binary", pipelines = list(chosen = sb_pipeline(learner = "logistic")),
    development_task = sb_binary_task(x$data, "id", "y", "gene"),
    comparison = list(comparison = data.frame(model = "chosen")),
    details = list(chosen = list(fit = list(status = "success", selection = list(state = NULL))))
  )
  local_mocked_bindings(
    sb_predict_final = function(...) x$prediction,
    sb_fit_final = function(...) stop("Unexpected refit"),
    sb_compare_cv = function(...) stop("Unexpected model search"))
  original_run <- run
  result <- sb_complete_run(run, "chosen", list(validation = x$data[8:1, ]), plots = FALSE)
  expect_identical(run, original_run)
  expect_identical(result$predictions$validation, x$prediction$predictions)
  expect_identical(result$direction_correction$validation$reversed, TRUE)
  evaluated <- result$evaluated_predictions$validation
  expect_equal(result$metrics$validation$roc_auc$estimate, 1)
  expect_equal(evaluated$probability, 1 - evaluated$probability_original)
  expect_identical(evaluated$patient_id, rev(x$prediction$predictions$patient_id))
  raw <- sb_complete_run(run, "chosen", list(validation = x$data), plots = FALSE, reverse = FALSE)
  expect_equal(raw$metrics$validation$roc_auc$estimate, 0)
  expect_null(raw$direction_correction$validation)
  path <- tempfile("merged-export-")
  sb_export_run(result, path)
  correction <- read.csv(file.path(path, "cohort-001-direction-correction.csv"))
  scored <- read.csv(file.path(path, "cohort-001-evaluated-predictions.csv"))
  saved <- read.csv(file.path(path, "cohort-001-predictions.csv"))
  original_metrics <- read.csv(file.path(path, "cohort-001-original-metrics.csv"))
  expect_identical(correction$reversed, TRUE)
  expect_equal(correction$auc_original, 0)
  expect_equal(correction$auc_processed, 1)
  expect_equal(scored$probability, evaluated$probability)
  expect_equal(saved$probability, x$prediction$predictions$probability)
  expect_equal(original_metrics$estimate[original_metrics$metric == "roc_auc"], 0)
})

test_that("reverse exports still expose the expanded APIs and complete menus", {
  expected <- c("bi_run", "surv_run", "multi_run", "sb_binary_menu_pipelines",
    "sb_survival_menu_pipelines", "sb_selected_features", "sb_algorithm_diagnostics",
    "sb_compare_multiclass", "sb_model_shap", "sb_complete_run", "sb_export_run")
  expect_identical(all(expected %in% getNamespaceExports("stormbreaker")), TRUE)
  expect_identical(formals(sb_complete_run)$reverse, TRUE)
})
