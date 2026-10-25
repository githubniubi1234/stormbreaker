test_that("external metric tables copy all supported saved scalar values", {
  for (type in c("binary", "survival")) {
    x <- make_sb_curve_view_fixture(type, comparison = TRUE)$external
    before <- x
    metrics <- if (type == "binary") c("roc_auc", "average_precision", "brier_score", "log_loss") else "harrell_c"
    for (metric in metrics) {
      tab <- sb_result_table(x, metric = metric, cohort_label = "Independent cohort")
      expect_identical(tab$pipeline_id, names(x$evaluations))
      expect_identical(tab$cohort_label, rep("Independent cohort", 2L))
      for (i in seq_along(x$evaluations)) {
        saved <- x$evaluations[[i]]
        expect_identical(tab$estimate[i], saved$metrics[[metric]]$estimate)
        expect_identical(tab$status[i], saved$metrics[[metric]]$status)
        expect_identical(tab$run_status[i], saved$status)
        expect_identical(tab$scope[i], saved$metrics[[metric]]$scope)
        expect_identical(tab$aggregation[i], saved$metrics[[metric]]$aggregation)
        expect_identical(tab$task_type[i], type)
      }
      single <- sb_result_table(x$evaluations[[1]], metric = metric)
      expect_identical(single$cohort_label, NA_character_)
      expect_identical(single$estimate, tab$estimate[1])
    }
    expect_identical(x, before)
  }
})

test_that("external predictions and evaluation runs preserve their distinct meaning", {
  for (type in c("binary", "survival")) {
    x <- make_sb_curve_view_fixture(type, comparison = TRUE)$external
    for (i in seq_along(x$evaluations)) {
      saved <- x$evaluations[[i]]
      tab <- sb_result_table(saved, "predictions", cohort_label = "Test cohort")
      expect_identical(tab[names(saved$scored_predictions)], saved$scored_predictions)
      expect_identical(tab$task_type, rep(type, nrow(tab)))
      expect_identical(tab$cohort_label, rep("Test cohort", nrow(tab)))
      runs <- sb_result_table(saved, "runs")
      expect_identical(nrow(runs), 1L)
      expect_identical(runs$status, saved$status)
      expect_identical(runs$patient_match_status, saved$patient_match$status)
      expect_identical(runs$error_stage, saved$error_stage)
      expect_identical(runs$error_class, saved$error_class)
      expect_identical(runs$error_message, saved$error_message)
      expect_identical(runs$n, saved$n)
      expect_identical("fold" %in% names(runs), FALSE)
    }
    expected <- do.call(rbind, lapply(unname(x$evaluations), `[[`, "scored_predictions"))
    rownames(expected) <- NULL
    tab <- sb_result_table(x, "predictions")
    expect_identical(tab[names(expected)], expected)
    expect_identical(sb_result_table(x, "runs")$pipeline_id, names(x$evaluations))
  }
})

test_that("real external evaluation failure remains named without invented metrics", {
  d <- sb_example_binary(120, seed = 1)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  ext <- sb_example_binary(100, seed = 2)
  pred <- sb_predict_final(fit, ext[c("patient_id", "gene_a", "gene_b")], id_col = "patient_id")
  ep <- sb_binary_endpoint(ext, "patient_id", "outcome")
  good <- sb_evaluate_external(pred, ep)
  pred$pipeline_id <- "failed_id_match"
  pred$predictions$pipeline_id <- "failed_id_match"
  pred$predictions$patient_id[1] <- "unknown_patient"
  failed <- sb_evaluate_external(pred, ep)
  expect_identical(failed$status, "failed")
  x <- structure(list(evaluations = list(failed, good)), class = "stormbreaker_external_comparison")
  tab <- sb_result_table(x, metric = "roc_auc")
  expect_identical(tab$pipeline_id, c("failed_id_match", good$pipeline_id))
  expect_identical(tab$estimate, c(NA_real_, good$metrics$roc_auc$estimate))
  expect_identical(tab$status, c("not_saved", good$metrics$roc_auc$status))
  expect_identical(tab$run_status, c("failed", good$status))
  expect_identical(tab$scope[1], NA_character_)
  expect_identical(tab$aggregation[1], NA_character_)
  expect_identical(tab$error_stage[1], failed$error_stage)
  expect_identical(tab$error_message[1], failed$error_message)
  expect_identical(nrow(sb_result_table(failed, "predictions")), 0L)
  expect_identical(nrow(sb_result_table(x, "predictions")), nrow(good$scored_predictions))
})

test_that("external display is read-only and rejects unsupported requests", {
  x <- make_sb_curve_view_fixture("binary")$external
  before <- x
  local_mocked_bindings(
    sb_fit_final = function(...) stop("must not fit"),
    sb_evaluate_external = function(...) stop("must not evaluate"),
    sb_external_binary_auc_metric = function(...) stop("must not compute")
  )
  expect_identical(sb_result_table(x, metric = "roc_auc")$estimate, x$metrics$roc_auc$estimate)
  expect_identical(x, before)
  expect_error(sb_result_table(x, metric = "harrell_c"), class = "stormbreaker_result_metric_not_saved")
  expect_error(sb_result_table(x, metric = "calibration"), class = "stormbreaker_bad_result_metric")
  expect_error(sb_result_table(x, "runs", metric = "roc_auc"), class = "stormbreaker_unused_result_metric")
  for (label in list("", "  ", NA_character_, c("A", "B"), 1)) {
    expect_error(sb_result_table(x, metric = "roc_auc", cohort_label = label), class = "stormbreaker_bad_cohort_label")
  }
  expect_error(sb_result_table(list(x), metric = "roc_auc"), class = "stormbreaker_bad_result_view")
  expect_error(sb_result_table(structure(list(evaluations = list()), class = "stormbreaker_external_comparison"), metric = "roc_auc"), class = "stormbreaker_bad_result_view")
})
