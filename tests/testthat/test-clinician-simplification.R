test_that("ordinary model comparison defaults to three declared models", {
  local_mocked_bindings(
    sb_compare_cv = function(task, pipelines, ...) list(comparison = data.frame(
      pipeline_id = vapply(pipelines, function(p) p$id, character(1)))),
    sb_attach_result_page = function(x, ...) x)
  for (type in c("binary", "survival")) {
    d <- if (type == "binary") sb_example_binary(60) else sb_example_survival(60)
    task <- if (type == "binary") sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b")) else
      sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month")
    expected <- if (type == "binary") c("logistic", "glmnet_logistic", "random_forest_binary") else
      c("coxph", "glmnet_cox", "ranger_survival")
    result <- sb_compare_models(task, view = FALSE)
    expect_identical(result$requested_models, expected)
    expect_identical(sb_compare_models(task, models = NULL, view = FALSE)$requested_models,
      sb_model_catalog(type)$model_id)
  }
  expect_identical(formals(sb_compare_multiclass)$methods,
    quote(c("lasso", "random_forest", "svm_radial")))
})

test_that("brief clinical views do not compute optional plots", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_curve_view_fixture("binary")$cv
  before <- serialize(x, NULL)
  full <- sb_clinician_views(x, detail = "full")
  local_mocked_bindings(sb_plot_pr = function(...) stop("Not requested"),
    sb_probability_distribution = function(...) stop("Not requested"))
  brief <- sb_clinician_views(x)
  expect_identical(brief$status$view, c("ROC", "Calibration", "Decision curve"))
  expect_identical(brief$plots$ROC$data, full$plots$ROC$data)
  expect_identical(serialize(x, NULL), before)
  expect_identical(all(c("Precision-recall", "Probability distribution") %in% full$status$view), TRUE)
})
