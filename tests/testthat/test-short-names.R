test_that("short names are exact aliases with unchanged arguments", {
  pairs <- c(
    models = "sb_models", bi_task = "sb_binary_task", surv_task = "sb_surv_task",
    multi_task = "sb_multiclass_task", bi_endpoint = "sb_binary_endpoint",
    surv_endpoint = "sb_surv_endpoint", run_cv = "sb_run_cv", compare_cv = "sb_compare_cv",
    fit_final = "sb_fit_final", predict_final = "sb_predict_final",
    eval_external = "sb_evaluate_external", compare_external = "sb_compare_external",
    finish_run = "sb_complete_run", clinical_plots = "sb_clinician_views", plots = "sb_plots",
    result_table = "sb_result_table", selected_genes = "sb_selected_features",
    export_run = "sb_export_run", clinical_report = "sb_render_clinician_report",
    model_shap = "sb_model_shap", model_coefs = "sb_model_coefficients",
    bi_example = "sb_example_binary", surv_example = "sb_example_survival",
    plot_roc = "sb_plot_roc", plot_pr = "sb_plot_pr", plot_calib = "sb_plot_calibration",
    plot_dca = "sb_plot_dca", plot_perf = "sb_plot_performance",
    plot_heatmap = "sb_plot_performance_heatmap", plot_confusion = "sb_plot_confusion",
    plot_prob = "sb_plot_probability_distribution", plot_survival = "sb_plot_survival_curves",
    plot_nomogram = "sb_plot_nomogram", plot_coefs = "sb_plot_coefficients",
    plot_shap = "sb_plot_shap_beeswarm", plot_shap_importance = "sb_plot_shap_importance",
    plot_shap_dependence = "sb_plot_shap_dependence", plot_shap_waterfall = "sb_plot_shap_waterfall",
    plot_rank = "sb_plot_model_ranking", plot_cohorts = "sb_plot_cohort_metrics",
    plot_panel = "sb_plot_clinician_panel", plot_pages = "sb_plot_pages",
    plot_forest = "sb_plot_meta_forest", plot_selector = "sb_plot_selector_diagnostic",
    plot_algorithm = "sb_plot_algorithm_diagnostic"
  )
  exports <- getNamespaceExports("stormbreaker")
  expect_length(pairs, 45L)
  expect_contains(exports, names(pairs))
  expect_contains(exports, unname(pairs))
  for (short in names(pairs)) {
    old <- getExportedValue("stormbreaker", pairs[[short]])
    new <- getExportedValue("stormbreaker", short)
    expect_identical(new, old, info = short)
    expect_identical(formals(new), formals(old), info = short)
    expect_lt(nchar(short), nchar(pairs[[short]]))
  }
})

test_that("short workflow preserves predictions metrics and clinical coordinates", {
  skip_if_not_installed("ggplot2")
  development <- bi_example(120, seed = 42)
  validation <- bi_example(80, seed = 17)
  task <- bi_task(development, "patient_id", "outcome", c("gene_a", "gene_b"))
  expect_identical(task, sb_binary_task(development, "patient_id", "outcome", c("gene_a", "gene_b")))
  run <- bi_run(task, models = "Logistic", folds = 3, plots = FALSE)
  old <- sb_complete_run(run, "Logistic", external = list(validation = validation), plots = FALSE)
  new <- finish_run(run, "Logistic", external = list(validation = validation), plots = FALSE)
  expect_identical(new$predictions, old$predictions)
  expect_identical(new$metrics, old$metrics)
  expect_identical(result_table(new$evaluations$validation, metric = "roc_auc"),
    sb_result_table(old$evaluations$validation, metric = "roc_auc"))
  expect_identical(new$direction_correction$validation$reverse, TRUE)
  expect_identical(clinical_plots(new)$plots$ROC$data, sb_clinician_views(old)$plots$ROC$data)
  expect_identical(plot_roc(new)$data, sb_plot_roc(old)$data)
})

test_that("short survival workflow preserves raw risk and direction defaults", {
  data <- surv_example(80, seed = 42)
  genes <- c("gene_a", "gene_b")
  old_task <- sb_surv_task(data, "patient_id", "time", "event", genes, "month")
  new_task <- surv_task(data, "patient_id", "time", "event", genes, "month")
  expect_identical(new_task, old_task)
  route <- surv_pipeline("Cox")
  old_fit <- sb_fit_final(old_task, route)
  new_fit <- fit_final(new_task, route)
  expect_identical(new_fit$status, "success")
  query <- data[c("patient_id", genes)]
  old_prediction <- sb_predict_final(old_fit, query, "patient_id")
  new_prediction <- predict_final(new_fit, query, "patient_id")
  expect_identical(new_prediction$predictions, old_prediction$predictions)
  endpoint <- surv_endpoint(data, "patient_id", "time", "event", "month")
  old_evaluation <- sb_evaluate_external(old_prediction, endpoint)
  new_evaluation <- eval_external(new_prediction, endpoint)
  expect_identical(new_evaluation$metrics, old_evaluation$metrics)
  expect_null(new_evaluation$direction_correction)
})

test_that("short task names preserve validation errors", {
  data <- sb_example_binary(30, seed = 42)
  data$outcome[1] <- 2
  old <- tryCatch(sb_binary_task(data, "patient_id", "outcome", c("gene_a", "gene_b")), error = identity)
  new <- tryCatch(bi_task(data, "patient_id", "outcome", c("gene_a", "gene_b")), error = identity)
  expect_s3_class(old, "error")
  expect_s3_class(new, "error")
  expect_identical(class(new), class(old))
  expect_identical(conditionMessage(new), conditionMessage(old))
})
