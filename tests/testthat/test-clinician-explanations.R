test_that("clinical reports reuse exact saved explanations and declared dependence features", {
  skip_if_not_installed("ggplot2")
  d <- sb_example_binary(n = 80L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  pipeline <- sb_pipeline(learner = "logistic")
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 41L)
  fit <- sb_fit_final(task, pipeline)
  shap <- sb_additive_shap(fit, d, "patient_id", d)
  nomogram <- sb_nomogram(fit, d)
  before <- serialize(list(fit, shap, nomogram), NULL, version = 2L)
  values <- d[rev(seq_len(nrow(d))), ]
  views <- sb_clinician_views(result, shap = shap, feature_values = values,
    id_col = "patient_id", dependence_features = "gene_b", nomogram = nomogram)
  expect_identical(serialize(list(fit, shap, nomogram), NULL, version = 2L), before)
  expect_identical(views$plots[["SHAP importance"]]$data, sb_plot_shap_importance(shap)$data)
  expect_identical(views$plots[["SHAP beeswarm"]]$data,
    sb_plot_shap_beeswarm(shap, feature_values = values, id_col = "patient_id")$data)
  expect_identical(views$plots[["SHAP dependence: gene_b"]]$data,
    sb_plot_shap_dependence(shap, "gene_b", values, "patient_id")$data)
  expect_identical(views$plots[["Original-model nomogram"]]$data, sb_plot_nomogram(nomogram)$data)
  expect_null(views$plots[["SHAP dependence: gene_a"]])
  expect_length(views$height, length(views$plots))
  expect_length(views$width, length(views$plots))
  expect_true(all(views$status$status[views$status$group == "Model explanations"] == "available"))
  file <- tempfile(fileext = ".html")
  sb_render_clinician_report(views, file, complete = TRUE)
  html <- paste(readLines(file), collapse = "\n")
  expect_true(grepl("Original-model nomogram", html, fixed = TRUE))
  expect_true(grepl("SHAP dependence: gene_b", html, fixed = TRUE))
  expect_false(grepl("patient card", html, fixed = TRUE))
  broken <- values
  broken$patient_id[1L] <- "wrong-id"
  failed <- sb_clinician_views(result, shap = shap, feature_values = broken,
    id_col = "patient_id", dependence_features = "gene_b")
  expect_identical(failed$status$status[failed$status$view == "SHAP dependence: gene_b"], "error")
  expect_null(failed$plots[["SHAP dependence: gene_b"]])
  expect_error(sb_clinician_views(result, shap = shap, dependence_features = "gene_b"),
    class = "stormbreaker_bad_report_input")
  expect_error(sb_clinician_views(result, nomogram = fit), class = "stormbreaker_bad_report_input")
})
