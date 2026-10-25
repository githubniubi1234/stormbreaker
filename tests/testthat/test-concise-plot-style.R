test_that("shared figure themes keep explanatory prose outside the drawing", {
  skip_if_not_installed("ggplot2")

  expect_s3_class(sb_viz_theme()$plot.caption, "element_blank")
  expect_s3_class(sb_viz_theme()$plot.subtitle, "element_blank")
  expect_s3_class(sb_shap_theme()$plot.caption, "element_blank")
  expect_s3_class(sb_shap_theme()$plot.subtitle, "element_blank")
  expect_equal(sb_viz_theme()$text$size, 11)
})

test_that("compacting keeps explanations and numeric plot content", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(x = c(1, 2), y = c(.2, .8))
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y)) + ggplot2::geom_point() +
    ggplot2::labs(title = "Saved model coefficients", x = "Time (month)",
      subtitle = "Reference population", caption = "Conditional association")
  attr(p, "stormbreaker_notes") <- "Original interpretation"
  before <- ggplot2::ggplot_build(p)$data

  q <- sb_compact_plot(p)

  expect_identical(q$labels$title, "Model coefficients")
  expect_null(q$labels$subtitle)
  expect_null(q$labels$caption)
  expect_setequal(attr(q, "stormbreaker_notes"),
    c("Original interpretation", "Reference population", "Conditional association"))
  expect_identical(q$data, d)
  expect_identical(q$labels$x, "Time (month)")
  expect_identical(ggplot2::ggplot_build(q)$data, before)
})

test_that("short clinical labels retain cohort and fold distinctions", {
  models <- data.frame(id = "model_id", model = "Logistic regression")
  clean <- sb_clinician_labeler(models)

  expect_identical(clean("model_id"), "Logistic")
  expect_identical(clean("Development CV\nOut-of-fold"), "CV")
  expect_identical(clean("Validation A\nAssessment fold 2"), "Validation A\nFold 2")
  expect_identical(clean("Validation B\nAssessment fold 2"), "Validation B\nFold 2")
})

test_that("ROC AUC labels use colour instead of repeating model names", {
  skip_if_not_installed("ggplot2")
  run <- bi_run(sb_example_binary(80, seed = 42), "outcome",
    c("gene_a", "gene_b"), "patient_id", models = c("Logistic", "NaiveBayes"),
    folds = 3, plots = FALSE)
  original <- lapply(run$comparison$results, `[[`, "predictions")
  local_mocked_bindings(sb_fit_final = function(...) stop("Unexpected refit"))

  p <- sb_plot_roc(run)

  statistics <- attr(p, "stormbreaker_roc_statistics")
  expect_identical(statistics$label, paste0("AUC ", sprintf("%.3f", statistics$auc)))
  expect_identical(all.vars(p$layers[[3]]$mapping$colour), "pipeline_id")
  expect_identical(lapply(run$comparison$results, `[[`, "predictions"), original)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  views <- sb_clinician_views(run)
  expect_match(views$captions[match("ROC", names(views$plots))],
    "AUC labels describe the displayed curve", fixed = TRUE)
})

test_that("performance axes name the metric without instruction text", {
  skip_if_not_installed("ggplot2")
  run <- bi_run(sb_example_binary(80, seed = 42), "outcome",
    c("gene_a", "gene_b"), "patient_id", models = "Logistic", folds = 3,
    plots = FALSE)

  p <- sb_plot_performance(run$comparison, "roc_auc")

  expect_identical(p$labels$y, "ROC AUC")
  expect_identical(p$scales$get_scales("x")$get_labels(), "Logistic")
  expect_s3_class(p$theme$plot.caption, "element_blank")
  expect_s3_class(p$theme$plot.subtitle, "element_blank")
  expect_match(p$labels$caption, "No confidence intervals", fixed = TRUE)
})
