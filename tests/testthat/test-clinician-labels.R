test_that("clinical labels hide internal status without changing model grouping", {
  skip_if_not_installed("ggplot2")
  fixture <- make_sb_curve_view_fixture("binary", comparison = TRUE)
  views <- sb_clinician_views(fixture$cv)
  clean <- sb_clinician_labeler(views$models)
  expect_identical(clean("Unnamed cohort\nRole: pipeline\nAssessment fold 2\ncurve: success run: success\nSaved entry: root/1"), "Fold 2")
  expect_identical(views$models$model, c("Logistic regression", "Gaussian naive Bayes"))
  expect_identical(clean(views$models$id), c("Logistic", "Gaussian NB"))
  expect_identical(grepl("smooth|none__", paste(clean(views$models$id), collapse = " ")), FALSE)
  expect_null(views$plots[["Confusion matrix"]])
  for (name in c("ROC", "Calibration", "Decision curve")) {
    p <- views$plots[[name]]
    rendered <- unlist(p$facet$params$labeller(data.frame(panel = unique(p$data$panel))))
    expect_identical(any(grepl("Unnamed cohort|Role:|curve: success|run: success|Saved entry:|none__", rendered)), FALSE)
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  }
  expect_identical(views$plots$ROC$data, sb_plot_roc(fixture$cv)$data)
})

test_that("same-family configurations get distinct readable labels", {
  leaf <- function(id) list(value = list(task_type = "binary",
    pipeline = list(id = id, learner = "logistic", selector = "none")))
  models <- sb_clinician_models(list(leaf("one"), leaf("two")))
  expect_length(unique(models$model), 2L)
  expect_identical(models$model, c("Logistic regression (1)", "Logistic regression (2)"))
  expect_identical(sb_clinician_labeler(models)(models$id), c("Logistic (1)", "Logistic (2)"))
})
