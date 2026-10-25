test_that("nomograms keep instructions outside the plotted figure", {
  d <- sb_example_binary(80, seed = 11)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  nomogram <- sb_nomogram(fit, d)

  plot <- sb_plot_nomogram(nomogram)

  expect_null(plot$labels$caption)
  expect_null(plot$labels$subtitle)
  expect_identical(plot$labels$title, "Nomogram")
  expect_match(paste(attr(plot, "stormbreaker_notes"), collapse = " "),
    "Sum feature points", fixed = TRUE)
  expect_equal(sb_nomogram_values(nomogram, d)$probability,
    sb_predict_final(fit, d, "patient_id")$predictions$probability)
})
