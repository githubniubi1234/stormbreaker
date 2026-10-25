test_that("dependence aligns IDs and preserves all contributions", {
  skip_if_not_installed("ggplot2")
  m <- matrix(c(1, 2, -1, -2), 2, dimnames = list(c("secret1", "secret2"), c("a", "b")))
  x <- structure(list(status = "success", pipeline_id = "logistic", score_scale = "log_odds",
    reference_value = 0, contributions = m, scores = data.frame(patient_id = rownames(m), score = rowSums(m)),
    background = list(n = 10L)), class = "stormbreaker_additive_shap")
  values <- data.frame(id = c("secret2", "secret1"), a = c(20, 10), b = c(4, 3))
  before <- x
  p <- sb_plot_shap_dependence(x, "a", values, "id", colour_by = "b")
  expect_equal(p$data$feature_value, c(10, 20))
  expect_equal(p$data$contribution, c(1, 2))
  expect_equal(p$data$colour_value, c(3, 4))
  expect_identical(x, before)
  expect_identical(any(grepl("secret", unlist(p$data))), FALSE)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})
