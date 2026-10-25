make_sb_shap_plot_fixture <- function() {
  m <- rbind(c(a = 1, b = -0.3, c = 0.2), c(a = -0.5, b = 0.4, c = -0.1))
  rownames(m) <- c("patient-secret-1", "patient-secret-2")
  structure(list(status = "success", pipeline_id = "logistic", score_scale = "log_odds",
    reference_value = -0.4, contributions = m,
    scores = data.frame(patient_id = rownames(m), score = -0.4 + rowSums(m)),
    background = list(n = 10L)), class = "stormbreaker_additive_shap")
}

test_that("importance carries numeric labels without altering the ranked values", {
  x <- make_sb_shap_plot_fixture()
  p <- sb_plot_shap_importance(x)
  expect_identical(p$data$value_label, c("0.75", "0.35", "0.15"))
  expect_equal(p$data$mean_absolute_shap, unname(colMeans(abs(x$contributions))))
})

test_that("waterfall has a readable separate value column and unchanged endpoints", {
  x <- make_sb_shap_plot_fixture()
  p <- sb_plot_shap_waterfall(x, top_n = 1)
  expect_identical(p$data$value_label, c("-0.4", "+1", "-0.1", "0.5"))
  expect_equal(length(unique(p$data$label_x)), 1L)
  expect_gt(p$data$label_x[1], max(p$data$end))
  expect_equal(tail(p$data$end, 1), x$scores$score[1])
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

test_that("zero and tiny contributions remain finite and reconstruct exactly", {
  x <- make_sb_shap_plot_fixture()
  x$contributions[,] <- 0
  x$reference_value <- 0
  x$scores$score[] <- 0
  p <- sb_plot_shap_waterfall(x)
  expect_identical(all(is.finite(p$data$label_x)), TRUE)
  expect_identical(p$data$value_label, c("0", "+0", "+0", "+0", "0"))
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  x$contributions[1,1:2] <- c(1e-9, -1e-9)
  q <- sb_plot_shap_waterfall(x)
  expect_identical(q$data$value_label[2:3], c("+1e-09", "-1e-09"))
  expect_equal(tail(q$data$end, 1), 0)
})
