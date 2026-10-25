test_that("Cox forest matches native Cox and retains unsuccessful genes", {
  d <- sb_example_survival(100, 42)
  d$constant <- 1
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "constant"), "month")
  result <- sb_cox_forest(task)
  native <- survival::coxph(survival::Surv(time, event) ~ gene_a, data = d)
  row <- result$table[1, ]
  expect_equal(row$hr, unname(exp(stats::coef(native))))
  expect_equal(c(row$lower, row$upper), unname(as.numeric(exp(stats::confint(native)))))
  expect_equal(result$table$p_adjusted, stats::p.adjust(result$table$p_value, "BH"))
  expect_identical(result$table$status[3], "failed")
  expect_identical(result$table$plot_available[3], FALSE)
  expect_s3_class(result$plot, "ggplot")
  expect_equal(nrow(ggplot2::ggplot_build(result$plot)$data[[3]]), 2)
})
