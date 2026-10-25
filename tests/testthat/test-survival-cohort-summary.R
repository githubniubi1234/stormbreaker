test_that("cohort bars preserve supplied metric values", {
  d <- data.frame(cohort = c("A", "B"), c_index = c(.68, .75))
  p <- sb_plot_cohort_metrics(d)
  expect_equal(p$data$value, d$c_index)
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})

test_that("risk meta agrees with native Cox and inverse-variance estimates", {
  withr::local_seed(102)
  cohorts <- lapply(1:3, function(i) {
    risk <- rnorm(120)
    d <- data.frame(patient_id = paste(i, 1:120), risk = risk,
      time = rexp(120, exp(.5 * risk)), event = rbinom(120, 1, .8),
      group = ifelse(risk > 0, "High", "Low"))
    list(patients = d, cutoff = 0)
  })
  names(cohorts) <- c("A", "B", "C")
  for (type in c("group", "continuous")) {
    m <- surv_meta_risk(cohorts, type)
    expect_identical(m$status, "success")
    e <- m$estimates
    expect_equal(m$pooled$log_hr[1], weighted.mean(e$log_hr, 1/e$se^2))
    expect_equal(sum(e$weight_random), 100)
    expect_s3_class(ggplot2::ggplotGrob(m$plot), "gtable")
  }
  d <- cohorts[[1]]$patients
  native <- survival::coxph(survival::Surv(time, event) ~ risk, d)
  expect_equal(sb_risk_cox_estimate(d, "continuous")$log_hr, unname(coef(native)))
  cohorts[[1]]$patients$event <- 0
  expect_identical(surv_meta_risk(cohorts)$status, "not_evaluable")
  cohorts[[1]] <- cohorts[[2]]
  expect_error(surv_meta_risk(cohorts), "overlap")
})
