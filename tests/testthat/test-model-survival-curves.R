test_that("fixed Efron baseline agrees with native Cox and original predictions", {
  d <- make_sb_fixture(n = 72L, seed = 9701L)
  d$time <- ceiling(d$time)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b"), "month",
    prediction_horizon = 6
  )
  fit <- sb_fit_final(task, sb_pipeline())
  augmented <- sb_add_survival_curve(fit, task)
  expect_identical(augmented$learner_state, fit$learner_state)
  expect_identical(
    sb_predict_final(augmented, d, "patient_id"),
    sb_predict_final(fit, d, "patient_id")
  )
  native <- survival::coxph(survival::Surv(time, event) ~ gene_a + gene_b, d)
  base <- survival::basehaz(native, centered = FALSE)
  at <- augmented$survival_curve$table
  expect_equal(at$cumulative_hazard,
    c(0, base$hazard)[findInterval(at$time, base$time) + 1L],
    tolerance = 1e-7
  )
  curves <- sb_survival_curves(augmented, d[1:3, ], "patient_id")
  expect_equal(dim(curves$survival), c(nrow(at), 3L))
  expect_equal(curves$survival[1, ], rep(1, 3))
  expect_lte(max(apply(curves$survival, 2, diff)), 0)
  expect_equal(
    sb_survival_curves(augmented, d[1, ], "patient_id", times = 6)$survival,
    curves$survival[at$time == 6, 1, drop = FALSE]
  )
  skip_if_not_installed("ggplot2")
  expect_s3_class(sb_plot_survival_curves(curves), "ggplot")
  expect_s3_class(sb_plot_survival_curves(curves, rows = integer()), "ggplot")
  sparse <- sb_survival_curves(augmented, d[1:3, ], "patient_id", times = 6)
  expect_equal(sb_plot_survival_curves(sparse)$data, sb_plot_survival_curves(curves)$data)
})

test_that("Breslow recovery equals an independent risk-set calculation", {
  d <- make_sb_fixture(n = 55L, seed = 9715L)
  d$time <- seq_len(nrow(d))
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 20
  )
  fit <- sb_fit_final(task, sb_pipeline())
  # With unique event times the native Efron and Breslow estimates coincide.
  fit$learner_state$baseline_hazard$ties <- "breslow"
  augmented <- sb_add_survival_curve(fit, task)
  risk <- drop(as.matrix(d[fit$required_features]) %*% fit$learner_state$coefficient)
  t <- d$time[d$event == 1]
  jumps <- vapply(t, function(point) 1 / sum(exp(risk[d$time >= point])), numeric(1))
  expected <- vapply(
    augmented$survival_curve$table$time,
    function(point) sum(jumps[t <= point]), numeric(1)
  )
  expect_equal(augmented$survival_curve$table$cumulative_hazard, expected, tolerance = 1e-12)
})

test_that("wrong development and extrapolation are refused", {
  d <- make_sb_fixture(n = 60L, seed = 105L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 6
  )
  fit <- sb_fit_final(task, sb_pipeline())
  other <- d
  other$patient_id <- paste0("other", other$patient_id)
  wrong <- sb_surv_task(other, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 6
  )
  expect_error(sb_add_survival_curve(fit, wrong), class = "stormbreaker_curve_development_mismatch")
  augmented <- sb_add_survival_curve(fit, task)
  expect_error(sb_survival_curves(augmented, d, "patient_id", max(d$time) + 1),
    class = "stormbreaker_curve_time_unsupported"
  )
  changed <- d
  changed$gene_a <- changed$gene_a + 10
  changed_task <- sb_surv_task(changed, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 6
  )
  expect_error(sb_add_survival_curve(fit, changed_task), class = "stormbreaker_curve_baseline_mismatch")
})
