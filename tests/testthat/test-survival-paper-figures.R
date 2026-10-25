paper_survival_task <- function() {
  withr::local_seed(32)
  d <- as.data.frame(matrix(rnorm(180 * 6), 180, 6))
  names(d) <- paste0("g", 1:6)
  d$id <- paste0("p", seq_len(nrow(d)))
  d$time <- rexp(nrow(d), exp(d$g1 - d$g2 + .7 * d$g3))
  d$event <- rbinom(nrow(d), 1, .8)
  sb_surv_task(d, "id", "time", "event", paste0("g", 1:6), "month")
}

test_that("C-index ranking reuses complete-cohort means without AUC labels", {
  d <- data.frame(model = rep(c("A", "B", "C"), each = 2), cohort = rep(c("V1", "V2"), 3),
    c_index = c(.8, .7, .9, .5, .99, NA))
  z <- sb_plot_model_ranking(d, c("V1", "V2"), metric = "c_index")
  expect_equal(z$ranking$mean_c_index, c(.75, .7, NA))
  expect_identical(z$ranking$status, c("ranked", "ranked", "incomplete"))
  expect_identical(z$heatmap$labels$fill, "C-index")
  expect_identical(z$bar$labels$x, "Mean C-index")
  expect_s3_class(ggplot2::ggplot_build(z$bar), "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(z$heatmap), "ggplot_built")
})

test_that("survival Enet matches native coefficients and keeps patient groups", {
  task <- paper_survival_task()
  x <- task$data[task$predictors]
  groups <- rep(1:90, each = 2)
  adapter <- sb_resolve_selector_adapter("enet_cox_nonzero", list(inner_folds = 3L, seed = 1L))
  s <- adapter$fit(x, task$data$time, task$data$event, group_ids = groups)
  native <- glmnet::glmnet(as.matrix(x), survival::Surv(task$data$time, task$data$event),
    family = "cox", cox.ties = "breslow", alpha = .5, standardize = TRUE)
  expect_equal(s$scores$coefficient, as.numeric(stats::coef(native, s = s$lambda)))
  expect_equal(s$alpha, .5)
  expect_identical(all(vapply(split(s$fold_id, groups), function(g) length(unique(g)) == 1L, logical(1))), TRUE)
  expect_identical(surv_selector("Enet"), surv_selector("ElasticNet"))
  for (p in list(surv_pipeline("Cox", surv_selector("Enet")),
    surv_pipeline("ElasticNet", surv_selector("Enet")),
    surv_pipeline("Cox", surv_selector(c("Enet", "StepCox"))))) {
    expect_identical(sb_lint_pipeline(p)$status, "admitted")
    fit <- sb_fit_final(task, p)
    expect_identical(fit$status, "success")
    expect_identical(sb_run_cv(task, p, folds = 3)$status, "success")
  }
})

test_that("risk panels share order and never recalibrate on validation outcomes", {
  task <- paper_survival_task()
  fit <- sb_fit_final(task, surv_pipeline("Cox", surv_selector("Enet")))
  ref <- surv_risk_reference(fit, task)
  z <- surv_risk_plots(ref, task)
  sm <- summary(z$km, censored = TRUE)
  expect_equal(z$plots$km$data$lower[z$plots$km$data$time > 0], sm$lower)
  expect_equal(z$plots$km$data$upper[z$plots$km$data$time > 0], sm$upper)
  expect_equal(z$cox_group$log_hr, unname(stats::coef(survival::coxph(
    survival::Surv(time, event) ~ group, z$patients))))
  scores <- sb_risk_predict(fit, task)
  expect_equal(ref$cutoff, stats::median(scores))
  expect_identical(z$patients$patient_id, task$data$id[order(scores)])
  expect_identical(rownames(z$expression), z$patients$patient_id)
  expected <- scale(as.matrix(task$data[order(scores), ref$genes, drop = FALSE]), center = ref$center, scale = ref$scale)
  expect_equal(unname(z$expression), unname(expected), ignore_attr = TRUE)
  test <- survival::survdiff(survival::Surv(time, event) ~ group, data = z$patients)
  expect_equal(z$logrank_p, stats::pchisq(test$chisq, 1, lower.tail = FALSE))
  expect_equal(z$at_risk$n_risk, mapply(function(t, g) sum(z$patients$time >= t & z$patients$group == g),
    z$at_risk$time, z$at_risk$group))
  altered <- task$data
  altered$time <- rev(altered$time)
  altered$event <- rev(altered$event)
  validation <- sb_surv_task(altered, "id", "time", "event", task$predictors, "month")
  changed <- surv_risk_plots(ref, validation)
  expect_identical(changed$patients[c("patient_id", "risk", "group")], z$patients[c("patient_id", "risk", "group")])
  expect_identical(changed$expression, z$expression)
  for (p in z$plots) expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  single <- ref; single$cutoff <- Inf
  expect_identical(surv_risk_plots(single, task)$status, "single_group")
})

test_that("survival ROC handles ties, censoring, ID order and unavailable times", {
  d <- data.frame(id = letters[1:5], time = c(1, 2, 3, 5, 6), event = c(1, 0, 1, 1, 0), x = 1:5)
  task <- sb_surv_task(d, "id", "time", "event", "x", "month")
  risk <- c(3, 0, 1, 1, 2)
  z <- surv_plot_roc(task, risk, c(.5, 4, 7))
  expect_identical(z$statistics$status, c("not_evaluable", "success", "not_evaluable"))
  expect_equal(z$statistics$n_early_censored[2], 1)
  # At t=4 case weights are 1 and 4/3; the second case earns 1/4 pair credit.
  expect_equal(z$statistics$auc[2], (1 + (4/3) * .25) / (1 + 4/3))
  shuffled <- setNames(rev(risk), rev(d$id))
  expect_equal(surv_plot_roc(task, shuffled, 4)$statistics$auc, z$statistics$auc[2])
  tied <- surv_plot_roc(task, rep(1, 5), 4)
  expect_equal(tied$statistics$auc, .5)
  q <- sb_time_auc_stratum(d$time, d$event, plogis(-risk), rep(TRUE, 5),
    sb_fit_censoring_km(d$time, d$event), 4)
  expect_equal(z$statistics$auc[2], q$estimate)
  expect_s3_class(ggplot2::ggplot_build(z$plot), "ggplot_built")
  expect_error(surv_plot_roc(task, risk, c(4, 4)), "unique positive")
})
