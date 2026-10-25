## ---- scenario-survival ----
library(stormbreaker)

development <- read.csv(system.file("extdata", "demo-survival-development.csv",
  package = "stormbreaker"), check.names = FALSE)
validation <- read.csv(system.file("extdata", "demo-survival-external.csv",
  package = "stormbreaker"), check.names = FALSE)
genes <- c("gene_a", "gene_b", "gene_c")
task <- surv_task(development, "patient_id", "time", "event", genes,
  time_unit = "month", prediction_horizon = 6,
  decision_thresholds = c(.1, .2, .3, .5, .7))
route <- surv_pipeline("Cox")
cv <- run_cv(task, route, folds = 3, seed = 42, progress = FALSE)
fit <- fit_final(task, route)
prediction <- predict_final(fit, validation[c("patient_id", genes)], "patient_id")
endpoint <- surv_endpoint(validation, "patient_id", "time", "event",
  time_unit = "month", prediction_horizon = 6)
evaluation <- eval_external(prediction, endpoint)
metrics <- result_table(evaluation, metric = "harrell_c", cohort_label = "Validation")
horizon_metrics <- do.call(rbind, lapply(c("time_auc", "ipcw_brier"), function(metric) {
  result_table(evaluation, metric = metric, cohort_label = "Validation")
}))
figures <- list(Calibration = plot_calib(evaluation), DCA = plot_dca(evaluation))
example_result <- list(input = development, external = validation, task = task,
  cv = cv, fit = fit, prediction = prediction, endpoint = endpoint,
  evaluation = evaluation, metrics = metrics, horizon_metrics = horizon_metrics, figures = figures)

if (interactive()) print(figures$Calibration)
