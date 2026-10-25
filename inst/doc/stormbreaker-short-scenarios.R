## ----setup, include=FALSE-----------------------------------------------------
library(stormbreaker)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message = FALSE,
  error = FALSE, fig.width = 7, fig.height = 4.8, dpi = 130, dev = "png",
  dev.args = list(type = "cairo"))
options(markdown.HTML.options = c("toc", "base64_images"))
scripts <- c("01-binary.R", "02-compare.R", "03-external.R", "04-survival.R",
  "05-diagnostics.R", "06-shap.R")
for (file in scripts) {
  knitr::read_chunk(system.file("examples", "scenarios", file,
    package = "stormbreaker", mustWork = TRUE))
}
knitr::read_chunk(system.file("examples", "specialty-short-names.R",
  package = "stormbreaker", mustWork = TRUE))

## ----input-binary, echo=FALSE-------------------------------------------------
binary_input <- read.csv(system.file("extdata", "demo-binary-development.csv",
  package = "stormbreaker"), check.names = FALSE)
knitr::kable(head(binary_input, 3), digits = 3)

## ----scenario-binary----------------------------------------------------------
library(stormbreaker)

development <- read.csv(system.file("extdata", "demo-binary-development.csv",
  package = "stormbreaker"), check.names = FALSE)
genes <- c("gene_a", "gene_b", "gene_c")

run <- bi_run(development, outcome = "outcome", id = "patient_id",
  predictors = genes, models = "Logistic", folds = 3, seed = 42,
  decision_thresholds = c(.1, .2, .3, .5, .7), plots = FALSE)
comparison <- run$comparison
metrics <- result_table(comparison, metric = "roc_auc")
figures <- list(
  ROC = plot_roc(comparison),
  Calibration = plot_calib(comparison),
  DCA = plot_dca(comparison)
)
example_result <- list(input = development, run = run,
  comparison = comparison, metrics = metrics, figures = figures)

if (interactive()) print(figures$ROC)

## ----binary-output, echo=FALSE------------------------------------------------
knitr::kable(metrics[c("pipeline_id", "metric", "estimate", "status", "scope")], digits = 3)
figures$ROC

## ----scenario-compare---------------------------------------------------------
library(stormbreaker)

development <- read.csv(system.file("extdata", "demo-binary-development.csv",
  package = "stormbreaker"), check.names = FALSE)
genes <- c("gene_a", "gene_b", "gene_c")
task <- bi_task(development, "patient_id", "outcome", genes)
routes <- list(
  Logistic = bi_pipeline("Logistic"),
  NaiveBayes = bi_pipeline("NaiveBayes")
)
comparison <- compare_cv(task, routes, folds = 3, seed = 42, progress = FALSE)
metrics <- result_table(comparison, metric = "roc_auc")
figures <- list(
  Performance = plot_perf(comparison, metric = "roc_auc"),
  ROC = plot_roc(comparison)
)
example_result <- list(input = development, task = task, routes = routes,
  comparison = comparison, metrics = metrics, figures = figures)

if (interactive()) print(figures$Performance)

## ----compare-output, echo=FALSE-----------------------------------------------
knitr::kable(metrics[c("pipeline_id", "metric", "estimate", "status")], digits = 3)
figures$Performance

## ----scenario-external--------------------------------------------------------
library(stormbreaker)

development <- read.csv(system.file("extdata", "demo-binary-development.csv",
  package = "stormbreaker"), check.names = FALSE)
validation_a <- read.csv(system.file("extdata", "demo-binary-external-a.csv",
  package = "stormbreaker"), check.names = FALSE)
validation_b <- read.csv(system.file("extdata", "demo-binary-external-b.csv",
  package = "stormbreaker"), check.names = FALSE)
genes <- c("gene_a", "gene_b", "gene_c")

run <- bi_run(development, outcome = "outcome", id = "patient_id",
  predictors = genes, models = "Logistic", folds = 3, seed = 42,
  decision_thresholds = c(.1, .2, .3, .5, .7), plots = FALSE)
finished <- finish_run(run, "Logistic",
  external = list(Validation_A = validation_a, Validation_B = validation_b),
  plots = TRUE)
views <- clinical_plots(finished)
metrics <- do.call(rbind, lapply(names(finished$evaluations), function(cohort) {
  result_table(finished$evaluations[[cohort]], metric = "roc_auc", cohort_label = cohort)
}))
figures <- list(ROC = views$plots$ROC)
example_result <- list(input = development,
  external = list(Validation_A = validation_a, Validation_B = validation_b),
  run = run, finished = finished, views = views, metrics = metrics, figures = figures)

if (interactive()) print(figures$ROC)

## ----external-output, echo=FALSE----------------------------------------------
knitr::kable(head(finished$predictions$Validation_A[c("patient_id", "probability")], 3), digits = 3)
knitr::kable(metrics[c("cohort_label", "metric", "estimate", "status")], digits = 3)
figures$ROC

## ----input-survival, echo=FALSE-----------------------------------------------
survival_input <- read.csv(system.file("extdata", "demo-survival-development.csv",
  package = "stormbreaker"), check.names = FALSE)
knitr::kable(head(survival_input, 3), digits = 3)

## ----scenario-survival--------------------------------------------------------
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

## ----survival-output, echo=FALSE----------------------------------------------
knitr::kable(head(prediction$predictions[c("patient_id", "risk", "survival_probability",
  "prediction_horizon", "time_unit")], 3), digits = 3)
knitr::kable(metrics[c("metric", "estimate", "status")], digits = 3)
knitr::kable(horizon_metrics[c("metric", "estimate", "status", "prediction_horizon", "time_unit")], digits = 3)
figures$Calibration

## ----scenario-diagnostics, eval=requireNamespace("glmnet", quietly=TRUE)------
library(stormbreaker)

development <- read.csv(system.file("extdata", "demo-binary-development.csv",
  package = "stormbreaker"), check.names = FALSE)
genes <- c("gene_a", "gene_b", "gene_c")
task <- bi_task(development, "patient_id", "outcome", genes)
route <- bi_pipeline(c("Lasso", "Logistic"))
fit <- fit_final(task, route, capture_diagnostics = TRUE)
diagnostics <- plots(fit)
retained_genes <- data.frame(gene = fit$required_features)
coefficients <- model_coefs(fit, reference_data = development[genes])
figures <- c(diagnostics$plots, list(Coefficients = plot_coefs(coefficients)))
example_result <- list(input = development, task = task, route = route, fit = fit,
  genes = retained_genes, diagnostics = diagnostics,
  coefficients = coefficients, figures = figures)

if (interactive()) print(figures$Coefficients)

## ----diagnostics-output, echo=FALSE, eval=requireNamespace("glmnet", quietly=TRUE)----
knitr::kable(retained_genes)
knitr::kable(data.frame(plot = names(figures)))
figures$selector_1_cv

## ----specialty-lasso, eval=requireNamespace("glmnet", quietly=TRUE)-----------
library(stormbreaker)

development <- read.csv(system.file("extdata", "demo-binary-development.csv",
  package = "stormbreaker"), check.names = FALSE)
genes <- c("gene_a", "gene_b", "gene_c")
task <- bi_task(development, "patient_id", "outcome", genes)
fit <- fit_final(task, bi_pipeline("Lasso"), capture_diagnostics = TRUE)
cv <- path_diag(fit, type = "cv_curve")
path <- path_diag(fit, type = "coefficient_path")
native <- native_diag(fit, complete = TRUE)
figures <- list(CV = plot_diag(cv), Path = plot_diag(path))
example_result <- list(input = development, task = task, fit = fit,
  cv = cv, path = path, native = native, figures = figures)

if (interactive()) print(figures$CV)

## ----specialty-output, echo=FALSE, eval=requireNamespace("glmnet", quietly=TRUE)----
knitr::kable(head(cv$table, 3), digits = 4)
figures$CV

## ----scenario-shap------------------------------------------------------------
library(stormbreaker)

development <- read.csv(system.file("extdata", "demo-binary-development.csv",
  package = "stormbreaker"), check.names = FALSE)
validation <- read.csv(system.file("extdata", "demo-binary-external-a.csv",
  package = "stormbreaker"), check.names = FALSE)
genes <- c("gene_a", "gene_b", "gene_c")

task <- bi_task(development, "patient_id", "outcome", genes)
fit <- fit_final(task, bi_pipeline("Logistic"))
query <- validation[1:8, c("patient_id", genes)]
background <- development[genes]
shap <- model_shap(fit, query, "patient_id", background_data = background,
  permutations = 16, background_n = 12, seed = 42)
figures <- list(
  Importance = plot_shap_imp(shap),
  Beeswarm = plot_shap(shap, feature_values = query, id_col = "patient_id"),
  Dependence = plot_shap_dep(shap, "gene_a", query, "patient_id"),
  Waterfall = plot_shap_one(shap, row = 1)
)
example_result <- list(input = development, task = task, fit = fit,
  query = query, background = background, shap = shap, figures = figures)

if (interactive()) print(figures$Beeswarm)

## ----shap-output, echo=FALSE--------------------------------------------------
knitr::kable(head(shap$contributions, 3), digits = 4)
figures$Beeswarm
figures$Waterfall

