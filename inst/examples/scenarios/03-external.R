## ---- scenario-external ----
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
