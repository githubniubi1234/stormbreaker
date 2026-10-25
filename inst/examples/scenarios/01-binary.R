## ---- scenario-binary ----
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
