## ---- scenario-compare ----
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
