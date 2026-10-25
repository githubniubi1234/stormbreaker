## ---- scenario-shap ----
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
