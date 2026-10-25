## ---- scenario-diagnostics ----
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
