## ---- specialty-lasso ----
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
