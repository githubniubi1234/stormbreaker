# Historical filename retained for compatibility. Use the unified modeling entry.
library(stormbreaker)
data <- bi_example(n = 100, seed = 42)
genes <- c("gene_a", "gene_b", "gene_c")
r <- bi_run(data, outcome = "outcome", id = "patient_id", predictors = genes,
            models = c("Lasso", "Logistic"), stages = 2, folds = 3, plots = FALSE)
finished <- finish_run(r, names(r$pipelines)[[1L]], plots = FALSE)
finished$genes
