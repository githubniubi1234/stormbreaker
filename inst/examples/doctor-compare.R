library(stormbreaker)

data <- bi_example(n = 100, seed = 42)
genes <- c("gene_a", "gene_b", "gene_c")
task <- bi_task(data, "patient_id", "outcome", predictors = genes)


run <- bi_run(task, models = c("Logistic", "Lasso", "RF"), folds = 3, plots = FALSE)
comparison <- run$comparison
comparison$comparison[c("status", "roc_auc", "brier_score")]

views <- clinical_plots(comparison)
print(views$plots$ROC)


route <- bi_pipeline(c("Lasso", "RF", "GBM"))
# bi_run(task, pipelines = list(Lasso_RF_GBM = route), folds = 3)
# bi_run(task, models = c("Lasso", "RF", "GBM"), stages = c(1, 2, 3))
# finished <- finish_run(run, "RF", plots = FALSE)
# finished$genes
