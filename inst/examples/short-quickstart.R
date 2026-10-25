library(stormbreaker)

development <- bi_example(120, seed = 42)
validation <- bi_example(80, seed = 17)
genes <- c("gene_a", "gene_b")

run <- bi_run(development, outcome = "outcome", predictors = genes,
              id = "patient_id", models = "Logistic", plots = FALSE)
finished <- finish_run(run, "Logistic", external = list(validation = validation))

views <- clinical_plots(finished)
roc <- views$plots$ROC
metrics <- result_table(finished$evaluations$validation, metric = "roc_auc")

if (interactive()) print(roc)
