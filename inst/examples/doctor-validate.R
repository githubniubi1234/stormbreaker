library(stormbreaker)


development <- bi_example(n = 100, seed = 42)
genes <- c("gene_a", "gene_b", "gene_c")
task <- bi_task(development, "patient_id", "outcome", predictors = genes)
fit <- fit_final(task, bi_pipeline("Logistic"))


external <- bi_example(n = 80, seed = 43)
external$patient_id <- paste0("external_", external$patient_id)
prediction <- predict_final(fit, external[c("patient_id", genes)], "patient_id")
endpoint <- bi_endpoint(external, "patient_id", "outcome")
evaluation <- eval_external(prediction, endpoint)
evaluation$metrics


