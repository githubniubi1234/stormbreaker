library(stormbreaker)
genes <- c("gene_a", "gene_b", "gene_c")

binary_data <- bi_example(n = 100, seed = 42)
binary_task <- bi_task(
  binary_data, "patient_id", "outcome", predictors = genes
)
binary_result <- sb_compare_models(
  binary_task, models = c("logistic", "gaussian_naive_bayes"),
  folds = 5, seed = 42
)
binary_result$comparison[c("model_id", "status", "roc_auc", "brier_score")]

survival_data <- surv_example(n = 100, seed = 42)
survival_task <- surv_task(
  survival_data, "patient_id", "time", "event",
  predictors = genes, time_unit = "month"
)
survival_result <- sb_compare_models(
  survival_task, models = c("coxph", "survreg_weibull"),
  folds = 5, seed = 42
)
survival_result$comparison[c("model_id", "status", "harrell_c")]
