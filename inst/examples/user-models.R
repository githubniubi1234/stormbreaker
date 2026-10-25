# User-controlled models and genes. These functions do not infer the endpoint.
# Replace paths, endpoint columns and genes with YOUR predeclared study inputs.
# Lines below are illustrative, not executed automatically by the package.

library(stormbreaker)

# Read your own gene list: one exact predictor column name per line.
genes <- readLines("D:/my_study/my_genes.txt", encoding = "UTF-8")
development <- read.csv("D:/my_study/development.csv", check.names = FALSE,
  colClasses = c(patient_id = "character"))

# Choose ONE task definition yourself:
task <- bi_task(development, "patient_id", "outcome", predictors = genes)
# For survival use this INSTEAD:
# task <- surv_task(development, "patient_id", "time", "event",
#   predictors = genes, time_unit = "month")

sb_model_catalog("binary") # use "survival" for a survival task

# Specific requirements: run just these models, on the same user gene set.
compared <- sb_compare_models(task,
  models = c("logistic", "ridge_logistic", "ranger_binary"),
  folds = 5L, seed = 2026L)

# No model preference: use this INSTEAD to run the complete task-specific menu.
# compared <- sb_compare_models(task, models = NULL, folds = 5L, seed = 2026L)

compared$comparison
compared$elapsed_seconds
# Use new file paths for each run; do not replace past patient results.
# write.csv(compared$comparison, "D:/my_study/run_01/comparison.csv", row.names = FALSE)
# saveRDS(compared, "D:/my_study/run_01/comparison.rds")

# After development-only review, explicitly choose a compared specification.
# selected <- "logistic"  # example, NOT an automatically selected winner
# final_fit <- fit_final(task, compared$model_pipelines[[selected]])
# saveRDS(final_fit, "D:/my_study/run_01/model.rds")
# Predict external predictors first, then supply outcomes for external evaluation;
# see guides/workflow.md for the outcome-blind sequence.
