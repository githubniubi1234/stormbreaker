# Three-class example. Replace iris with a numeric gene/clinical feature table.
library(stormbreaker)
data <- iris
data$patient_id <- seq_len(nrow(data))
genes <- names(data)[1:4]
task <- multi_task(data, "patient_id", "Species", predictors = genes)

# Independent full-data feature exploration, not held-out performance.
analysis <- analyze_methods(task, c("Boruta", "varSelRF", "VSURF", "ABESS"))
analysis$status
selected_genes(analysis, "Boruta")
selected_genes(analysis, "ABESS")
analysis$results$abess$plots$support_criterion

# Each learner receives the original candidate features inside its training fold.
# Do not substitute genes selected on the full data above.
comparison <- sb_compare_multiclass(task,
  c("svm_radial", "ranger", "xgboost", "lightgbm", "catboost"), folds = 3)
comparison$comparison

# Fit one chosen model and predict from its saved native fit.
model <- analyze_methods(task, "ranger")
probability <- sb_predict_multiclass(model$results$ranger$fit, data[1:3, ])
probability
# This last prediction uses training rows only to demonstrate the API.
# For clinical use, pass genuinely new patients with the same predictor columns.
