# Iris is a public demonstration dataset, not clinical validation evidence.
library(stormbreaker)
d <- datasets::iris
d$id <- seq_len(nrow(d))
task <- multi_task(d, "id", "Species", names(d)[1:4],
  classes = c("virginica", "setosa", "versicolor"))
methods <- c("lda", "qda", "gaussian_naive_bayes", "decision_tree")
analysis <- analyze_methods(task, methods)
analysis$status
sb_method_features(analysis, "gaussian_naive_bayes", class = "setosa")
subset(sb_method_features(analysis, "decision_tree"), used_for_split)
analysis$results$gaussian_naive_bayes$plots$class_statistics
sb_predict_multiclass(analysis$results$qda$fit, d[1:3, ])
comparison <- sb_compare_multiclass(task, methods, folds = 3)
comparison$comparison
# No automatic files, reports or confusion plots.
