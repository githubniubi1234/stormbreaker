# Public iris data demonstrate the interface; these are not clinical results.
library(stormbreaker)
d <- datasets::iris
d$sample_id <- seq_len(nrow(d))
task <- multi_task(d, "sample_id", "Species", names(d)[1:4],
  classes = c("virginica", "setosa", "versicolor"))

analysis <- analyze_methods(task, c("lasso", "elastic_net", "random_forest"),
  inner_folds = 3)
analysis$status
selected_genes(analysis, "lasso", class = "setosa")
analysis$results$lasso$plots$coefficient_path
sb_predict_multiclass(analysis$results$lasso$fit, d[1:3, ])

comparison <- sb_compare_multiclass(task, c("lasso", "elastic_net", "random_forest"),
  folds = 3, inner_folds = 3)
comparison$comparison
comparison$results$lasso$metrics$per_class
# Optional only: plot_multi_cm(comparison$results$lasso$metrics)
