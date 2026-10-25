# Multiclass analysis

Use multi_task() to declare the outcome classes and predictors, analyze_methods()
for independent method fits, and sb_compare_multiclass() for predictive validation.
These are multiclass interfaces rather than a Cartesian extension of the binary
menu. Class labels in the examples are illustrative and must match your data.

selected_genes() can summarize supported nonzero-feature sets across classes;
sb_method_features() retains class-specific coefficients or evidence. Forest
importance provides rankings rather than an automatically chosen gene set.
sb_predict_multiclass() returns per-class probabilities for supported saved fits.

The supported method catalog includes task-specific regularized regression,
forest, SVM, discriminant, naive Bayes, tree and boosting implementations.
Inspect availability and dependency status rather than assuming all registered
methods can fit a given dataset. LDA/QDA assumptions and small-class restrictions
remain relevant. Native diagnostics differ across algorithms.

Mutual-information extractors such as mRMR, JMI, CMIM, DISR and JMIM are separate
selection methods. Feature filtering and native tuning must be repeated within
the relevant analysis partitions when evaluating a predictive pipeline.
Independent full-data analysis is not a replacement for that validation.

Use confusion matrices only when class predictions and outcomes are available;
keep multiclass metrics separate from binary discrimination measures.
The fragments below retain the individual supported interfaces and method IDs.

## Code examples

The following fragments use the input and object names described above. Supply your prepared data and the required engines before running them.

### Example 1

```r
task <- multi_task(
  data, id_col = "patient_id", outcome_col = "subtype",
  predictors = genes,
  classes = c("Control", "Subtype_A", "Subtype_B")
)
```

### Example 2

```r
analysis <- analyze_methods(task, c("LASSO", "Elastic Net", "random_forest"))
analysis$status


selected_genes(analysis, "LASSO")


selected_genes(analysis, "LASSO", class = "Subtype_A")
sb_method_features(analysis, "LASSO", class = "Subtype_A")


sb_method_features(analysis, "random_forest", class = "Subtype_A")

analysis$results$lasso$plots$coefficient_path
analysis$results$lasso$plots$cv_curve
analysis$results$random_forest$plots$class_importance
```

### Example 3

```r
p <- sb_predict_multiclass(analysis$results$lasso$fit, new_patients)
```

### Example 4

```r
comparison <- sb_compare_multiclass(
  task, c("lasso", "elastic_net", "random_forest"),
  folds = 5, inner_folds = 5
)
comparison$comparison
comparison$results$lasso$metrics$per_class
```

### Example 5

```r
plot_multi_cm(comparison$results$lasso$metrics)
```

### Example 6

```r
methods <- c("lda", "qda", "gaussian_naive_bayes", "decision_tree")
analysis <- analyze_methods(task, methods)
comparison <- sb_compare_multiclass(task, methods)

sb_method_features(analysis, "lda", class = "Subtype_A")
analysis$results$gaussian_naive_bayes$plots$class_statistics
subset(sb_method_features(analysis, "decision_tree"), used_for_split)
```

### Example 7

```r
analysis <- analyze_methods(task, c("boruta", "varselrf", "vsurf", "abess"))
analysis$status
selected_genes(analysis, "boruta")
selected_genes(analysis, "abess", class = "Subtype_A")
analysis$results$boruta$plots$boruta_history
analysis$results$varselrf$plots$selection_history
analysis$results$vsurf$plots$vsurf_stages
analysis$results$abess$plots$support_criterion
analysis$results$abess$plots$class_coefficients


comparison <- sb_compare_multiclass(task, "abess", folds = 5)
```

### Example 8

```r
methods <- c("ridge", "svm_linear", "svm_radial", "svm_radial_weighted",
  "svm_polynomial", "svm_sigmoid", "kernel_naive_bayes", "rda", "nnet", "knn",
  "ranger", "ranger_extratrees", "c50", "c50_rules", "ctree", "cforest",
  "bagging", "rrf", "tree", "evtree")
analysis <- analyze_methods(task, methods)
comparison <- sb_compare_multiclass(task, methods, folds = 5, inner_folds = 5)
```

### Example 9

```r
methods <- c("fda", "mda", "pda", "slda", "sda", "loclda", "pamr",
  "gaussian_process", "kernel_weighted_knn", "liblinear_l2_logistic",
  "liblinear_l1_logistic", "xgboost", "xgboost_dart", "xgboost_linear",
  "lightgbm", "lightgbm_dart", "lightgbm_goss", "lightgbm_rf", "catboost")
analysis <- analyze_methods(task, methods)
p <- sb_predict_multiclass(analysis$results$xgboost$fit, new_patients)
comparison <- sb_compare_multiclass(task, methods, folds = 5)
```

### Example 10

```r
analysis <- analyze_methods(task, c("mrmr", "jmi", "cmim", "disr", "jmim"))
selected_genes(analysis, "cmim")
sb_method_features(analysis, "jmi")
analysis$results$jmi$plots$selection_scores
```
