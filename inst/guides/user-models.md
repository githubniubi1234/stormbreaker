# Choose models and candidate genes

Choose the task explicitly with bi_task(), surv_task() or multi_task(), and supply
the candidate feature names. The task is not inferred merely from a 0/1 column.
Gene names must match the prepared data exactly.

sb_model_catalog() lists the registered model IDs and dependency status, including
entries whose optional engines are absent. Lasso, Ridge and Elastic Net use the
documented task-specific configurations; their common glmnet alpha settings are
1, 0 and 0.5 respectively. This catalog is not every possible hyperparameter
combination. bi_run()/surv_run() request their complete menus when models is omitted;
the older sb_compare_models() interface retains its own smaller default set.

Declare the comparison metric and candidate list before inspecting test outcomes.
sb_compare_models() keeps failed models and unavailable metrics in its results.
Choose one route explicitly for final fitting; the largest development score
alone does not establish external performance.

Risk-ranking capability and fixed-horizon survival-probability capability differ.
Do not substitute a risk score for a probability when a calibration or decision
curve requires probabilities. Inspect the model's saved prediction capabilities.

## Code examples

The following fragments use the input and object names described above. Supply your prepared data and the required engines before running them.

### Example 1

```r
genes <- readLines("D:/my_study/my_genes.txt", encoding = "UTF-8")
task <- bi_task(data, "patient_id", "outcome", predictors = genes)

# task <- surv_task(data, "patient_id", "time", "event",
#                     predictors = genes, time_unit = "month")
```

### Example 2

```r
menu <- sb_model_catalog("binary") 
menu[c("model_id", "learner", "glmnet_alpha", "dependency_status")]
```

### Example 3

```r
result <- sb_compare_models(task,
  models = c("logistic", "ridge_logistic", "ranger_binary"),
  folds = 5L, seed = 2026L)
```

### Example 4

```r
result <- sb_compare_models(task, folds = 5L, seed = 2026L)


```

### Example 5

```r
result$comparison           
result$elapsed_seconds      
result$input_genes          
result$requested_models     
result$results             
```

### Example 6

```r
selected <- "logistic" 
fit <- fit_final(task, result$model_pipelines[[selected]])
```
