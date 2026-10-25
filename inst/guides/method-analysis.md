# Analyze methods independently

analyze_methods() analyzes each chosen method on the same candidate features.
Unlike an ordered pipeline, it does not pass one method's selected genes to
the next. Use bi_task(), surv_task() or multi_task() to declare the endpoint and
column roles; inspect sb_method_catalog() for task support.

selected_genes() returns a method's selected set when such a rule is defined.
sb_method_features() returns the full feature evidence, including retained and
unretained genes. Random forests provide importance and rankings, not an
automatically selected set. Request the ranking table for those results.

| Method | Evidence | Selection rule |
|---|---|---|
| Boruta | Importance summaries, hit rates and history | Confirmed features after RoughFix |
| varSelRF | Subset size and OOB error | Smallest subset within one standard error of the minimum |
| VSURF | Importance and stage sizes | Prediction-stage subset |
| Lasso / Elastic Net | Coefficient path and CV | Nonzero coefficients at the selected lambda |
| MCP / SCAD | Coefficient path and BIC | Nonzero coefficients of the minimum-BIC fit |
| ABESS | Subset size and GIC | GIC-selected subset |
| Random forest | Importance and tree size; binary OOB where available | Ranking only |

Lasso defaults to lambda.1se; Elastic Net uses alpha=0.5. Their internal CV uses
five folds and inner_seed=1701 in this independent-analysis interface.
Method-specific settings, evidence and diagnostic tables remain in the result.
Boruta summaries show medians and ranges, not boxplots. Survival forest diagnostics
use the supported ranger quantities rather than an invented binary OOB trajectory.

A failure in one method is recorded without discarding subsequent methods.
Plotting failures retain fitted results and appear in plot_status.
These are exploratory fits: internal CV and OOB curves are diagnostics, not
out-of-sample validation of an entire preselected gene workflow. Use pipeline CV
when evaluating feature selection followed by prediction.

## Code examples

The following fragments use the input and object names described above. Supply your prepared data and the required engines before running them.

### Example 1

```r
task <- bi_task(
  data = expression_data,
  id_col = "patient_id",
  outcome_col = "outcome",
  predictors = gene_names
)

analysis <- analyze_methods(
  task,
  methods = c("Boruta", "varSelRF", "VSURF", "LASSO", "Elastic Net",
              "MCP", "SCAD", "ABESS", "random_forest")
)
analysis$status
```

### Example 2

```r
selected_genes(analysis, "LASSO")
selected_genes(analysis, "Boruta")
sb_method_features(analysis, "LASSO")
sb_method_features(analysis, "random_forest")
```

### Example 3

```r
analysis$results$lasso$plots$coefficient_path
analysis$results$lasso$plots$cv_curve
analysis$results$boruta$plots$boruta_history
analysis$results$random_forest$plots$oob_error

ggplot2::ggsave(
  "lasso-path.pdf",
  analysis$results$lasso$plots$coefficient_path,
  width = 7, height = 5
)
```
