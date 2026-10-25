# Workflow short names

Short names are aliases of the original public functions: they do not define
additional algorithms or change model settings. Use a qualified name such as
stormbreaker::fit_final() if another package has an overlapping name.

The complete mapping is in [preferred-name-map.csv](preferred-name-map.csv). The CSV records each short
name, original name, function group and English description. Native diagnostic
availability still depends on the fitted model and retained state.

Use plot_diag() for a diagnostic object, plots(fit) for the overall algorithm
view, and clinical_plots() for clinical evaluation views. A selector's diagnostic
must be requested from that stage rather than confused with the final learner.

## Code examples

The following fragments use the input and object names described above. Supply your prepared data and the required engines before running them.

### Example 1

```r
library(stormbreaker)
d <- bi_example(120, seed = 42)
r <- bi_run(d, outcome = "outcome", id = "patient_id",
            predictors = c("gene_a", "gene_b"), models = "Logistic")
clinical_plots(r)$plots$ROC
```

### Example 2

```r
task <- bi_task(development, "patient_id", "outcome", genes)
route <- bi_pipeline(c("mRMR", "xgboost_linear"))
fit <- fit_final(task, route, capture_diagnostics = TRUE)
prediction <- predict_final(fit, validation[c("patient_id", genes)], "patient_id")
evaluation <- eval_external(prediction, bi_endpoint(validation, "patient_id", "outcome"))
clinical_plots(evaluation)$plots$ROC
```
