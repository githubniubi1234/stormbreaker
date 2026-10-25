# Getting started

Start with a small, explicitly named model set. STORMBREAKER compares binary or
survival pipelines, then continues a chosen pipeline through final fitting,
prediction, evaluation and plotting. Multiclass analysis has its own task interface.

## Prepare the data

Use one row per patient, a unique identifier and numeric predictor columns.
For binary prediction, code the outcome as 0/1 with 1 as the positive class.
For survival prediction, provide follow-up time, an event indicator (1 for an
observed event, 0 for censoring), the time unit and, when needed, a prediction horizon.
Align expression and clinical tables by identifier rather than row position.
Declare patient groups in a task when observations are repeated.

## Compare, choose and continue

Use bi_run() or surv_run() with a short model list for the first analysis.
The run's comparison contains cross-validation results; names(run$pipelines)
contains the exact names accepted by finish_run(). Select a route using
development evidence, then supply independently prepared cohorts to finish_run().
No winner is selected automatically by these entrypoints.

The completed run keeps the final fit, stage-wise feature tables, cohort
predictions, evaluations, plots and issues. fit_final() produces a reusable model;
predict_final() applies it without fitting again. If outcomes are unavailable,
predict only and evaluate later when the corresponding outcomes are available.

## Inspect results

Use result_table() on a comparison or evaluation object, clinical_plots() for
clinical evaluation views, and plots(fit) for model-specific diagnostics.
Read the availability/status tables when a requested plot is absent.
export_run() writes into a new directory; saveRDS()/readRDS() retain the final fit.

The [six runnable tutorials](scenarios-en.md) and the installed
stormbreaker-short-scenarios vignette demonstrate the complete object sequence.
See [data import](data-import.md), [model choices](user-models.md) and
[visualization](visualization-map.md) for details.

## Code examples

The following fragments use the input and object names described above. Supply your prepared data and the required engines before running them.

### Example 1

```r
r <- bi_run(data, outcome = "diagnosis", id = "patient_id",
            predictors = genes, models = c("Logistic", "Lasso", "RF"))
r$comparison$comparison
```

### Example 2

```r
route <- bi_pipeline(c("Lasso", "RF", "GBM"))
r <- bi_run(task, pipelines = list(Lasso_RF_GBM = route), folds = 5)
finished <- finish_run(r, "Lasso_RF_GBM", plots = FALSE)
finished$genes
export_run(finished, "results-new")
```

### Example 3

```r
r <- bi_run(task, models = c("Lasso", "RF", "GBM"),
            stages = c(1, 2, 3), folds = 5, plots = FALSE)
names(r$pipelines)
r$unavailable
```
