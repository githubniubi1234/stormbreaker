# Tutorials

These six scripts demonstrate the main STORMBREAKER workflows using the packaged
synthetic CSV files. Each script returns `example_result` and `figures` without
automatically exporting files or opening a browser.

## Run a tutorial

```r
library(stormbreaker)
source(system.file("examples", "scenarios", "01-binary.R",
  package = "stormbreaker", mustWork = TRUE))
metrics
figures$ROC
```

Install `ggplot2` and `patchwork` for plotting and `glmnet` for the Lasso diagnostic
example. The scripts do not install dependencies.

## Choose a workflow

| Script | Workflow | Inputs | Main outputs |
|---|---|---|---|
| `01-binary.R` | Binary classification | 120-row development cohort | Cross-validation results, AUC, ROC, calibration and decision curves |
| `02-compare.R` | Model comparison | Development cohort and specified pipelines | Shared-fold comparison table, performance plot and ROC |
| `03-external.R` | External evaluation | Development cohort and two validation cohorts | One final model, cohort predictions, metrics and ROC |
| `04-survival.R` | Survival prediction | Development and validation cohorts; time in months | Risk scores, survival probabilities, concordance, time-dependent AUC, Brier score, calibration and decision curves |
| `05-diagnostics.R` | Feature selection and diagnostics | Development cohort; Lasso followed by logistic regression | Retained genes, selection diagnostics and model coefficients |
| `06-shap.R` | Prediction interpretation | A fitted model, eight query samples and development background data | SHAP contributions, importance, beeswarm, dependence and waterfall plots |

Replace the script name in the call above to run another tutorial.

## Input data

- Binary: `patient_id`, `outcome`, `gene_a`, `gene_b`, `gene_c`. The outcome is 0/1;
  1 is the positive class.
- Survival: `patient_id`, `time`, `event`, `gene_a`, `gene_b`, `gene_c`. An event is
  coded as 1 and censoring as 0. Tutorial times are in months.
- Prediction only: the identifier and the features required by the saved model;
  outcomes are not needed.

Each row represents a patient. Transpose a gene-by-sample expression matrix and
align it to the clinical identifiers before constructing a task. For repeated
observations, use unique row identifiers and declare `group_col` in the task.

Data generation, seeds and units are described in
[the example-data notes](../extdata/demo-data.txt).

## Working with results

`run$comparison` holds cross-validation results. A completed workflow contains
`finished$fit`, cohort-specific raw `finished$predictions` and
`finished$evaluations`. Pass a comparison or evaluation object to `result_table()`:

```r
result_table(finished$evaluations$Validation_A, metric = "roc_auc")
```

`clinical_plots()` returns clinical evaluation plots; `plots(fit)` returns
model-specific diagnostics. The final model's required features are available in
`fit$required_features`.

