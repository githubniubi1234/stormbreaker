# STORMBREAKER

STORMBREAKER is an R package for transcriptomic prediction with binary
classification and right-censored survival outcomes. It connects a diverse model
library and built-in parameter configurations with feature selection, model
comparison, final fitting, external evaluation and interpretation.

## Features

- **Model diversity:** choose from 92 binary and 88 survival prediction
  configurations.
- **Built-in settings:** compare statistical and machine-learning models using
  preset parameters and shared cross-validation folds.
- **Ordered pipelines:** combine feature selection and prediction through a
  common modeling interface.
- **Connected workflow:** fit a selected pipeline, save the fitted model and
  predict new cohorts.
- **Evaluation and interpretation:** assess discrimination, calibration and
  clinical utility; explore selected features, model-specific diagnostics and
  SHAP explanations.
- **Results and reports:** export predictions, feature lists, evaluation tables
  and plots, and create HTML reports.

Use `models("binary")` or `models("survival")` to browse the prediction menus.

## Installation

STORMBREAKER requires R 4.3 or later. To install from GitHub:

```r
# install.packages("remotes")
remotes::install_github("githubniubi1234/stormbreaker",
  dependencies = NA, build_vignettes = FALSE)
```

Required dependencies are installed with the package. Install additional model
engines and plotting packages for the methods you choose. The six tutorial
scripts use `ggplot2`, `patchwork` and `glmnet`.

## Quick start

Start with a sample-by-feature table containing numeric predictors and clinical
outcomes.

### Binary classification

The package includes synthetic development and validation cohorts. Each row
represents a patient; `outcome` is coded as 0/1, with 1 as the positive class.

```r
library(stormbreaker)

development <- read.csv(system.file(
  "extdata", "demo-binary-development.csv", package = "stormbreaker"))
genes <- c("gene_a", "gene_b", "gene_c")

binary <- bi_run(
  development, outcome = "outcome", id = "patient_id", predictors = genes,
  models = c("Logistic", "NaiveBayes"), folds = 3, seed = 42, plots = FALSE
)
result_table(binary$comparison, metric = "roc_auc")
```

Set `models` to select the candidates for each comparison.

### External evaluation

Choose a pipeline using the development results, then fit it and evaluate an
independent cohort:

```r
validation <- read.csv(system.file(
  "extdata", "demo-binary-external-a.csv", package = "stormbreaker"))

finished <- finish_run(
  binary, "Logistic", external = list(Validation = validation),
  plots = FALSE
)
result_table(finished$evaluations$Validation, metric = "roc_auc")
```

### Survival prediction

Survival data include a follow-up time and an event indicator: 1 for an observed
event and 0 for censoring. Specify the time unit explicitly.

```r
survival_data <- read.csv(system.file(
  "extdata", "demo-survival-development.csv", package = "stormbreaker"))

survival_run <- surv_run(
  survival_data, time = "time", event = "event", id = "patient_id",
  predictors = genes, time_unit = "month", models = "Cox",
  folds = 3, seed = 42, plots = FALSE
)
result_table(survival_run$comparison, metric = "harrell_c")
```

### Save and export

Save the final fitted model for later prediction and export the completed run:

```r
model_file <- tempfile(fileext = ".rds")
saveRDS(finished$fit, model_file)
fit <- readRDS(model_file)
prediction <- predict_final(fit, validation[c("patient_id", genes)], "patient_id")

export_run(finished, "stormbreaker-results")
```

`export_run()` writes CSV tables and available PDF plots to a new directory.
Create HTML reports with `sb_render_report()`.

## Documentation

- [Six runnable tutorials](inst/guides/scenarios-en.md): binary classification,
  model comparison, external evaluation, survival prediction, diagnostics and SHAP.
- The installed vignette: `vignette("stormbreaker-short-scenarios")`.
- [Getting started](inst/guides/getting-started.md) and
  [advanced workflows](inst/guides/advanced-workflows-en.md).
- [Example data](inst/extdata/demo-data.txt): input columns, seeds and time units.
- Function help in R, for example `?bi_run`, `?surv_run` and `?finish_run`.
- [Changelog](NEWS.md) and [release notes](RELEASE_NOTES.md).

## Citation

To cite the package, run:

```r
citation("stormbreaker")
```

## License

MIT. See [LICENSE.md](LICENSE.md) and the third-party notices in
[inst/COPYRIGHTS](inst/COPYRIGHTS).
