# Installation

## GitHub or release archives

Install the development version from GitHub, or use a source installation
archive or a Windows R 4.5 binary from the package Release attachments.

For installation from GitHub:

```r
install.packages("remotes") # only if missing
remotes::install_github("githubniubi1234/stormbreaker",
  dependencies = NA, build_vignettes = FALSE)
```

For local archive installation:

```r
# Run only if these required packages are missing:
install.packages(c("lpSolve", "rlang", "scales", "survival", "withr"))
# Select the supplied stormbreaker_0.3.47.zip on Windows R 4.5:
install.packages(file.choose(), repos = NULL, type = "win.binary")
# On other supported platforms select stormbreaker_0.3.47.tar.gz instead:
# install.packages(file.choose(), repos = NULL, type = "source")
```

R 4.3 or newer is required. stormbreaker has no compiled source. Built archives
include offline HTML tutorials, so readers do not need Pandoc. Tutorial builders
need knitr, rmarkdown, markdown and Pandoc. Missing optional engines are not installed automatically.

The six scenarios need ggplot2 for plots, patchwork for linked panels and glmnet
for scenario 5. Install only the extras you choose to use:

```r
install.packages(c("ggplot2", "patchwork", "glmnet")) # only if missing
```

Restart R after updating
a loaded package; check packageVersion("stormbreaker") and find.package("stormbreaker").
For a trial installation, pass lib=your_project_library and leave the global library unchanged.

## Add model engines

Start with `sb_start()` and install only the dependencies for the chosen workflow.
Do not use `dependencies = TRUE` to install every optional engine. For a small selected
binary comparison, add `glmnet`, `randomForest` and `ggplot2`; for a selected survival set,
add `glmnet`, `ranger` and `ggplot2`; for default multiclass, add `glmnet`,
`randomForest`, `e1071` and `ggplot2`. The screening example additionally uses
`Boruta` and `ranger`. Advanced engine installation remains optional below.

Use the catalog to see the packages for the models you want:

```r
catalog <- sb_model_catalog("binary")
catalog[c("model_id", "r_packages", "dependency_status")]
```

| Models / features | Packages | Installation |
| --- | --- | --- |
| Lasso, Ridge, elastic net | glmnet | `install.packages("glmnet")` |
| ranger forests | ranger | `install.packages("ranger")` |
| Boruta selection | Boruta, ranger | `install.packages(c("Boruta", "ranger"))` |
| varSelRF selection | varSelRF, randomForest | `install.packages(c("varSelRF", "randomForest"))` |
| FDA, penalized FDA, mixture discriminant analysis | mda | `install.packages("mda")` |
| XGBoost | xgboost | `install.packages("xgboost")` |
| LightGBM | lightgbm | `install.packages("lightgbm")` |
| Plots and cohort synthesis | ggplot2, metafor | `install.packages(c("ggplot2", "metafor"))` |

For Bioconductor engines such as `mixOmics`, `cancerclass` and `switchBox`, use
[BiocManager](https://www.bioconductor.org/install/), which selects the release
for your R version:

```r
install.packages("BiocManager")
BiocManager::install(c("mixOmics", "cancerclass", "switchBox"))
```

CatBoost provides [official R installation instructions](https://catboost.ai/docs/en/installation/r-installation-binary-installation).
The historical `fastAdaboost` engine is available through the
[CRAN archive](https://cran.r-project.org/src/contrib/Archive/fastAdaboost/).
The catalog and model-specific help identify additional engines.

Akritas also needs `distr6` and `survivalmodels`. distr6 is archived on CRAN;
its maintained R-universe distribution was used for local verification:

```r
install.packages("distr6", repos = c(
  "https://mlr-org.r-universe.dev", "https://cloud.r-project.org"))
install.packages("survivalmodels")
```

The distr6 installer also resolves required packages such as set6 and param6.
See the [distr6 repository](https://github.com/xoopR/distr6) and
[R-universe package page](https://mlr-org.r-universe.dev/distr6) for current platform builds.

## Saved models and engine updates

Save final models with `saveRDS()` and keep the engine version recorded in the
result provenance. Version 0.3.5 writes new XGBoost states in JSON, and supports
both the 1.x and current loading interfaces. Existing states written by older
stormbreaker versions retain their original XGBoost storage format; replay
those with the engine version used to create them before migrating a project.

## Input format

Each row represents a patient. Supply a unique patient ID and numeric predictor
columns. Binary outcomes use 0/1; survival tasks specify time, event and time
unit. The data-preparation and external-validation examples are in the
[workflow guide](workflow.md).
