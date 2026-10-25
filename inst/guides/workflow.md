# From a cohort to a saved model

The development task defines the endpoint, candidate features, patient groups,
time units and optional horizons or decision thresholds. A pipeline defines
the ordered feature-selection steps and final learner; it is not a fitted model.

## Plan before comparing

Specify the primary metric, resampling scheme, candidate routes and final
selection rule before examining test outcomes. Feature selection belongs inside
the resampling procedure. Stability summaries describe the saved fold-specific
selections; they do not make a feature causal or guarantee future performance.

## Development and final evaluation

compare_cv() compares routes on shared folds. bi_run()/surv_run() provide a
short interface around the task-specific workflows. Omitting the model list
requests the full supported menu and can create many routes.

After reviewing development results, fit a declared pipeline and save the final
model. Predict external samples from their features, then evaluate against the
matching endpoint. Displayed CV, apparent training and external results have
different scopes and should remain separately identified.

## Objects and outputs

- run$comparison: cross-validation evidence.
- finished$fit: the final fitted pipeline.
- finished$genes: features retained at each stage.
- finished$predictions and finished$evaluations: cohort-specific outputs.
- run$issues and finished$issues: failures and unavailable components.
- fit$required_features: the columns needed for subsequent prediction.

Output directories and report paths must be new. A saved model may contain
patient-derived information; review its contents before sharing it.
The [six scenarios](scenarios-en.md) are the recommended complete runnable examples.

## Code examples

The following fragments use the input and object names described above. Supply your prepared data and the required engines before running them.

### Example 1

```r
library(stormbreaker)
source(system.file("examples", "user-workflow.R", package = "stormbreaker"))


demo <- workflow_demo("D:/my_stormbreaker/synthetic_demo_01")
```

### Example 2

```r
dev <- workflow_develop(
  development_csv = "D:/my_stormbreaker/input/binary_development.csv",
  output_dir = "D:/my_stormbreaker/binary_development_01",
  kind = "binary",
  predictors = c("gene_a", "gene_b", "gene_c"),
  folds = 5L,
  seed = 5201314L,
  near_tie = 0.01
)
dev$comparison
dev$suggested_candidate
```

### Example 3

```r
final <- workflow_finalize(
  development_rds = "D:/my_stormbreaker/binary_development_01/development.rds",
  selected_candidate = "logistic",  
  selection_reason = "Selected by the prespecified rule with calibration and stability review.",
  external_predictors_csv = "D:/my_stormbreaker/input/binary_external_predictors.csv",
  external_outcomes_csv = "D:/my_stormbreaker/input/binary_external_outcomes.csv",
  output_dir = "D:/my_stormbreaker/binary_final_01"
)
```

### Example 4

```r
fit <- readRDS("D:/my_stormbreaker/binary_final_01/final_model.rds")
new_data <- workflow_read_csv("D:/my_stormbreaker/input/new_patients.csv")
prediction <- predict_final(fit, new_data, id_col = "patient_id")
prediction$predictions
```
