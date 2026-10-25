# Import your own cohort

Each expression table has one sample per row, a patient_id column and numeric
gene columns. A separate clinical table may be supplied with the same identifiers
in a different order: sb_read_cohort() matches by ID and preserves textual IDs
such as 001. A combined expression/clinical table can be read directly.

Binary outcomes are coded as 0/1. Survival outcomes require time, event and an
explicit, consistent unit. Do not merge distinct endpoints such as overall and
disease-free survival. The packaged demo CSVs show the supported format.

Specify candidate genes before task construction. Prepare missing values,
expression scales, platform mapping and repeated-patient decisions upstream.
Missing genes are not replaced with zeros.

## Development, validation and test cohorts

sb_compare_cohorts() can compare routes on development CV and designated
validation cohorts. sb_select_model() explicitly selects an eligible route by
the declared validation metric; multiple validation cohorts receive equal
weight and ties follow the original declaration order. Test cohorts must not
participate in selection. Validation outcomes are not added to the training fit.

When a plan is fixed, save the model and predict new samples from their IDs and
features. Read outcomes separately to construct an endpoint and evaluate the
saved predictions. Keep patient-level models and predictions in appropriate storage.

## Common input failures

- Missing or extra clinical IDs: inspect the mismatches rather than forcing a
  row-wise join or silently dropping samples.
- Duplicate IDs: determine whether records are repeated measurements or technical
  replicates; do not disguise one patient as several independent patients.
- Missing predictors: check gene names and platform mapping.
- No eligible model: inspect the recorded fitting and metric failures rather
  than choosing a replacement from the final test set.
- Long runs: request a small pool first; forest selection and nested tuning are
  real computations whose cost depends on the data and settings.

## Code examples

The following fragments use the input and object names described above. Supply your prepared data and the required engines before running them.

### Example 1

```r
library(stormbreaker)
genes <- c("EGFR", "TP53") 

train <- sb_read_cohort("Training_expression.csv", "Training_clinical.csv")
valid <- sb_read_cohort("Validation_expression.csv", "Validation_clinical.csv")
train_task <- bi_task(train, "patient_id", "outcome", genes)
valid_task <- bi_task(valid, "patient_id", "outcome", genes)

compared <- sb_compare_cohorts(
  train_task, list(Validation = valid_task),
  roles = c(Validation = "validation"),
  models = c("logistic", "glmnet_logistic", "ridge_logistic"),
  folds = 5, seed = 919, inner_folds = 3, inner_seed = 917,
  progress = TRUE
)
model <- sb_select_model(compared)
model$model_selection$ranking
saveRDS(model, "final_model.rds") 
```

### Example 2

```r
train_task <- surv_task(train, "patient_id", "time", "event", genes,
  time_unit = "month", prediction_horizon = 60)
valid_task <- surv_task(valid, "patient_id", "time", "event", genes,
  time_unit = "month", prediction_horizon = 60)
```

### Example 3

```r
library(stormbreaker)
model <- readRDS("final_model.rds")
new_data <- sb_read_cohort("Test_expression.csv") 
prediction <- predict_final(model, new_data, "patient_id")
if (!identical(prediction$status, "success")) stop(prediction$error_message)
write.csv(prediction$predictions, "test_predictions.csv", row.names = FALSE)
saveRDS(prediction, "test_prediction_object.rds")
```
