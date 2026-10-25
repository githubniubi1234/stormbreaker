# Interpret a saved model

Coefficient views, SHAP contributions and predicted survival curves answer
different questions. They use an existing fit rather than selecting a new gene
set or fitting a surrogate model.

## Coefficients and SHAP

model_coefs() retains original-scale coefficients. The reference_sd view multiplies
them by the original reference standard deviations for display; it does not change
the prediction model. Zero coefficients remain in the table. Logistic intercepts
use log-odds units; absolute Cox probabilities depend on baseline information.

additive_shap() provides exact additive-model contributions on the supported
linear score scale: log-odds for Logistic and linear log-risk for Cox.
model_shap() instead supplies the documented approximate permutation explanation
for supported saved prediction states. Check score_scale; contributions are not
automatically probability percentage-point effects. Use a declared development
background and reuse the computed explanations for plots.

## Survival curves

add_surv_curve() requires the original unchanged development task to recover a
supported complete baseline while retaining the existing coefficients. It checks
the saved identifiers, time unit and baseline point. New-data prediction then
reads only identifiers and features. Validation outcomes must not be used to
reconstruct a training baseline.

The baseline corresponds to x beta=0, not an average patient. Curves are model
predictions, not Kaplan-Meier estimates; no unsupported confidence bands are
invented and the display does not extrapolate beyond observed training follow-up.

## Reports

sb_render_model_report() writes an offline report using the saved model,
reference features and requested samples. Larger cohorts generate larger reports.
svglite is optional. Use a new output path and review patient-derived contents
before sharing; display-row numbers are not an anonymization procedure.

## Code examples

The following fragments use the input and object names described above. Supply your prepared data and the required engines before running them.

### Example 1

```r
parameters <- model_coefs(model, reference_data = training_features)
parameters$table
parameters$intercept
plot_coefs(parameters)
plot_coefs(parameters, scale = "reference_sd")
```

### Example 2

```r
shap <- additive_shap(model, new_features, "patient_id",
  background_data = training_features)
plot_shap(shap, feature_values = new_features, id_col = "patient_id")
plot_shap_imp(shap)
plot_shap_one(shap, row = 1)
```

### Example 3

```r
model_with_curve <- add_surv_curve(model, original_development_task)
saveRDS(model_with_curve, "new_model_with_curve.rds") 
curves <- surv_curves(model_with_curve, new_features, "patient_id")
plot_surv(curves, rows = 1)
plot_surv(curves, rows = integer(0)) 
```

### Example 4

```r
sb_render_model_report(model, new_features, "patient_id", training_features,
  file = "new_binary_model_explanation.html")

sb_render_model_report(model_with_curve, new_features, "patient_id", training_features,
  file = "new_survival_model_explanation.html")
```
