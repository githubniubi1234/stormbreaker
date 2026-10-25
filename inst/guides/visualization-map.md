# Choose a result view

Use clinical_plots() for saved clinical evaluation results and plots(fit) for
algorithm diagnostics. Neither changes the fitted learner.

| Question | Functions | Required evidence |
|---|---|---|
| Binary discrimination | plot_roc(), plot_pr() | Saved predictions and outcomes |
| Calibration | plot_calib() | Saved calibration coordinates; a horizon for survival |
| Clinical utility | plot_dca() | Prespecified decision thresholds |
| Fixed-threshold errors | plot_cm(), plot_thresholds() | A declared threshold |
| Model performance | plot_perf(), plot_heatmap() | Comparable saved metrics |
| Cohort uncertainty | plot_ext_ci(), plot_forest() | Requested bootstrap or meta-analysis results |
| Survival risk groups | surv_risk() | Development risk reference and follow-up outcomes |
| Horizon-specific discrimination | surv_roc() | Survival task, scores and requested times |
| Feature stability | plot_freq(), plot_jaccard(), plot_upset() | Saved resampling selections |
| Native model evidence | plots(fit), native_diag() | Supported retained diagnostics |
| Prediction contributions | model_shap(), plot_shap(), plot_shap_one() | A supported fit and declared background |
| Linear score display | model_coefs(), nomogram(), plot_nomogram() | Original model coefficients and reference ranges |
| Multiclass errors | plot_multi_cm() | Multiclass results |

Check returned status/availability tables rather than assuming every plot exists
for every model. CV, apparent-training and external views have different scopes.
Do not pool different horizons or cohort roles into a single performance value.

## Code examples

The following fragments use the input and object names described above. Supply your prepared data and the required engines before running them.

### Example 1

```r
r <- bi_run(data, "outcome", predictors = genes, id = "patient_id")
v <- clinical_plots(r)
```

### Example 2

```r
library(stormbreaker)
d <- bi_example(120, seed = 42)
genes <- c("gene_a", "gene_b")

r <- bi_run(d, "outcome", predictors = genes, id = "patient_id",
            models = c("Logistic", "NaiveBayes"),
            decision_thresholds = c(.1, .3, .5))
v <- clinical_plots(r)
v$plots$ROC
v$plots$Calibration
v$plots[["Decision curve"]]
```

### Example 3

```r
v <- clinical_plots(r)  


finished <- finish_run(r, "Logistic", external = list(validation = external))
v_external <- clinical_plots(finished)
v_external$plots$ROC
v_external$plots$Calibration
```

### Example 4

```r
r_surv <- surv_run(data, time = "time", event = "event",
  predictors = genes, id = "patient_id", models = "Cox",
  time_unit = "month", prediction_horizon = 36,
  decision_thresholds = c(.1, .3, .5))
v_surv <- clinical_plots(r_surv)
v_surv$plots$Calibration
v_surv$plots[["Decision curve"]]
```
