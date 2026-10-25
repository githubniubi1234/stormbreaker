# Prepare figures from saved results

Use saved result objects for ROC, calibration, decision curves, model ranking
heatmaps, coefficient displays and nomograms. Keep apparent training, development
CV and external evaluation clearly distinguished in the figure caption.

plot_roc() uses saved out-of-fold predictions for CV results; by_fold is an
optional diagnostic view. External uncertainty intervals require the corresponding
resampling calculation. A heatmap uses one compatible value per model/cohort;
unavailable values remain NA rather than zero.

Nomograms represent supported original linear models, not surrogate fits.
Bootstrap calibration and decision curves require their declared reference,
horizon and thresholds. Use plotting metadata and returned notes for captions.
The functions generate new figures from results; they do not reproduce another
paper's hidden data, undocumented model settings or image pixels.

## Code examples

The following fragments use the input and object names described above. Supply your prepared data and the required engines before running them.

### Example 1

```r

pred <- predict_final(fit, task$data[c(task$id_col, task$predictors)], task$id_col)
p_train <- plot_pred_roc(
  task$data[[task$outcome_col]], pred$predictions$probability,
  kind = "apparent", cohort = "Training", ci_reps = 1000
)


p_external <- plot_roc(external_evaluation, ci_reps = 1000)


p_oof <- plot_roc(cv_result)

```

### Example 2

```r


ranking <- plot_rank(
  performance, ranking_cohorts = c("Validation1", "Validation2"), top_n = 20
)
ranking$plot
ranking$ranking
```

### Example 3

```r
nom <- nomogram(fit, task$data)
p_nom <- plot_nomogram(nom, reference_data = task$data,
  patient = task$data[1, , drop = FALSE])

cal <- sb_bootstrap_calibration(task, pipeline, reps = 1000, seed = 42)
p_cal <- plot_boot_calib(cal)
cal$status                 
cal$curve$effective_reps    


p_dca <- plot_dca(cv_result, standardized = TRUE, cost_benefit_axis = TRUE)
```
