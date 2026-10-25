# Diagnostic and plotting short names

Short names are aliases of the original public functions: they do not define
additional algorithms or change model settings. Use a qualified name such as
stormbreaker::fit_final() if another package has an overlapping name.

The complete mapping is in [specialty-name-map.csv](specialty-name-map.csv). The CSV records each short
name, original name, function group and English description. Native diagnostic
availability still depends on the fitted model and retained state.

Use plot_diag() for a diagnostic object, plots(fit) for the overall algorithm
view, and clinical_plots() for clinical evaluation views. A selector's diagnostic
must be requested from that stage rather than confused with the final learner.

## Code examples

The following fragments use the input and object names described above. Supply your prepared data and the required engines before running them.

### Example 1

```r
g <- plots(fit)
names(g$plots)
g$unavailable

d <- forest_diag(fit, type = "importance")
plot_diag(d)
```

### Example 2

```r
fit <- fit_final(task, bi_pipeline("Lasso"), capture_diagnostics = TRUE)
cv <- path_diag(fit, type = "cv_curve")
plot_diag(cv)
```
