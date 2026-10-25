# Inspect algorithm diagnostics

plots(fit) is the general entrypoint for available native diagnostics. Capture
diagnostics during fitting when a method requires them; plotting does not refit
a missing coefficient path or generate unsupported evidence.

Use path_diag() for saved regularization paths and CV, forest_diag() for supported
importance/tree/OOB evidence, pamr_diag() for shrunken centroids, mclust_diag()
for component parameters, and svm_diag()/pls_diag()/nb_diag() for their native
structures. native_diag(fit, complete=TRUE)$status lists availability.

Selector evidence and final-learner evidence are separate. A Logistic final
predictor after Lasso selection is not itself a Lasso model. Request diagnostics
from the corresponding stage. SHAP is an additional explanation calculation,
not a plot automatically supplied by every algorithm.

The full short-name table is in [specialty-short-names](specialty-short-names.md).

## Code examples

The following fragments use the input and object names described above. Supply your prepared data and the required engines before running them.

### Example 1

```r
g <- plots(fit)
names(g$plots)
g$plots[[1]]
g$unavailable
```

### Example 2

```r
g <- plots(fit, new_data = query, id_col = "patient_id",
              background_data = development_features)
```
