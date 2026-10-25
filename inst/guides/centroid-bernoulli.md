# Centroid and Bernoulli classifiers

NearestCentroid uses training class means, Euclidean distance, equal priors and
no shrinkage or automatic standardization. Feature units therefore affect
distance. The two squared distances are mapped through a softmax score;
ties produce 0.5. This score is not a separately calibrated clinical risk.

BernoulliNB uses the fixed rule x > 0 (zero belongs to the absent category),
Laplace smoothing alpha=1 and priors estimated from the training partition.
Both feature presence and absence contribute. Positive-only expression data
can lose almost all quantitative information under this rule; no threshold
search is performed automatically.

Both classifiers use the existing pipeline, CV, final-fit and saved-prediction
interfaces. Saved states contain centroids or probability tables, not a new
data-augmentation procedure.

## Code examples

The following fragments use the input and object names described above. Supply your prepared data and the required engines before running them.

### Example 1

```r
p1 <- sb_pipeline(learner = "nearest_centroid_binary")
p2 <- sb_pipeline(learner = "bernoulli_naive_bayes_binary")
cv <- run_cv(task, p1, folds = 5, seed = 42)
plot_roc(cv)
fit <- fit_final(task, p2)
prediction <- predict_final(fit, new_data, "patient_id")
```
