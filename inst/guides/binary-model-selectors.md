# Binary feature-selection pipelines

Supported extractors can be composed with registered binary prediction methods.
Selection uses the current training partition; supported nested tuning or
calibration repeats the relevant selection inside its analysis partitions.
External prediction applies the saved feature set without selecting again.

| Extractor | Rule | Configuration |
|---|---|---|
| lasso_binary_nonzero | Nonzero Lasso coefficients | inner_folds=5, seed=1, lambda_rule=lambda.1se |
| stepglm_binary_aic | Variables retained by stepwise logistic AIC | direction=both |
| glmboost_binary_nonzero | Nonzero glmBoost coefficients | iterations=100, shrinkage=0.1 |
| rf_binary_importance | Importance strictly above the declared threshold | trees=500, seed=1, permutation importance, threshold=0 |

Lasso also accepts lambda.min; stepwise selection supports forward/backward;
forest importance may use Gini. Declare the rule before evaluation.
An empty selection is reported as failure rather than replaced by an arbitrary
top feature. Some native learners restrict compatible extractors; a registered
method does not imply every possible pairing is supported.

## Code examples

The following fragments use the input and object names described above. Supply your prepared data and the required engines before running them.

### Example 1

```r
library(stormbreaker)

task <- bi_task(data, "patient_id", "outcome", genes)
pipeline <- sb_pipeline("lasso_binary_nonzero", "gbm_binary",
  selector_options = list(lambda_rule = "lambda.min"))

cv <- run_cv(task, pipeline, folds = 5, seed = 42)
fit <- fit_final(task, pipeline)
fit$selection$selected_features  
fit$selection$scores             
cv$selections                    
pred <- predict_final(fit, new_data, "patient_id")


p2 <- sb_pipeline("stepglm_binary_aic", "gaussian_naive_bayes",
  selector_options = list(direction = "forward"))
p3 <- sb_pipeline("rf_binary_importance", "logistic",
  selector_options = list(threshold = 0))
```
