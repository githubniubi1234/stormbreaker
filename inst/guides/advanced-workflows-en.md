# Advanced workflows

For serial execution of `ranger_cindex_survival`, opt into up to four
native fitting threads with
`withr::with_options(list(stormbreaker.ranger_cindex_threads = 4L), your_call)`.
The default remains one thread and statistical settings/IDs do not change.
Saved model parameters record the effective thread count. Predictions matched
exactly in paired tests; native importance sums may differ at rounding precision.
Built-in parallel repeated CV forces one native thread per worker. When using
your own outer parallel code, also leave this option at one to avoid nested
parallelism. Other learners and prediction threading are unaffected.

**Build prediction models from candidate genes and evaluate them in independent cohorts.**

For repeated CV, `sb_compare_repeated_cv(..., workers = 2)` runs independent
repetitions in two local processes; the default `workers = 1` remains serial.
Each repetition retains its complete pipeline group and shared computations.
Results stay in input order, without averaging or selecting a winner. Allow
memory for two workers; this does not speed up a single repetition. The native
BlackBoost-Cox, RSF and AdaBoost learners currently require serial execution
because their caller-RNG stream must not be silently changed. This restriction
does not apply to the locally seeded RSF feature selector.

For expensive selectors, advanced callers can declare Boruta
`selector_options = list(num_trees = 100, maximum_runs = 20)` in
`sb_pipeline()`. These are smaller
experimental specifications, not equivalent replacements for the defaults
(500 trees / 100 maximum Boruta runs).
Choose budgets before evaluation, not from external-cohort scores. Existing
selector reuse shares eligible identical fits. Different
training folds, feature orders, endpoints or resampling settings never share fits.


stormbreaker is an R package for binary, multiclass and survival prediction with
transcriptomic data. Start with a patient table and candidate genes, compare
model specifications, then fit a selected model and evaluate it in new cohorts.

## Start here: binary and survival prediction

Akritas survival prediction is available as `akritas_survival`
(`Akritas` in the model menu). It uses the native `distr6` empirical CDF and
`survivalmodels` conditional product-limit engine with fixed bandwidth 0.5.
The scalar risk is negative step-RMST through maximum analysis follow-up;
horizon probabilities are direct native predictions, not Cox calibration.
Its saved state necessarily contains the training predictor basis and paired
follow-up/event values, but no patient IDs or row names. Handle it as sensitive
training data. `distr6` is available from the author's r-universe repository.
The 0.1.191 survivalmodels release retains the native engine but not its former
public R wrapper; this adapter calls that native engine directly.

High-dimensional Akritas inputs can produce identical native conditional curves.
Final predictions and comparison tables explicitly label exactly constant
risk/probability outputs; valid metrics such as C = 0.5 remain visible.
`sb_select_model()` excludes candidates constant in any validation cohort or
within every development assessment fold. It never consults test/external
predictions for eligibility. This does not change bandwidths, add jitter, or
reject a legitimate single-patient prediction request.

The current native LPC automatic threshold search can return an empty threshold
with fewer than 50 predictors. Such an observed failure is reported as
`stormbreaker_lpc_dimension_unsupported`; no features or thresholds are
fabricated. This is a limitation of the native automatic search, not a universal
lower-dimensional restriction on LPC theory. CV-tuned PLS-Cox retains its
component search and warning checks while using tighter native Cox likelihood
convergence (eps = 1e-11, at most 50 iterations); remaining numerical warnings
still fail the affected route rather than silently lowering its component cap.

`bi_pipeline(c("PLSR", "Logistic"))` uses native logistic PLS (at most two scaled components)
and ranks absolute coefficients multiplied by analysis-fold predictor SD.
It retains the top half, with at least two positive weights when available,
matching the declared coefficient-ranking convention of survival PLSR.
This ranking/subset rule is a STORMBREAKER selector, not native PLS sparsity or
an assertion of a historical paper's exact selection specification. It is
refitted inside glmnet analysis folds and never pads zero-score genes.

`bi_pipeline(c("XGBoost", "Logistic"))` and `surv_selector("XGBoost")` retain features with
positive native tree gain from the corresponding binary-logistic or Cox model
(100 rounds, eta 0.1, depth 3, full sampling, one thread, local seed 1).
They refit in analysis folds, including glmnet inner folds, and may retain all
features. No target count, importance threshold search or fallback is used.

Survival `surv_selector("GamBoost")` selects active native Cox spline terms
(`bbs`, dfbase 4, 100 steps, shrinkage 0.1). It is distinct from linear
`surv_selector("glmBoost")`. Original selected genes, not spline-expanded
columns, are passed downstream; selection is refitted in analysis folds,
including learner-inner folds for glmnet. Diagnostic scores are spline
coefficient-vector norms, not linear effects or calibrated gene importance.
`bi_pipeline(c("GamBoost", "Logistic"))` uses the corresponding native Binomial-logit spline
model with the same 100-step/dfbase-4/shrinkage-0.1 specification and local seed 1.

For routine multi-model comparisons, start with roughly 100 or fewer candidates
with a stated biological or independent-data rationale to limit computational
cost and simplify diagnostics. This is not a software limit or a guarantee of
adequate sample size: feasibility depends on patient and event/class counts,
collinearity and the chosen model. Do not tune the candidate set on external
performance. Outcome-based screening, such as univariable Cox screening, must be
repeated within the relevant training folds for cross-validation. High-dimensional
analyses can use suitable regularization or fold-local selection; the package
does not automatically truncate the input to 100 genes.

```r
library(stormbreaker)
d <- bi_example(100, seed = 42)
r <- bi_run(d, outcome = "outcome", id = "patient_id",
            predictors = c("gene_a", "gene_b"),
            models = c("Logistic", "NaiveBayes"))
finished <- finish_run(r, pipeline = "Logistic")
export_run(finished, "results-new")
```

Use `surv_run()` for survival data; completion and export are shared.
Use `models()` to inspect available algorithms. Omitting `models` requests
all default candidates, including preset screening chains, and can take time.
With automatic plotting (`plots = NULL`), multiple candidates produce comparison
plots only; a single candidate also gets its algorithm-specific plots.
`finish_run()` fits the chosen candidate without repeating comparison.
Explicit `plots = TRUE` requests all candidate plots; `FALSE` skips plotting.

Computations do not open a browser. Existing HTML helpers remain optional.
See the [getting-started guide](getting-started.md).
Multiclass is supplementary and has a separate interface; arbitrary multiclass
selector-plus-learner composition is not provided.

## Compare models on shared folds

Choose models from a task-specific menu or combine a feature selector with a
learner. Models use the same cross-validation splits, with selection and tuning
performed within the training folds. Patient-level out-of-fold predictions are
retained for discrimination, calibration and decision-curve analysis.

A standalone results page shows selectors and learners as rows, with development,
validation and test cohorts in separate metric columns. Opening it is optional.

### Execution design: reuse computation, preserve the analysis

Multi-pipeline comparison separates a statistical specification from its execution.
Eligible pipelines can share identical fold-local selection steps and supported
selection-chain prefixes, rather than fitting the same selector for every learner.
Reuse requires matching ordered predictors, sample identities, endpoints, method
settings and random-number configuration; matching a selector name is not enough.
Outer analysis folds, inner analysis folds and final development data remain separate.

Within a scoped reuse context, eligible identical CV and final learner fits can
also be shared, including their complete analysis-local searches. Each pipeline retains its own identity, diagnostics, provenance and
prediction outputs. Validated prediction paths avoid recomputing a native survival
curve separately for risk and horizon survival probability. Pipeline construction
also passes one local registry table to its helpers rather than rebuilding it.

These optimizations do not reduce candidates, trees, iterations or CV folds. Errors
are not reused as successful fits, and native warnings remain visible. Reuse is
in-memory by default and ends with its scope; it is not a persistent cross-run model cache.
Bounded research controllers can additionally hand off completed Boruta selections
between workers. These sensitive job-local checkpoints
contain matching analysis inputs and endpoints, not just portable prediction state;
they must not be published as model artifacts. Consumers still match complete inputs.
CV and final learner entries share a 256 MiB retention budget. Least-recently-used
entries are released as needed; a model larger than the budget is returned but
not retained. Released models are fitted normally if needed again. This limits
retained learner entries, not total process memory, active training, returned
results or selector/preparation storage. Reuse diagnostics report bytes,
evictions and oversized entries; eligibility does not guarantee a cache hit.
The design is validated against unshared execution using selected features, fold
assignments, patient-level predictions, metrics, diagnostics and failure behavior.
Elapsed-time improvements are reported for the measured workload, not assumed for
every method. Smaller experimental training budgets are a separate specification,
not an equivalent acceleration.

See the [execution-design guide](selector-reuse-and-ranking.md)
for supported reuse boundaries, controls and measured examples.

## Installation

See [the current installation guide](installation.md).

## Advanced: explicit task interface

These examples use synthetic expression data bundled with the package.
Replace the data and `genes` with your own patient table and gene set.

```r
library(stormbreaker)
genes <- c("gene_a", "gene_b", "gene_c")

binary_data <- bi_example(n = 100, seed = 42)
binary_task <- bi_task(
  binary_data, "patient_id", "outcome", predictors = genes
)
binary_result <- sb_compare_models(
  binary_task, models = c("logistic", "gaussian_naive_bayes"),
  folds = 5, seed = 42
)
binary_result$comparison[c("model_id", "status", "roc_auc", "brier_score")]

survival_data <- surv_example(n = 100, seed = 42)
survival_task <- surv_task(
  survival_data, "patient_id", "time", "event",
  predictors = genes, time_unit = "month"
)
survival_result <- sb_compare_models(
  survival_task, models = c("coxph", "survreg_weibull"),
  folds = 5, seed = 42
)
survival_result$comparison[c("model_id", "status", "harrell_c")]

```

Browse the model menu with `sb_model_catalog("binary")` or
`sb_model_catalog("survival")`. In this advanced `sb_compare_models()` interface, omitted `models` runs three fixed starting
models. Explicitly passing `models = NULL` to
`sb_compare_models()` runs the full task-specific menu. The catalog lists
engine availability; the comparison retains a status for each requested model.

## Explore the comparison

Create an editable performance plot:

```r
# install.packages("ggplot2")
plot_perf(binary_result, metric = "roc_auc")
```

Browser opening is off by default. Request `view = TRUE` explicitly,
or open a saved result:

```r
sb_view_results(binary_result)
```

The page uses standard algorithm names, with settings and internal IDs in details.
It supports metric selection, search, dataset visibility, CSV and SVG export.

## Compare several cohorts

Supply one development task and named independent cohort tasks:

```r
# development_task, validation_task and test_task are task objects for your data.
result <- sb_compare_cohorts(
  development_task,
  cohorts = list(Validation = validation_task, Test = test_task),
  roles = c(Validation = "validation", Test = "test"),
  models = c("logistic", "gaussian_naive_bayes"),
  folds = 5
)
```

Every pipeline is fitted on development patients only, then evaluated in every
supplied cohort. Development columns show out-of-fold performance, not training
fit performance. Default sorting uses development CV; optional validation means
exclude test cohorts. Patient IDs must identify distinct people across cohorts.
See the [runnable example](../examples/multi-cohort.R) and
[multi-cohort guide](multi-cohort-results.md).

## Fit a final model

Choose the specification for your study, fit it on the development data,
and supply the new patients' predictor columns:

```r
final_model <- fit_final(binary_task, sb_pipeline(learner = "logistic"))
new_patients <- bi_example(n = 40, seed = 43)
predictions <- predict_final(
  final_model, new_patients[c("patient_id", genes)], id_col = "patient_id"
)
predictions$predictions
```

Use `bi_endpoint()` or `surv_endpoint()` and
`eval_external()` to evaluate saved predictions against independent
outcomes. The [complete workflow](../examples/minimum-core.R) covers both tasks.

## Learn more

- [Saved coefficients, existing SHAP and survival-curve explanations](model-explanation.md)
- [Choosing models and gene sets](user-models.md)
- [Development and external-validation workflow](workflow.md)
- [Shared-execution design and the leaderboard](selector-reuse-and-ranking.md)
- Package tutorials: `browseVignettes("stormbreaker")` after installing a built Release package.
- [Changes in this release](../../NEWS.md)

## Help and citation

Share a small synthetic reproducer with the maintainer. After publication, use the issue tracker listed in DESCRIPTION.
Run `citation("stormbreaker")` for the software citation.
Contributions are welcome; see [CONTRIBUTING.md](../../CONTRIBUTING.md).

## License

MIT. See [LICENSE.md](../../LICENSE.md).
# PAMR feature selection

`bi_pipeline(c("PAMR", "Logistic"))` and `surv_selector("PAMR")` retain the native
nonzero shrunken-centroid genes at threshold 1. Binary selection uses empirical
class priors; survival selection uses PAMR's two censoring-aware groups.
No threshold search or gene padding is performed. Selection is repeated within
each analysis fold, including glmnet inner analysis folds. Saved selector state
contains feature-level scores and selected names, not patient rows or the engine.


## Unified binary modeling

The current binary menu has 104 method/configuration entries: 92 can predict,
23 can extract features, and 11 have both capabilities. These are not 104
independent algorithm families. Removed REML GAM and binary SLOPE learners are
not restored by the unified menu. Inspect `can_predict`, `can_extract_features`
and `stage_selector` for actual capabilities.

```r
r <- bi_run(task, models = c("Lasso", "RF", "GBM"),
            stages = c(1, 2, 3), folds = 5, plots = FALSE)
finished <- finish_run(r, pipeline = "Lasso + RF + GBM",
                            external = list(Validation = validation))
export_run(finished, "results-new")
```

Within the chosen pool, one stage fits a predictor, two stages extract features
then predict, and three stages use two ordered extractors then a predictor.
Repetitions are allowed. The three methods above construct 39 routes; fitting
success still depends on data and engines. The whole menu defaults to stages
1 and 2, never implicit three-stage expansion. An explicit method list defaults
to single-stage models. Feature-only methods cannot terminate a route.

`bi_selector()` is internal compatibility code, not an exported entrypoint.
`analyze_methods()` remains an advanced analysis interface. Survival retains
its existing `surv_pipeline()` / `surv_selector()` syntax and fixed three-stage
presets; multiclass remains separate. This release does not extend them.

For the six binary learners with a formal-run constant-sum repair (brglm2, LDA,
HDA, FDA, MDA and local LDA), training-only Helmert coordinates remove only the
known constant-sum direction. Fitted states retain the transform and check its
assessment-time applicability. Ordinary input representations and old native
states remain unchanged. This is not arbitrary PCA or new feature selection.
See [source and compatibility notes](../../NEWS.md).
