# Shared execution and result ranking

STORMBREAKER separates a declared statistical pipeline from reusable execution
steps. Equivalent fold-local selections, supported CV/final fits and selected
survival-curve predictions can be shared without merging pipeline identities.

## Reuse boundaries

Reuse requires matching ordered inputs, patient identities, outcomes, settings,
seeds and relevant RNG/diagnostic context. A shared method name alone is not
enough. Different outer/inner analysis partitions remain separate. Inner tuning
must not reuse a selection learned from an outer or complete development set.

Supported ordered selection chains can share identical prefixes while distinct
downstream stages still execute separately. Registered learner reuse additionally
matches tuning context and, where required, the starting and resulting random
state. Failed or warning-producing fits are not treated as successful reusable
entries. Saved diagnostics and provenance remain pipeline-specific.

Some verified survival paths share one native prediction calculation between
risk and horizon probabilities; this is not assumed for every survival learner.
The context is in memory and released at scope exit. CV/final learner entries
use a 256 MiB retention budget with least-recently-used eviction; this is not a
limit on total process memory, active training or returned results.

## Use and inspect

compare_cv() defaults to reuse_selectors=TRUE. Set it to FALSE for a paired
comparison. sb_with_selector_reuse() gives an explicitly scoped batch and can
return reuse diagnostics. Separate processes do not automatically share their
in-memory contexts. Keep scope sizes appropriate for available memory.

Reuse does not reduce candidate routes, trees, bootstrap draws, iterations or
fold counts. Verify numerical equivalence using selected features/order, splits,
settings, patient-level predictions, metrics, diagnostics and random behavior.
Elapsed-time gains depend on how many computations are genuinely shared; no
universal speedup is implied. Historical project-specific benchmarks are not
part of this user guide.

The result viewer ranks saved development/validation summaries. Sorting does
not select, refit or merge models. See [multi-cohort results](multi-cohort-results.md).

```r
result <- compare_cv(task, pipelines, folds = 5, reuse_selectors = FALSE)
measured <- sb_with_selector_reuse({
  cv_results <- lapply(pipelines, function(p) run_cv(task, p, folds = 5))
  final_fits <- lapply(pipelines, function(p) fit_final(task, p))
  list(cv = cv_results, fits = final_fits)
}, diagnostics = TRUE)
measured$reuse
```
